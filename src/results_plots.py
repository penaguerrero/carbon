import numpy
import os
import string
import emcee
import corner
import matplotlib
import scipy.optimize as op
from scipy import stats
from matplotlib import pyplot as plt
from copy import deepcopy
from scipy.stats._stats_mstats_common import linregress

############################################################################################################################################
'''
This scripts takes the found metallicities and produces:
 * plot of C/O vs Otot
 * plot of N/O vs Otot
 * plot of Ne/O vs Otot
 
 NOTE: The data for N/O and Ne/O was taken all from Lopez-Sanchez & Esteban 2010
'''

use_our_sample_ONLY = True
save_images = False
# Do you want to make the C/O vs O/H with points from literature with respect to solar: [C/O] vs [O/H]?
solar = True
# Do you want to correct values?
correct_values = False
# Type of image to be saved?
typeofimage = '.jpg'
# Show fitted line?
show_fitted_line = False
# Use results from measurements with C_Hbeta?
C_Hbeta = True

############################################################################################################################################

img_name1 = 'COvsOH'
img_name1a = 'COvsOH_literature'
img_name1b = 'COvsOH_Cldy'
img_name2 = 'NOvsOH'
img_name3 = 'NeOvsOH'
img_name4 = 'CNvsOH'
img_name5 = 'C+2N+2vsOH'
img_name6 = 'NCvsCH_v2'
img_name7 = 'NOvsCN_v2'
otherrefs = '_plusOtherRefs'


#### FUNCTIONS

def fit_line(arrx, arry):
    '''
    This function fits a line to the given array.
    arrx, arry = numpy arrays 
    RETURNS:
    - constants of the line
    - fitted line
    '''
    order = 1
    polyinfo = numpy.polyfit(arrx, arry, order, full=True)
    print 'polyinfo: ', polyinfo
    # polyinfo = coefficients array, residuals, rank, singular values array, and conditioning threshold
    coefficients = polyinfo[0]
    polynomial = numpy.poly1d(coefficients)
    f_pol = polynomial(arrx)
    #fitted_line = np.array([arry, f_pol])
    #print 'this is x and y of the fitted_line = ', fitted_line
    # Errors -- taken from http://mathworld.wolfram.com/LeastSquaresFitting.html
    avgx = sum(arrx)/float(len(arrx))
    avgy = sum(arry)/float(len(arry))
    ssx, ssy, ssxy = [], [], []
    for xi, yi in zip(arrx, arry):
        ssxi = (xi- avgx)**2 
        ssyi = (yi- avgy)**2 
        ssxyi = (xi - avgx)*(yi - avgy)
        ssx.append(ssxi)
        ssy.append(ssyi)
        ssxy.append(ssxyi)
    ssx = sum(ssx)
    ssy = sum(ssy)
    ssxy = sum(ssxy)
    n = float(len(arrx))
    s = numpy.sqrt((ssy - (ssxy**2 / ssx))/(n-2.0))
    b_err = s * numpy.sqrt(1.0/n + avgx**2/ssx)
    m_err = s / numpy.sqrt(ssx)
    errs = [m_err, b_err]
    return coefficients, errs, f_pol

# the true line
def line_eq(theta, x):
    m, b = theta
    return m * x + b

# likelihood function
def lnlike(theta, xobs, yobs, yerrobs):
    #print 'len(theta), theta', len(theta), theta
    #print 'len(xobs), len(yobs), len(yerrobs), xobs, yobs, yerrobs', len(xobs), len(yobs), len(yerrobs), xobs, yobs, yerrobs
    model = line_eq(theta, xobs)
    #print 'type(yobs), type(model), type(yerrobs)', type(yobs), type(model), type(yerrobs)
    chi2 = (yobs - model)**2 / yerrobs**2
    chi2 = chi2.sum()
    return - chi2/2.0   # we are returning the log of the likelihood function

def lntophat(x, a, b):
    # b has to be grater than a
    if a > b:
        bb = a
        aa = b
    elif a == b:
        print 'Boundaries are the same in lntophat function... change them please.'
        exit()
    else:
        aa = a
        bb = b
    if (aa < x) and (x < bb):
        return numpy.log(1/(bb-aa))
    else:
        return -numpy.inf

# define the priors
def lnprior(theta):
    m, b = theta
    pm = lntophat(m, -0.5, -0.01)     
    pb = lntophat(b, 0.1, 3.5)     
    if pm != -numpy.inf and pb != -numpy.inf:
        return pm + pb
    return -numpy.inf

# then the probability function will be
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not numpy.isfinite(lp):
        return -numpy.inf
    return lp + lnlike(theta, x, y, yerr)

def approxC(N, O):
    ''' This function uses equation of Carbon paper to determine C from N and O arrays.'''
    #return 0.5 - 0.18 * (N - O) + N   # previous
    #return 0.35 - 0.30 * (N - O) + N   # MCMC
    #return 0.30 - 0.18 * (N - O) + N   # MCMC without Sun, IZw18, Orion, and 30 Dor, deriving C/N
    #return N - 0.297 - 0.223*(N - O)  # MCMC without Sun, IZw18, Orion, and 30 Dor, deriving N/C
    #return 0.31*(N - O) + 0.38 + N  # MCMC without Sun, IZw18, Orion, and 30 Dor, average of all lines
    #return 0.8*O - 5.96 + N  # derived from the C/N vs O/H diagram
    return 0.31 - 0.22*(N - O) + N  # MCMC without Sun, IZw18, Orion, and 30 Dor, deriving C/N

def err_approxC(N, Nerr, O, Oerr):
    ''' This function uses equation of Carbon paper to determine the error in C from N and O arrays.'''
    e1 = numpy.sqrt(0.18 * Nerr**2 + 0.18 * Oerr**2)
    err = numpy.sqrt( e1**2 + Nerr**2 )
    return err

def convert2solar(X, Xsun):
    ''' This function returns values in terms of solar values.'''
    return X - Xsun


############################################################################################################################################

#### CODE

# Go into the directory of each object. Files assumed to be in: /Users/home_direcotry/Documents/AptanaStudio3/src/
results_path = "../results"
# but just to make sure we are in the right place, get and work with full paths
full_results_path = os.path.abspath(results_path)
# read the results file
if use_our_sample_ONLY:
    resfile = os.path.join(full_results_path, 'abunds_OCNNe.txt')
    if C_Hbeta:
        resfile = os.path.join(full_results_path, 'abunds_OCNNe_CHbeta.txt')
    rows2skip = 5
    rf = numpy.loadtxt(resfile, skiprows=rows2skip, usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18), unpack=True)
    OH, OHerr, CO, COerr, NO, NOerr, NeO, NeOerr, C2O2, C2O2err, XO2, O2O1, O2O1err, N2O2, N2O2err, Z, ifc, ifcerr = rf
    solarZ = 0.020
    Z_respect2solar = []
    #print 'Values with respect to solar metallicity:'
    for eachZ in Z:
        Z2solar = eachZ / solarZ
        Z_respect2solar.append(Z2solar)
        #print '%1.3f  %1.3f' % (eachZ, Z2solar)
else:
    resfile = os.path.join(full_results_path, 'abunds_OCNNe_plusotherrefs.txt')
    if C_Hbeta:
        resfile = os.path.join(full_results_path, 'abunds_OCNNe_plusotherrefs_CHbeta.txt')
    rows2skip = 10
    rf = numpy.loadtxt(resfile, skiprows=rows2skip, usecols=(1,2,3,4,5,6,7,8), unpack=True)
    OH, OHerr, CO, COerr, NO, NOerr, NeO, NeOerr = rf
        
# These are the results are for the sample objects in the following order:
ordered_sample = ['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
                  'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
# now divide the sample in objects for in which we found Te and those in which we did not
objects_with_Te = ['iiizw107', 'mrk5', 'pox4', 'sbs0948', 'sbs1054', 'tol1457', 'iras08208', 'mrk1087', 
                   'mrk1199', 'mrk960', 'ngc1741', 'sbs0218', 'sbs0926', 'sbs1319', 'arp252', 'sbs1415']
objects_Te_literature = ['iras08339', 'tol9']

# objects not in the sample
objects_ref = ['30Dor', 'Orion', 'izw18', 'Sun'] # these are useful points of reference
# additional data points not in the sample
objects_not_in_sample = ['ngc346', 'ngc456', 'ngc6822']

CN = (CO+OH) - (NO+OH)   # though OH cancels out, I left it explicit to better trace the error propagation 
CNerr = []
for coe, noe, oe in zip(COerr, NOerr, OHerr):
    cne = numpy.sqrt(coe**2 + noe**2 + 2*oe**2)
    CNerr.append(cne)
CNerr = numpy.array(CNerr)

if correct_values:
    O_dust_correction = 0.1
    C_dust_correction = 0.2
    img_name1 = 'COvsOH_corr'
    img_name2 = 'NOvsOH_corr'
    img_name3 = 'NeOvsOH_corr'
    img_name4 = 'CNvsOH_corr'
    OH_non_corr = deepcopy(OH)
    CO_non_corr = deepcopy(CO)
    NO_non_corr = deepcopy(NO)
    NeO_non_corr = deepcopy(NeO)
    #NeOor_non_corr = deepcopy(NeOor) 
    print '* The format is: Object    Non-corrected values   Corrected values'
    print '{:<10} {:>10} {:>12} {:>9} {:>6}'.format('', 'O/H', 'C/O', 'Z', '[Z]')
    for i in range(len(ordered_sample)):
        OH[i] = OH[i] + O_dust_correction
        CO[i] = CO[i] + C_dust_correction - O_dust_correction
        NO[i] = NO[i] + C_dust_correction - O_dust_correction
        NeO[i] = NeO[i] + C_dust_correction - O_dust_correction
        #NeOor[i] = NeOor[i] + C_dust_correction - O_dust_correction
        if use_our_sample_ONLY:
            print '{:<10} {:>8.2f} {:>5.2f} {:>6.2f} {:>5.2f} {:>6.3f} {:>5.3f}'.format(ordered_sample[i], OH_non_corr[i], OH[i], CO_non_corr[i], CO[i], Z[i], Z_respect2solar[i])
        else:
            print '{:<10} {:>8.2f} {:>5.2f} {:>6.2f} {:>5.2f}'.format(ordered_sample[i], OH_non_corr[i], OH[i], CO_non_corr[i], CO[i])
    # Now add the correction to the reference objects in the same order as in the text file EXCEPT for the Sun
    for j in range(len(objects_ref)-1):
        OH[i+j+1] = OH[i+j+1] + O_dust_correction
        CO[i+j+1] = CO[i+j+1] + C_dust_correction - O_dust_correction
        NO[i+j+1] = NO[i+j+1] + C_dust_correction - O_dust_correction
        NeO[i+j+1] = NeO[i+j+1] + C_dust_correction - O_dust_correction
        #NeOor[i+j+1] = NeOor[i+j+1] + C_dust_correction - O_dust_correction
    # Do not add anything for the Sun
    OH[i+j+1] = OH[i+j+1]
    CO[i+j+1] = CO[i+j+1]
    NO[i+j+1] = NO[i+j+1]
    NeO[i+j+1] = NeO[i+j+1]
    #NeOor[i+j+1] = NeOor[i+j+1] 
    # Now add the correction to the objects not in the sample in the same order as in the text file
    for k in range(len(objects_not_in_sample)-1):
        OH[i+j+k+1] = OH[i+j+k+1] + O_dust_correction
        CO[i+j+k+1] = CO[i+j+k+1] + C_dust_correction - O_dust_correction
        NO[i+j+k+1] = NO[i+j+k+1] + C_dust_correction - O_dust_correction
        NeO[i+j+k+1] = NeO[i+j+k+1] + C_dust_correction - O_dust_correction
        #NeOor[i+j+k+1] = NeOor[i+j+k+1] + C_dust_correction - O_dust_correction

# Compile all the names for the labeling in the same order
objects_list = []
indeces_list = []
i = 0
for obj in ordered_sample:
    objects_list.append(obj)
    indeces_list.append(i)
    #print i, obj
    i = i+1
for obj in objects_ref:
    objects_list.append(obj)
    indeces_list.append(i)
    #print i, obj
    i = i+1
if use_our_sample_ONLY == False:
    for obj in objects_not_in_sample:
        objects_list.append(obj)
        indeces_list.append(i)
        #print i, obj
        i = i+1

# Abundances obtained with Cloudy simulations
no_errs = False
if no_errs:
    Cldy_best_abunds = 'carbon/results/Cloudy_best_abunds.txt'
else:
    Cldy_best_abunds = 'carbon/results/Cloudy_best_abunds_witherrs.txt'    
abunds_path_list = string.split(os.getcwd(), sep='carbon')
fname = os.path.join(abunds_path_list[0], Cldy_best_abunds)
Clhe, errClheu, errClhed, Clo, errClou, errClod, Clc, errClcu, errClcd, Cln, errClnu, errClnd, Clne, errClneu, errClned, Cls, errClsu, errClsd, Clto3, errClto3u, errClto3d, Clto2, errClto2u, errClto2d = numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([])
if no_errs:
    Clhe, Clo, Clc, Cln, Clne, Cls, Clto3, Clto2 = numpy.loadtxt(fname, skiprows=1, usecols=(2,3,4,5,6,7,8,9), unpack=True)
else:
    Clhe, errClheu, errClhed, Clo, errClou, errClod, Clc, errClcu, errClcd, Cln, errClnu, errClnd, Clne, errClneu, errClned, Cls, errClsu, errClsd, Clto3, errClto3u, errClto3d, Clto2, errClto2u, errClto2d = numpy.loadtxt(fname, skiprows=1, usecols=(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25), unpack=True)


# PLOTS
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}
matplotlib.rc('font', **font)

# C/O vs O/H from OBSERVATIONS
fig1 = plt.figure(1, figsize=(12, 10))
#plt.errorbar(OH, CO, xerr=OHerr, yerr=COerr, fmt='ko')    # this plots ALL points with the same symbol
for obj, i in zip(objects_list, indeces_list):
    if obj in objects_with_Te:
        fmt='ko'
        markersize=8
    elif obj in objects_Te_literature:
        fmt='k*'
        markersize=13
    elif (obj in objects_not_in_sample) or (obj in objects_ref):
        fmt='b^'
        markersize=9
    plt.errorbar(OH[i], CO[i], xerr=OHerr[i], yerr=COerr[i], fmt=fmt, ecolor='k', markersize=markersize)
plt.xlim(6.8, 9.0)
yup = 1.1
ylo = -1.9
plt.ylim(ylo, yup)
plt.xticks(numpy.arange(6.8, 9.1, 0.2))
plt.yticks(numpy.arange(-1.8, yup, 0.2))
plt.minorticks_on()
for x, xe, y, ye, z in zip(OH, OHerr, CO, COerr, objects_list):
    subxcoord = -3
    subycoord = 5
    side = 'right'
    if (z == 'sbs0948') or (z == 'mrk5') or (z == 'Sun'):
        subxcoord = -4
        subycoord = -12
    if (z == 'sbs0218') or (z == 'iras08208') or (z == 'mrk1087') or (z == 'sbs1319') or (z == 'ngc456') or (z == 'ngc1741') or (z == 'ngc6822') or (z == 'Orion') or (z == 'izw18'):
        subxcoord = 5
        side = 'left'
    if (z == 'ngc346') or (z == 'iiizw107') or (z == 'sbs1054') or (z == 'pox4') or (z == '30Dor'):
        subxcoord = 4
        subycoord = -11
        side = 'left'
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
plt.title('C/O vs O/H')
plt.xlabel('12 + log (O/H)')
plt.ylabel('log (C/O)')
# Linear fit to linear part of diagram -- 2 different methods for this
cooh_coeffs, errs, cooh_line_fit = fit_line(OH[:-4], CO[:-4])   # 1
print '* Coefficients of linear part of the plot: C/O vs O/H'
m = cooh_coeffs[0]
b = cooh_coeffs[1]
print '  m = %0.2f     b = %0.2f      errs: ' % (m, b), errs
m, b, r_value, p_value, err = stats.linregress(OH[:-4], CO[:-4])  # 2
print 'm, b, r_value, p_value, err: ', m, b, r_value, p_value, err
plt.plot(OH[:-4], cooh_line_fit, 'g', lw=2)   # fit to linear part of diagram
plt.plot(OH[:-4], cooh_line_fit+0.35, 'g', lw=1, alpha=0.4)   # fit to linear plus error
plt.plot(OH[:-4], cooh_line_fit-0.35, 'g', lw=1, alpha=0.4)   # fit to linear minus error
# Garnett (1995) line
OHgar = numpy.array([7.3, 7.5, 7.7, 7.9, 8.1, 8.3, 8.5, 8.7])
COgar = 1.01 + 0.43 * (OHgar-12)
plt.plot(OHgar, COgar, 'm--', lw=1)   # fit of Garnett (1995)
if save_images:
    img_name = img_name1 + typeofimage
    if use_our_sample_ONLY == False:
        img_name = img_name1 + otherrefs + typeofimage
    destination = os.path.join(full_results_path+'/plots', img_name)
    plt.savefig(destination)
    print('Plot %s was saved!' % destination)
plt.show()

# Re-plot our sample PLUS the objects from literature
fig1 = plt.figure(1, figsize=(12, 10))
ax1 = fig1.add_subplot(111)
if solar:
    yup, ylo = 1.0, -2.0
    xmin, xmax = -3.0, 0.6
    # plot guide lines
    plt.vlines(0.0, ylo, yup, colors='k', linestyles='dashed')
    plt.hlines(0.0, xmin, xmax, colors='k', linestyles='dashed')
    plt.title('[C/O] vs [O/H]')
    plt.xlabel('[O/H]')
    plt.ylabel('[C/O]')
    img_name = img_name1a + '_inSolar' + typeofimage
else:
    yup, ylo = 0.8, -1.8
    xmin, xmax = 6.0, 8.8
    plt.title('C/O vs O/H')
    plt.xlabel('12+log(O/H)')
    plt.ylabel('log(C/O)')
    img_name = img_name1a + typeofimage
    plt.xticks(numpy.arange(xmin, xmax, 0.2))
    plt.yticks(numpy.arange(ylo, yup+0.01, 0.2))
plt.ylim(ylo, yup)
plt.xlim(xmin, xmax)
plt.minorticks_on()
# PROTOSOLAR values from Asplund et al. (2009)
COsun = -0.26 # because C/H = 8.47+-0.05
COsunerr = 0.07
OHsun = 8.73 
OHsunerr = 0.05
# OUR STIS OBSERVATIONS
COsol_arr, OHsol_arr = [], []
for obj, i in zip(objects_list, indeces_list):
    if obj in objects_with_Te:
        fmt='ko'
        markersize=8
    elif obj in objects_Te_literature:
        fmt='k*'
        markersize=13
    elif (obj in objects_not_in_sample) or (obj in objects_ref):
        fmt='b^'
        markersize=9
    #plt.errorbar(OH[i], CO[i], xerr=OHerr[i], yerr=COerr[i], fmt=fmt, ecolor='k')   # with errors
    if solar:
        COsol = convert2solar(CO[i], COsun)
        OHsol = convert2solar(OH[i], OHsun)
        COsol_arr.append(COsol)
        OHsol_arr.append(OHsol)
        plt.plot(OHsol, COsol, fmt, ms=7)   # without errors AND in terms of solar values
        if obj=='izw18':
            OHsolerr = numpy.sqrt(OHerr[i]**2 + OHsunerr**2)
            COsolerr = numpy.sqrt(COerr[i]**2 + COsunerr**2)
            plt.errorbar(OHsol, COsol, xerr=OHsolerr, yerr=COsolerr, fmt=fmt, ecolor='k', markersize=markersize)   # with errors
    else:
        plt.plot(OH[i], CO[i], fmt, ms=7)   # without errors
# Linear fit to linear part of diagram -- 2 different methods for this
coeffs, errs, line_fit = fit_line(OHsol_arr[:-4], COsol_arr[:-4])   # 1
plt.plot(OHsol_arr[:-4], line_fit, 'k', lw=1)   # fit to linear part of diagram
plt.plot(OHsol_arr[:-4], line_fit+0.35, 'k', lw=1, alpha=0.4)   # fit to linear plus error
plt.plot(OHsol_arr[:-4], line_fit-0.35, 'k', lw=1, alpha=0.4)   # fit to linear minus error
print 'In solar units, m, b, errs: ', m, b, errs
m, b, r_value, p_value, err = stats.linregress(OH[:-4], CO[:-4])  # 2
print 'm, b, r_value, p_value, err: ', m, b, r_value, p_value, err
# Garnett (1995) linear fit
COgar = convert2solar(COgar, COsun)
OHgar = convert2solar(OHgar, OHsun)
plt.plot(OHgar, COgar, 'k-.', lw=1)   # fit of Garnett (1995)
# Insert labels of type of object responsible for C
#plt.annotate('Massive stars and/or additional source', xy=(1.21, -0.21), xycoords='axes fraction', xytext=(0.1, 0.75))
#plt.annotate('Massive stars', xy=(1.21, -0.21), xycoords='axes fraction', xytext=(0.43, 0.07))
#plt.annotate('Low- and intermediate-mass stars', xy=(1.21, -0.21), xycoords='axes fraction', xytext=(0.6, 0.85))
# Data from Cunha & Lambert (1994)
_, clam, clamerr, nlam, nlamerr, olam, olamerr = numpy.loadtxt('../results/Cunha94.txt', comments='#', unpack=True)
colam = clam-olam
# Data from Nieva & Simon-Diaz (2011)
_, nievc, nievcerr, nievn, nievnerr, nievo = numpy.loadtxt('../results/Nieva2011.txt', comments='#', unpack=True)
nievco = nievc-nievo
# Data from Kilian (1992)
kilc, kilcerr, kiln, kilnerr, kilo, kiloerr = numpy.loadtxt('../results/Kilian1992.txt', comments='#', usecols=(1,2,3,4,5,6),unpack=True)
kilco = kilc-kilo
# Data from Daflon
_, ref, dafc, dafcerr, dafn, dafnerr, dafo, dafoerr = numpy.loadtxt('../results/Daflon.txt', comments='#', unpack=True)
dafco = dafc - dafo
# Data from Morel et al. (2008)
_, morc, morcerr, morn, mornerr, moro, moroerr = numpy.loadtxt('../results/Morel2008.txt', comments='#', unpack=True)
morco = morc - moro
if solar:
    scolam = convert2solar(colam, COsun)
    solam = convert2solar(olam, OHsun)
    plt.plot(scolam, solam, 'gx')
    snievco = convert2solar(nievco, COsun)
    snievo = convert2solar(nievo, OHsun)
    plt.plot(snievco, snievo, 'g+')
    skilco = convert2solar(kilco, COsun)
    skilo = convert2solar(kilo+12, OHsun)
    plt.plot(skilco, skilo, 'g.', alpha=0.4)
    sdafco = convert2solar(dafco, COsun)
    sdafo = convert2solar(dafo, OHsun)
    plt.plot(sdafco, sdafo, 'g<', alpha=0.4)
    smorco = convert2solar(morco, COsun)
    smoro = convert2solar(moro, OHsun)
    plt.plot(smorco, smoro, 'gs', alpha=0.4)
else:
    plt.plot(colam, olam, 'gx')
    plt.plot(nievco, nievo, 'g+')
    plt.plot(kilco, kilo, 'g.', alpha=0.4)
    plt.plot(dafco, dafo, 'g<', alpha=0.4)
    plt.plot(morco, moro, 'gs', alpha=0.4)
# James et al. (2015)
jam15 = '../results/James15.txt'
jamOH, jamOHerr, jamNO, jamNOerr = numpy.loadtxt(jam15, skiprows=3, usecols=(2,3,4,5), unpack=True)
jamN = jamNO + jamOH
jamCH = approxC(jamN, jamOH)
jamCO = jamCH - jamOH
jamNerr = numpy.sqrt(jamNOerr**2 + jamOHerr**2)
jamCOerr = numpy.sqrt(jamNerr**2 + jamOHerr**2)
#plt.errorbar(jamOH, jamCO, xerr=jamOHerr, yerr=jamCOerr, fmt='g*')   # with errors
if solar:
    jamCOsol = convert2solar(jamCO, COsun)
    jamOHsol = convert2solar(jamOH, OHsun)
    plt.plot(jamOHsol, jamCOsol, 'g*', ms=10)   # in solar terms
else:
    plt.plot(jamOH, jamCO, 'g*', ms=10)   
# Other literature points
lit = '../results/literatureCandO.txt'
ref, litCH, litCHerr, litOH, litOHerr = numpy.loadtxt(lit, comments='#', usecols=(1,4,5,6,7), unpack=True)
litCO = litCH - litOH
litnames = ['IZw18', 'SBS0335-052', 'SBS1415+437', 'NGC5253-1', 'NGC5253-2', 'NGC4670', 'SDSSJ0035-0918', 'SDSSJ1111+1332', 
           'Q0913+072', 'SDSSJ1016+4040', 'SDSSJ1558+4053', 'Q2206-199', 'SDSSJ0137-4224', 'SDSSJ2155+1358']
for ln, lo, lc, lco in zip(litnames, litOH, litCH, litCO):
    print '{:<20} {:<6} {:>6} {:>6}'.format(ln, lo, lc, lco)
litCOerr = numpy.sqrt(litCHerr**2 + litOHerr**2)
#plt.errorbar(litOH, litCO, xerr=litOHerr, yerr=litCOerr, fmt='rs')   # with errors
for idx, _ in enumerate(litCO):
    if idx <= 6:
        fmt='r>'
        litCO[idx] = litCH[idx]+0.5 - litOH[idx]
    elif idx == 7 or idx == 8:
        fmt='ys'
    elif idx >= 9:
        fmt='cd'
    if solar:
        litCOsol = convert2solar(litCO[idx], COsun)
        litOHsol = convert2solar(litOH[idx], OHsun)
        #print idx, fmt, litOHsol, litCOsol
        if ref[idx] != 1:   # do not plot results from James et al. (2014)
            plt.plot(litOHsol, litCOsol, fmt)   # in solar terms
            if idx==0:
                OHsolerr = numpy.sqrt(litOHerr[idx]**2 + OHsunerr**2)
                COsolerr = numpy.sqrt(litCOerr[idx]**2 + COsunerr**2)
                plt.errorbar(litOHsol, litCOsol, xerr=OHsolerr, yerr=COsolerr, fmt=fmt, ecolor='k')   # with errors
    else:
        plt.plot(litOH[idx], litCO[idx], fmt)    
# Pettini et al (2008)
pet08 = '../results/Pettini08.txt'
petOH, petNO = numpy.loadtxt(pet08, skiprows=5, usecols=(1,2), unpack=True)
petNH = petNO + petOH
petCH = approxC(petNH, petOH)
petCO = petCH - petOH
petnames = ['Q0000-2620', 'Q0100+1300', 'Q0112-306', 'Q0201+1120', 'J0307-4945', 'Q0347-3819', 'Q0407-4410', 'Q0407-4410',
            'Q0528-2505',  'HS0741+4741', 'J0841+1256', 'J0841+1256', 'J0900+4215', 'Q0913+072', 'Q0930+2858', 'Q1108-077',
            'Q1210+17', 'Q1232+0815',  'Q1331+170', 'Q1337+113', 'Q1409+095', 'J1435+5359', 'J1443+2724', 'SDSSJ1558-0031',
            'SDSSJ1558+4053', 'GB1759+7539', 'Q2059-360', 'Q2206-199', 'Q2230+02', 'Q2231-002', 'HE2243-6031', 'Q2332-094',
            'Q2342+342', 'Q2343+1232', 'Q2348-1444']
print ''
for pn, po, pc, pco in zip(petnames, petOH, petCH, petCO):
    if pn in litnames:
        print '{:<20} {:<6.2f} {:>6.2f} {:>6.2f}'.format(pn, po, pc, pco)
if solar:
    petCOsol = convert2solar(petCO, COsun)
    petOHsol = convert2solar(petOH, OHsun)
    plt.plot(petOHsol, petCOsol, 'mD')   # in solar terms
else:
    plt.plot(petOH, petCO, 'mD')  
# Data from Nava et al. (2006)
nav06 = '../results/Nava06.txt'
navOH, navOHerr, navNO, navNOerr = numpy.loadtxt(nav06, skiprows=2, usecols=(0,1,2,3), unpack=True)
navNH = navNO + navOH
navNHerr = numpy.sqrt(navNOerr**2 + navOHerr**2)
navCH = approxC(navNH, navOH)
navCHerr = err_approxC(navNH, navNHerr, navOH, navOHerr)
navCO = navCH - navOH
navCOerr = numpy.sqrt(navCHerr**2 + navOHerr**2)
#for navidx,navc in enumerate(navCO):
#    print navidx, navCH[navidx], navc 
if solar:
    navCOsol = convert2solar(navCO, COsun)
    navOHsol = convert2solar(navOH, OHsun)
    plt.plot(navOHsol, navCOsol, 'r>')   # in solar terms
    for idx, _ in enumerate(navOH):
        if idx==0 or idx ==1:
            OHsolerr = numpy.sqrt(navOHerr[idx]**2 + OHsunerr**2)
            COsolerr = numpy.sqrt(navCOerr[idx]**2 + COsunerr**2)
            plt.errorbar(navOHsol[idx], navCOsol[idx], xerr=OHsolerr, yerr=COsolerr, fmt='r+', ecolor='k')   # with errors
else:
    plt.plot(navOH, navCO, 'r+')  
# Galactic stellar data from Akerman et al.(2004) and Gustafsson et al.(1999) - but the last only for solar values
ake04 = '../results/Akerman04.txt'
akeOH, akeOHsol, akeCH, akeCOsol = numpy.loadtxt(ake04, skiprows=2, usecols=(8,10,11,12), unpack=True)
akeCO = akeCH - akeOH
if solar:
    gus99 = '../results/Gustafsson99.txt'
    gusCHsol, gusOHsol, gusage = numpy.loadtxt(gus99, skiprows=3, usecols=(1,2,3), unpack=True)
    gusCOsol = gusCHsol - gusOHsol
    plt.plot(akeOHsol, akeCOsol, 'wo', ms=5)   # in solar terms
    ax1.plot(gusOHsol, gusCOsol, 'ws', ms=5)   # in solar terms
    # Find correlation between metallicity and log age
    #coeffs = numpy.polyfit(gusOHsol, gusage, 1.0)
    #logAge_poly = numpy.poly1d(coeffs)
    #agefit = logAge_poly(gusOHsol)
    #ax2 = ax1.twiny()
    # get the age for xmin and xmax metallicities
    #ax2.set_xlim( logAge_poly(xmin), logAge_poly(xmax) )
    #ax2.set_xlabel("log (age of disk stars) [Gyr]")
    # Find correlation between C/O and log age
    #coeffs = numpy.polyfit(gusCOsol, gusage, 1.0)
    #COlogAge_poly = numpy.poly1d(coeffs)
    #agefit = COlogAge_poly(gusCOsol)
    #ax2 = ax1.twinx()
    # get the age for xmin and xmax C/O values
    #ax2.set_ylim( COlogAge_poly(xmin), COlogAge_poly(xmax) )
    #ax2.set_ylabel("log (age of disk stars) [Gyr]")
else:
    plt.plot(akeOH, akeCO, 'wo')  
# Save plot
if save_images:
    destination = os.path.join(full_results_path+'/plots', img_name)
    plt.savefig(destination)
    print('Plot %s was saved!' % destination)
plt.show()

# C/O vs O/H from SIMULATIONS
fig1 = plt.figure(1, figsize=(12, 10))
Clco = Clc-Clo
errClcou = numpy.sqrt( errClcu**2 + errClou**2 )
errClcod = numpy.sqrt( errClcd**2 + errClod**2 )
# add the reference objects
'''
for i in range(18, len(indeces_list)):
    Clo = numpy.append(Clo, OH[i])
    Clco = numpy.append(Clco, CO[i])
    errClcou = numpy.append(errClcou, COerr[i])
    errClcod = numpy.append(errClcod, COerr[i])
    errClou = numpy.append(errClou, OHerr[i])
    errClod = numpy.append(errClod, OHerr[i])
'''
errClcou = numpy.round(errClcou, decimals=2) 
errClcod = numpy.round(errClcod, decimals=2) 
#if no_errs:
#plt.plot(Clo, Clco, 'ro')
#else:
plt.errorbar(Clo, Clco, xerr=[errClou, errClod], yerr=[errClcou, errClcod], fmt='ro', ecolor='k')
'''
    for obj, i in zip(objects_list, indeces_list):
        #if obj in objects_with_Te:
        #    fmt='ro'
        #elif obj in objects_Te_literature:
        #    fmt='wo'
        if (obj in objects_not_in_sample) or (obj in objects_ref):
            fmt='b^'
        print i, Clo[i], Clco[i], [errClou[i], errClod[i]], [errClcou[i], errClcod[i]]
        plt.errorbar(Clo[i], Clco[i], xerr=[errClou[i], errClod[i]], yerr=[errClcou[i], errClcod[i]], fmt=fmt, ecolor='k')
        #plt.errorbar(Clo[i], Clco[i], xerr=errClou[i], yerr=errClcou[i], fmt=fmt, ecolor='r')
    '''
plt.xlim(6.8, 9.0)
yup = 1.1
ylo = -1.9
plt.ylim(ylo, yup)
plt.xticks(numpy.arange(6.8, 9.1, 0.2))
plt.yticks(numpy.arange(-1.8, yup, 0.2))
plt.minorticks_on()
for x, y, z in zip(Clo, Clco, objects_list):
    subxcoord = -1
    subycoord = 5
    side = 'right'
    if (z == 'pox4') or (z == 'sbs0948'):
        subxcoord = 5
        subycoord = -6
        side = 'left'
    if (z == 'sbs0926'):
        subxcoord = 25
    if (z == 'sbs0218'):
        subycoord = -13
    if (z == 'mrk5') or (z == 'sbs1319') or (z == 'sbs1054') or (z == 'mrk1087'):
        subxcoord = 5
        subycoord = -2
        side = 'left'
    if (z == 'iras08339') or (z == 'Sun'):
        subycoord = -13
        side = 'left'
    if (z == 'iiizw107') or (z == 'ngc1741') or (z == 'ngc456') or (z == 'iras08208') or (z == 'Orion'):
        subxcoord = 3
        side = 'left'
    if (z == 'ngc346'):
        subxcoord = 4
        subycoord = -11
        side = 'left'
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
plt.title('MCMC modeled C/O vs O/H')
plt.xlabel('12 + log (O/H)')
plt.ylabel('log (C/O)')
if save_images:
    img_name = img_name1b + typeofimage
    if use_our_sample_ONLY == False:
        img_name = img_name1b + otherrefs + typeofimage
    destination = os.path.join(full_results_path+'/plots', img_name)
    plt.savefig(destination)
    print('Plot %s was saved!' % destination)
plt.show()

'''
# C2/O2 vs O/H
fig1 = plt.figure(1, figsize=(12, 10))
for obj, i in zip(objects_list, indeces_list):
    if obj in objects_with_Te:
        fmt='ko'
    elif obj in objects_Te_literature:
        fmt='wo'
    elif (obj in objects_not_in_sample) or (obj in objects_ref):
        fmt='b^'
    plt.errorbar(OH[i], C2O2[i], xerr=OHerr[i], yerr=C2O2err[i], fmt=fmt, ecolor='k')
plt.xlim(7.0, 9.0)
#yup = 1.02
#ylo = -1.8
#plt.ylim(ylo, yup)
plt.xticks(numpy.arange(7.0, 9.0, 0.1))
#plt.yticks(numpy.arange(ylo, yup, 0.2))
for x, xe, y, ye, z in zip(OH, OHerr, CO, COerr, objects_list):
    subxcoord = -4
    subycoord = 5
    side = 'right'
    if (z == 'pox4'):
        subxcoord = -8
        subycoord = -7
    if (z == 'mrk960') or (z == 'sbs1319') or (z == 'ngc1741') or (z == 'Sun'):
        subycoord = -15
    if (z == 'iiizw107') or (z == 'ngc456') or (z == 'sbs0948') or (z == 'iras08208') or (z == 'iras08339') or (z == 'arp252') or (z == 'Orion'):
        subxcoord = 5
        side = 'left'
    if (z == 'ngc346'):
        subxcoord = 4
        subycoord = -11
        side = 'left'
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
plt.title('C2/O2 vs O/H')
plt.xlabel('12 + log (O/H)')
plt.ylabel('log (C2/O2)')
if save_images:
    img_name = img_name1 + typeofimage
    if use_our_sample_ONLY == False:
        img_name = img_name1 + otherrefs + typeofimage
    destination = os.path.join(full_results_path+'/plots', img_name)
    plt.savefig(destination)
    print('Plot %s was saved!' % destination)
plt.show()
'''

# N/O vs O/H
fig1 = plt.figure(1, figsize=(12, 10))
#plt.errorbar(OH, NO, xerr=OHerr, yerr=NOerr, fmt='ko')
for obj, i in zip(objects_list, indeces_list):
    if obj in objects_with_Te:
        fmt='ko'
        markersize=8
    elif obj in objects_Te_literature:
        fmt='k*'
        markersize=13
    elif (obj in objects_not_in_sample) or (obj in objects_ref):
        fmt='b^'
        markersize=9
    plt.errorbar(OH[i], NO[i], xerr=OHerr[i], yerr=NOerr[i], fmt=fmt, ecolor='k', markersize=markersize)
plt.xlim(6.8, 9.0)
yup = -0.0
ylo = -2.1
plt.ylim(ylo, yup)
plt.xticks(numpy.arange(6.8, 9.1, 0.2))
plt.yticks(numpy.arange(ylo, yup+0.1, 0.2))
plt.minorticks_on()
# Linear fit to linear part of diagram -- 2 different methods for this
coeffs, errs, line_fit = fit_line(OH[:-4], NO[:-4])   # 1
#plt.plot(OH[:-4], line_fit, 'g', lw=2)   # fit to linear part of diagram
print '\n * Coefficients of linear part of the plot: N/O vs O/H'
m = coeffs[0]
b = coeffs[1]
print '  m = %0.2f     b = %0.2f      errs: ' % (m, b), errs
m, b, r_value, p_value, err = stats.linregress(OH[:-4], NO[:-4])  # 2
print 'm, b, r_value, p_value, err: ', m, b, r_value, p_value, err
for x, y, z in zip(OH, NO, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    #print z, x, y
    subxcoord = -2
    subycoord = 5
    side = 'right'
    if (z == 'mrk1087') or (z == 'mrk1199'):
        subxcoord = 7
        side = 'left'
    if (z == 'ngc6822') or (z == 'ngc456') or (z == 'ngc346'):
        subxcoord = 4
        subycoord = 4
        side = 'left'
    if (z == 'sbs0926') or (z == 'iras08208') or (z == 'pox4') or (z == 'tol1457') or (z == '30Dor'):
        subycoord = -12
    if (z == 'sbs1319') or (z == 'sbs0948') or (z == 'Sun') or (z == 'Orion') or (z == 'izw18') or (z == 'ngc1741'):
        subxcoord = 4
        subycoord = -14
        side = 'left'
    if (z == 'sbs1054'):
        subxcoord = 33
        subycoord = -1
        side = 'left'
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
plt.title('N/O vs O/H')
plt.xlabel('12 + log (O/H)')
plt.ylabel('log (N/O)')
if save_images:
    img_name = img_name2 + typeofimage
    if use_our_sample_ONLY == False:
        img_name = img_name2 + otherrefs + typeofimage
    destination = os.path.join(full_results_path+'/plots', img_name)
    plt.savefig(destination)
    print('Plot %s was saved!' % destination)
plt.show()


# Ne/O vs O/H
fig1 = plt.figure(1, figsize=(12, 10))
#plt.errorbar(OH, NeO, xerr=OHerr, yerr=NeOerr, fmt='ko')
for obj, i in zip(objects_list, indeces_list):
    if obj in objects_with_Te:
        fmt='ko'
        markersize=8
    elif obj in objects_Te_literature:
        fmt='k*'
        markersize=13
    elif (obj in objects_not_in_sample) or (obj in objects_ref):
        fmt='b^'
        markersize=9
    plt.errorbar(OH[i], NeO[i], xerr=OHerr[i], yerr=NeOerr[i], fmt=fmt, ecolor='k', markersize=markersize)
plt.xlim(7.6, 9.0)
yup = 0.1
ylo = -1.2
plt.ylim(ylo, yup)
plt.xticks(numpy.arange(7.4, 9.0, 0.2))
plt.yticks(numpy.arange(ylo, yup, 0.1))
plt.minorticks_on()
for x, y, z in zip(OH, NeO, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    #print z, x, y
    subxcoord = -5
    subycoord = 5
    side = 'right'
    if (z == 'mrk5') or (z == 'ngc346'):
        subycoord = -11
    if (z == 'tol9') or (z == 'sbs1319'):
        subxcoord = 5
        subycoord = -11
        side = 'left'
    if (z == 'iras08208') or (z == 'pox4') or (z == 'arp252') or (z == 'iras08339') or (z == 'sbs1054'):
        subxcoord = 4
        subycoord = 4
        side = 'left'
    if (z == 'ngc6822') or (z == 'sbs0948'):
        subxcoord = 4
        subycoord = -14
        side = 'left'
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
plt.title('Ne/O vs O/H')
plt.xlabel('12 + log (O/H)')
plt.ylabel('log (Ne/O)')
if save_images:
    img_name = img_name3 + typeofimage
    if use_our_sample_ONLY == False:
        img_name = img_name3 + otherrefs + typeofimage
    destination = os.path.join(full_results_path+'/plots', img_name)
    plt.savefig(destination)
    print('Plot %s was saved!' % destination)
plt.show()


# C/N vs O/H
fig1 = plt.figure(1, figsize=(12, 10))
for obj, i in zip(objects_list, indeces_list):
    if obj in objects_with_Te:
        fmt='ko'
        markersize=8
    elif obj in objects_Te_literature:
        fmt='k*'
        markersize=13
    elif (obj in objects_not_in_sample) or (obj in objects_ref):
        fmt='b^'
        markersize=9
    plt.errorbar(OH[i], CN[i], xerr=OHerr[i], yerr=CNerr[i], fmt=fmt, ecolor='k', markersize=markersize)
plt.xlim(6.8, 9.0)
yup = 1.6
ylo = -0.7
plt.ylim(ylo, yup)
plt.xticks(numpy.arange(6.9, 9.0, 0.2))
plt.yticks(numpy.arange(ylo, yup, 0.2))
plt.minorticks_on()
# Linear fit to linear part of diagram -- 2 different methods for this
coeffs, errs, line_fit = fit_line(OH[:-4], CN[:-4])   # 1
plt.plot(OH[:-4], line_fit, 'g', lw=2)   # fit to linear part of diagram
plt.plot(OH[:-4], line_fit+err, 'g', lw=1, alpha=0.4)   # fit to linear plus error
plt.plot(OH[:-4], line_fit-err, 'g', lw=1, alpha=0.4)   # fit to linear minus error
print '\n * Coefficients of linear part of the plot: C/N vs O/H'
m = coeffs[0]
b = coeffs[1]
print '  m = %0.2f     b = %0.2f      errs: ' % (m, b), errs
m, b, r_value, p_value, err = stats.linregress(OH[:-4], CN[:-4])  # 2
print 'm, b, r_value, p_value, err: ', m, b, r_value, p_value, err
for x, y, z in zip(OH, CN, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    subxcoord = -5
    subycoord = 5
    side = 'right'
    if (z == 'ngc346'):
        subycoord = -10
    if (z == 'sbs0948') or (z == 'pox4') or (z == 'iras08339') or (z == '30Dor') or (z == 'ngc1741') or (z == 'sbs0948'):
        subxcoord = 5
        subycoord = -12
        side = 'left'
    if (z == 'mrk960') or (z == 'arp252') or (z == 'izw18') or (z == 'sbs1319') or (z == 'Orion'):
        subxcoord = 5
        subycoord = 4
        side = 'left'
    if (z == 'ngc6822'):
        subxcoord = 4
        subycoord = -14
        side = 'left'
    if (z == 'sbs0926') or (z == 'iiizw107') or (z == 'tol9'):
        subycoord = -12
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
plt.title('C/N vs O/H')
plt.xlabel('12 + log (O/H)')
plt.ylabel('log (C/N)')
if save_images:
    img_name = img_name4 + typeofimage
    if use_our_sample_ONLY == False:
        img_name = img_name4 + otherrefs + typeofimage
    destination = os.path.join(full_results_path+'/plots', img_name)
    plt.savefig(destination)
    print('Plot %s was saved!' % destination)
plt.show()

'''
# C2/N2 vs O/H
fig1 = plt.figure(1, figsize=(12, 10))
for obj, i in zip(objects_list, indeces_list):
    if obj in objects_with_Te:
        fmt='ko'
    elif obj in objects_Te_literature:
        fmt='k*'
    elif (obj in objects_not_in_sample) or (obj in objects_ref):
        fmt='b^'    
    errC2N2 = numpy.sqrt(C2O2err[i]**2 + N2O2err[i]**2)
    plt.errorbar(OH[i], C2O2[i]-N2O2[i], xerr=OHerr[i], yerr=errC2N2, fmt=fmt, ecolor='k')
plt.xlim(7.0, 9.0)
#yup = -0.2
#ylo = -1.2
#plt.ylim(ylo, yup)
plt.xticks(numpy.arange(7.0, 9.0, 0.2))
#plt.yticks(numpy.arange(ylo, yup, 0.2))
C2N2 = C2O2 - N2O2
for x, y, z in zip(OH, C2N2, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    #print z, x, y
    subxcoord = -5
    subycoord = 5
    side = 'right'
    if (z == 'ngc1741') or (z == 'pox4') or (z == 'sbs0948') or (z == 'ngc346') or (z == '30Dor'):
        subycoord = -10
    if(z == 'sbs1319') or (z == 'iiizw107') or (z == 'mrk5') or (z == 'sbs0948') or (z == 'arp252') or (z == 'Orion'):
        subxcoord = 5
        subycoord = 4
        side = 'left'
    if (z == 'ngc6822'):
        subxcoord = 4
        subycoord = -14
        side = 'left'
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
plt.title('C2/N2 vs O/H')
plt.xlabel('12 + log (O/H)')
plt.ylabel('log (C2+/N2+)')
if save_images:
    img_name = img_name5 + typeofimage
    if use_our_sample_ONLY == False:
        img_name = img_name5 + otherrefs + typeofimage
    destination = os.path.join(full_results_path+'/plots', img_name)
    plt.savefig(destination)
    print('Plot %s was saved!' % destination)
plt.show()
'''

'''# N/C vs C/H
# Calculate errors for C/N and C/H
CH = CO+OH
CH2fit = CH[:-4]
errCH = numpy.sqrt( COerr**2 + OHerr**2 )
NC = (NO+OH) - (CO+OH)
NC2fit = NC[:-4]
errNC = numpy.sqrt( NOerr**2 + COerr**2 + 2*OHerr**2)
# Adjust a linear fit to the plot
coeffs, errs, line_fit = fit_line(CH2fit, NC2fit)
print '\n Coefficients of initial guess to the plot of: N/C vs C/H'
m = coeffs[0]
b = coeffs[1]
print 'm = %0.3f     b = %0.3f    errs: ' % (m, b), errs
m, b, r_value, p_value, err = stats.linregress(CH[:-4], NC[:-4])  # 2
print 'm, b, r_value, p_value, err: ', m, b, r_value, p_value, err
fig1 = plt.figure(1, figsize=(12, 10))
for obj, i in zip(objects_list, indeces_list):
    if obj in objects_with_Te:
        fmt='ko'
    elif obj in objects_Te_literature:
        fmt='k*'
    elif (obj in objects_not_in_sample) or (obj in objects_ref):
        fmt='b^'
    #plt.plot(CH[i], NC[i], 'ko')
    plt.errorbar(CH[i], NC[i], xerr=errCH[i], yerr=errNC[i], fmt=fmt, ecolor='k')
    #plt.errorbar(CH[i], CN[i], xerr=errCH[i], yerr=CNerr[i], fmt=fmt, ecolor='k')
plt.xlim(5.6, 9.5)
yup = 0.7
ylo = -2.6
plt.ylim(0, yup)
plt.xticks(numpy.arange(5.6, 9.5, 0.2))
plt.yticks(numpy.arange(ylo, yup+0.2, 0.2))
plt.minorticks_on()
for x, y, z in zip(CH, NC, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    subxcoord = -5
    subycoord = 5
    side = 'right'
    #print z, x, y
    if (z == 'mrk960'):
        subxcoord = -60
        subycoord = 0
    if (z == 'pox4'):
        subxcoord = -22
        subycoord = -2
    if (z == 'sbs0926'):
        subxcoord = -38
        subycoord = -5
    if (z == 'sbs1054'):
        subxcoord = -5
        subycoord = -12
    if (z == 'mrk1087') or (z == 'sbs1319'):
        subxcoord = 4
        subycoord = -14
        side = 'left'
    if (z == 'iras08339') or (z == 'mrk1199') or (z == 'tol9'):
        subxcoord = 4
        side = 'left'
    if (z == 'Orion'):
        subxcoord = -3
        subycoord = -14 
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
plt.title('N/C vs C/H')
plt.xlabel('12 + log (C/H)')
plt.ylabel('log (N/C)')
# Show line presented in paper:
mp = -0.36
bp = 2.20
yp = mp*CH + bp
plt.plot(CH, yp, 'g', lw=2)
if save_images:
    img_name = img_name6 + typeofimage
    if use_our_sample_ONLY == False:
        img_name = img_name6 + otherrefs + typeofimage
    destination = os.path.join(full_results_path+'/plots', img_name)
    plt.savefig(destination)
    print('Plot %s was saved!' % destination)
if show_fitted_line:
    plt.plot(CH2fit, line_fit)   # plot fitted line
# Initialize the chain with a Gaussian ball around the maximum likelihood result, for which we use optimize
x = CH2fit
y = NC2fit
yerr = NCerr[:-4]
ndim, nwalkers, nruns = 2, 100, 300
randadd2point = lambda x: x+numpy.random.rand(1)*-1
p0 = numpy.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[x, y, yerr])
pos, prob, state = sampler.run_mcmc(p0, nruns)   
# best model
wh = numpy.where( prob == prob.max() )[0][0]
p = pos[ wh, : ]
line1 ='values of the %i dimensions that best fit the data according to Chi^2 in %i runs are the following:' % (ndim, nruns)
line2 = 'm = %0.3f   b = %0.3f' % (p[0], p[1])
print line1
print line2
allm, allb = [], []
subsample = numpy.array([]).reshape(0, 2)
positives, negatives = 0, 0
for theta in pos:
    #print theta
    allm.append(theta[0])
    allb.append(theta[1])
    subsample = numpy.vstack((subsample, theta))
    y = line_eq(theta, x)
    #plt.plot( x, y, "b", alpha=0.1 )
    if theta[0] < 0.0:
        negatives += 1
    else:
        positives += 1
avgm = sum(allm)/len(allm)
avgb = sum(allb)/len(allb)
print 'Average values:  m = %0.3f   b= %0.3f' % (avgm, avgb)
theta = [avgm, avgb]
if show_fitted_line:
    plt.plot( x, line_eq( p, x ), "b", lw=5, alpha=0.4 )
samples = sampler.chain[:, nruns*0.2, :].reshape((-1, ndim))
print ('positives, negatives: ', positives, negatives)
#fig = corner.corner(samples, labels=["$m$", "$b$"], truths=[m, b])
#fig.show()
# Calculate the uncertainties based on the 25, 50, and 75th percentiles
m_mcmc, b_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*numpy.percentile(subsample, [25, 50, 75], axis=0)))
print 'm_mcmc, b_mcmc:', m_mcmc, b_mcmc
m_mcmc2, b_mcmc2 = map(lambda v: (v), zip(*numpy.percentile(subsample, [25, 50, 75], axis=0)))
print 'm_mcmc2, b_mcmc2:', m_mcmc2, b_mcmc2

# Cloudy data
errClc = (errClcu+errClcd)/2.0
Clnc = Cln-Clc
errClncu = numpy.sqrt( errClnu**2 + errClcu**2 )
errClncd = numpy.sqrt( errClnd**2 + errClcd**2 )
errClncu = numpy.round(errClncu, decimals=2) 
errClncd = numpy.round(errClncd, decimals=2) 
errClnc = (errClncu+errClncd)/2.0
# Adjust a linear fit to the plot
clcoeffs, errs, clline_fit = fit_line(Clnc, Clc)
print '* Coefficients of CLOUDY initial guess to the plot of: N/O vs C/N'
clm = clcoeffs[0]
clb = clcoeffs[1]
print '  Cloudy_m = %0.3f     Cloudy_b = %0.3f   errs: ' % (clm, clb), errs
m, b, r_value, p_value, err = stats.linregress(Clnc, Clc)  # 2
print 'm, b, r_value, p_value, err: ', m, b, r_value, p_value, err
ypcl = clm*Clnc + clb
plt.plot(Clnc, ypcl, 'm', lw=2)
# Initialize the chain with a Gaussian ball around the maximum likelihood result, for which we use optimize
x = Clnc
y = Clc
yerr = errClc
ndim, nwalkers, nruns = 2, 100, 300
randadd2point = lambda x: x+numpy.random.rand(1)*-1
p0 = numpy.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[x, y, yerr])
pos, prob, state = sampler.run_mcmc(p0, nruns)   
# best model
wh = numpy.where( prob == prob.max() )[0][0]
p = pos[ wh, : ]
line1 =' CLOUDY values of the %i dimensions that best fit the data according to Chi^2 in %i runs are the following:' % (ndim, nruns)
line2 = ' m = %0.3f   b = %0.3f' % (p[0], p[1])
print line1
print line2
allm, allb = [], []
subsample = numpy.array([]).reshape(0, 2)
positives, negatives = 0, 0
for theta in pos:
    if theta[0] < 0.0:
        negatives += 1
    else:
        positives += 1
    #plt.plot( x, y, "r", alpha=0.1 )
    allm.append(theta[0])
    allb.append(theta[1])
    subsample = numpy.vstack((subsample, theta))
avgm = sum(allm)/len(allm)
avgb = sum(allb)/len(allb)
print 'CLOUDY Average values:  m = %0.3f   b= %0.3f' % (avgm, avgb)
theta = [avgm, avgb]
if show_fitted_line:
    plt.plot( x, line_eq( p, x ), "r", lw=5, alpha=0.4 )
samples = sampler.chain[:, nruns*0.2, :].reshape((-1, ndim))
print (' positives, negatives: ', positives, negatives)
#fig = corner.corner(samples, labels=["$m$", "$b$"], truths=[m, b])
#fig.show()
# Calculate the uncertainties based on the 25, 50, and 75th percentiles
m_mcmc, b_mcmc = map(lambda v: (v[1], numpy.abs(v[2])-numpy.abs(v[1]), numpy.abs(v[1])-numpy.abs(v[0])), zip(*numpy.percentile(subsample, [25, 50, 75], axis=0)))
print ' m_mcmc, b_mcmc:', m_mcmc, b_mcmc
m_mcmc2, b_mcmc2 = map(lambda v: (v), zip(*numpy.percentile(subsample, [25, 50, 75], axis=0)))
print ' m_mcmc2, b_mcmc2:', m_mcmc2, b_mcmc2

plt.show()
'''

# C/N vs C/H
# Calculate errors for C/N and C/H
CH = CO+OH
CH2fit = CH[:-4]
errCH = numpy.sqrt( COerr**2 + OHerr**2 )
CN2fit = CN[:-4]
# Adjust a linear fit to the plot
coeffs, errs, line_fit = fit_line(CH2fit, CN2fit)
print '\n Coefficients of initial guess to the plot of: C/N vs C/H'
m = coeffs[0]
b = coeffs[1]
print 'm = %0.3f     b = %0.3f    errs: ' % (m, b), errs
m, b, r_value, p_value, err = stats.linregress(CH[:-4], CN[:-4])  # 2
print 'm, b, r_value, p_value, err: ', m, b, r_value, p_value, err
fig1 = plt.figure(1, figsize=(12, 10))
for obj, i in zip(objects_list, indeces_list):
    if obj in objects_with_Te:
        fmt='ko'
        markersize=8
    elif obj in objects_Te_literature:
        fmt='k*'
        markersize=13
    elif (obj in objects_not_in_sample) or (obj in objects_ref):
        fmt='b^'
        markersize=9
    #plt.plot(CH[i], CN[i], 'ko')
    plt.errorbar(CH[i], CN[i], xerr=errCH[i], yerr=CNerr[i], fmt=fmt, ecolor='k', markersize=markersize)
plt.xlim(5.6, 9.5)
yup = 1.6
ylo = -0.5 
plt.ylim(0, yup)
plt.xticks(numpy.arange(5.6, 9.5, 0.2))
plt.yticks(numpy.arange(ylo, yup+0.2, 0.2))
plt.minorticks_on()
for x, y, z in zip(CH, CN, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    subxcoord = -4
    subycoord = 5
    side = 'right'
    #print z, x, y
    if (z == 'mrk960'):
        subxcoord = -60
        subycoord = -2
    if (z == 'pox4'):
        subxcoord = -22
        subycoord = -2
    if (z == 'sbs0926') or (z == 'sbs0948'):
        subxcoord = -38
        subycoord = -2
    if (z == 'iras08339') or (z == 'ngc1741') or (z == '30Dor'):
        subxcoord = -5
        subycoord = -12
    if (z == 'sbs1054') or (z == 'sbs1415'):
        subxcoord = 4
        subycoord = -14
        side = 'left'
    if (z == 'mrk1199') or (z == 'tol9') or (z == 'mrk1087') or (z == 'sbs1319'):
        subxcoord = 4
        side = 'left'
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
plt.title('C/N vs C/H')
plt.xlabel('12 + log (C/H)')
plt.ylabel('log (C/N)')
# Show line presented in paper:
mp = 0.36
bp = -2.20
yp = mp*CH + bp
plt.plot(CH, yp, 'g', lw=2)
plt.plot(CH, yp+0.3, 'g', lw=1, alpha=0.4)   # fit to linear plus error
plt.plot(CH, yp-0.3, 'g', lw=1, alpha=0.4)   # fit to linear minus error
if save_images:
    img_name = img_name6 + typeofimage
    if use_our_sample_ONLY == False:
        img_name = img_name6 + otherrefs + typeofimage
    destination = os.path.join(full_results_path+'/plots', img_name)
    plt.savefig(destination)
    print('Plot %s was saved!' % destination)
if show_fitted_line:
    plt.plot(CH2fit, line_fit)   # plot fitted line
# Initialize the chain with a Gaussian ball around the maximum likelihood result, for which we use optimize
x = CH2fit
y = CN2fit
yerr = CNerr[:-4]
ndim, nwalkers, nruns = 2, 100, 300
randadd2point = lambda x: x+numpy.random.rand(1)*-1
p0 = numpy.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[x, y, yerr])
pos, prob, state = sampler.run_mcmc(p0, nruns)   
# best model
wh = numpy.where( prob == prob.max() )[0][0]
p = pos[ wh, : ]
line1 ='values of the %i dimensions that best fit the data according to Chi^2 in %i runs are the following:' % (ndim, nruns)
line2 = 'm = %0.3f   b = %0.3f' % (p[0], p[1])
print line1
print line2
allm, allb = [], []
subsample = numpy.array([]).reshape(0, 2)
positives, negatives = 0, 0
for theta in pos:
    #print theta
    allm.append(theta[0])
    allb.append(theta[1])
    subsample = numpy.vstack((subsample, theta))
    y = line_eq(theta, x)
    #plt.plot( x, y, "b", alpha=0.1 )
    if theta[0] < 0.0:
        negatives += 1
    else:
        positives += 1
avgm = sum(allm)/len(allm)
avgb = sum(allb)/len(allb)
print 'Average values:  m = %0.3f   b= %0.3f' % (avgm, avgb)
theta = [avgm, avgb]
if show_fitted_line:
    plt.plot( x, line_eq( p, x ), "b", lw=5, alpha=0.4 )
samples = sampler.chain[:, nruns*0.2, :].reshape((-1, ndim))
print ('positives, negatives: ', positives, negatives)
#fig = corner.corner(samples, labels=["$m$", "$b$"], truths=[m, b])
#fig.show()
# Calculate the uncertainties based on the 25, 50, and 75th percentiles
m_mcmc, b_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*numpy.percentile(subsample, [25, 50, 75], axis=0)))
print 'm_mcmc, b_mcmc:', m_mcmc, b_mcmc
m_mcmc2, b_mcmc2 = map(lambda v: (v), zip(*numpy.percentile(subsample, [25, 50, 75], axis=0)))
print 'm_mcmc2, b_mcmc2:', m_mcmc2, b_mcmc2

# Cloudy data
errClc = (errClcu+errClcd)/2.0
Clcn = Clc-Cln
errClcnu = numpy.sqrt( errClcu**2 + errClnu**2 )
errClcnd = numpy.sqrt( errClcd**2 + errClnd**2 )
errClcnu = numpy.round(errClcnu, decimals=2) 
errClcnd = numpy.round(errClcnd, decimals=2) 
errClcn = (errClcnu+errClcnd)/2.0
print '  min(Clcn)=', min(Clcn), '  max(Clcn)=', max(Clcn)
print '  min(Clc)=', min(Clc), '  max(Clc)=', max(Clc)
# Adjust a linear fit to the plot
clcoeffs, errs, clline_fit = fit_line(Clcn, Clc)
print '* Coefficients of CLOUDY initial guess to the plot of: N/O vs C/N'
clm = clcoeffs[0]
clb = clcoeffs[1]
print '  Cloudy_m = %0.3f     Cloudy_b = %0.3f   errs: ' % (clm, clb), errs
m, b, r_value, p_value, err = stats.linregress(Clcn, Clc)  # 2
print 'm, b, r_value, p_value, err: ', m, b, r_value, p_value, err
ypcl = clm*Clcn + clb
plt.plot(Clcn, ypcl, 'm', lw=2)
# Initialize the chain with a Gaussian ball around the maximum likelihood result, for which we use optimize
x = Clcn
y = Clc
yerr = errClc
ndim, nwalkers, nruns = 2, 100, 300
randadd2point = lambda x: x+numpy.random.rand(1)*-1
p0 = numpy.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[x, y, yerr])
pos, prob, state = sampler.run_mcmc(p0, nruns)   
# best model
wh = numpy.where( prob == prob.max() )[0][0]
p = pos[ wh, : ]
line1 =' CLOUDY values of the %i dimensions that best fit the data according to Chi^2 in %i runs are the following:' % (ndim, nruns)
line2 = ' m = %0.3f   b = %0.3f' % (p[0], p[1])
print line1
print line2
allm, allb = [], []
subsample = numpy.array([]).reshape(0, 2)
positives, negatives = 0, 0
for theta in pos:
    if theta[0] < 0.0:
        negatives += 1
    else:
        positives += 1
    #plt.plot( x, y, "r", alpha=0.1 )
    allm.append(theta[0])
    allb.append(theta[1])
    subsample = numpy.vstack((subsample, theta))
avgm = sum(allm)/len(allm)
avgb = sum(allb)/len(allb)
print 'CLOUDY Average values:  m = %0.3f   b= %0.3f' % (avgm, avgb)
theta = [avgm, avgb]
if show_fitted_line:
    plt.plot( x, line_eq( p, x ), "r", lw=5, alpha=0.4 )
samples = sampler.chain[:, nruns*0.2, :].reshape((-1, ndim))
print (' positives, negatives: ', positives, negatives)
#fig = corner.corner(samples, labels=["$m$", "$b$"], truths=[m, b])
#fig.show()
# Calculate the uncertainties based on the 25, 50, and 75th percentiles
m_mcmc, b_mcmc = map(lambda v: (v[1], numpy.abs(v[2])-numpy.abs(v[1]), numpy.abs(v[1])-numpy.abs(v[0])), zip(*numpy.percentile(subsample, [25, 50, 75], axis=0)))
print ' m_mcmc, b_mcmc:', m_mcmc, b_mcmc
m_mcmc2, b_mcmc2 = map(lambda v: (v), zip(*numpy.percentile(subsample, [25, 50, 75], axis=0)))
print ' m_mcmc2, b_mcmc2:', m_mcmc2, b_mcmc2

plt.show()


# N/O vs C/N
# Adjust a linear fit to the plot
NO2fit = NO[:-4]
CN2fit = CN[:-4]
coeffs, errs, line_fit = fit_line(NO, CN)
#coeffs, line_fit = fit_line(NO2fit, CN2fit)
#print 'Linear fit: y = mx + b'
#print line_fit
print '\n Coefficients of initial guess to the plot of: N/O vs C/N'
m = coeffs[0]
b = coeffs[1]
print 'm = %0.3f     b = %0.3f    errs: ' % (m, b), errs
m, b, r_value, p_value, err = stats.linregress(NO2fit, CN2fit)  # 2
print 'm, b, r_value, p_value, err: ', m, b, r_value, p_value, err
fig1 = plt.figure(1, figsize=(12, 10))
for obj, i in zip(objects_list, indeces_list):
    if obj in objects_with_Te:
        fmt='ko'
        markersize=8
    elif obj in objects_Te_literature:
        fmt='k*'
        markersize=13
    elif (obj in objects_not_in_sample) or (obj in objects_ref):
        fmt='b^'
        markersize=9
    #plt.plot(NO[i], CN[i], 'ko')   # plot points without errors
    plt.errorbar(NO[i], CN[i], xerr=NOerr[i], yerr=CNerr[i], fmt=fmt, ecolor='k', markersize=markersize)
if show_fitted_line:
    plt.plot(NO, line_fit)   # plot fitted line
    #plt.plot(NO2fit, line_fit)   # plot fitted line
# Show line previously presented in paper:
mp = -0.18
bp = 0.50
yp = mp*NO + bp
#plt.plot(NO, yp, 'm')

bp2 = 0.31#+0.10-0.10   #0.35 
mp2 = -0.21#+0.09-0.06  #-0.30
yp2 = mp2*NO + bp2
plt.plot(NO, yp2, 'g', lw=2)
plt.plot(NO, yp2+err, 'g', lw=1, alpha=0.4)   # fit to linear plus error
plt.plot(NO, yp2-err, 'g', lw=1, alpha=0.4)   # fit to linear minus error

# Set limits
plt.xlim(-2.1, 0.0)
yup = 1.6
ylo = -0.8
plt.ylim(ylo, yup)
plt.xticks(numpy.arange(-2.1, 0.1, 0.2))
plt.yticks(numpy.arange(ylo, yup+0.2, 0.2))
plt.minorticks_on()
for x, y, z in zip(NO, CN, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    #print z, x, y
    subxcoord = -2
    subycoord = 5
    side = 'right'
    if (z == 'sbs0948') or (z == 'Sun') or (z == 'pox4') or (z == 'sbs1415'):
        subycoord = -13
    if (z == 'mrk960') or (z == 'iiizw107') or (z == '30Dor') or (z == 'mrk1199') or (z == 'ngc1741'):
        subxcoord = 4
        subycoord = -13
        side = 'left'
    if (z == 'arp252') or (z == 'iras08339') or (z == 'sbs1054') or (z == 'sbs0926') or (z == 'sbs1319'):
        subxcoord = 4
        side = 'left'
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
plt.title('N/O vs C/N')
plt.xlabel('log (N/O)')
plt.ylabel('log (C/N)')
if save_images:
    img_name = img_name7 + typeofimage
    if use_our_sample_ONLY == False:
        img_name = img_name7 + otherrefs + typeofimage
    destination = os.path.join(full_results_path+'/plots', img_name)
    #plt.savefig(destination)
    #print('Plot %s was saved!' % destination)
# Data from Cunha & Lambert (1994)
plt.plot(nlam-olam, clam-nlam, 'md')
# Data from Nieva & Simon-Diaz (2011)
plt.plot(nievn-nievo, nievc-nievn, 'cs')
# Data from Kilian (1992)
plt.plot(kiln-kilo, kilc-kiln, 'y>')
# Data from Daflon
plt.plot(dafn-dafo, dafc-dafn, 'rD')
# Data from Morel 2008
plt.plot(morn-moro, morc-morn, 'w*')
# Initialize the chain with a Gaussian ball around the maximum likelihood result, for which we use optimize
x = NO2fit
y = CN2fit 
yerr = CNerr[:-4]
ndim, nwalkers, nruns = 2, 100, 500
randadd2point = lambda x: x+numpy.random.rand(1)
p0 = numpy.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[x, y, yerr])
pos, prob, state = sampler.run_mcmc(p0, nruns)   
# best model
wh = numpy.where( prob == prob.max() )[0][0]
p = pos[ wh, : ]
#theta = [avgm, avgb]
#theta = [m, b]
line1 ='values of the %i dimensions that best fit the data according to Chi^2 in %i runs are the following:' % (ndim, nruns)
line2 = 'm = %0.3f   b = %0.3f' % (p[0], p[1])
print line1
print line2
allm, allb = [], []
subsample = numpy.array([]).reshape(0, 2)
for theta in pos:
    if theta[0] < 0.0:
        allm.append(theta[0])
        allb.append(theta[1])
        subsample = numpy.vstack((subsample, theta))
avgm = sum(allm)/len(allm)
avgb = sum(allb)/len(allb)
print 'Average values:  m = %0.3f   b= %0.3f' % (avgm, avgb)
if show_fitted_line:
    plt.plot( x, line_eq( p, x ), "b", lw=5, alpha=0.4 )   # best fit according to data (min Chi2)
    avgp = [avgm, avgb]
    plt.plot( x, line_eq( avgp, x ), "c", lw=5, alpha=0.4 )   # average line fit
positives, negatives = 0, 0
for p in pos:
    y = line_eq(p, x)
    extrap_xmin = numpy.interp(NO[-2], x, y)
    #x = numpy.append(x, NO[-2])
    #y = numpy.append(y, extrap_xmin)
    if p[0] < 0.0:
        negatives += 1
        plt.plot( x, y, "r", alpha=0.1 )
    else:
        positives += 1
#samples = sampler.chain[:, nruns*0.2, :].reshape((-1, ndim))
#fig = corner.corner(samples, labels=["$m$", "$b$"], truths=[m, b])
#fig.show()
# to save with mcmc lines uncomment following line
if save_images:
    plt.savefig(destination)
    print('Plot %s was saved!' % destination)
# Calculate the uncertainties based on the 25, 50, and 75 percentiles
m_mcmc, b_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*numpy.percentile(subsample, [25, 50, 75], axis=0)))
print 'm_mcmc, b_mcmc:', m_mcmc, b_mcmc
m_mcmc2, b_mcmc2 = map(lambda v: (v), zip(*numpy.percentile(subsample, [25, 50, 75], axis=0)))
print 'm_mcmc2, b_mcmc2:', m_mcmc2, b_mcmc2

# Cloudy data
Clno = Cln-Clo
errClnou = numpy.sqrt( errClnu**2 + errClou**2 )
errClnod = numpy.sqrt( errClnd**2 + errClod**2 )
errClnou = numpy.round(errClnou, decimals=2) 
errClnod = numpy.round(errClnod, decimals=2) 
errClno = (errClnou+errClnod)/2.0
Clcn = Clc-Cln
errClcnu = numpy.sqrt( errClcu**2 + errClnu**2 )
errClcnd = numpy.sqrt( errClcd**2 + errClnd**2 )
errClcnu = numpy.round(errClcnu, decimals=2) 
errClcnd = numpy.round(errClcnd, decimals=2) 
errClcn = (errClcnu+errClcnd)/2.0
print '  min(Clcn)=', min(Clcn), '  max(Clcn)=', max(Clcn)
print '  min(Clno)=', min(Clno), '  max(Clno)=', max(Clno)
# Adjust a linear fit to the plot
clcoeffs, errs, clline_fit = fit_line(Clno, Clcn)
print '* Coefficients of CLOUDY initial guess to the plot of: N/O vs C/N'
clm = clcoeffs[0]
clb = clcoeffs[1]
print '  Cloudy_m = %0.3f     Cloudy_b = %0.3f    errs' % (clm, clb), errs
print ' fitted cloudy line   y= ', clline_fit
print ' fitted cloudy line   x= ', Clno
#plt.plot(Clno, clline_fit, 'm', lw=2)   # initial guess for cloudy line fit
m, b, r_value, p_value, err = stats.linregress(Clno, Clcn)  # 2
print 'm, b, r_value, p_value, err: ', m, b, r_value, p_value, err
'''
plt.plot(Clno, Clcn, 'ro')
for x, y, z in zip(Clno, Clcn, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    #print z, x, y
    subxcoord = -2
    subycoord = 5
    side = 'right'
    if (z == 'sbs1319') or (z == 'pox4') or (z == 'sbs0948') or (z == 'Sun'):
        subycoord = -9
    if (z == 'sbs0926') or (z == 'iiizw107') or (z == '30Dor') or (z == 'sbs1415') or (z == 'iras08339') or (z == 'mrk1199'):
        subxcoord = 4
        subycoord = -12
        side = 'left'
    if (z == 'arp252'):
        subxcoord = 4
        side = 'left'
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
'''
# Initialize the chain with a Gaussian ball around the maximum likelihood result, for which we use optimize
x = Clno
y = Clcn 
yerr = errClcn
ndim, nwalkers, nruns = 2, 100, 500
randadd2point = lambda x: x+numpy.random.rand(1)
p0 = numpy.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[x, y, yerr])
pos, prob, state = sampler.run_mcmc(p0, nruns)   
# best model
wh = numpy.where( prob == prob.max() )[0][0]
p = pos[ wh, : ]
line1 =' CLOUDY values of the %i dimensions that best fit the data according to Chi^2 in %i runs are the following:' % (ndim, nruns)
line2 = ' m = %0.3f   b = %0.3f' % (p[0], p[1])
print line1
print line2
allm, allb = [], []
subsample = numpy.array([]).reshape(0, 2)
for theta in pos:
    if theta[0] < 0.0:
        allm.append(theta[0])
        allb.append(theta[1])
        subsample = numpy.vstack((subsample, theta))
avgm = sum(allm)/len(allm)
avgb = sum(allb)/len(allb)
print ' CLOUDY Average values:  m = %0.3f   b= %0.3f' % (avgm, avgb)
if show_fitted_line:
    plt.plot( x, line_eq( p, x ), "y", lw=5, alpha=0.4 )   # line fit of min Chi2
    avgp = [avgm, avgb]
    plt.plot( x, line_eq( avgp, x ), "r", lw=5, alpha=0.4 )   # avg line fit
positives, negatives = 0, 0
for p in pos:
    y = line_eq(p, x)
    extrap_xmin = numpy.interp(NO[-2], x, y)
    if p[0] < 0.0:
        negatives += 1
        #plt.plot( x, y, "r", alpha=0.1 )
    else:
        positives += 1
#img_name = img_name7 + typeofimage
#destination = os.path.join(full_results_path+'/plots', img_name)
#plt.savefig(destination)
#print('Plot %s was saved!' % destination)
#samples = sampler.chain[:, nruns*0.2, :].reshape((-1, ndim))
#fig = corner.corner(samples, labels=["$m$", "$b$"], truths=[m, b])
#fig.show()
# Calculate the uncertainties based on the 25, 50, and 75 percentiles
m_mcmc, b_mcmc = map(lambda v: (v[1], numpy.abs(v[2])-numpy.abs(v[1]), numpy.abs(v[1])-numpy.abs(v[0])), zip(*numpy.percentile(subsample, [25, 50, 75], axis=0)))
print '  m_mcmc, b_mcmc:', m_mcmc, b_mcmc
m_mcmc2, b_mcmc2 = map(lambda v: (v), zip(*numpy.percentile(subsample, [25, 50, 75], axis=0)))
print '  m_mcmc2, b_mcmc2:', m_mcmc2, b_mcmc2

plt.show()


print ' Code finished!'
