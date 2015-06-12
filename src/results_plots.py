import numpy
import os
import string
import emcee
import triangle
import scipy.optimize as op
from matplotlib import pyplot as plt
from copy import deepcopy

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

############################################################################################################################################

img_name1 = 'COvsOH'
img_name1a = 'COvsOH_literature'
img_name1b = 'COvsOH_Cldy'
img_name2 = 'NOvsOH'
img_name3 = 'NeOvsOH'
img_name4 = 'CNvsOH'
img_name5 = 'C+2N+2vsOH'
img_name6 = 'NCvsCH'
img_name7 = 'NOvsCN'
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
    coefficients = numpy.polyfit(arrx, arry, order)
    polynomial = numpy.poly1d(coefficients)
    f_pol = polynomial(arrx)
    #fitted_line = np.array([arry, f_pol])
    #print 'this is x and y of the fitted_line = ', fitted_line
    return coefficients, f_pol

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
    return 0.5 - 0.18 * (N - O) + N   

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
    rows2skip = 10
    rf = numpy.loadtxt(resfile, skiprows=rows2skip, usecols=(1,2,3,4,5,6,7,8), unpack=True)
    OH, OHerr, CO, COerr, NO, NOerr, NeO, NeOerr = rf
        
# These are the results are for the sample objects in the following order:
ordered_sample = ['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
                  'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
# now divide the sample in objects for in which we found Te and those in which we did not
objects_with_Te = ['mrk5', 'pox4', 'sbs0948', 'sbs1054', 'tol1457', 'iras08208']
objects_Te_literature = ['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk960', 'ngc1741', 'sbs0218',
                         'sbs0926', 'sbs1319', 'tol9', 'arp252', 'sbs1415']

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
no_errs = True
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

# C/O vs O/H from OBSERVATIONS
fig1 = plt.figure(1, figsize=(12, 10))
#plt.errorbar(OH, CO, xerr=OHerr, yerr=COerr, fmt='ko')    # this plots ALL points with the same symbol
for obj, i in zip(objects_list, indeces_list):
    if obj in objects_with_Te:
        fmt='ko'
    elif obj in objects_Te_literature:
        fmt='wo'
    elif (obj in objects_not_in_sample) or (obj in objects_ref):
        fmt='b^'
    plt.errorbar(OH[i], CO[i], xerr=OHerr[i], yerr=COerr[i], fmt=fmt, ecolor='k')
plt.xlim(6.8, 9.0)
yup = 1.1
ylo = -1.9
plt.ylim(ylo, yup)
plt.xticks(numpy.arange(6.8, 9.1, 0.1))
plt.yticks(numpy.arange(-1.8, yup, 0.2))
for x, xe, y, ye, z in zip(OH, OHerr, CO, COerr, objects_list):
    subxcoord = -3
    subycoord = 5
    side = 'right'
    if (z == 'mrk960') or (z == 'mrk5') or (z == 'ngc1741') or (z == 'Sun'):
        subxcoord = -4
        subycoord = -12
    if (z == 'sbs0218') or (z == 'sbs0948') or (z == 'iras08208') or (z == 'arp252') or (z == 'pox4') or (z == 'ngc456') or (z == 'ngc6822') or (z == 'Orion') or (z == 'izw18'):
        subxcoord = 5
        side = 'left'
    if (z == 'ngc346') or (z == 'sbs1319') or (z == 'iiizw107'):
        subxcoord = 4
        subycoord = -11
        side = 'left'
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
plt.title('C/O vs O/H')
plt.xlabel('12 + log (O/H)')
plt.ylabel('log (C/O)')
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
# PROTOSOLAR values from Asplund et al. (2009)
COsun = -0.26 # because C/H = 8.47+-0.05
COsunerr = 0.07
OHsun = 8.73 
OHsunerr = 0.05
# OUR STIS OBSERVATIONS
for obj, i in zip(objects_list, indeces_list):
    if obj in objects_with_Te:
        fmt='ko'
    elif obj in objects_Te_literature:
        fmt='ko'
    elif (obj in objects_not_in_sample) or (obj in objects_ref):
        fmt='b^'
    #plt.errorbar(OH[i], CO[i], xerr=OHerr[i], yerr=COerr[i], fmt=fmt, ecolor='k')   # with errors
    if solar:
        COsol = convert2solar(CO[i], COsun)
        OHsol = convert2solar(OH[i], OHsun)
        plt.plot(OHsol, COsol, fmt, ms=7)   # without errors AND in terms of solar values
        if obj=='izw18':
            OHsolerr = numpy.sqrt(OHerr[i]**2 + OHsunerr**2)
            COsolerr = numpy.sqrt(COerr[i]**2 + COsunerr**2)
            plt.errorbar(OHsol, COsol, xerr=OHsolerr, yerr=COsolerr, fmt=fmt, ecolor='k')   # with errors
    else:
        plt.plot(OH[i], CO[i], fmt, ms=7)   # without errors
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
litCH, litCHerr, litOH, litOHerr = numpy.loadtxt(lit, skiprows=8, usecols=(4,5,6,7), unpack=True)
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
# Galactic stellar data from Akerman et al.(2004) and Gustafsson et al.(1999) - but the last only for solar values
ake04 = '../results/Akerman04.txt'
akeOH, akeOHsol, akeCH, akeCOsol = numpy.loadtxt(ake04, skiprows=2, usecols=(8,10,11,12), unpack=True)
akeCO = akeCH - akeOH
if solar:
    gus99 = '../results/Gustafsson99.txt'
    gusCHsol, gusOHsol = numpy.loadtxt(gus99, skiprows=3, usecols=(1,2), unpack=True)
    gusCOsol = gusCHsol - gusOHsol
    plt.plot(akeOHsol, akeCOsol, 'wo', ms=5)   # in solar terms
    plt.plot(gusOHsol, gusCOsol, 'ws', ms=5)   # in solar terms
else:
    plt.plot(akeOH, akeCO, 'wo')  
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
    plt.plot(navOHsol, navCOsol, 'r+')   # in solar terms
    for idx, _ in enumerate(navOH):
        if idx==0 or idx ==1:
            OHsolerr = numpy.sqrt(navOHerr[idx]**2 + OHsunerr**2)
            COsolerr = numpy.sqrt(navCOerr[idx]**2 + COsunerr**2)
            plt.errorbar(navOHsol[idx], navCOsol[idx], xerr=OHsolerr, yerr=COsolerr, fmt='r+', ecolor='k')   # with errors
else:
    plt.plot(navOH, navCO, 'r+')  
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
if no_errs:
    plt.plot(Clo, Clco, 'ro')
else:
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
plt.xticks(numpy.arange(6.8, 9.1, 0.1))
plt.yticks(numpy.arange(-1.8, yup, 0.2))
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
    if (z == 'mrk5') or (z == 'sbs1319') or (z == 'sbs1054'):
        subxcoord = 5
        subycoord = -2
        side = 'left'
    if (z == 'iras08339') or (z == 'Sun'):
        subycoord = -13
        side = 'left'
    if (z == 'iiizw107') or (z == 'ngc1741') or (z == 'ngc456') or (z == 'iras08208') or (z == 'arp252') or (z == 'Orion'):
        subxcoord = 3
        side = 'left'
    if (z == 'ngc346'):
        subxcoord = 4
        subycoord = -11
        side = 'left'
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
plt.title('C/O vs O/H from Cloudy')
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
    elif obj in objects_Te_literature:
        fmt='wo'
    elif (obj in objects_not_in_sample) or (obj in objects_ref):
        fmt='b^'
    plt.errorbar(OH[i], NO[i], xerr=OHerr[i], yerr=NOerr[i], fmt=fmt, ecolor='k')
plt.xlim(6.8, 9.0)
yup = -0.0
ylo = -2.1
plt.ylim(ylo, yup)
plt.xticks(numpy.arange(6.8, 9.1, 0.1))
plt.yticks(numpy.arange(ylo, yup+0.1, 0.1))
for x, y, z in zip(OH, NO, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    #print z, x, y
    subxcoord = -2
    subycoord = 5
    side = 'right'
    if (z == 'mrk5') or (z == 'mrk1087') or (z == 'arp252'):
        subxcoord = 7
        side = 'left'
    if (z == 'ngc6822') or (z == 'ngc456') or (z == 'ngc346'):
        subxcoord = 4
        subycoord = 4
        side = 'left'
    if (z == 'sbs0926') or (z == 'sbs0218') or (z == 'sbs1054') or (z == 'iras08208') or (z == 'pox4') or (z == 'tol1457') or (z == 'mrk1199') or (z == '30Dor'):
        subycoord = -12
    if (z == 'iiizw107') or (z == 'sbs1319') or (z == 'sbs0948') or (z == 'Sun') or (z == 'Orion') or (z == 'izw18'):
        subxcoord = 4
        subycoord = -14
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
    elif obj in objects_Te_literature:
        fmt='wo'
    elif (obj in objects_not_in_sample) or (obj in objects_ref):
        fmt='b^'
    plt.errorbar(OH[i], NeO[i], xerr=OHerr[i], yerr=NeOerr[i], fmt=fmt, ecolor='k')
plt.xlim(7.6, 9.0)
yup = -0.1
ylo = -1.2
plt.ylim(ylo, yup)
plt.xticks(numpy.arange(7.4, 9.0, 0.1))
plt.yticks(numpy.arange(ylo, yup, 0.1))
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
    if (z == 'iras08208') or (z == 'pox4') or (z == 'arp252') or (z == 'iras08339'):
        subxcoord = 4
        subycoord = 4
        side = 'left'
    if (z == 'sbs1054') or (z == 'ngc6822') or (z == 'sbs0948'):
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
    elif obj in objects_Te_literature:
        fmt='wo'
    elif (obj in objects_not_in_sample) or (obj in objects_ref):
        fmt='b^'
    plt.errorbar(OH[i], CN[i], xerr=OHerr[i], yerr=CNerr[i], fmt=fmt, ecolor='k')
plt.xlim(6.8, 9.0)
yup = 2.2
ylo = -0.7
plt.ylim(ylo, yup)
plt.xticks(numpy.arange(6.9, 9.0, 0.1))
plt.yticks(numpy.arange(ylo, yup, 0.2))
for x, y, z in zip(OH, CN, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    subxcoord = -5
    subycoord = 5
    side = 'right'
    if (z == 'sbs0948') or (z == 'ngc346'):
        subycoord = -10
    if (z == 'mrk960') or (z == 'pox4') or (z == 'sbs1319') or (z == 'iras08339') or (z == '30Dor'):
        subxcoord = 5
        subycoord = -12
        side = 'left'
    if (z == 'sbs0948') or (z == 'arp252') or (z == 'tol9') or (z == 'izw18'):
        subxcoord = 5
        subycoord = 4
        side = 'left'
    if (z == 'ngc6822'):
        subxcoord = 4
        subycoord = -14
        side = 'left'
    if (z == 'sbs0926') or (z == 'iiizw107'):
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
        fmt='wo'
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

# N/C vs C/H
# Calculate errors for N/C and C/H
CH = CO+OH
errCH = numpy.sqrt( COerr**2 + OHerr**2 )
NC = (NO+OH) - (CO+OH)
errNC = numpy.sqrt( NOerr**2 + COerr**2 + 2*OHerr**2)
# Adjust a linear fit to the plot
coeffs, line_fit = fit_line(CH, NC)
print 'Coefficients of initial guess to the plot of: N/C vs C/H'
m = coeffs[0]
b = coeffs[1]
print 'm = %0.3f     b = %0.3f' % (m, b)
fig1 = plt.figure(1, figsize=(12, 10))
for obj, i in zip(objects_list, indeces_list):
    if obj in objects_with_Te:
        fmt='ko'
    elif obj in objects_Te_literature:
        fmt='wo'
    elif (obj in objects_not_in_sample) or (obj in objects_ref):
        fmt='b^'
    #plt.plot(CH[i], NC[i], 'ko')
    plt.errorbar(CH[i], NC[i], xerr=errCH[i], yerr=errNC[i], fmt=fmt, ecolor='k')
plt.xlim(5.6, 9.5)
yup = 0.7
ylo = -2.6
plt.ylim(ylo, yup)
plt.xticks(numpy.arange(5.6, 9.5, 0.2))
plt.yticks(numpy.arange(ylo, yup+0.2, 0.2))
for x, y, z in zip(CH, NC, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    subxcoord = -5
    subycoord = 5
    side = 'right'
    #print z, x, y
    if (z == 'sbs1319'):
        subxcoord = -42
        subycoord = 0
    if (z == 'sbs0926'):
        subxcoord = -42
        subycoord = -5
    if (z == 'sbs1054'):
        subxcoord = -5
        subycoord = -12
    if (z == 'mrk1199'):
        subxcoord = 4
        subycoord = -14
        side = 'left'
    if (z == 'iras08339'):
        subxcoord = 4
        side = 'left'
    if (z == 'ngc1741') or (z == 'tol9') or (z == 'Orion'):
        subycoord = -14 
    if  (z == 'sbs0948'):
        subxcoord = 5
        subycoord = 5
        side = 'left'
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
plt.title('N/C vs C/H')
plt.xlabel('12 + log (C/H)')
plt.ylabel('log (N/C)')
if save_images:
    img_name = img_name6 + typeofimage
    if use_our_sample_ONLY == False:
        img_name = img_name6 + otherrefs + typeofimage
    destination = os.path.join(full_results_path+'/plots', img_name)
    plt.savefig(destination)
    print('Plot %s was saved!' % destination)
if show_fitted_line:
    plt.plot(CH, line_fit)   # plot fitted line
# Initialize the chain with a Gaussian ball around the maximum likelihood result, for which we use optimize
x = CH
y = NC
yerr = CNerr
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
for theta in pos:
    if theta[0] < 0.0:
        allm.append(theta[0])
        allb.append(theta[1])
        subsample = numpy.vstack((subsample, theta))
avgm = sum(allm)/len(allm)
avgb = sum(allb)/len(allb)
print 'Average values:  m = %0.3f   b= %0.3f' % (avgm, avgb)
theta = [avgm, avgb]
if show_fitted_line:
    plt.plot( x, line_eq( p, x ), "r", lw=5, alpha=0.4 )
samples = sampler.chain[:, nruns*0.2, :].reshape((-1, ndim))
fig = triangle.corner(samples, labels=["$m$", "$b$"], truths=[m, b])
fig.show()
# Calculate the uncertainties based on the 25, 50, and 75th percentiles
m_mcmc, b_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*numpy.percentile(subsample, [25, 50, 75], axis=0)))
print 'm_mcmc, b_mcmc:', m_mcmc, b_mcmc
m_mcmc2, b_mcmc2 = map(lambda v: (v), zip(*numpy.percentile(subsample, [25, 50, 75], axis=0)))
print 'm_mcmc2, b_mcmc2:', m_mcmc2, b_mcmc2
plt.show()


# N/O vs C/N
# Adjust a linear fit to the plot
coeffs, line_fit = fit_line(NO, CN)
#print 'Linear fit: y = mx + b'
#print line_fit
print 'Coefficients of initial guess to the plot of: N/O vs C/N'
m = coeffs[0]
b = coeffs[1]
print 'm = %0.3f     b = %0.3f' % (m, b)
fig1 = plt.figure(1, figsize=(12, 10))
for obj, i in zip(objects_list, indeces_list):
    if obj in objects_with_Te:
        fmt='ko'
    elif obj in objects_Te_literature:
        fmt='wo'
    elif (obj in objects_not_in_sample) or (obj in objects_ref):
        fmt='b^'
    #plt.plot(NO[i], CN[i], 'ko')   # plot points without errors
    plt.errorbar(NO[i], CN[i], xerr=NOerr[i], yerr=CNerr[i], fmt=fmt, ecolor='k')
if show_fitted_line:
    plt.plot(NO, line_fit)   # plot fitted line
plt.xlim(-2.1, 0.0)
yup = 2.60
ylo = -0.8
plt.ylim(ylo, yup)
plt.xticks(numpy.arange(-2.1, 0.1, 0.2))
plt.yticks(numpy.arange(ylo, yup+0.2, 0.2))
for x, y, z in zip(NO, CN, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    #print z, x, y
    subxcoord = -2
    subycoord = 5
    side = 'right'
    if (z == 'sbs1319') or (z == 'pox4') or (z == 'Sun'):
        subycoord = -12
    if (z == 'sbs0926') or (z == 'iiizw107') or (z == '30Dor') or (z == 'sbs1415') or (z == 'iras08339'):
        subxcoord = 4
        subycoord = -12
        side = 'left'
    if (z == 'arp252') or (z == 'sbs1054'):
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
    plt.savefig(destination)
    print('Plot %s was saved!' % destination)
# Initialize the chain with a Gaussian ball around the maximum likelihood result, for which we use optimize
ndim, nwalkers, nruns = 2, 100, 500
x = NO
y = CN 
yerr = CNerr
randadd2point = lambda x: x+numpy.random.rand(1)
p0 = numpy.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[x, y, yerr])
pos, prob, state = sampler.run_mcmc(p0, nruns)   
# best model
wh = numpy.where( prob == prob.max() )[0][0]
p = pos[ wh, : ]
#theta = [avgm, avgb]
#theta = [m, b]
if show_fitted_line:
    plt.plot( x, line_eq( p, x ), "r", lw=5, alpha=0.4 )
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
samples = sampler.chain[:, nruns*0.2, :].reshape((-1, ndim))
fig = triangle.corner(samples, labels=["$m$", "$b$"], truths=[m, b])
fig.show()
# Calculate the uncertainties based on the 25, 50, and 75 percentiles
m_mcmc, b_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*numpy.percentile(subsample, [25, 50, 75], axis=0)))
print 'm_mcmc, b_mcmc:', m_mcmc, b_mcmc
m_mcmc2, b_mcmc2 = map(lambda v: (v), zip(*numpy.percentile(subsample, [25, 50, 75], axis=0)))
print 'm_mcmc2, b_mcmc2:', m_mcmc2, b_mcmc2
plt.show()


print ' Code finished!'
