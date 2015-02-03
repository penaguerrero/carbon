import numpy
import os
import string
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

use_our_sample_ONLY = False
save_images = False
use_our_sample_ONLY = True
save_images = True
# Do you want to correct values?
correct_values = False
# Type of image to be saved?
typeofimage = '.eps'

############################################################################################################################################

img_name1 = 'COvsOH'
img_name1b = 'COvsOH_Cldy'
img_name2 = 'NOvsOH'
img_name2b = 'NCvsCH'
img_name3 = 'NeOvsOH'
img_name4 = 'CNvsOH'
img_name5 = 'C+2N+2vsOH'
otherrefs = '_plusOtherRefs'

# Go into the directory of each object. Files assumed to be in: /Users/home_direcotry/Documents/AptanaStudio3/src/
results_path = "../results"
# but just to make sure we are in the right place, get and work with full paths
full_results_path = os.path.abspath(results_path)
# read the results file
if use_our_sample_ONLY:
    resfile = os.path.join(full_results_path, 'abunds_OCNNe.txt')
    rows2skip = 5
    rf = numpy.loadtxt(resfile, skiprows=rows2skip, usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), unpack=True)
    OH, OHerr, CO, COerr, NO, NOerr, NeO, NeOerr, NeOor, NeOorerr, C2O2, C2O2err, XO2, O2O1, O2O1err, N2O2, N2O2err, Z, ifc, ifcerr = rf
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

CN = CO - NO
CNerr = []
for coe, noe in zip(COerr, NOerr):
    cne = numpy.sqrt(coe**2 + noe**2)
    CNerr.append(cne)

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
    i = i+1
for obj in objects_ref:
    objects_list.append(obj)
    indeces_list.append(i)
    i = i+1
if use_our_sample_ONLY == False:
    for obj in objects_not_in_sample:
        objects_list.append(obj)
        indeces_list.append(i)
        i = i+1

# Abundances obtained with Cloudy simulations
Cldy_best_abunds = 'carbon/results/Cloudy_best_abunds.txt'
abunds_path_list = string.split(os.getcwd(), sep='carbon')
fname = os.path.join(abunds_path_list[0], Cldy_best_abunds)
Clhe, Clo, Clc, Cln, Clne, Cls, Clto3, Clto2 = numpy.loadtxt(fname, skiprows=1, usecols=(2,3,4,5,6,7,8,9), unpack=True)


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
plt.xlim(7.0, 9.0)
yup = 1.02
ylo = -1.8
plt.ylim(ylo, yup)
plt.xticks(numpy.arange(7.0, 9.0, 0.1))
plt.yticks(numpy.arange(ylo, yup, 0.2))
for x, xe, y, ye, z in zip(OH, OHerr, CO, COerr, objects_list):
    subxcoord = -3
    subycoord = 5
    side = 'right'
    if (z == 'mrk960') or (z == 'mrk5') or (z == 'sbs1319') or (z == 'ngc1741') or (z == 'Sun'):
        subxcoord = -4
        subycoord = -12
    if (z == 'iiizw107') or (z == 'sbs0218') or (z == 'sbs0948') or (z == 'iras08208') or (z == 'iras08339') or (z == 'arp252') or (z == 'pox4') or (z == 'ngc456') or (z == 'ngc6822') or (z == 'Orion'):
        subxcoord = 5
        side = 'left'
    if (z == 'ngc346'):
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

# C/O vs O/H from SIMULATIONS
fig1 = plt.figure(1, figsize=(12, 10))
for obj, i in zip(objects_list, indeces_list):
    if obj in objects_with_Te:
        fmt='ko'
    elif obj in objects_Te_literature:
        fmt='wo'
    elif (obj in objects_not_in_sample) or (obj in objects_ref):
        fmt='b^'
    #plt.errorbar(OH[i], CO[i], xerr=OHerr[i], yerr=COerr[i], fmt=fmt, ecolor='k')
    plt.plot(Clo, Clc-Clo, 'ro')
plt.xlim(7.0, 9.0)
yup = 1.02
ylo = -1.8
plt.ylim(ylo, yup)
plt.xticks(numpy.arange(7.0, 9.0, 0.1))
plt.yticks(numpy.arange(ylo, yup, 0.2))
for x, xe, y, ye, z in zip(Clo, OHerr, Clc-Clo, COerr, objects_list):
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
        img_name = img_name1 + otherrefs + typeofimage
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
plt.xlim(7.0, 9.0)
yup = -0.5
ylo = -1.8
plt.ylim(ylo, yup)
plt.xticks(numpy.arange(7.0, 9.0, 0.2))
plt.yticks(numpy.arange(ylo, yup, 0.2))
for x, y, z in zip(OH, NO, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    #print z, x, y
    subxcoord = -2
    subycoord = 5
    side = 'right'
    if (z == 'mrk5'):
        subxcoord = 7
        side = 'left'
    if (z == 'Sun') or (z == 'ngc6822') or (z == 'ngc456') or (z == 'ngc346'):
        subxcoord = 4
        subycoord = 4
        side = 'left'
    if (z == 'sbs0926') or (z == 'sbs0218') or (z == 'pox4') or (z == 'tol1457') or (z == 'mrk1199') or (z == '30Dor'):
        subycoord = -12
    if (z == 'iiizw107') or (z == 'sbs1319') or (z == 'mrk1087') or (z == 'iras08339') or (z == 'sbs0948') or (z == 'Orion'):
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

# N/C vs C/H
# Calculate errors for N/C and C/H
CH = CO+OH
errCH = numpy.sqrt( COerr**2 + OHerr**2 )
NC = (NO+OH) - (CO+OH)
#nc_errup = 10**(NO+OH - 12.0) - 10**()
errNC = numpy.sqrt( NOerr**2 + OHerr**2 + COerr**2 + OHerr**2)
#            O2err_up = 10**(logO2+O2logerr - 12.0) - 10**(logO2 - 12.0)
#            O2err_down = 10**(logO2 - 12.0) - 10**(logO2-O2logerr - 12.0)
#            O2err = (O2err_up + O2err_down) / 2.0
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
plt.xlim(5.5, 9.0)
yup = 0.7
ylo = -2.6
plt.ylim(ylo, yup)
plt.xticks(numpy.arange(5.5, 9.0, 0.2))
plt.yticks(numpy.arange(ylo, yup, 0.2))
for x, y, z in zip(CH, NC, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    #print z, x, y
    subxcoord = -2
    subycoord = 5
    side = 'right'
    if (z == 'sbs0948') or (z == 'sbs1054') or (z == 'sbs1319') or (z == 'sbs0218') or (z == 'Orion'):
        subxcoord = 4
        subycoord = -14
        side = 'left'
    if (z == 'sbs0926'):
        subycoord = -14        
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
plt.title('N/C vs C/H')
plt.xlabel('12 + log (C/H)')
plt.ylabel('log (N/C)')
if save_images:
    img_name = img_name2 + typeofimage
    if use_our_sample_ONLY == False:
        img_name = img_name2b + otherrefs + typeofimage
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
plt.xlim(7.0, 9.0)
yup = -0.2
ylo = -1.0
plt.ylim(ylo, yup)
plt.xticks(numpy.arange(7.0, 9.0, 0.2))
plt.yticks(numpy.arange(ylo, yup, 0.2))
for x, y, z in zip(OH, NeO, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    #print z, x, y
    subxcoord = -5
    subycoord = 5
    side = 'right'
    if (z == 'ngc1741') or (z == 'pox4') or (z == 'sbs0948') or (z == 'ngc346'):
        subycoord = -10
    if (z == 'mrk960') or (z == 'sbs1319') or (z == 'iiizw107') or (z == 'mrk5'):
        subxcoord = 5
        subycoord = 4
        side = 'left'
    if (z == 'sbs1054') or (z == 'ngc6822'):
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
plt.xlim(7.0, 9.0)
#yup = -0.2
#ylo = -1.2
#plt.ylim(ylo, yup)
plt.xticks(numpy.arange(7.0, 9.0, 0.2))
#plt.yticks(numpy.arange(ylo, yup, 0.2))
for x, y, z in zip(OH, CN, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    #print z, x, y
    subxcoord = -5
    subycoord = 5
    side = 'right'
    if (z == 'ngc1741') or (z == 'pox4') or (z == 'sbs0948') or (z == 'ngc346') or (z == '30Dor'):
        subycoord = -10
    if(z == 'sbs1319') or (z == 'iiizw107') or (z == 'mrk5') or (z == 'sbs0948') or (z == 'arp252'):
        subxcoord = 5
        subycoord = 4
        side = 'left'
    if (z == 'ngc6822'):
        subxcoord = 4
        subycoord = -14
        side = 'left'
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

print ' Code finished!'
