import numpy
import os
from matplotlib import pyplot as plt


############################################################################################################################################
'''
This scripts takes the found metallicities and produces:
 * plot of C/O vs Otot
 * plot of N/O vs Otot
 * plot of Ne/O vs Otot
 
 NOTE: The data for N/O and Ne/O was taken all from Lopez-Sanches & Esteban 2010
'''

use_our_sample_ONLY = False
save_images = False

############################################################################################################################################

# Go into the directory of each object. Files assumed to be in: /Users/home_direcotry/Documents/AptanaStudio3/src/
results_path = "../results"
# but just to make sure we are in the right place, get and work with full paths
full_results_path = os.path.abspath(results_path)
# read the results file
if use_our_sample_ONLY:
    resfile = os.path.join(full_results_path, 'abunds_OCNNe.txt')
    rows2skip = 5
else:
    resfile = os.path.join(full_results_path, 'abunds_OCNNe_plusotherrefs.txt')
    rows2skip = 10
rf = numpy.loadtxt(resfile, skiprows=rows2skip, usecols=(1,2,3,4,5,6,7,8), unpack=True)
OH, OHerr, CO, COerr, NO, NOerr, NeO, NeOerr = rf
# These are the results are for the sample objects in the following order:
ordered_sample = ['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
                  'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
# now divide the sample in objects for in which we found Te and those in which we did not
objects_with_Te = ['pox4', 'sbs0948', 'sbs1054', 'tol1457', 'iras08208']
objects_Te_literature = ['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'sbs0218',
                         'sbs0926', 'sbs1319', 'tol9', 'arp252', 'sbs1415']

# objects not in the sample
objects_ref = ['30Dor', 'Orion', 'izw18', 'Sun'] # these are useful points of reference
# additional data points not in the sample
objects_not_in_sample = ['ngc346', 'ngc456', 'ngc6822']

# Complie all the names for the labeling in the same order
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

# PLOTS
# C/O vs O/H
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
    subxcoord = -4
    subycoord = 5
    side = 'right'
    if (z == 'pox4') or (z == 'sbs1054'):
        subxcoord = -8
        subycoord = -7
    if (z == 'mrk960') or (z == 'sbs1319') or (z == 'ngc1741') or (z == 'Sun'):
        subycoord = -15
    if (z == 'ngc456') or (z == 'sbs0948') or (z == 'iras08208') or (z == 'iras08339') or (z == 'arp252') or (z == 'Orion'):
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
    img_name = 'COvsOH.jpg'
    if use_our_sample_ONLY == False:
        img_name = 'COvsOH_plusOtherRefs.jpg'
    destination = os.path.join(full_results_path+'/plots', img_name)
    plt.savefig(destination)
    print('Plot %s was saved!' % destination)
plt.show()

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
    if (z == 'sbs0926') or (z == 'sbs0948') or (z == 'sbs0218') or (z == 'pox4') or (z == 'tol1457') or (z == 'mrk1199') or (z == '30Dor'):
        subycoord = -12
    if (z == 'iiizw107') or (z == 'sbs1319') or (z == 'mrk1087') or (z == 'iras08339') or (z == 'Orion'):
        subxcoord = 4
        subycoord = -14
        side = 'left'
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
plt.title('N/O vs O/H')
plt.xlabel('12 + log (O/H)')
plt.ylabel('log (N/O)')
if save_images:
    img_name = 'NOvsOH.jpg'
    if use_our_sample_ONLY == False:
        img_name = 'NOvsOH_plusOtherRefs.jpg'
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
ylo = -1.2
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
    if (z == 'mrk960') or (z == 'sbs1319') or (z == 'iiizw107') or (z == 'mrk5') or (z == 'Orion'):
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
    img_name = 'NeOvsOH.jpg'
    if use_our_sample_ONLY == False:
        img_name = 'NeOvsOH_plusOtherRefs.jpg'
    destination = os.path.join(full_results_path+'/plots', img_name)
    plt.savefig(destination)
    print('Plot %s was saved!' % destination)
plt.show()

print ' Code finished!'