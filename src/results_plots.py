import numpy
import os
from matplotlib import pyplot as plt


############################################################################################################################################
'''
This scripts takes the found metallicities and produces:
 * plot of C/O vs Otot
 * plot of N/O vs Otot
 * plot of Ne/O vs Otot
'''

# Go into the directory of each object. Files assumed to be in: /Users/home_direcotry/Documents/AptanaStudio3/src/
results_path = "../results/"
# but just to make sure we are in the right place, get and work with full paths
full_results_path = os.path.abspath(results_path)
# read the results file
resfile = os.path.join(full_results_path, 'abunds_OCNNe.txt')
rf = numpy.loadtxt(resfile, skiprows=1, usecols=(1,2,3,4,5,6,7,8), unpack=True)
OH, OHerr, CO, COerr, NO, NOerr, NeO, NeOerr = rf
# These are the results are for the sample objects in the following order:
#                objects in the sample 
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415',
               # objects not in the sample
               '30Dor', 'Orion', 'izw18', 'Sun']

# PLOTS
# C/O vs O/H
fig1 = plt.figure(1, figsize=(12, 10))
plt.errorbar(OH, CO, xerr=OHerr, yerr=COerr, fmt='ko')
plt.xlim(7.0, 9.0)
yup = 1.0
ylo = -2.0
plt.ylim(ylo, yup)
plt.xticks(numpy.arange(7.0, 9.0, 0.1))
plt.yticks(numpy.arange(ylo, yup, 0.5))
for x, y, z in zip(OH, CO, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    #print z, x, y
    subxcoord = -5
    subycoord = 5
    side = 'right'
    if (z == 'mrk960') or (z == 'sbs1054'):
        subycoord = -15
    if (z == 'iras08208'):
        subxcoord = 8
        side = 'left'
    if (z == 'sbs1319') or (z == 'iiizw107') or (z == 'Orion'):
        subxcoord = 8
        subycoord = -15
        side = 'left'
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
plt.title('C/O vs O/H')
plt.xlabel('12 + log (O/H)')
plt.ylabel('log (C/O)')
img_name = 'COvsOH.jpg'
destination = os.path.join(full_results_path, img_name)
plt.savefig(destination)
print('Plot %s was saved!' % destination)
plt.show()

# N/O vs O/H
fig1 = plt.figure(1, figsize=(12, 10))
plt.errorbar(OH, NO, xerr=OHerr, yerr=NOerr, fmt='ko')
plt.xlim(7.0, 9.0)
plt.ylim(-2.0, 0.1)
plt.xticks(numpy.arange(7.0, 9.0, 0.2))
plt.yticks(numpy.arange(-2.0, 0.1, 0.5))
for x, y, z in zip(OH, NO, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    #print z, x, y
    subxcoord = -5
    subycoord = 5
    side = 'right'
    if (z == 'mrk5') or(z == 'iras08339'):
        subxcoord = 10
        side = 'left'
    if (z == 'iras08208'):
        subxcoord = 5
        subycoord = 10
        side = 'left'
    if (z == 'ngc1741') or (z == 'sbs0218') or (z == 'tol1457') or (z == '30Dor'):
        subycoord = -15
    if (z == 'iiizw107') or (z == 'pox4') or (z == 'mrk1199') or (z == 'sbs0948') or (z == 'Orion'):
        subxcoord = 8
        subycoord = -15
        side = 'left'
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
plt.title('N/O vs O/H')
plt.xlabel('12 + log (O/H)')
plt.ylabel('log (N/O)')
img_name = 'NOvsOH.jpg'
destination = os.path.join(full_results_path, img_name)
plt.savefig(destination)
print('Plot %s was saved!' % destination)
plt.show()

# Ne/O vs O/H
fig1 = plt.figure(1, figsize=(12, 10))
plt.errorbar(OH, NeO, xerr=OHerr, yerr=NeOerr, fmt='ko')
plt.xlim(7.0, 9.0)
plt.ylim(-4.0, 0.1)
plt.xticks(numpy.arange(7.0, 9.0, 0.2))
plt.yticks(numpy.arange(-4.0, 0.1, 0.5))
for x, y, z in zip(OH, NeO, objects_list):
    # Annotate the points 5 _points_ above and to the left of the vertex
    #print z, x, y
    subxcoord = -5
    subycoord = 5
    side = 'right'
    if (z == 'pox4') or (z == 'mrk960') or (z == 'sbs1319') or (z == 'Orion'):
        subxcoord = 5
        subycoord = 10
        side = 'left'
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
plt.title('Ne/O vs O/H')
plt.xlabel('12 + log (O/H)')
plt.ylabel('log (Ne/O)')
img_name = 'NeOvsOH.jpg'
destination = os.path.join(full_results_path, img_name)
plt.savefig(destination)
print('Plot %s was saved!' % destination)
plt.show()

print ' Code finished!'