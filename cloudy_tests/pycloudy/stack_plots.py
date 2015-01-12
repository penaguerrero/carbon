from __future__ import division
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import string
import emcee
import triangle
import argparse
import itertools
import time
from collections import OrderedDict


'''
This script creates the following output for ALL objects in the sample:
            - plot or contours of temperature of [O III] vs C/O
            - plot or contours of C/O vs O/H
            - plot or contours of C/H vs O/H
            - plot or contours of temperature of [O III] vs C/H 
            - triangle plot after the MCMC chain has been ran
'''

#######################################################################################################################


# Specify if you want to save the plots and the type of image to be saved
save_figs = True
img_format = '.jpg'

# Do you want to see the plots?
show_plots = False


#####################################################################################################################

# Sample objects
#                  0            1           2          3         4        5         6         7         8
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                 9           10         11         12         13       14        15         16         17

#### FUNCTIONS

def get_Te_percentiles(percentiles_arr, Te_arr):
    # Calculate percentiles at 16, 50, and 80
    Te_pencentiles = []
    theta_p16 = np.array([[percentiles_arr[0][0], percentiles_arr[1][0], percentiles_arr[2][0], percentiles_arr[3][0], percentiles_arr[4][0], percentiles_arr[5][0]]])
    theta_p50 = np.array([[percentiles_arr[0][1], percentiles_arr[1][1], percentiles_arr[2][1], percentiles_arr[3][1], percentiles_arr[4][1], percentiles_arr[5][1]]])
    theta_p80 = np.array([[percentiles_arr[0][2], percentiles_arr[1][2], percentiles_arr[2][2], percentiles_arr[3][2], percentiles_arr[4][2], percentiles_arr[5][2]]])
    wp16 = np.where(samples == theta_p16)[0][0]
    wp50 = np.where(samples == theta_p50)[0][0]
    wp80 = np.where(samples == theta_p80)[0][0]
    Tep16 = Te_arr[wp16]
    Tep50 = Te_arr[wp50]
    Tep80 = Te_arr[wp80]
    Te_pencentiles.append(Tep16)
    Te_pencentiles.append(Tep50)
    Te_pencentiles.append(Tep80)
    return Te_pencentiles

def get_percentiles(samples, TO3, TO2):
    # Calculate the uncertainties based on the 16th, 50th and 84th percentiles
    percentiles0 = 'MCMC values and uncertainties according to 16th, 50th, and 84th percentiles:'
    p_mcmc1 = map(lambda v: (v), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    p_mcmc2 = map(lambda v: (v[1], np.abs(v[2]-v[1]), np.abs(v[1]-v[0])), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    # find temperatures of these sets:
    TO3_perc = get_Te_percentiles(p_mcmc1, TO3)
    TO2_perc = get_Te_percentiles(p_mcmc1, TO2)
    perc = [p_mcmc1[0], p_mcmc1[1], p_mcmc1[2], p_mcmc1[3], p_mcmc1[4], p_mcmc1[5], TO3_perc, TO2_perc]
    percentiles1 = '{:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:>8} {:<9} \n {:>17.2f} {:>7.2f} {:>6.2f} {:>6.2f}'.format(
                     p_mcmc1[0][0], p_mcmc1[1][0], p_mcmc1[1][0], p_mcmc1[3][0], p_mcmc1[4][0], p_mcmc1[5][0], int(TO3_perc[0]), int(TO2_perc[0]),
                     p_mcmc1[2][0]+p_mcmc1[1][0], p_mcmc1[3][0]+p_mcmc1[1][0], p_mcmc1[4][0]+p_mcmc1[1][0], p_mcmc1[5][0]+p_mcmc1[1][0] )
    percentiles2 = '{:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:>8} {:<9} \n {:>17.2f} {:>7.2f} {:>6.2f} {:>6.2f}'.format(
                     perc[0][1], perc[1][1], perc[2][1], perc[3][1], perc[4][1], perc[5][1], int(perc[6][1]), int(perc[7][1]),
                     perc[2][1]+perc[1][1], perc[3][1]+perc[1][1], perc[4][1]+perc[1][1], perc[5][1]+perc[1][1] )
    percentiles3 = '{:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:>8} {:<9} \n {:>17.2f} {:>7.2f} {:>6.2f} {:>6.2f}'.format(
                     perc[0][2], perc[1][2], perc[2][2], perc[3][2], perc[4][2], perc[5][2], int(perc[6][2]), int(perc[7][2]),
                     perc[2][2]+perc[1][2], perc[3][2]+perc[1][2], perc[4][2]+perc[1][2], perc[5][2]+perc[1][2] )
    print '\n'
    print percentiles0
    print 'He      O       C      N      Ne      S      Te_O3   Te_O2'
    print '16th percentile:'
    print percentiles1
    print '50th percentile:'
    print percentiles2
    print '80th percentile:'
    print percentiles3
    print '\n DIFFERENCES of percentiles with respect to 50th percentile: '
    print p_mcmc2

def get_averages(He, O, CO, NO, NeO, SO, TO3, TO2):
    # Calculate the averages for all models
    avgHe = np.mean(He) 
    avgO = np.mean(O)
    avgCO = np.mean(CO+O)
    avgNO = np.mean(NO+O)
    avgNeO = np.mean(NeO+O)
    avgSO = np.mean(SO+O)
    avgTO3 = np.mean(TO3)
    avgTO2 = np.mean(TO2)
    avgabunds0 = 'Average abundances: '
    avgabunds1 = 'He = %0.2f   O = %0.2f   C/O = %0.2f   N/O = %0.2f   Ne/O = %0.2f   S/O = %0.2f   Te_O3 = %d   Te_O2 = %d' % (
                                   avgHe, avgO, avgCO-avgO, avgNO-avgO, avgNeO-avgO, avgSO-avgO, int(avgTO3), int(avgTO2))
    avgabunds2 = 'He = %0.2f   O = %0.2f     C = %0.2f      N = %0.2f      Ne = %0.2f      S = %0.2f ' % (
                                                            avgHe, avgO, avgCO, avgNO, avgNeO, avgSO )
    print '\n'+avgabunds0
    print avgabunds1
    print avgabunds2

def unique(a, arrange_TsProbs=None):
    order = np.lexsort(a.T)
    a = a[order]
    diff = np.diff(a, axis=0)
    ui = np.ones(len(a), 'bool')
    ui[1:] = (diff != 0).any(axis=1) 
    if arrange_TsProbs!=None:
        TO3, TO2, prob = arrange_TsProbs
        return a[ui], TO3[ui], TO2[ui], prob[ui]
    return a[ui]

def count_repetitions(samples):
    # create a list of samples in order to count the repetitions
    samples_list = []
    for item in samples:
        samples_list.append(item.tolist())  
    # count the repetitions and store this info in a dictionary
    repetitions = {}
    for item in samples_list:
        r = samples_list.count(item)
        repetitions[str(item)] = r
    ordered_repeats = sorted(repetitions.values(), reverse=True)   # get the list of the number of times each set is repeated
    ordered_repeats_sets = sorted(repetitions, key=repetitions.__getitem__, reverse=True)   # get the list of set repetitions
    ordered_repetitions = OrderedDict()
    for r, set in zip(ordered_repeats, ordered_repeats_sets):
        ordered_repetitions[str(set)] = r
    # now make an array of the same length as the samples with ONLY the number of times each set is repeated
    frequency_arr = np.array([])
    for set in samples:
        set = str(set.tolist())
        if set in ordered_repeats_sets:
            i = ordered_repeats_sets.index(set)
            value = ordered_repeats[i]
            #print set, 'is repeated ', value, 'times'
            frequency_arr = np.append(frequency_arr, value)
    return ordered_repetitions, frequency_arr

def get_zarr(x, y):
    # calculate the z array for the contour plots    
    z = np.array([]).reshape(0, 2)
    for xi, yi in zip(x, y):
        zi = [xi, yi]
        z = np.vstack((z, zi))
    return x, y, z


#### CODE

# start the timer to compute the whole running time
start_time = time.time()

# Load the chain and create the samples
not_final = True
study_NO_repetitions = False

allhe, allo, allco, allno, allneo, allso, allto3, allto2 = [], [], [], [], [], [], [], []
for object_name in objects_list:
    chain_file = os.path.abspath('mcmc_'+object_name+'_chain.dat')
    if not_final:
        chain_file = os.path.abspath('mcmc_'+object_name+'_chain0.dat')
    print 'Now reading file: ', chain_file
    He, O, CO, NO, NeO, SO, TO3, TO2, prob = np.loadtxt(chain_file, skiprows=17, unpack=True)
    for i, _ in enumerate(He):
        allhe.append(He[i])
        allo.append(O[i])
        allco.append(CO[i])
        allno.append(NO[i])
        allneo.append(NeO[i])
        allso.append(SO[i])
        allto3.append(TO3[i])
        allto2.append(TO2[i])
        
    if study_NO_repetitions:
        # count the repetitions
        ordered_repetitions, frequency_arr = count_repetitions(samples)
        # eliminate the possible repetitions of models for personal curiosity
        arrange_TsProbs = [TO3, TO2, prob]
        subsamples, subTO3, subTO2, subprob = unique(samples, arrange_TsProbs)
        imax = 20
        print 'The top ', imax, ' repeated sets are the following:'
        x = itertools.islice(ordered_repetitions.items(), 0, imax)
        for key, value in x:
            print key, 'is repeated ', value, 'times'
        print 'The SUBsamples* shape is: ', np.shape(subsamples)
        print '* SUBsamples is the set of models without repeated abundance and temperature sets.'
        
# Create the samples for ALL the sample    
print 'lengths of the arrays to use: '
print 'He, O, CO, NO, NeO, SO, TO3, TO2: ', len(allhe), len(allo), len(allco), len(allno), len(allneo), len(allso), len(allto3), len(allto2)
print '\nCreating the arrays for the plots...'
allhe = np.array(allhe)
allo = np.array(allo)
allco = np.array(allco)
allno = np.array(allno)
allneo = np.array(allneo)
allso = np.array(allso)
allto3 = np.array(allto3)
allto2 = np.array(allto2)
samples = np.array([]).reshape(0, 6)
for he, o, co, no, neo, so in zip(allhe, allo, allco, allno, allneo, allso):
    theta = [he, o, co, no, neo, so]
    samples = np.vstack((samples, theta))
if study_NO_repetitions:
    print 'ok, got the samples WITHOUT repetitions! The shape of array is: ', np.shape(samples)
else:
    print 'ok, got the samples! shape of array is: ', np.shape(samples)

# STATISTICS
get_averages(allhe, allo, allco, allno, allneo, allso, allto3, allto2)
get_percentiles(samples, allto3, allto2)

# PLOTS
print '\n READY TO PLOT!'    

# plots without the benchmark values
object_name = 'WholeSample'

fig = plt.figure(1, figsize=(12, 10))
plt.title('C/O vs O3_Temps')
xlab = 'Te_[OIII]'
ylab = 'log (C/O)'
plt.xlabel(xlab)
plt.ylabel(ylab)
x, y, z = get_zarr(allto3, allco)
fig = triangle.corner(z, labels=[xlab, ylab], extents=[(7000.0, 18000.0), (-1.6, 1.7)])
if save_figs:
    fn = object_name+'_tempsVsCO.jpg'
    fig.savefig(os.path.abspath(fn))
    print 'Figure ', fn, 'saved!'
if show_plots:
    plt.plot(x, y, 'k.')
    plt.show()

fig = plt.figure(1, figsize=(12, 10))
plt.title('C/O vs O2_Temps')
xlab = 'Te_[OII]'
ylab = 'log (C/O)'
plt.xlabel(xlab)
plt.ylabel(ylab)
x, y, z = get_zarr(allto2, allco)
fig = triangle.corner(z, labels=[xlab, ylab], extents=[(7000.0, 18000.0), (-1.6, 1.7)])
if save_figs:
    fn = object_name+'_tempO2VsCO.jpg'
    fig.savefig(os.path.abspath(fn))
    print 'Figure ', fn, 'saved!'
if show_plots:
    plt.plot(x, y, 'k.')
    plt.show()

fig = plt.figure(1, figsize=(12, 10))
plt.title('C/O vs O/H')
xlab = '12 + log (O/H)'
ylab ='log (C/O)'
plt.xlabel(xlab)
plt.ylabel(ylab)
x, y, z = get_zarr(allo, allco)
fig = triangle.corner(z, labels=[xlab, ylab], extents=[(7.5, 8.7), (-1.6, 1.7)])
if save_figs:
    fn = object_name+'_COvsOH.jpg'
    fig.savefig(os.path.abspath(fn))
    print 'Figure ', fn, 'saved!' 
if show_plots:
    plt.plot(x, y, 'k.')
    plt.show()

fig = plt.figure(1, figsize=(12, 10))
plt.title('C/H vs O/H')
xlab= '12 + log (O/H)'
ylab ='12 + log (C/H)'
plt.xlabel(xlab)
plt.ylabel(ylab)
ymin = 6.0
ymax = 8.7
plt.ylim(ymin, ymax)
x, y, z = get_zarr(allo, allco+allo)
fig = triangle.corner(z, labels=[xlab, ylab], extents=[(7.5, 8.7), (ymin, ymax)])
if save_figs:
    fn = object_name+'_CHvsOH.jpg'
    fig.savefig(os.path.abspath(fn))
    print 'Figure ', CHOH, 'saved!' 
if show_plots:
    plt.plot(x, y, 'k.')
    plt.show()

fig = plt.figure(1, figsize=(12, 10))
plt.title('C/H vs O3_Temps')
xlab = 'Te_[OIII]'
ylab = '12 + log (C/H)'
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.ylim(ymin, ymax)
x, y, z = get_zarr(allto3, allco+allo)
fig = triangle.corner(z, labels=[xlab, ylab], extents=[(7000.0, 18000.0), (ymin, ymax)])
if save_figs:
    fn = object_name+'_tempsVsCH.jpg'
    fig.savefig(os.path.abspath(fn))
    print 'Figure ', fn, 'saved!' 
if show_plots:
    plt.plot(x, y, 'k.')
    plt.show()

fig = plt.figure(1, figsize=(12, 10))
plt.title('C/H vs O2_Temps')
xlab = 'Te_[OII]'
ylab = '12 + log (C/H)'
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.ylim(ymin, ymax)
x, y, z = get_zarr(allto2, allco+allo)
fig = triangle.corner(z, labels=[xlab, ylab], extents=[(7000.0, 18000.0), (ymin, ymax)])
if save_figs:
    fn = object_name+'_tempO2VsCH.jpg'
    fig.savefig(os.path.abspath(fn))
    print 'Figure ', fn, 'saved!' 
if show_plots:
    plt.plot(x, y, 'k.')
    plt.show()

fig = plt.figure(1, figsize=(12, 10))
plt.title('N/H vs O/H')
xlab= '12 + log (O/H)'
ylab ='12 + log (N/H)'
plt.xlabel(xlab)
plt.ylabel(ylab)
ymin = 6.0
ymax = 8.7
plt.ylim(ymin, ymax)
x, y, z = get_zarr(allo, allno+allo)
fig = triangle.corner(z, labels=[xlab, ylab], extents=[(7.5, 8.7), (ymin, ymax)])
if save_figs:
    fn = object_name+'_NHvsOH.jpg'
    fig.savefig(os.path.abspath(fn))
    print 'Figure ', fn, 'saved!' 
if show_plots:
    plt.plot(x, y, 'k.')
    plt.show()

fig = plt.figure(1, figsize=(12, 10))
plt.title('Ne/H vs O/H')
xlab= '12 + log (O/H)'
ylab ='12 + log (Ne/H)'
plt.xlabel(xlab)
plt.ylabel(ylab)
ymin = 6.0
ymax = 8.7
plt.ylim(ymin, ymax)
x, y, z = get_zarr(allo, allneo+allo)
fig = triangle.corner(z, labels=[xlab, ylab], extents=[(7.5, 8.7), (ymin, ymax)])
if save_figs:
    fn = object_name+'_NeHvsOH.jpg'
    fig.savefig(os.path.abspath(fn))
    print 'Figure ', fn, 'saved!' 
if show_plots:
    plt.plot(x, y, 'k.')
    plt.show()

fig = plt.figure(1, figsize=(12, 10))
plt.title('S/H vs O/H')
xlab= '12 + log (O/H)'
ylab ='12 + log (S/H)'
plt.xlabel(xlab)
plt.ylabel(ylab)
ymin = 5.0
ymax = 7.0
plt.ylim(ymin, ymax)
x, y, z = get_zarr(allo, allso+allo)
fig = triangle.corner(z, labels=[xlab, ylab], extents=[(7.5, 8.7), (ymin, ymax)])
if save_figs:
    fn = object_name+'_SHvsOH.jpg'
    fig.savefig(os.path.abspath(fn))
    print 'Figure ', fn, 'saved!' 
if show_plots:
    plt.plot(x, y, 'k.')
    plt.show()


nwalkers = 100
nruns = 100
# plot of the ratios without the benchmark abundances
fig = triangle.corner(samples, 
                      labels=["$He$", "$O$", "$C/O$", "$N/O$", "$Ne/O$", "$S/O$"],
                      # Limits:    He           O            C/O          N/O           Ne/O         S/O
                      extents=[(9.5, 11.9), (7.5, 8.7), (-1.6, 1.7), (-1.7, -0.4), (-1.0, 0.01), (-2.3, -1.4)]                      
                      )
pltwithoutbench = 'mcmc_'+object_name+"_ratios2_"+repr(nwalkers)+"w"+repr(nruns)+"r"+"_initb.jpg"
if save_figs:
    fig.savefig(os.path.abspath(pltwithoutbench))   
    print 'Figure ', pltwithoutbench, 'saved!'  

time2run = 'Took  %s  minutes to finish reading and processing the data, and making the plots.' % ( 
                                                               (time.time() - start_time) / 60.0 )
print time2run
    
print '\nCode finished! \n'
