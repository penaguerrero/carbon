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
from collections import OrderedDict


'''
This script creates the following output:
            - determines statistics of MCMC: averages, best fit with Chi2 minimization, and temperature subset
            - plot or contours of temperature of [O III] vs C/O
            - plot or contours of C/O vs O/H
            - plot or contours of C/H vs O/H
            - plot or contours of temperature of [O III] vs C/H 
            - triangle plot after the MCMC chain has been ran
            - write all statistics with chain into a text file
'''

#######################################################################################################################

parser = argparse.ArgumentParser(description='Count number of models ran so far.')
parser.add_argument("object_name",
                    action='store',
                    default=None,
                    help='The abbreviated name of the object, i.e. sbs1319.')
parser.add_argument("-n",
                    action='store_true',
                    dest="not_final",
                    default=False,
                    help='If the chain is not finished, use -n = not_final to read the correct chain file.')
parser.add_argument("-s",
                    action='store_true',
                    dest="use_subset",
                    default=False,
                    help='Use the subset of samples WITHOUT repetitions.')
parser.add_argument("-p",
                    action='store_true',
                    dest="mk_plots",
                    default=False,
                    help='Create and save the plots.')
parser.add_argument("-c",
                    action='store_true',
                    dest="contours",
                    default=False,
                    help='Create and save contours.')
parser.add_argument("-f",
                    action='store_true',
                    dest="wt_file",
                    default=False,
                    help='Write file with calculated statistics.')
args = parser.parse_args()

# name of the object
#                  0            1           2          3         4        5         6         7         8
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                 9           10         11         12         13       14        15         16         17

#### FUNCTIONS

def get_true_abunds(object_name):
    object_number = objects_list.index(object_name)
    print object_name, 'is object number ', object_number, ' in the sample list.'
    #fname = os.path.abspath('../../results/Final_meas_abunds.txt')
    Final_meas_abunds = 'carbon/results/Final_meas_abunds.txt'
    Final_meas_abunds_path_list = string.split(os.getcwd(), sep='carbon')
    fname = os.path.join(Final_meas_abunds_path_list[0], Final_meas_abunds)
    print '-->  looking to read ', fname, '...'
    he, o, c, n, ne, s, to3, to2 = np.loadtxt(fname, skiprows=1, usecols=(2,3,4,5,6,7,8,9), unpack=True)
    abunds = [he[object_number], o[object_number], c[object_number], n[object_number], ne[object_number], s[object_number]]
    true_temps = [to3[object_number], to2[object_number]]
    # set the ratios with respect to O, theta
    true_abunds = []
    true_abunds.append(abunds[0])
    true_abunds.append(abunds[1])
    for i in range(2, len(abunds)):
        true_abunds.append(abunds[i]-abunds[1])
    #print 'these are the benchmark abundances: \n', true_abunds
    return true_abunds, true_temps

def get_Te_percentiles(percentiles_arr, Te_arr):
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
    #print theta_p16, '\n', theta_p50, '\n', theta_p80
    #print Tep16, Tep50, Tep80
    Te_pencentiles.append(Tep16)
    Te_pencentiles.append(Tep50)
    Te_pencentiles.append(Tep80)
    return Te_pencentiles

def get_bestmodel(He, O, CO, NO, NeO, SO, TO3, TO2, prob, samples, true_abunds, true_temps):
    trueabs0 = '** Values of the BENCHMARK abundances: '
    trueabs1 = 'He = %0.2f   O = %0.2f   C/O = %0.2f   N/O = %0.2f   Ne/O = %0.2f   S/O = %0.2f' % (true_abunds[0], true_abunds[1], true_abunds[2], 
                                                                                                    true_abunds[3], true_abunds[4], true_abunds[5])
    trueabs2 = 'He = %0.2f   O = %0.2f     C = %0.2f      N = %0.2f      Ne = %0.2f      S = %0.2f' % (true_abunds[0], true_abunds[1], true_abunds[2]+true_abunds[1],
                                                                                           true_abunds[3]+true_abunds[1], true_abunds[4]+true_abunds[1], 
                                                                                           true_abunds[5]+true_abunds[1])
    temps0 = '** Benchmark temperatures: '
    temps1 = 'Te_O3 = %s     Te_O2 = %s' % (true_temps[0], true_temps[1])
    print '\n'+trueabs0
    print trueabs1
    print trueabs2
    print '\n'+temps0
    print temps1

    # best model from average
    #print He, '\n', O, '\n', CO, '\n', NO, '\n', NeO, '\n', SO
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
    
    # now find averages of those models with temperatures within +-1500 of the benchmark values
    subsamples = np.array([]).reshape(0, len(true_abunds))
    He_nearby, O_nearby, CO_nearby, NO_nearby, NeO_nearby, SO_nearby, TO3_nearby, TO2_nearby, prob_nearby = [],[],[],[],[],[],[],[], []
    nearby = 2500.0
    for he, o, co, no, neo, so, t3, t2, p in zip(He, O, CO, NO, NeO, SO, TO3, TO2, prob):
        diff_t3 = np.abs(t3 - true_temps[0])
        diff_t2 = np.abs(t2 - true_temps[1])
        if (diff_t3 <= nearby) and (diff_t2 <= nearby):
            theta = [he, o, co, no, neo, so]
            subsamples = np.vstack((subsamples, theta))
            #print 'these are the differences: ', t3, '-', true_temps[0], '=', diff_t3, '   ', t2, '-', true_temps[1], '=', diff_t2
            #print 'for ', theta, ' differences are: ', diff_t3, diff_t2
            #print ' differences are: ', diff_t3, diff_t2
            He_nearby.append(he)
            O_nearby.append(o)
            CO_nearby.append(co)
            NO_nearby.append(no)
            NeO_nearby.append(neo)
            SO_nearby.append(so)
            TO3_nearby.append(t3)
            TO2_nearby.append(t2)
            prob_nearby.append(p)
    avgHe = sum(He_nearby)/len(He_nearby) 
    avgO = sum(O_nearby)/len(He_nearby) 
    avgCO = sum(CO_nearby+O_nearby)/len(He_nearby) 
    avgNO = sum(NO_nearby+O_nearby)/len(He_nearby) 
    avgNeO = sum(NeO_nearby+O_nearby)/len(He_nearby) 
    avgSO = sum(SO_nearby+O_nearby)/len(He_nearby) 
    avgTO3 = sum(TO3_nearby)/len(He_nearby) 
    avgTO2 = sum(TO2_nearby)/len(He_nearby) 
    subset = str(len(O_nearby))
    subavgabunds0 = 'Average abundances with O temperatures within '+str(nearby)+' from the benchmark values (i.e. averages of '+subset+' models): '
    subavgabunds1 = 'He = %0.2f   O = %0.2f   C/O = %0.2f   N/O = %0.2f   Ne/O = %0.2f   S/O = %0.2f   Te_O3 = %d   Te_O2 = %d' % (
                                   avgHe, avgO, avgCO-avgO, avgNO-avgO, avgNeO-avgO, avgSO-avgO, int(avgTO3), int(avgTO2))
    subavgabunds2 = 'He = %0.2f   O = %0.2f     C = %0.2f      N = %0.2f      Ne = %0.2f      S = %0.2f ' % (
                                                                    avgHe, avgO, avgCO, avgNO, avgNeO, avgSO )
    print '\n'
    print subavgabunds0
    print subavgabunds1 
    print subavgabunds2
    
    # Chi^2 minimization within subsets
    whsub = np.where( prob_nearby == max(prob_nearby) )[0][0]
    psub = subsamples[ whsub, : ]
    TO3whsub = TO3_nearby[whsub]
    TO2whsub = TO2_nearby[whsub]
    subchi0 = 'Values that best fit the data  from Chi^2 minimization within then '+str(nearby)+'K temperature tolerance are the following:'
    subchi1 = 'He = %0.2f   O = %0.2f   C/O = %0.2f   N/O = %0.2f   Ne/O = %0.2f   S/O = %0.2f  Te_O3 = %d   Te_O2 = %d' % (
                                        psub[0], psub[1], psub[2],  psub[3], psub[4], psub[5], int(TO3whsub), int(TO2whsub) )
    subchi2 = 'He = %0.2f   O = %0.2f     C = %0.2f      N = %0.2f      Ne = %0.2f      S = %0.2f' % (psub[0], psub[1], 
                      psub[2]+true_abunds[1], psub[3]+true_abunds[1], psub[4]+true_abunds[1], psub[5]+true_abunds[1] )
    print '\n'
    print subchi0
    print subchi1
    print subchi2
    
    # best model from Chi^2 minimization
    nwalkers, ndim = np.shape(samples)
    wh = np.where( prob == prob.max() )[0][0]
    p = samples[ wh, : ]
    TO3wh = TO3[wh]
    TO2wh = TO2[wh]
    #nruns = 100
    #line1 = 'Values of the %i dimensions that best fit the data in %i runs with %i walkers, are the following:' % (ndim, nruns, nwalkers)
    minch0 = 'Values that best fit the data from Chi^2 minimization are the following:'
    minch1 = 'He = %0.2f   O = %0.2f   C/O = %0.2f   N/O = %0.2f   Ne/O = %0.2f   S/O = %0.2f  Te_O3 = %d   Te_O2 = %d' % (
                                                               p[0], p[1], p[2], p[3], p[4], p[5], int(TO3wh), int(TO2wh) )
    minch2 = 'He = %0.2f   O = %0.2f     C = %0.2f      N = %0.2f      Ne = %0.2f      S = %0.2f' % (p[0], p[1], p[2]+true_abunds[1], 
                                                               p[3]+true_abunds[1], p[4]+true_abunds[1], p[5]+true_abunds[1])
    # Calculate the uncertainties based on the 16th, 50th and 84th percentiles
    #samples[:, ndim-1] = np.exp(samples[:, ndim-1])
    percentiles0 = 'MCMC values and uncertainties according to 16th, 50th, and 84th percentiles:'
    p_mcmc1 = map(lambda v: (v), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    p_mcmc2 = map(lambda v: (v[1], np.abs(v[2]-v[1]), np.abs(v[1]-v[0])), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    # find temperatures of these sets:
    TO3_perc = get_Te_percentiles(p_mcmc1, TO3)
    TO2_perc = get_Te_percentiles(p_mcmc1, TO2)
    perc = [p_mcmc1[0], p_mcmc1[1], p_mcmc1[2], p_mcmc1[3], p_mcmc1[4], p_mcmc1[5], TO3_perc, TO2_perc]
    #print p_mcmc1[5][0]#perc[6], perc[7]
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
    print minch0
    print minch1
    print minch2
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
    
    true_abstemps = [trueabs0, trueabs1, trueabs2, temps0, temps1]
    avgabunds = [avgabunds0, avgabunds1, avgabunds2]
    subavgabunds = [subavgabunds0, subavgabunds1, subavgabunds2]
    subchi = [subchi0, subchi1, subchi2]
    minchi_percentiles = [minch0, minch1, minch2, p, percentiles0, percentiles1, percentiles2, percentiles3, p_mcmc2]
    stat_info = [true_abstemps, avgabunds, subavgabunds, subchi, minchi_percentiles]
    return stat_info


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

# Set the the variables
object_name = args.object_name
not_final = args.not_final
use_subset = args.use_subset
mk_plots = args.mk_plots
contours = args.contours
wt_file = args.wt_file

# Obtain the benchmark or "true" abundances
true_abunds, true_temps = get_true_abunds(object_name)

# Load the chain and create the samples
chain_file = os.path.abspath('mcmc_'+object_name+'_chain.dat')
if not_final:
    chain_file = os.path.abspath('mcmc_'+object_name+'_chain0.dat')   # this is equivalent to cSFRg40hZsolarUV
He, O, CO, NO, NeO, SO, TO3, TO2, prob = np.loadtxt(chain_file, skiprows=17, unpack=True)
samples = np.array([]).reshape(0, len(true_abunds))
for he, o, co, no, neo, so in zip(He, O, CO, NO, NeO, SO):
    #carbon = co + o
    #if carbon >=7.7 and carbon <= 8.6:
    theta = [he, o, co, no, neo, so]
        #print o, co, carbon, no+o
    samples = np.vstack((samples, theta))
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
print 'ok, got the samples! shape of array is: ', np.shape(samples)
print 'The SUBsamples* shape is: ', np.shape(subsamples)
print '* SUBsamples is the set of models without repeated abundance and temperature sets.'
    
# PLOTS
if mk_plots or contours:
    # plot of the temperatures without the benchmark values
    fig1 = plt.figure(1, figsize=(12, 10))
    plt.title('C/O vs O3_Temps')
    xlab = 'Te_[OIII]'
    ylab = 'log (C/O)'
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    if mk_plots:
        plt.plot(TO3, CO, 'k.')
    elif contours:
        x, y, z = get_zarr(TO3, CO)
        plt.plot(x, y, 'k.')
        fig1 = triangle.corner(z, labels=[xlab, ylab])
    COTemp = object_name+'_tempsVsCO.jpg'
    fig1.savefig(os.path.abspath(COTemp))
    #plt.show()
    
    fig2 = plt.figure(1, figsize=(12, 10))
    plt.title('C/O vs O/H')
    xlab = '12 + log (O/H)'
    ylab ='log (C/O)'
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    if mk_plots:
        plt.plot(O, CO, 'k.')
    elif contours:
        x, y, z = get_zarr(O, CO)
        plt.plot(x, y, 'k.')
        fig2 = triangle.corner(z, labels=[xlab, ylab])
    COHO = object_name+'_COvsOH.jpg'
    fig2.savefig(os.path.abspath(COHO))
    #plt.show()

    fig3 = plt.figure(1, figsize=(12, 10))
    plt.title('C/H vs O/H')
    xlab= '12 + log (O/H)'
    ylab ='12 + log (C/H)'
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    ymin = 6.0
    ymax = 8.8
    plt.ylim(ymin, ymax)
    if mk_plots:
        plt.plot(O, CO+O, 'k.')
    elif contours:
        x, y, z = get_zarr(O, CO+O)
        plt.plot(x, y, 'k.')
        fig3 = triangle.corner(z, labels=[xlab, ylab], extents=[(7.5, 9.0), (ymin, ymax)])
    CHHO = object_name+'_CHvsOH.jpg'
    fig3.savefig(os.path.abspath(CHHO))
    #plt.show()

    fig4 = plt.figure(1, figsize=(12, 10))
    plt.title('C/H vs O3_Temps')
    xlab = 'Te_[OIII]'
    ylab = '12 + log (C/H)'
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.ylim(ymin, ymax)
    if mk_plots:
        plt.plot(TO3, CO+O, 'k.')
    elif contours:
        x, y, z = get_zarr(TO3, CO+O)
        plt.plot(x, y, 'k.')
        fig4 = triangle.corner(z, labels=[xlab, ylab], extents=[(8000.0, 18000.0), (ymin, ymax)])
    CTemp = object_name+'_tempsVsCH.jpg'
    fig4.savefig(os.path.abspath(CTemp))
    #plt.show()

    nwalkers = 100
    nruns = 100
    # plot of abundance ratios including the benchmark abundances
    fig = triangle.corner(samples, 
                          labels=["$He$", "$O$", "$C/O$", "$N/O$", "$Ne/O$", "$S/O$"], 
                          truths=[true_abunds[0], true_abunds[1], true_abunds[2], true_abunds[3],
                                  true_abunds[4], true_abunds[5]] ,
                          # Limits:    He           O            C/O          N/O           Ne/O         S/O
                          extents=[(9.5, 11.7), (7.5, 8.6), (-1.6, 1.6), (-1.7, -0.4), (-1.0, 0.01), (-2.3, -1.3)]
                          )
    pltbench = 'mcmc_'+object_name+"_ratios_"+repr(nwalkers)+"w"+repr(nruns)+"r"+"_initb.jpg"
    fig.savefig(os.path.abspath(pltbench))
        
    # plot of the ratios without the benchmark abundances
    fig = triangle.corner(samples, labels=["$He$", "$O$", "$C/O$", "$N/O$", "$Ne/O$", "$S/O$"])
    pltwithoutbench = 'mcmc_'+object_name+"_ratios2_"+repr(nwalkers)+"w"+repr(nruns)+"r"+"_initb.jpg"
    fig.savefig(os.path.abspath(pltwithoutbench))    

if use_subset:
    He = subsamples[:, 0]
    O = subsamples[:, 1]
    CO = subsamples[:, 2]
    NO = subsamples[:, 3]
    NeO = subsamples[:, 4]
    SO = subsamples[:, 5]
    stat_info = get_bestmodel(He, O, CO, NO, NeO, SO, subTO3, subTO2, subprob, subsamples, 
                                                      true_abunds, true_temps)
else:
    stat_info = get_bestmodel(He, O, CO, NO, NeO, SO, TO3, TO2, prob, samples, true_abunds, true_temps)

if wt_file:
    # write statistics in final chain file...
    final_chain_file = os.path.abspath('mcmc_'+object_name+"_chain.txt")   
    true_abstemps, avgabunds, subavgabunds, subchi, minchi_percentiles = stat_info 
    trueabs0, trueabs1, trueabs2, temps0, temps1 = true_abstemps
    avgabunds0, avgabunds1, avgabunds2 = avgabunds
    subavgabunds0, subavgabunds1, subavgabunds2 = subavgabunds
    subchi0, subchi1, subchi2 = subchi
    minch0, minch1, minch2, p, percentiles0, percentiles1, percentiles2, percentiles3, p_mcmc2 = minchi_percentiles
    stats = [trueabs0, trueabs1, trueabs2, temps0, temps1, avgabunds0, avgabunds1, avgabunds2,
             subavgabunds0, subavgabunds1, subavgabunds2, subchi0, subchi1, subchi2,
             minch0, minch1, minch2, p, percentiles0, percentiles1, p_mcmc2]
    He, O, CO, NO, NeO, SO, TO3, TO2, Chi2 = np.loadtxt(chain_file, skiprows=16, unpack=True)
    f1 = open(final_chain_file, "w")
    for item in stats:
        print >> f1, item
    print >> f1, "\n{:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<15} {:<15} {:<20}".format('He', 'O', 'C/O', 'N/O', 'Ne/O', 'S/O', 'TeO3', 'TeO2', 'ln (-0.5*Chi2)')
    for he, o, co, no, neo, so, to3, to2, chi2 in zip(He, O, CO, NO, NeO, SO, TO3, TO2, Chi2):
        print >> f1, "{:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<15} {:<15} {:<20.3f}".format(he, o, co, no, neo, so, to3, to2, chi2)
    f1.close()
    print '\nPlots made: '
    print pltbench+'\n', pltwithoutbench+'\n', COTemp+'\n', COHO+'\n', CHHO+'\n'
    
    
print '\nCode finished! \n'
