from __future__ import division
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import string
import emcee
import corner
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

Run from terminal as:  >> python mktriplot.py sbs1415 -n -c -CHb
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
parser.add_argument("-full",
                    action='store_true',
                    dest="FULL_chain_file",
                    default=False,
                    help='Use -full to read FULL_chain_file.')
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
parser.add_argument("-CHb",
                    action='store_true',
                    dest="C_Hbeta",
                    default=False,
                    help='Used CHbeta [True] or E(B-V) [False].')
args = parser.parse_args()

# name of the object
#                  0            1           2          3         4        5         6         7         8
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                 9           10         11         12         13       14        15         16         17

#### FUNCTIONS

def get_true_abunds(object_name, C_Hbeta):
    object_number = objects_list.index(object_name)
    print object_name, 'is object number ', object_number, ' in the sample list.'
    #fname = os.path.abspath('../../results/Final_meas_abunds.txt')
    Final_meas_abunds = 'carbon/results/Final_meas_abunds.txt'
    if C_Hbeta:
        Final_meas_abunds = 'carbon/results/Final_meas_abunds_CHbeta.txt'
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
    #print '  and benchmark temperatures: \n', true_temps
    #raw_input()
    return true_abunds, true_temps

def get_Te_percentiles(percentiles_arr, Te_arr, samples2use):
    Te_pencentiles = []
    theta_p1 = np.array([[percentiles_arr[0][0], percentiles_arr[1][0], percentiles_arr[2][0], percentiles_arr[3][0], percentiles_arr[4][0], percentiles_arr[5][0]]])
    theta_p2 = np.array([[percentiles_arr[0][1], percentiles_arr[1][1], percentiles_arr[2][1], percentiles_arr[3][1], percentiles_arr[4][1], percentiles_arr[5][1]]])
    theta_p3 = np.array([[percentiles_arr[0][2], percentiles_arr[1][2], percentiles_arr[2][2], percentiles_arr[3][2], percentiles_arr[4][2], percentiles_arr[5][2]]])
    wp1 = np.where(samples2use == theta_p1)[0][0]
    wp2 = np.where(samples2use == theta_p2)[0][0]
    wp3 = np.where(samples2use == theta_p3)[0][0]
    Tep1 = Te_arr[wp1]
    Tep2 = Te_arr[wp2]
    Tep3 = Te_arr[wp3]
    Te_pencentiles.append(Tep1)
    Te_pencentiles.append(Tep2)
    Te_pencentiles.append(Tep3)
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
    
    medHe = np.median(He) 
    medO = np.median(O)
    medCO = np.median(CO+O)
    medNO = np.median(NO+O)
    medNeO = np.median(NeO+O)
    medSO = np.median(SO+O)
    medTO3 = np.median(TO3)
    medTO2 = np.median(TO2)
    medabunds0 = 'Median abundances: '
    medabunds1 = 'He = %0.2f   O = %0.2f   C/O = %0.2f   N/O = %0.2f   Ne/O = %0.2f   S/O = %0.2f   Te_O3 = %d   Te_O2 = %d' % (
                                   medHe, medO, medCO-medO, medNO-medO, medNeO-medO, medSO-medO, int(medTO3), int(medTO2))
    medabunds2 = 'He = %0.2f   O = %0.2f     C = %0.2f      N = %0.2f      Ne = %0.2f      S = %0.2f ' % (
                                                            medHe, medO, medCO, medNO, medNeO, medSO )
    print '\n'+medabunds0
    print medabunds1
    print medabunds2
    
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
    avg = [avgHe, avgO, avgCO-avgO, avgNO-avgO, avgNeO-avgO, avgSO-avgO, avgTO3, avgTO2]
    subavgabunds0 = 'Average abundances with O temperatures within '+str(nearby)+' from the benchmark values (i.e. averages of '+subset+' models): '
    subavgabunds1 = 'He = %0.2f   O = %0.2f   C/O = %0.2f   N/O = %0.2f   Ne/O = %0.2f   S/O = %0.2f   Te_O3 = %d   Te_O2 = %d' % (
                                   avg[0], avg[1], avg[2], avg[3], avg[4], avg[5], int(avg[6]), int(avg[7]))
    subavgabunds2 = 'He = %0.2f   O = %0.2f     C = %0.2f      N = %0.2f      Ne = %0.2f      S = %0.2f ' % (
                                                                    avgHe, avgO, avgCO, avgNO, avgNeO, avgSO )
    # find temperatures of these sets:
    p_subavgabunds0 = 'MCMC values and uncertainties according to 25th, 50th, and 75th percentiles:'
    percentiles = [25, 50, 75]
    p_subavgabunds1 = map(lambda v: (v), zip(*np.percentile(subsamples, percentiles, axis=0)))
    p_subavgabunds1a = map(lambda v: (v[1], np.abs(v[2]-v[1]), np.abs(v[1]-v[0])), zip(*np.percentile(subsamples, percentiles, axis=0)))
    TO3_p = np.percentile(np.array(TO3_nearby), percentiles)
    TO2_p = np.percentile(np.array(TO2_nearby), percentiles)
    print 'TO3_p = ', TO3_p
    print 'errors with respect to average value: ', int(avgTO3-TO3_p[0]), int(avgTO3-TO3_p[1]), int(TO3_p[2]-avgTO3)
    print 'TO2_p = ', TO2_p
    print 'errors with respect to average value: ', int(avgTO2-TO2_p[0]), int(avgTO2-TO2_p[1]), int(TO2_p[2]-avgTO2)
    TO3_perc = get_Te_percentiles(p_subavgabunds1, np.array(TO3_nearby), subsamples)
    TO2_perc = get_Te_percentiles(p_subavgabunds1, np.array(TO2_nearby), subsamples)
    perc = [p_subavgabunds1[0], p_subavgabunds1[1], p_subavgabunds1[2], p_subavgabunds1[3], p_subavgabunds1[4], p_subavgabunds1[5], TO3_perc, TO2_perc]
    p_subavgabunds2 = '{:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:>8} {:<9} \n {:>17.2f} {:>7.2f} {:>6.2f} {:>6.2f} \n {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:>8} {:<9}'.format(
                     perc[0][0], perc[1][0], perc[2][0], perc[3][0], perc[4][0], perc[5][0], int(perc[6][0]), int(perc[7][0]),
                     perc[2][0]+perc[1][0], perc[3][0]+perc[1][0], perc[4][0]+perc[1][0], perc[5][0]+perc[1][0], 
                     avg[0]-perc[0][0], avg[1]-perc[1][0], avg[2]-perc[2][0], avg[3]-perc[3][0], avg[4]-perc[4][0], avg[5]-perc[5][0], int(avg[6]-perc[6][0]), int(avg[7]-perc[7][0]) )
    p_subavgabunds3 = '{:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:>8} {:<9} \n {:>17.2f} {:>7.2f} {:>6.2f} {:>6.2f}'.format(
                     perc[0][1], perc[1][1], perc[2][1], perc[3][1], perc[4][1], perc[5][1], int(perc[6][1]), int(perc[7][1]),
                     perc[2][1]+perc[1][1], perc[3][1]+perc[1][1], perc[4][1]+perc[1][1], perc[5][1]+perc[1][1] )
    p_subavgabunds4 = '{:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:>8} {:<9} \n {:>17.2f} {:>7.2f} {:>6.2f} {:>6.2f} \n {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:<6.2f} {:>8} {:<9}'.format(
                     perc[0][2], perc[1][2], perc[2][2], perc[3][2], perc[4][2], perc[5][2], int(perc[6][2]), int(perc[7][2]),
                     perc[2][2]+perc[1][2], perc[3][2]+perc[1][2], perc[4][2]+perc[1][2], perc[5][2]+perc[1][2], 
                     perc[0][2]-avg[0], perc[1][2]-avg[1], perc[2][2]-avg[2], perc[3][2]-avg[3], perc[4][2]-avg[4], perc[5][2]-avg[5], int(perc[6][2]-avg[6]), int(perc[7][2]-avg[7]) )
    print '\n'
    print subavgabunds0
    print subavgabunds1 
    print subavgabunds2
    print p_subavgabunds0
    print 'He      O       C      N      Ne      S      Te_O3   Te_O2'
    print '25th percentile:'
    print p_subavgabunds2
    print '50th percentile:'
    print p_subavgabunds3
    print '75th percentile:'
    print p_subavgabunds4
    print '\n DIFFERENCES of percentiles with respect to 50th percentile: '
    print p_subavgabunds1a
    
    # Chi^2 minimization within subsets
    whsub = np.where( prob_nearby == max(prob_nearby) )[0][0]
    psub = subsamples[ whsub, : ]
    TO3whsub = TO3_nearby[whsub]
    TO2whsub = TO2_nearby[whsub]
    subchi0 = 'Values of subset that best fit the data  from Chi^2 minimization within then '+str(nearby)+'K temperature tolerance are the following:'
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
    minch0 = 'Values whole set that best fit the data from Chi^2 minimization are the following:'
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
    TO3_perc = get_Te_percentiles(p_mcmc1, TO3, samples)
    TO2_perc = get_Te_percentiles(p_mcmc1, TO2, samples)
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

def fit_line(arrx, arry):
    '''
    This function fits a line to the given array.
    arrx, arry = numpy arrays 
    RETURNS:
    - constants of the line
    - fitted line
    '''
    order = 1
    coefficients = np.polyfit(arrx, arry, order)
    polynomial = np.poly1d(coefficients)
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

# define the priors
def lnprior(theta):
    m, b = theta
    mmod = lntophat(m, -5.5, 0.0) 
    bmod = lntophat(b,-0.7, 10.0) 
    if mmod != -np.inf and bmod != -np.inf:
        return mmod + bmod
    return -np.inf

# then the probability function will be
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)

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
        return np.log(1/(bb-aa))
    else:
        return -np.inf


#### CODE 

# Set the the variables
object_name = args.object_name
not_final = args.not_final
FULL_chain_file = args.FULL_chain_file
use_subset = args.use_subset
mk_plots = args.mk_plots
contours = args.contours
wt_file = args.wt_file
C_Hbeta = args.C_Hbeta

# Format of images to save
img_type = '.jpg'

# Obtain the benchmark or "true" abundances
true_abunds, true_temps = get_true_abunds(object_name, C_Hbeta)

# Load the chain and create the samples
chain_file = os.path.abspath('mcmc_'+object_name+'_chain.dat')
if not_final:
    chain_file = os.path.abspath('mcmc_'+object_name+'_chain0.dat')   # this is equivalent to cSFRg40hZsolarUV
if FULL_chain_file:
    chain_file = os.path.abspath('mcmc_'+object_name+'_FULL_chain0.dat')   # this is equivalent to cSFRg40hZsolarUV
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
quantiles = [0.25, 0.5, 0.75] #[0.16, 0.5, 0.84]
if mk_plots or contours:
    # plot of the temperatures without the benchmark values
    fig1 = plt.figure(1, figsize=(12, 10))
    plt.title('C/O vs O3_Temps')
    xlab = 'Te_[OIII]'
    ylab = 'log (C/O)'
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    x = TO3
    y = CO
    if mk_plots:
        plt.plot(x, y, 'k.')
    elif contours:
        x, y, z = get_zarr(x, y)
        plt.plot(x, y, 'k.')
        print 'These are the quantiles for Te_OIII vs C/O'
        fig1 = corner.corner(z, labels=[xlab, ylab], quantiles=quantiles)
    COTemp = object_name+'_tempsVsCO'+img_type
    fig1.savefig(os.path.abspath(COTemp))
    #plt.show()
    
    fig2 = plt.figure(1, figsize=(12, 10))
    plt.title('C/O vs O/H')
    xlab = '12 + log (O/H)'
    ylab ='log (C/O)'
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    x = O
    y = CO
    plt.xlim(7.5, 8.9)
    if mk_plots:
        plt.plot(x, y, 'k.')
    elif contours:
        x, y, z = get_zarr(x, y)
        plt.plot(x, y, 'k.')
        fig2 = corner.corner(z, labels=[xlab, ylab])
    COHO = object_name+'_COvsOH'+img_type
    fig2.savefig(os.path.abspath(COHO))
    #plt.show()

    fig3 = plt.figure(1, figsize=(12, 10))
    # Polynomial fit
    x = O
    y = CO+O
    coeffs, line_fit = fit_line(x, y)
    print 'Linear fit: y = mx + b'
    print line_fit
    print 'Coefficients: '
    print coeffs
    plt.title('C/H vs O/H')
    xlab= '12 + log (O/H)'
    ylab ='12 + log (C/H)'
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    ymin = 6.0
    ymax = 8.8
    plt.ylim(ymin, ymax)
    plt.xlim(7.5, 8.9)
    if mk_plots:
        plt.plot(x, y, 'k.')
    elif contours:
        x, y, z = get_zarr(x, y)
        plt.plot(x, y, 'k.')
        plt.plot(x, line_fit, 'b:')
        print 'These are the quantiles for C/H to O/H'
        fig3 = corner.corner(z, labels=[xlab, ylab], #quantiles=quantiles, 
                               range=[(7.5, 9.0), (ymin, ymax)])
    CHHO = object_name+'_CHvsOH'+img_type
    fig3.savefig(os.path.abspath(CHHO))
    #plt.show()

    fig4 = plt.figure(1, figsize=(12, 10))
    plt.title('C/H vs O3_Temps')
    xlab = 'Te_[OIII]'
    ylab = '12 + log (C/H)'
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.ylim(ymin, ymax)
    x = TO3
    y = CO+O
    if mk_plots:
        plt.plot(x, y, 'k.')
    elif contours:
        x, y, z = get_zarr(x, y)
        plt.plot(x, y, 'k.')
        fig4 = corner.corner(z, labels=[xlab, ylab], range=[(8000.0, 18000.0), (ymin, ymax)])
    CTemp = object_name+'_tempsVsCH'+img_type
    fig4.savefig(os.path.abspath(CTemp))
    #plt.show()

    fig5 = plt.figure(1, figsize=(12, 10))
    plt.title('N/C vs C/H')
    xlab = '12 + log (C/H)'
    ylab = 'log (N/C)'
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.ylim(6.0, 8.8)
    x = CO+O
    y = NO-CO
    if mk_plots:
        plt.plot(x, y, 'k.')
    elif contours:
        x, y, z = get_zarr(x, y)
        plt.plot(x, y, 'k.')
        fig5 = corner.corner(z, labels=[xlab, ylab], #quantiles=quantiles, 
                               range=[(6.0, 8.8), (-2.6, 0.7)])
    NC = object_name+'_NCvsCH'+img_type
    fig5.savefig(os.path.abspath(NC))
    #plt.show()

    fig = plt.figure(1, figsize=(12, 10))
    plt.title('N/O vs C/N')
    xlab= 'log (N/O)'
    ylab ='log (C/N)'
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    ymin = -0.60
    ymax = 2.60
    plt.ylim(ymin, ymax)
    CN = CO - NO
    x, y, z = get_zarr(NO, CN)
    fig = corner.corner(z, labels=[xlab, ylab], range=[(-1.7, -0.6), (ymin, ymax)])
    fn = object_name+'_NOvsCN.jpg'
    fig.savefig(os.path.abspath(fn))
    print 'Figure ', fn, 'saved!' 
    # Adjust a linear fit to the plot
    coeffs, line_fit = fit_line(x, y)
    print 'Coefficients of initial guess to the plot of: ', fn
    m = coeffs[0]
    b = coeffs[1]
    print 'm = %0.3f     b = %0.3f' % (m, b)
    # Initialize the chain with a Gaussian ball around the maximum likelihood result, for which we use optimize
    ndim, nwalkers, nruns = 2, 100, 100
    yerr = []
    for yi in y:
        ye = yi * 0.1
        yerr.append(ye)
    yerr = np.array(yerr)
    #randadd2point = lambda x: x+np.random.rand(1)
    #p0 = np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))
    p0 = [[np.random.uniform(-9.5, 0.01), np.random.uniform(-2.0, 0.0)] for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[x, y, yerr])
    pos, prob, state = sampler.run_mcmc(p0, nruns)   # this allows for the first 50 steps to be "burn-in" type
    # best model
    wh = np.where( prob == prob.max() )[0][0]
    p = pos[ wh, : ]
    plt.plot( x, line_eq( p, x ), "r", lw=5, alpha=0.4 )
    line1 ='values of the %i dimensions that best fit the data according to Chi^2 in %i runs are the following:' % (ndim, nruns)
    line2 = 'm = %0.3f   b = %0.3f' % (p[0], p[1])
    print line1
    print line2
    allm, allb = [], []
    subsample = np.array([]).reshape(0, 2)
    #print 'pos :', pos
    for theta in pos:
        if theta[0] < 0.0:
            allm.append(theta[0])
            allb.append(theta[1])
            subsample = np.vstack((subsample, theta))
    avgm = sum(allm)/len(allm)
    avgb = sum(allb)/len(allb)
    print 'Average values:  m = %0.3f   b= %0.3f' % (avgm, avgb)
    #linesamples = sampler.chain[:, nruns*0.2, :].reshape((-1, ndim))
    #fig = corner.corner(linesamples, labels=["$m$", "$b$"], truths=[m, b])
    #fig.show()

    nwalkers = 100
    nruns = 100
    # plot of abundance ratios including the benchmark abundances
    print 'These are the quantiles for triangle plot: '
    labels = ["$He$", "$O$", "$C/O$", "$N/O$", "$Ne/O$", "$S/O$"]
    fig = corner.corner(samples, 
                          labels=labels, 
                          truths=[true_abunds[0], true_abunds[1], true_abunds[2], true_abunds[3],
                                  true_abunds[4], true_abunds[5]],
                          #quantiles=quantiles,
                          # Limits:    He           O            C/O          N/O           Ne/O         S/O
                          range=[(9.5, 11.7), (7.55, 8.6), (-1.6, 1.6), (-1.7, -0.4), (-1.0, 0.01), (-2.3, -1.4)]
                          )
    pltbench = 'mcmc_'+object_name+"_ratios_"+repr(nwalkers)+"w"+repr(nruns)+"r"+"_initb"+img_type
    if FULL_chain_file:
        pltbench = 'mcmc_'+object_name+"_ratios_"+repr(nwalkers)+"w"+repr(nruns)+"r"+"_FULL_initb"+img_type
    fig.savefig(os.path.abspath(pltbench))
        
    # plot of the ratios without the benchmark abundances
    fig = corner.corner(samples, labels=labels)
    pltwithoutbench = 'mcmc_'+object_name+"_ratios2_"+repr(nwalkers)+"w"+repr(nruns)+"r"+"_initb"+img_type
    if FULL_chain_file:
        pltwithoutbench = 'mcmc_'+object_name+"_ratios2_"+repr(nwalkers)+"w"+repr(nruns)+"r"+"_FULL_initb"+img_type
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
