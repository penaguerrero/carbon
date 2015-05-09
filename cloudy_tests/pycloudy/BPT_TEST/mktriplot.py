from __future__ import division
import os
import numpy as np
import matplotlib
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
            - BPT diagram 
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
parser.add_argument("-g0l",
                    action='store_true',
                    dest="g0l",
                    default=False,
                    help='Used the Geneva tracks with v00, Solar UV radiation, and low metallicity.')
parser.add_argument("-g0h",
                    action='store_true',
                    dest="g0h",
                    default=False,
                    help='Used the Geneva tracks with v00, Solar UV radiation, and high metallicity.')
parser.add_argument("-g4l",
                    action='store_true',
                    dest="g4l",
                    default=False,
                    help='Used the Geneva tracks with v40, Solar UV radiation, and low metallicity.')
parser.add_argument("-g4h",
                    action='store_true',
                    dest="g4h",
                    default=False,
                    help='Used the Geneva tracks with v40, Solar UV radiation, and high metallicity.')
parser.add_argument("-g4ha",
                    action='store_true',
                    dest="g4ha",
                    default=False,
                    help='Used the Geneva tracks with v40, Solar UV radiation, and high metallicity, but when the total abundances instead of the ratios.')
args = parser.parse_args()

# name of the object
#                  0            1           2          3         4        5         6         7         8
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                 9           10         11         12         13       14        15         16         17

full_names_list = ['IIIZw 107', 'IRAS 08339+6517', 'MRK 1087', 'MRK 1199', 'MRK 5', 'MRK 960', 'NGC 1741', 'POX 4', 'SBS 0218+003',
               'SBS 0948+532', 'SBS 0926+606', 'SBS 1054+365', 'SBS 1319+579', 'TOL 1457-262', 'TOL 9', 'ARP 252', 'IRAS 08208+2816',
               'SBS 1415+437']

#### FUNCTIONS

def get_true_abunds(object_name):
    object_number = objects_list.index(object_name)
    print object_name, 'is object number ', object_number, ' in the sample list.'
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

def get_Te_percentiles(percentiles_arr, Te_arr, samples2use):
    Te_pencentiles = []
    theta_p16 = np.array([[percentiles_arr[0][0], percentiles_arr[1][0], percentiles_arr[2][0], percentiles_arr[3][0], percentiles_arr[4][0], percentiles_arr[5][0]]])
    theta_p50 = np.array([[percentiles_arr[0][1], percentiles_arr[1][1], percentiles_arr[2][1], percentiles_arr[3][1], percentiles_arr[4][1], percentiles_arr[5][1]]])
    theta_p80 = np.array([[percentiles_arr[0][2], percentiles_arr[1][2], percentiles_arr[2][2], percentiles_arr[3][2], percentiles_arr[4][2], percentiles_arr[5][2]]])
    wp16 = np.where(samples2use == theta_p16)[0][0]
    wp50 = np.where(samples2use == theta_p50)[0][0]
    wp80 = np.where(samples2use == theta_p80)[0][0]
    Tep16 = Te_arr[wp16]
    Tep50 = Te_arr[wp50]
    Tep80 = Te_arr[wp80]
    Te_pencentiles.append(Tep16)
    Te_pencentiles.append(Tep50)
    Te_pencentiles.append(Tep80)
    return Te_pencentiles

def get_bestmodel(He, O, CO, NO, NeO, SO, TO3, TO2, prob, samples, true_abunds, true_temps, bpt_data):
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

    # BPT data for all accepted models
    linesample = mklinesample(bpt_data)
    I4861, I5007, I6548, I6563, I6584 = bpt_data

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
    avg4861 = np.mean(I4861)
    avg5007 = np.mean(I5007)
    avg6548 = np.mean(I6548)
    avg6563 = np.mean(I6563)
    avg6584 = np.mean(I6584)
    avgIs0 = '* Average line intensities for: \n4861      5007      6548     6563     6584'
    avgIs1 = '%0.2f    %0.2f    %0.2f    %0.2f    %0.2f' % (avg4861, avg5007, avg6548, avg6563, avg6584)
    
    print '\n'+avgabunds0
    print avgabunds1
    print avgabunds2
    print avgIs0
    print avgIs1
    
    # now find averages of those models with temperatures within +-2500 of the benchmark values 
    # AND of the corresponding line intensities for the BPT diagrams
    subsamples = np.array([]).reshape(0, len(true_abunds))
    He_nearby, O_nearby, CO_nearby, NO_nearby, NeO_nearby, SO_nearby, TO3_nearby, TO2_nearby, prob_nearby = [],[],[],[],[],[],[],[],[]
    I4861_nearby, I5007_nearby, I6548_nearby, I6563_nearby, I6584_nearby = [], [], [], [], []
    nearby = 2500.0
    idx = 0
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
            I4861_nearby.append(I4861[idx])
            I5007_nearby.append(I5007[idx])
            I6548_nearby.append(I6548[idx])
            I6563_nearby.append(I6563[idx])
            I6584_nearby.append(I6584[idx])
        idx = idx + 1
    avgHe = sum(He_nearby)/len(He_nearby) 
    avgO = sum(O_nearby)/len(He_nearby) 
    avgCO = sum(CO_nearby+O_nearby)/len(He_nearby) 
    avgNO = sum(NO_nearby+O_nearby)/len(He_nearby) 
    avgNeO = sum(NeO_nearby+O_nearby)/len(He_nearby) 
    avgSO = sum(SO_nearby+O_nearby)/len(He_nearby) 
    avgTO3 = sum(TO3_nearby)/len(He_nearby) 
    avgTO2 = sum(TO2_nearby)/len(He_nearby) 
    avg = [avgHe, avgO, avgCO-avgO, avgNO-avgO, avgNeO-avgO, avgSO-avgO, avgTO3, avgTO2]
    subset = str(len(O_nearby))
    subavgabunds0 = 'Average abundances with O temperatures within '+str(nearby)+' from the benchmark values (i.e. averages of '+subset+' models): '
    subavgabunds1 = 'He = %0.2f   O = %0.2f   C/O = %0.2f   N/O = %0.2f   Ne/O = %0.2f   S/O = %0.2f   Te_O3 = %d   Te_O2 = %d' % (
                                   avgHe, avgO, avgCO-avgO, avgNO-avgO, avgNeO-avgO, avgSO-avgO, int(avgTO3), int(avgTO2))
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
    
    # BPT lines
    avg4861 = sum(I4861_nearby)/len(I4861_nearby)
    avg5007 = sum(I5007_nearby)/len(I5007_nearby)
    avg6548 = sum(I6548_nearby)/len(I6548_nearby)
    avg6563 = sum(I6563_nearby)/len(I6563_nearby)
    avg6584 = sum(I6584_nearby)/len(I6584_nearby)
    avgIs = [avg4861, avg5007, avg6548, avg6563, avg6584]
    bpt_data_nearby = [I4861_nearby, I5007_nearby, I6548_nearby, I6563_nearby, I6584_nearby]
    linesubsample = mklinesample(bpt_data_nearby)
    percs4861 = np.percentile(np.array(I4861_nearby), percentiles) 
    percs5007 = np.percentile(np.array(I5007_nearby), percentiles) 
    percs6548 = np.percentile(np.array(I6548_nearby), percentiles) 
    percs6563 = np.percentile(np.array(I6563_nearby), percentiles) 
    percs6584 = np.percentile(np.array(I6584_nearby), percentiles) 
    print '* Average line intensities for: \n4861      5007      6548     6563     6584'
    print '%0.2f    %0.2f    %0.2f    %0.2f    %0.2f' % (avg4861, avg5007, avg6548, avg6563, avg6584)
    print '4861 25th, 50th, and 75th percentiles:'
    print '%0.2f  %0.2f  %0.2f' % (percs4861[0],percs4861[1],percs4861[2]), '  and with respect to average: %0.2f  %0.2f %0.2f' % (
                                         np.abs(percs4861[0]-avgIs[0]), np.abs(percs4861[1]-avgIs[0]), np.abs(percs4861[2]-avgIs[0]))
    print '5007  25th, 50th, and 75th percentiles:'
    print '%0.2f  %0.2f  %0.2f' % (percs5007[0], percs5007[1], percs5007[2]), '  and with respect to average: %0.2f  %0.2f %0.2f' % (
                                         np.abs(percs5007[0]-avgIs[1]), np.abs(percs5007[1]-avgIs[1]), np.abs(percs5007[2]-avgIs[1]))
    print '6548  25th, 50th, and 75th percentiles:'
    print '%0.2f  %0.2f  %0.2f' % (percs6548[0], percs6548[1], percs6548[2]), '  and with respect to average: %0.2f  %0.2f %0.2f' % (
                                         np.abs(percs6548[0]-avgIs[2]), np.abs(percs6548[1]-avgIs[2]), np.abs(percs6548[2]-avgIs[2]))
    print '6563  25th, 50th, and 75th percentiles:'
    print '%0.2f  %0.2f  %0.2f' % (percs6563[0], percs6563[1], percs6563[2]), '  and with respect to average: %0.2f  %0.2f %0.2f' % (
                                         np.abs(percs6563[0]-avgIs[3]), np.abs(percs6563[1]-avgIs[3]), np.abs(percs6563[2]-avgIs[3]))
    print '6584  25th, 50th, and 75th percentiles:'
    print '%0.2f  %0.2f  %0.2f' % (percs6584[0], percs6584[1], percs6584[2]), '  and with respect to average: %0.2f  %0.2f %0.2f' % (
                                         np.abs(percs6584[0]-avgIs[4]), np.abs(percs6584[1]-avgIs[4]), np.abs(percs6584[2]-avgIs[4]))
    
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
    percentiles0 = 'MCMC values and uncertainties according to 16th, 50th, and 84th percentiles FOR ALL models:'
    p_mcmc1 = map(lambda v: (v), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    p_mcmc2 = map(lambda v: (v[1], np.abs(v[2]-v[1]), np.abs(v[1]-v[0])), zip(*np.percentile(samples, [25, 50, 75], axis=0)))
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
    print '25th percentile:'
    print percentiles1
    print '50th percentile:'
    print percentiles2
    print '75th percentile:'
    print percentiles3
    print '\n DIFFERENCES of percentiles with respect to 50th percentile for ALL the sample: '
    print p_mcmc2
    
    true_abstemps = [trueabs0, trueabs1, trueabs2, temps0, temps1]
    avgabunds = [avgabunds0, avgabunds1, avgabunds2]
    subavgabunds = [subavgabunds0, subavgabunds1, subavgabunds2]
    subchi = [subchi0, subchi1, subchi2]
    minchi_percentiles = [minch0, minch1, minch2, p, percentiles0, percentiles1, percentiles2, percentiles3, p_mcmc2]
    stat_info = [true_abstemps, avgabunds, subavgabunds, subchi, minchi_percentiles]
    return stat_info, bpt_data_nearby


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

def get_bptdata(object_name, gen_tracks):
    '''
    This function reads the line intensities from the BPT file. 
    '''
    bpt_file = 'mcmc_'+object_name+gen_tracks+'_BPT.dat'
    bpt_data = np.loadtxt(bpt_file, skiprows=2, unpack=True)
    # bpt_data = I4861, I5007, I6548, I6563, I6584
    return bpt_data

def mklinesample(bpt_data):
    I4861, I5007, I6548, I6563, I6584 = bpt_data
    samples = np.array([]).reshape(0, len(bpt_data))
    for i4861, i5007, i6548, i6563, i6584 in zip(I4861, I5007, I6548, I6563, I6584):
        theta = [i4861, i5007, i6548, i6563, i6584]
        samples = np.vstack((samples, theta))
    return samples
    


#### CODE 

# Set the the variables
object_name = args.object_name
not_final = args.not_final
FULL_chain_file = args.FULL_chain_file
use_subset = args.use_subset
mk_plots = args.mk_plots
contours = args.contours
wt_file = args.wt_file
g0l = args.g0l
g0h = args.g0h
g4l = args.g4l
g4h = args.g4h
g4ha = args.g4ha

# Format of images to save
img_type = '.jpg'

# Obtain the benchmark or "true" abundances
true_abunds, true_temps = get_true_abunds(object_name)
if g4ha:
    true_abunds[2] = true_abunds[2] + true_abunds[1]
    true_abunds[3] = true_abunds[3] + true_abunds[1]
    true_abunds[4] = true_abunds[4] + true_abunds[1]
    true_abunds[5] = true_abunds[5] + true_abunds[1]

# Load the chain and create the samples
chain_file = os.path.abspath('mcmc_'+object_name+'_chain.dat')
gen_tracks =  ''
if not_final:
    chain_file = os.path.abspath('mcmc_'+object_name+'_chain0.dat')   # this is equivalent to cSFRg40hZsolarUV
if FULL_chain_file:
    chain_file = os.path.abspath('mcmc_'+object_name+'FULL_chain0.dat')   # this is equivalent to cSFRg40hZsolarUV
if g0h:
    chain_file = os.path.abspath('mcmc_'+object_name+'_cSFRg00hZsolarUV_chain0.dat')
    gen_tracks = '_cSFRg00hZsolarUV'
if g0h and FULL_chain_file:
    chain_file = os.path.abspath('mcmc_'+object_name+'_cSFRg00hZsolarUV_FULL_chain0.dat')
    gen_tracks = 'cSFRg00hZsolarUV_FULL'
if g0l:
    chain_file = os.path.abspath('mcmc_'+object_name+'_cSFRg00lZsolarUV_chain0.dat')
    gen_tracks = '_cSFRg00lZsolarUV'
if g0l and FULL_chain_file:
    chain_file = os.path.abspath('mcmc_'+object_name+'_cSFRg00lZsolarUV_FULL_chain0.dat')
    gen_tracks = '_cSFRg00lZsolarUV_FULL'
if g4h:
    chain_file = os.path.abspath('mcmc_'+object_name+'_cSFRg40hZsolarUV_chain0.dat')   # equivalent to runs of all models
    gen_tracks = '_cSFRg40hZsolarUV'
if g4h and FULL_chain_file:
    chain_file = os.path.abspath('mcmc_'+object_name+'_cSFRg40hZsolarUV_FULL_chain0.dat')   # equivalent to runs of all models
    gen_tracks = '_cSFRg40hZsolarUV_FULL'
if g4ha:
    chain_file = os.path.abspath('mcmc_'+object_name+'_cSFRg40_hZ_sUV_abundsXH_chain0.dat')   # equivalent to runs of all models
    gen_tracks = '_cSFRg40_hZ_sUV_abundsXH'
if g4ha and FULL_chain_file:
    chain_file = os.path.abspath('mcmc_'+object_name+'_cSFRg40_hZ_sUV_abundsXH_FULL_chain0.dat')   # equivalent to runs of all models
    gen_tracks = '_cSFRg40_hZ_sUV_abundsXH_FULL'
if g4l:
    chain_file = os.path.abspath('mcmc_'+object_name+'_cSFRg40lZsolarUV_chain0.dat')
    gen_tracks = '_cSFRg40lZsolarUV'
if g4l and FULL_chain_file:
    chain_file = os.path.abspath('mcmc_'+object_name+'_cSFRg40lZsolarUV_FULL_chain0.dat')
    gen_tracks = '_cSFRg40lZsolarUV_FULL'
He, O, CO, NO, NeO, SO, TO3, TO2, prob = np.loadtxt(chain_file, skiprows=17, unpack=True)
samples = np.array([]).reshape(0, len(true_abunds))
for he, o, co, no, neo, so in zip(He, O, CO, NO, NeO, SO):
    theta = [he, o, co, no, neo, so]
    samples = np.vstack((samples, theta))
'''
# select only even indexes: start from 0, use all, skip 2
samples = samples[0::2]
He = He[0::2]
O = O[0::2] 
CO = CO[0::2]
NO = NO[0::2] 
NeO = NeO[0::2] 
SO = SO[0::2] 
TO3 = TO3[0::2]
TO2 = TO2[0::2] 
prob = prob[0::2]
'''
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

# get the BPT data
bpt_data = get_bptdata(object_name, gen_tracks)

    
# PLOTS
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
    if g4ha:
        y = CO-O
    if mk_plots:
        plt.plot(x, y, 'k.')
    elif contours:
        x, y, z = get_zarr(x, y)
        plt.plot(x, y, 'k.')
        print 'These are the quantiles for Te_OIII vs C/O'
        fig1 = triangle.corner(z, labels=[xlab, ylab], quantiles=[0.16, 0.5, 0.84])
    COTemp = object_name+gen_tracks+'_tempsVsCO'+img_type
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
    if g4ha:
        y = CO-O 
    plt.xlim(7.5, 8.9)
    if mk_plots:
        plt.plot(x, y, 'k.')
    elif contours:
        x, y, z = get_zarr(x, y)
        plt.plot(x, y, 'k.')
        fig2 = triangle.corner(z, labels=[xlab, ylab])
    COHO = object_name+gen_tracks+'_COvsOH'+img_type
    fig2.savefig(os.path.abspath(COHO))
    #plt.show()

    fig3 = plt.figure(1, figsize=(12, 10))
    # Polynomial fit
    x = O
    y = CO+O
    if g4ha:
        y = CO
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
        fig3 = triangle.corner(z, labels=[xlab, ylab], quantiles=[0.16, 0.5, 0.84], extents=[(7.5, 9.0), (ymin, ymax)])
    CHHO = object_name+gen_tracks+'_CHvsOH'+img_type
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
    if g4ha:
        y = CO
    if mk_plots:
        plt.plot(x, y, 'k.')
    elif contours:
        x, y, z = get_zarr(x, y)
        plt.plot(x, y, 'k.')
        fig4 = triangle.corner(z, labels=[xlab, ylab], extents=[(8000.0, 18000.0), (ymin, ymax)])
    CTemp = object_name+gen_tracks+'_tempsVsCH'+img_type
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
    if g4ha:
        x = CO
    if mk_plots:
        plt.plot(x, y, 'k.')
    elif contours:
        x, y, z = get_zarr(x, y)
        plt.plot(x, y, 'k.')
        fig5 = triangle.corner(z, labels=[xlab, ylab], quantiles=[0.16, 0.5, 0.84], extents=[(6.0, 8.8), (-2.6, 0.7)])
    NC = object_name+gen_tracks+'_NCvsCH'+img_type
    fig5.savefig(os.path.abspath(NC))
    #plt.show()

    nwalkers = 100
    nruns = 100
    # plot of abundance ratios including the benchmark abundances
    print 'These are the quantiles for triangle plot: '
    labels = ["$He$", "$O$", "$C/O$", "$N/O$", "$Ne/O$", "$S/O$"]
    if g4ha:
        labels = ["$He$", "$O$", "$C$", "$N$", "$Ne$", "$S$"]
    fig = triangle.corner(samples, 
                          labels=labels, 
                          truths=[true_abunds[0], true_abunds[1], true_abunds[2], true_abunds[3],
                                  true_abunds[4], true_abunds[5]],
                          quantiles=[0.16, 0.5, 0.84]
                          # Limits:    He           O            C/O          N/O           Ne/O         S/O
                          #extents=[(9.5, 11.7), (7.5, 8.6), (-1.6, 1.6), (-1.7, -0.4), (-1.0, 0.01), (-2.3, -1.3)]
                          )
    pltbench = 'mcmc_'+object_name+gen_tracks+"_ratios_"+repr(nwalkers)+"w"+repr(nruns)+"r"+"_initb"+img_type
    fig.savefig(os.path.abspath(pltbench))
        
    # plot of the ratios without the benchmark abundances
    fig = triangle.corner(samples, labels=labels)
    pltwithoutbench = 'mcmc_'+object_name+gen_tracks+"_ratios2_"+repr(nwalkers)+"w"+repr(nruns)+"r"+"_initb"+img_type
    fig.savefig(os.path.abspath(pltwithoutbench))    

if use_subset:
    He = subsamples[:, 0]
    O = subsamples[:, 1]
    CO = subsamples[:, 2]
    NO = subsamples[:, 3]
    NeO = subsamples[:, 4]
    SO = subsamples[:, 5]
    stat_info, bpt_data_nearby = get_bestmodel(He, O, CO, NO, NeO, SO, subTO3, subTO2, subprob, subsamples, 
                                                      true_abunds, true_temps, bpt_data)
else:
    stat_info, bpt_data_nearby = get_bestmodel(He, O, CO, NO, NeO, SO, TO3, TO2, prob, samples, true_abunds, true_temps, bpt_data)

# Now make the BPT diagrams
I4861, I5007, I6548, I6563, I6584 = bpt_data
I4861_nearby, I5007_nearby, I6548_nearby, I6563_nearby, I6584_nearby = bpt_data_nearby
print 'len(I4861), len(I4861_nearby)', len(I4861), len(I4861_nearby)
fig = plt.figure(1, figsize=(12, 10))
#plt.title('BPT diagram')
xlab = r'log ([NII]6584/H$\alpha$)'
ylab = r'log ([OIII]5007/H$\beta$)'
plt.xlabel(xlab)
plt.ylabel(ylab)
x0 = np.log10(I6584/I6563)
y0 = np.log10(I5007/100.0)
print 'log ([NII]6584/Ha): ', x0
print 'log ([OIII]5007/Hb): ', y0
# SDSS outlines, columns are:
# log([NII]/Halpha), median of log([OIII]/Halpha) in each bin of log([NII]/Halpha), and 3*sigma deviation around the median, 
# where sigma is the standard deviation of log([OIII]/Halpha) values in each bin of log([NII]/Halpha) in following file:
bpt_contours = 'bpt_contours.txt'
logNHa, medOHaminus3sig, medOHa, medOHaplus3sig = np.loadtxt(bpt_contours, skiprows=1, unpack=True)
plt.xlim(-2.2, 0.5)
plt.ylim(-1.5, 1.5)
plt.plot(x0, y0, 'ko')
plt.plot(logNHa, medOHa, 'b')
plt.plot(logNHa, medOHaminus3sig, 'b--')
plt.plot(logNHa, medOHaplus3sig, 'b--')
# overplot models with T_benchmark+-2500 
x = np.log10(np.array(I6584_nearby)/np.array(I6563_nearby))
y = np.log10(np.array(I5007_nearby)/100.0)
plt.plot(x, y, 'ro')
# Insert name of object in the plot
idx = objects_list.index(object_name)
# Find the corresponding observed points
samp_bpt_txt = 'sampleBPTlines.txt'
obs4861, obs5007, obs6548, obs6563, obs6584 = np.loadtxt(samp_bpt_txt, skiprows=1, usecols=(1,3,5,7,9), unpack=True)
obsO3 = obs5007[idx]
obsN21 = obs6548[idx]
obsHa = obs6563[idx]
obsN22 = obs6584[idx]
obsx = np.log10(obsN22/obsHa)
obsy = np.log10(obsO3/100.0)
print 'observed point: ', obsx, obsy
plt.plot(obsx, obsy, 'g*', ms=10)
plt.text(0.0, -1.3, full_names_list[idx], size='large')
bptfig = object_name+'_BPTdiag.jpg'
#fig.savefig(bptfig)
#print 'Figure ', bptfig, ' saved!'
plt.show()
 

if wt_file:
    # write statistics in final chain file...
    final_chain_file = os.path.abspath('mcmc_'+object_name+gen_tracks+"_chain.txt")   
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
    if g4ha:
        print >> f1, "\n{:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<15} {:<15} {:<20}".format('He', 'O', 'C', 'N', 'Ne', 'S', 'TeO3', 'TeO2', 'ln (-0.5*Chi2)')
    else:
        print >> f1, "\n{:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<15} {:<15} {:<20}".format('He', 'O', 'C/O', 'N/O', 'Ne/O', 'S/O', 'TeO3', 'TeO2', 'ln (-0.5*Chi2)')
    for he, o, co, no, neo, so, to3, to2, chi2 in zip(He, O, CO, NO, NeO, SO, TO3, TO2, Chi2):
        print >> f1, "{:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<15} {:<15} {:<20.3f}".format(he, o, co, no, neo, so, to3, to2, chi2)
    f1.close()
    print '\nPlots made: '
    print pltbench+'\n', pltwithoutbench+'\n', COTemp+'\n', COHO+'\n', CHHO+'\n'
    
    
print '\nCode finished! \n'
