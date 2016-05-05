import os
import numpy as np
import matplotlib.pyplot as plt
import pyCloudy as pc
import time
import copy
import mcmc_infrastructure4 as mcmcis
from sys import argv

##################################################################################################################
# name of the object
#                  0            1           2          3         4        5         6         7         8
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                 9           10         11         12         13       14        15         16         17
object_number = int(argv[1]) 
recover = True   # set to True if you want to recover a previous chain and start a new one from there
debug = True

object_name = objects_list[object_number] 


# Define initial conditions for Cloudy model

#stb99_table = 'table star "constSFR.mod"'  # --> used for ALL models of ALL objects = cSFRg40_hZ_sUV.mod
#stb99_table = 'table star "cSFRg00_hZ_sUV.mod"'   # already ran for tol9 and sbs1415
#stb99_table = 'table star "cSFRg00_lZ_sUV.mod"'
stb99_table = 'table star "cSFRg40_hZ_sUV.mod"'  # equivalent to constSFR.mod
#stb99_table = 'table star "cSFRg40_lZ_sUV.mod"'

emis_tab = ['H  1  4861',   # must keep this format (as lines appear on Cloudy output)
            'H  1  6563',
            'H  1  4340',
            'He 1  5876',
            'He 2  4686',   
            'N  2  6548',    # could not separate N2 from Halpha
            'N  2  6584',    # could not separate N2 from Halpha
#            'Ne 3  3869',    # did not trust our Ne abundances
#            'O  1  6300',    # generally weak line
#            'O II  3726',
#            'O II  3729',
            'TOTL  3727',
            'O  3  5007',
#            'TOTL  4363',
#            'S II  6716',    # generally weak line
#            'S II  6731',
            'S  3  9532',
#            'Cl 3  5518',
#            'Cl 3  5538',
#            'C  2  4267',
            'C  3  1907',
            'C  3  1910',
            'TOTL  1909']    # Currently using 11 line to find Chi2

options = ('no molecules',
            'no level2 lines',
            'no fine opacities',
            'atom h-like levels small',
            'atom he-like levels small',
            'COSMIC RAY BACKGROUND',
            'element limit off -8',
            'print line optical depth', )

age = 'age=8.0'
dens = 2.7 #log cm-3

##################################################################################################################

model_name = 'mcmc_'+object_name

# Abundances of:
#                  He     O      C      N      Ne     S
abunds_lists = [[10.94,  8.22,  7.50,  7.07,  7.30,  6.42],   # 0 = iiizw107  ** S from Lop-San09
                [10.91,  8.42,  8.27,  7.51,  8.25,  7.13],   # 1 = iras08339  ** S from Lop-San09 = --
                [10.84,  8.35,  8.81,  7.76,  8.12,  6.79],   # 2 = mrk1087 ** does not have a "true" initial S value
                [10.79,  8.29,  8.88,  7.90,  7.66,  6.75],   # 3 = mrk1199  ** S from Lop-San09
                [10.91,  7.96,  6.95,  6.71,  7.54,  6.46],   # 4 = mrk5  ** S from Lop-San09
                [10.97,  7.86,  6.90,  7.34,  7.29,  6.71],   # 5 = mrk960  ** S from Lop-San09
                [10.94,  8.25,  7.83,  7.10,  7.68,  6.00],   # 6 = ngc1741  ** S from Lop-San09 = 6.36
                [10.91,  8.05,  7.36,  6.50,  7.29,  6.08],   # 7 = pox4  ** S from Lop-San09 = 6.24
                [10.88,  7.99,  6.82,  6.84,  7.27,  6.36],   # 8 = sbs0218  ** S from Lop-San09 = 6.29                
                [10.88,  8.03,  7.33,  6.16,  7.44,  6.34],   # 9 = sbs0948  ** S from Lop-San09 
                [10.94,  7.98,  7.29,  6.48,  6.81,  6.49],   # 10 = sbs0926  ** S from Lop-San09 = 6.29
                [10.88,  8.06,  7.41,  6.59,  7.62,  6.21],   # 11 = sbs1054  ** S from Lop-San09
                [10.94,  8.19,  7.66,  6.52,  7.46,  6.05],   # 12 = sbs1319  ** S from Lop-San09 = 6.29
                [10.95,  7.83,  6.36,  6.27,  7.09,  6.18],   # 13 = tol1457  ** S from Lop-San09
                [10.93,  8.58,  8.54,  7.80,  8.10,  6.96],   # 14 = tol9  ** S from Lop-San09 = 6.97
                [11.09,  8.13,  7.98,  7.71,  7.89,  6.89],   # 15 = arp252  ** S from Lop-San09 = --
                [10.97,  8.43,  7.54,  7.44,  7.82,  6.75],   # 16 = iras08208  ** S from Lop-San09 = 6.39
                [10.77,  7.56,  5.92,  6.04,  6.66,  5.89]]   # 17 = sbs1415  ** S from Lop-San09 
abunds = abunds_lists[object_number]

theta = []
theta.append(abunds[0])
theta.append(abunds[1])
for i in range(2, len(abunds)):
    theta.append(abunds[i]-abunds[1])

print 'abundances of: He =', abunds[0], ' O =', abunds[1], ' C =', abunds[2], 'N =', abunds[3], 'Ne =', abunds[4], 'S =', abunds[5]
print 'INITIAL values of the 6 dimensions used in the MCMC:'
print '    He =', theta[0], '   O =', theta[1], '   C/O = ', theta[2], '   N/O = ', theta[3], '   Ne/O = ', theta[4], '   S/O = ', theta[5]

# Was the measurement of the line intensities manual or with the code (this is for E(B-V))?
#                            0     1     2     3     4       5      6      7     8
manual_measurement_lists = [False, True, True, True, False, False, False, False, False, 
                            True, False, True, False, False, True, False, True, False]
#                            9     10     11     12    13     14     15    16    17
'''
# For CHbeta
#                            0     1     2     3      4       5      6      7     8
manual_measurement_lists = [False, True, False, False, False, False, False, False, False, 
                           True, False, False, False, False, True, False, False, False]
#                            9     10     11     12    13     14     15    16    17

'''
manual_measurement = manual_measurement_lists[object_number]

# start the timer to compute the whole running time
start_time = time.time()


dir = './'
verb = 0
iterations = 3
keep_files = None
initial_Cloudy_conditions = [model_name, dens, emis_tab, theta, stb99_table, age, dir, verb, options, iterations, keep_files] 

modeled_temperatures_list = []
nwalkers = 100 # if grater than 50, please use multiples of 50 
nruns = 100
threads = 25
true_abunds = copy.deepcopy(theta)

# With regular the probability function that allows to use HTCondor
init_Cldy_conds_file = os.path.abspath(object_name+'_initial_Cloudy_conditions.txt')
f = open(init_Cldy_conds_file, 'w+')
for item in initial_Cloudy_conditions:
    print >> f, item
f.close()

if recover == True:
    keep_recover = raw_input('Code set to recover from last run, keep this mode?  [y]  n   \n')
    if keep_recover == 'n':
        recover = False
mchammer = mcmcis.run_chain_and_plot(model_name, dir, true_abunds, theta, nwalkers, nruns, object_name, manual_measurement, 
                                     init_Cldy_conds_file, modeled_temperatures_list, threads=threads, recover=recover, debug=debug)

print '\nCode finished! Took  %s  days to finish.' % ((time.time() - start_time) / 86400.0)


