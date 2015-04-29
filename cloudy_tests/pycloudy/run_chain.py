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
abunds_lists = [[10.94,  8.20,  7.32,  7.07,  7.40,  6.42],   # 0 = iiizw107  ** S from Lop-San09
                [10.91,  8.31,  8.30,  7.51,  7.54,  6.55],   # 1 = iras08339  ** S from Lop-San09 = --
                [10.97,  8.20,  8.52,  6.77,  7.72,  6.50],   # 2 = mrk1087 ** does not have a "true" initial S value
                [10.79,  7.92,  7.84,  7.26,  7.81,  6.92],   # 3 = mrk1199  ** S from Lop-San09
                [10.91,  8.07,  6.82,  6.71,  7.35,  6.44],   # 4 = mrk5  ** S from Lop-San09
                [10.97,  8.01,  7.12,  7.34,  7.72,  6.65],   # 5 = mrk960  ** S from Lop-San09
                [10.94,  8.00,  6.94,  6.60,  7.23,  5.81],   # 6 = ngc1741  ** S from Lop-San09 = 6.36
                [10.88,  8.00,  7.16,  6.50,  7.26,  6.02],   # 7 = pox4  ** S from Lop-San09 = 6.24
                [10.94,  7.92,  6.59,  6.84,  7.24,  6.45],   # 8 = sbs0218  ** S from Lop-San09 = 6.29                
                [10.88,  8.07,  7.38,  6.16,  7.30,  6.34],   # 9 = sbs0948  ** S from Lop-San09 
                [10.94,  7.93,  7.20,  6.48,  6.38,  6.34],   # 10 = sbs0926  ** S from Lop-San09 = 6.29
                [10.88,  7.97,  7.33,  6.59,  7.33,  6.21],   # 11 = sbs1054  ** S from Lop-San09
                [10.94,  8.17,  7.21,  6.52,  7.35,  6.17],   # 12 = sbs1319  ** S from Lop-San09 = 6.29
                [10.95,  7.83,  6.32,  6.27,  7.00,  6.18],   # 13 = tol1457  ** S from Lop-San09
                [10.93,  8.64,  8.75,  7.80,  7.84,  7.10],   # 14 = tol9  ** S from Lop-San09 = 6.97
                [10.97,  8.32,  8.11,  7.71,  7.90,  7.00],   # 15 = arp252  ** S from Lop-San09 = --
                [11.09,  8.55,  7.81,  7.44,  7.66,  6.28],   # 16 = iras08208  ** S from Lop-San09 = 6.39
                [10.77,  7.56,  6.00,  6.04,  6.67,  5.89]]   # 17 = sbs1415  ** S from Lop-San09 
abunds = abunds_lists[object_number]

theta = []
theta.append(abunds[0])
theta.append(abunds[1])
for i in range(2, len(abunds)):
    theta.append(abunds[i]-abunds[1])

print 'abundances of: He =', abunds[0], ' O =', abunds[1], ' C =', abunds[2], 'N =', abunds[3], 'Ne =', abunds[4], 'S =', abunds[5]
print 'INITIAL values of the 6 dimensions used in the MCMC:'
print '    He =', theta[0], '   O =', theta[1], '   C/O = ', theta[2], '   N/O = ', theta[3], '   Ne/O = ', theta[4], '   S/O = ', theta[5]

# Was the measurement of the line intensities manual or with the code?
#                            0     1     2     3     4       5      6      7     8
manual_measurement_lists = [False, True, True, True, False, False, False, False, False, 
                            True, False, True, False, False, True, False, True, False]
#                            9     10     11     12    13     14     15    16    17
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


