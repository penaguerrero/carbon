import os
import numpy as np
import matplotlib.pyplot as plt
import pyCloudy as pc
import time
import mcmc_infrastructure as mcmcis

##################################################################################################################

# Define initial conditions for Cloudy model
model_name = 'mcmc_test_pox4'
stb99_table = 'table star "constSFR.mod"'

emis_tab = ['H  1  4861',   # must keep this format (as lines appear on Cloudy output)
            'H  1  6563',
            'H  1  4340',
            'He 1  5876',
            'He 2  4686',   
#            'N  2  6584',    # could not separate N2 from Halpha
#            'Ne 3  3869',    # did not trust our Ne abundances
#            'O  1  6300',    # generally weak line
#            'O II  3726',
#            'O II  3729',
            'TOTL  3727',
            'O  3  5007',
#            'TOTL  4363',
            'S II  6716',
#            'S II  6731',
            'S  3  9532',
#            'Cl 3  5518',
#            'Cl 3  5538',
#            'C  2  4267',
            'C  3  1907',
            'C  3  1910',
            'TOTL  1909']

options = ('no molecules',
            'no level2 lines',
            'no fine opacities',
            'atom h-like levels small',
            'atom he-like levels small',
            'COSMIC RAY BACKGROUND',
            'element limit off -8',
            'print line optical depth', )

age = 'age=5.0'
dens = 4. #log cm-3
#abunds = {'He' : 10.97 - 12, 'C' : 7.12 - 12, 'N' : 7.34 - 12, 'O' : 8.01 - 12, 'Ne' : 7.72 - 12, 'S' : 6.65 - 12}
# Abundances of:
#           He     O      C      N      Ne     S
#abunds = [10.97,  8.01, 7.12,  7.34,  7.72,  6.65]   # mrk960
abunds = [10.88,  8.00, 7.16,  6.50,  7.26,  6.24]   # pox4

##################################################################################################################
theta = []
theta.append(abunds[0])
theta.append(abunds[1])
for i in range(2, len(abunds)):
    theta.append(abunds[i]-abunds[1])

print 'abundances of: He =', abunds[0], ' O =', abunds[1], ' C =', abunds[2], 'N =', abunds[3], 'Ne =', abunds[4], 'S =', abunds[5]
print 'INITIAL values of the 6 dimensions used in the MCMC:'
#print '   density of H =', dens
#print '   abundances of: He =', abunds['He'], ' C =', abunds['C'], ' N =', abunds['N'], 'O =', abunds['O'], 'Ne =', abunds['Ne'], 'S =', abunds['S']
print '    He =', theta[0], '   O =', theta[1], '   C/O = ', theta[2], '   N/O = ', theta[3], '   Ne/O = ', theta[4], '   S/O = ', theta[5]

# start the timer to compute the whole running time
start_time = time.time()

# Changing the location and version of the cloudy executable. Will turn the relative path into an absolute path.
cloudyexe_path = '../../../../Addons/cloudy/c13.03all/c13.03/source/cloudy.exe'
mcmcis.find_cloudyexe(cloudyexe_path)

# Define path where to save the plots
pypics_path = os.path.abspath('pypics')

object_name = 'pox4' 
manual_measurement = False
dir = './'
verb = 3
iterations = 2
keep_files = None
initial_Cloudy_conditions = [model_name, dens, emis_tab, theta, stb99_table, age, dir, verb, options, iterations, keep_files] 
mchammer = mcmcis.MCMC(object_name, manual_measurement, initial_Cloudy_conditions)
mchammer.run_chain()

print '\nCode finished! Took  %s  seconds to finish.' % (time.time() - start_time)

