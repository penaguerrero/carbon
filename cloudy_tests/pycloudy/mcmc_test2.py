import os
import numpy as np
import matplotlib.pyplot as plt
import pyCloudy as pc
import time
#import mcmc_infrastructure as mcmcis
import mcmc_infrastruc_multipross as mcmcis

##################################################################################################################

# Define initial conditions for Cloudy model
model_name = 'mcmc_test3'
dens = 4. #log cm-3
stb99_table = 'table star "starburst99.mod"'
age = 'age=5.0'

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

abunds = {'He' : 10.97 - 12, 'C' : 7.12 - 12, 'N' : 7.34 - 12, 'O' : 8.01 - 12, 'Ne' : 7.72 - 12, 
         'S' : 6.65 - 12}

options = ('no molecules',
            'no level2 lines',
            'no fine opacities',
            'atom h-like levels small',
            'atom he-like levels small',
            'COSMIC RAY BACKGROUND',
            'element limit off -8',
            'print line optical depth', )

##################################################################################################################

# start the timer to compute the whole running time
start_time = time.time()

# Changing the location and version of the cloudy executable. Will turn the relative path into an absolute path.
cloudyexe_path = '../../../../Addons/cloudy/c13.03all/c13.03/source/cloudy.exe'
mcmcis.find_cloudyexe(cloudyexe_path)

# Define path where to save the plots
pypics_path = os.path.abspath('pypics')
verbosity = 3
iterations = 2
keep_files = None
dir = './'
initial_conditions = mcmcis.PyCloudy_model(model_name, dens, emis_tab, abunds, stb99_table, age, dir, verbosity, options, 
                                           iterations, keep_files)
cldymod = initial_conditions.mk_model()

modeled_lines = os.path.abspath(cldymod)
object_name = 'mrk960' 
type_of_measurement = '_RedCor_Ebv.txt'  # my code's line measurements
used_fiven_lineinfo = False
if used_fiven_lineinfo:
    type_of_measurement = '_measuredLI_RedCor_Ebv.txt'   # my manual line measurements
measured_lines = os.path.abspath('../../results/'+object_name+'/'+object_name+type_of_measurement)
files4chi2 = mcmcis.Get_Chi2_info(measured_lines, modeled_lines)
Chi2 = files4chi2.getChi2()


print '\nCode finished! Took  %s  seconds to finish.' % (time.time() - start_time)

