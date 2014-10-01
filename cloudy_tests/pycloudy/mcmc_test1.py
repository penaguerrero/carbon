import os
import numpy as np
import matplotlib.pyplot as plt
import pyCloudy as pc
import time
import pickle

# start the timer to compute the whole running time
start_time = time.time()

# Changing the location and version of the cloudy executable.
cloudyexe_path = '../../../../Addons/cloudy/c13.03all/c13.03/source/cloudy.exe'
full_cloudyexe_path = os.path.abspath(cloudyexe_path)
pc.config.cloudy_exe = full_cloudyexe_path
# Define plots path
pypics_path = os.path.abspath('pypics')
# Define the path for the cloudy data
cloudydata_path = '../../../../Addons/cloudy/c13.03all/c13.03/data'
full_cloudydata_path = os.path.abspath(cloudydata_path)

# Define verbosity to high level (will print errors, warnings and messages)
pc.log_.level = 0#3
# The directory in which we will have the model
# You may want to change this to a different place so that the current directory
# will not receive all the Cloudy files.
dir_ = './'

# Define some parameters of the model:
model_name = 'mcmc_test1'
full_model_name = '{0}{1}'.format(dir_, model_name)
dens = 2. #log cm-3
options = ('no molecules',
            'no level2 lines',
            'no fine opacities',
            'atom h-like levels small',
            'atom he-like levels small',
            'COSMIC RAY BACKGROUND',
            'element limit off -8',
            'print line optical depth', )

emis_tab = ['H  1  4861',
            'H  1  6563',
            'He 1  5876',
            'He 2  4686',
            'N  2  6584',
            'O  1  6300',
            'O II  3726',
            'O II  3729',
            'O  3  5007',
            'TOTL  4363',
            'S II  6716',
            'S II 6731',
            'Cl 3 5518',
            'Cl 3 5538',
            'O  1 6300',
            'C  2 4267',
            'C  3 1907']

abund = {'He' : 10.97 - 12, 'C' : 7.12 - 12, 'N' : 7.34 - 12, 'O' : 8.01 - 12, 'Ne' : 7.72 - 12, 
         'S' : 6.65 - 12}#, 'Ar' : -5.80, 'Fe' : -7.4, 'Cl' : -7.00}

# Defining the object that will manage the input file for Cloudy
c_input = pc.CloudyInput(full_model_name)
# Filling the object with the parameters
# Defining the ionizing SED: Effective temperature and luminosity.
# The lumi_unit is one of the Cloudy options, like "luminosity solar", "q(H)", "ionization parameter", etc... 
c_input.set_star(SED = 'table star "starburst99.mod"', SED_params = 'age=4.0', lumi_unit='f(nu)', lumi_value=-12.3316)
# Defining the density. You may also use set_dlaw(parameters) if you have a density law defined in dense_fabden.cpp.
c_input.set_cste_density(dens)
c_input.set_abund(ab_dict = abund, nograins = True)
c_input.set_other(options)
c_input.set_emis_tab(emis_tab) # better use read_emis_file(file) for long list of lines, where file is an external file.
c_input.set_iterate(3) # (0) for no iteration, () for one iteration, (N) for N iterations.
# Writing the Cloudy inputs. to_file for writing to a file (named by full_model_name). verbose to print on the screen.
c_input.print_input(to_file = True, verbose = False)
# Printing some message to the screen
pc.log_.message('Running {0}'.format(model_name), calling = 'stb99_test1')

# Running Cloudy with a timer. Here we reset it to 0.
#pc.log_.timer('Starting Cloudy', quiet = True, calling = 'stb99_test1')
c_input.run_cloudy()
#pc.log_.timer('Cloudy ended after seconds:', calling = 'stb99_test1')

# Reading the Cloudy outputs in the Mod CloudyModel object
a = pc.CloudyModel(full_model_name)
file_Name = "testfile"
# open the file for writing
fileObject = open(file_Name,'wb') 
# this writes the object a to the
# file named 'testfile'
pickle.dump(a,fileObject)   
# here we close the fileObject
fileObject.close()
# we open the file for reading
fileObject = open(file_Name,'r')  
# load the object from the file into var b
Mod = pickle.load(fileObject)  


dir(Mod) # This is the online answering way
Mod.print_stats()
Mod.print_lines()
Mod.get_ab_ion_vol_ne('O',2)
Mod.get_T0_ion_vol_ne('O', 2)
Mod.log_U_mean
Mod.log_U_mean_ne
print('T0 = {0:7.1f}K, t2 = {1:6.4f}'.format(Mod.T0, Mod.t2))
print('Hbeta Equivalent width = {0:6.1f}, Hbeta Surface Brightness = {1:4.2e}'.format(Mod.get_Hb_EW(), Mod.get_Hb_SB()))
# printing line intensities
for line in Mod.emis_labels:
    print('{0} {1:10.3e} {2:7.2f}'.format(line, Mod.get_emis_vol(line), Mod.get_emis_vol(line) / Mod.get_emis_vol('H__1__4861A') * 100.))


print '\nCode finished! Took  %s  seconds to finish.' % (time.time() - start_time)

