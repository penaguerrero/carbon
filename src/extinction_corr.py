# Sample extinction in PyNeb
# shows how to display available extinction laws, select one or define a new one,
# and do some simple dereddening calculations
# Further examples can be found in other sample scripts

import numpy as np
import matplotlib.pyplot as plt
import pyneb as pn
import os
import science

############################################################################################################################################

'''  Choose parameters to run script  '''

# 1) Select a number from objects_list, i = :
objects_list =['arp252', 'iiizw107', 'iras08208', 'iras08339', 'mrk5', 'mrk960', 'mrk1087', 'mrk1199', 'ngc1741', 
               'pox4', 'sbs0218', 'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol9', 'tol1457']
#       arp252 = 0,  iiizw107 = 1,  iras08208 = 2,  iras08339 = 3,  mrk5 = 4,  mrk960 = 5, mrk1087 = 6,  mrk1199 = 7,  ngc1741 = 8,  
#       pox4 =9,  sbs0218 = 10,  sbs0948 = 11, sbs0926 = 12,  sbs1054 = 13,  sbs1319 = 14,  tol9 =15,  tol1457 = 16
object_number = 9
object_name = objects_list[object_number]

# 2) Do you want to create a unified lines text file?
create_txt = False

############################################################################################################################################

# Convert wavelength to x
def x(wave):
    return 10000. / wave

# Define an extinction law (to be used below)
def my_X(wave):
    x = 10000. / wave
    Rv = 3.1
    X_lin = x/2. # linear part of the extinction law
    X_bump = 0.5*x**2. -6*x + 20. # bump part of the extinction law
    return Rv*np.where(x<5., X_lin, X_bump)

# Define a reddening correction object
RC = pn.RedCorr()

# List the available laws
#RC.printLaws()

# Plot the available laws
#RC.plot(laws='all')
#plt.show()

# Choose the one we intend to use 
RC.law = 'CCM 89'
# or define a new one
#RC.UserFunction = my_X
#RC.law = 'user'

# Plot the selected law as a function of x
# Define an array in lambda to do the plot
'''
wave= np.logspace(2.5, 5, 100)
# Plot commands
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylim([0, 15])
ax.set_xlim([0, 10])
ax.plot(x(wave), my_X(wave), label='%s' % (RC.law))
plt.xlabel('1/$\lambda$ ($\mu^{-1}$)')
plt.ylabel('A$_V$/E(B-V)')
plt.legend(loc='upper left')
plt.show()
'''
#### Read the observed lines from the table of lines_info.txt and normalize to Hbeta
# Path of the text files of wavelengths and fluxes for the objects. 
# This works as long as the folder structure is always the same. This file is assumed to be in
#            /Users/home_direcotry/Documents/AptanaStudio3/src/
# so I want to go one back to find the results folder
results_path = "../results/"
# but just to make sure we are in the right place, get and work with full paths
full_results_path = os.path.abspath(results_path)
text_files_path = os.path.join(full_results_path, "1Dspecs/")
# Go into the object folder
results4object_path = os.path.join(full_results_path, object_name)
# Use nuv, opt, and nir files
specs = [0, 1, 2]
add_str = "_lineinfo"
object_file = os.path.join(results4object_path, object_name)
text_file_list, _ = science.spectrum.get_obj_files2use(object_file, add_str, specs)

# Define the file name for for all the lines
name_out_file = os.path.join(results4object_path, object_name+"_linesNUV2NIR.txt")

if create_txt == True:
    # Get the single textfile with all the lines and save it into 
    cols_in_file = science.spectrum.gather_specs(text_file_list, name_out_file, reject=50.0, start_w=None, create_txt=create_txt)
else:
    cols_in_file = science.spectrum.readlines_from_lineinfo(name_out_file)
catalog_wavelength, observed_wavelength, element, ion, forbidden, how_forbidden, width, flux, continuum, EW = cols_in_file

# Create an array of the numbers
data = np.array([catalog_wavelength, observed_wavelength, flux, continuum])
# data: 0=catalog_wavelength, 1=observed_wavelength, 2=flux, 3=continuum

# Normalize observed line ratios to Hbeta
_, Hb_idx = science.spectrum.find_nearest(data[0], 4861.0)
Hbeta = data[2][Hb_idx]
norm_fluxes = (data[2] / Hbeta) * 100
print norm_fluxes[Hb_idx]
# Recalculate the equivalent widths: Possitive = Emission, Negative = absorption
EWs = []
for f,c in zip(data[2], data[3]):
    new_ew = f / (c) #* (-1)
    EWs.append(new_ew)
for w, f, oe, e in zip(data[0], norm_fluxes, EW, EWs):
    print w, '    norm_flux=', f, '    old_EW=', oe, '    newEW=', e

# Step 1 of first iteration of reddening correction: Assume that there is no collisional excitation
# get the EW_abs of the H and He lines with respect to EW_abs(Hbeta)
Hline_and_EWs, Heline_and_EWs = science.spectrum.readlines_EWabsRelHbeta()
# if both H and He lines are going to be used to correct for underlying absorption
undabs_wav = []
undabs_EW = []
for i in range(len(Hline_and_EWs[0])):
    undabs_wav.append(Hline_and_EWs[0][i])
    undabs_EW.append(Hline_and_EWs[1][i])
for i in range(len(Heline_and_EWs[0])):
    undabs_wav.append(Heline_and_EWs[0][i])
    undabs_EW.append(Heline_and_EWs[1][i])
lines_undabs_and_EW = [undabs_wav, undabs_EW]

# Set initial value of EWabsHbeta (this is a guessed value taken from HII regions)
EWabsHbeta = 0.280
# Remove underlying absorption for optical lines to get Intensities
intensities = []
for i in range(len(observed_wavelength)):
    if observed_wavelength[i] < 3000.0:
        I = flux[i]
    else:
        if observed_wavelength[i] in lines_undabs_and_EW[0]:
            EWabsLine = undabs_EW.index(observed_wavelength[i])
            I = EWabsLine * EWabsHbeta * continuum[i] * flux[i]
    intensities.append(I)

# Step 2 of first iteration of reddening correction: Using Seaton
# set theoretical Halpha/Hbeta ratio
I_theo_HaHb = 2.86 
# Round all catalog lines
rounded_catalog_wavelength = []
for item in catalog_wavelength:
    rw = int(item)
    rounded_catalog_wavelength.append(rw)
# Find observed Halpha/Hbeta ratio
Halpha_idx = rounded_catalog_wavelength.index(6563)
Halpha = norm_fluxes[Halpha_idx]
Hbeta_idx = rounded_catalog_wavelength.index(4861)
Hbeta = norm_fluxes[Hbeta_idx]
I_obs_HaHb = Halpha/Hbeta
print 'I_theo_HaHb = %0.2f      I_obs_HaHb = %0.2f' % (I_theo_HaHb, I_obs_HaHb)
# Correct based on the given law and the observed Ha/Hb ratio
RC = pn.RedCorr(law='CCM 89')

wave1 = 5007
_, idx = science.spectrum.find_nearest(data[0], 5007)
I_obs1 = 4.0
my_I_obs = norm_fluxes[idx]
wave2 = 4686
I_obs2 = 0.10
#RC.setCorr(I_obs_HaHb / I_theo_HaHb, 6563., 4861.)
RC.setCorr(I_obs_HaHb / I_theo_HaHb, Halpha, Hbeta)

print 'Correct based on the given law and the observed Ha/Hb ratio:'
print str(wave1) + ': I_obs =', I_obs1, ' I_dered =', I_obs1 * RC.getCorrHb(wave1), 'CHRISTOPHE'
print str(wave2) + ': I_obs =', I_obs2, ' I_dered =', I_obs2 * RC.getCorrHb(wave2), 'CHRISTOPHE'
print '5007' + ': I_obs =', my_I_obs, ' I_dered =', my_I_obs * RC.getCorrHb(wave1), 'YO'

# Correct based on the given law and c(Hb)
RC = pn.RedCorr(law='CCM 89', cHbeta=0.3)
print '\nCorrect based on the given law and the observed Ha/Hb ratio:'
print str(wave1) + ': I_obs =', I_obs1, ' I_dered =', I_obs1 * RC.getCorrHb(wave1)
print str(wave2) + ': I_obs =', I_obs2, ' I_dered =', I_obs2 * RC.getCorrHb(wave2)

