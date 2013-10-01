import numpy
import os
from science import spectrum

'''
This code finds the net fluxes (or raw intensities) and equivalent widths from the redshift-corrected spectra of the given object.
It writes a text file containing the line information it found for that wavelength range.
    * The possible lines the code looks for are all those in the line_catalog.txt file in the spectrum folder
'''

##################################################################################################################################

# name of the object
object_name = 'sbs1319'

# redshift of the object
z = 0.0077

# use all 3 files for NUV, optical, and NIR?
whole_spectrum = True
# if False, type which ones to use: nuv=0, opt=1, nir=2
specs = [1]

# write the text file with the line net fluxes and equivalent widths?
text_table = False

##################################################################################################################################

# Path of the text files of wavelengths and fluxes for the objects. 
# This works as long as the folder structure is always the same. This file is assumed to be in
#            /Users/home_direcotry/Documents/AptanaStudio3/src/
# so I want to go one back to find the results folder
results_path = "../results/"
# but just to make sure we are in the right place, get and work with full paths
full_results_path = os.path.abspath(results_path)
text_files_path = os.path.join(full_results_path, "1Dspecs/")
results4object_path = os.path.join(full_results_path, object_name)

nuv = object_name+"_selectedspecs_nuv.txt"
opt = object_name+"_selectedspecs_opt.txt"
nir = object_name+"_selectedspecs_nir.txt"
text_file_list = [nuv, opt, nir]

if whole_spectrum == False:
    new_text_file_list = []
    for item in specs:
        tf = text_file_list[item]
        new_text_file_list.append(tf)
    text_file_list = new_text_file_list
elif whole_spectrum == True:
    specs = [0, 1, 2]
    
# List of the data contained in each file
data = []
for i in range(len(text_file_list)):
    txtfile_path = os.path.join(text_files_path, text_file_list[i])
    f = open(txtfile_path, 'r')
    # The file is expected to be two columns without header: wavelengths, fluxes
    data_file = numpy.loadtxt(f, unpack=True)
    f.close()
    data.append(data_file)

# Terminations used for the lines text files
spectrum_region = ["_nuv", "_opt", "_nir"]
for d, s in zip(data, specs):
    print 'Working with:  %s' % text_file_list[s]
    # Correct spectra for redshift and calculate continuum with a polynomial of nth order
    object_spectra, fitted_continuum = spectrum.fit_continuum(d, z, nth=5, thresold_fraction=1, window_wdith=150, normalize=False)
    
    # Obtain the lines net fluxes and EWs
    new_file_name = object_name+spectrum_region[s]+"_lineinfo.txt"
    lineinfo_text_file = os.path.join(results4object_path, new_file_name)
    object_lines_info = spectrum.find_lines_info(object_spectra, fitted_continuum, lineinfo_text_file, text_table=text_table, vacuum=False)
    print ''
    
print 'Code finished!'
    
    
