import numpy
import os
from science import spectrum

'''
This code finds the net fluxes (or raw intensities) and equivalent widths from the redshift-corrected spectra of the given object.
It writes a text file containing the line information it found for that wavelength range.
    * The possible lines the code looks for are all those in the line_catalog.txt file in the spectrum folder
'''

##################################################################################################################################

''' OBJECTS OF THE SAMPLE '''

# name of the object
objects_list =['arp252', 'iiizw107', 'iras08208', 'iras08339', 'mrk5', 'mrk960', 
               'mrk1087', 'mrk1199', 'ngc1741', 'pox4', 'sbs0218', 'sbs0948', 
               'sbs0926', 'sbs1054', 'sbs1319', 'tol9', 'tol1457']
# corresponding redshifts
z_list = [0.032989, 0.01972, 0.04678, 0.19113, 0.002695, 0.021371,
          0.02877, 0.013454, 0.01348, 0.01997, 0.058420, 0.046240,
          0.013642, 0.002010, 0.006870, 0.011625, 0.01763]

##################################################################################################################################

# Choose parameters to run script

# 1) Select a number from objects_list, i = :
#              arp252 =0, iiizw107 =1, iras08208 =2, iras08339 =3, mrk5 =4, mrk960 =5, 
#              mrk1087 =6, mrk1199 =7, ngc1741 =8, pox4 =9, sbs0218 =10,  sbs0948 =11, 
#              sbs0926 =12, sbs1054 =13, sbs1319 =14, tol9 =15, tol1457 =16
object_number = 14
object_name = objects_list[object_number]
z = z_list[object_number]

# 2) use all 3 files for NUV, optical, and NIR?
whole_spectrum = False
# 2.a) if False, type which ones to use: nuv=0, opt=1, nir=2
specs = [1]

# 3) Do you want to use Vacuum wavelengths?
vacuum = False

# 4) Do you want to normalize the spectra to the continuum?
normalize = False

# 5) Choose the confidence interval to fit a good continuum
sigmas_away = 1.5

# 6) Order of the polynomial for the continuum fit
order = 11

# 6) Do you want to see the plots of the fitted continuum?
plot = True

# 7) write the text file with the line net fluxes and equivalent widths?
text_table = True


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
full_file_list = [nuv, opt, nir]

if whole_spectrum == False:
    new_text_file_list = []
    for item in specs:
        tf = full_file_list[item]
        new_text_file_list.append(tf)
    text_file_list = new_text_file_list
elif whole_spectrum == True:
    specs = [0, 1, 2]
    text_file_list = full_file_list
    
# List of the data contained in each file
data = []
for i in range(len(text_file_list)):
    txtfile_path = os.path.join(text_files_path, text_file_list[i])
    print 'Opening: ', txtfile_path
    # The file is expected to be two columns without header: wavelengths, fluxes
    data_file = numpy.loadtxt(txtfile_path, unpack=True)
    #print 'LIMITS', data_file[0][0], max(data_file[0])
    data.append(data_file)

# Terminations used for the lines text files
spectrum_region = ["_nuv", "_opt", "_nir"]
for d, s in zip(data, specs):
    print 'Working with:  %s' % full_file_list[s]
    # Correct spectra for redshift and calculate continuum with a polynomial of nth order
    object_spectra, fitted_continuum = spectrum.fit_continuum(object_name, d, z, order=order, sigmas_away=sigmas_away, window=350, plot=plot, z_correct=True, normalize=normalize)

    # Obtain the lines net fluxes and EWs
    new_file_name = object_name+spectrum_region[s]+"_lineinfo.txt"
    lineinfo_text_file = os.path.join(results4object_path, new_file_name)
    object_lines_info = spectrum.find_lines_info(object_spectra, fitted_continuum, lineinfo_text_file, text_table=text_table, vacuum=vacuum)
    print ''
    
print 'Code finished!'
    
    
