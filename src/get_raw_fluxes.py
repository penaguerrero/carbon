import numpy
import os
from science import spectrum

'''
This code finds the net fluxes (or raw intensities) and equivalent widths from the redshift-corrected spectra of the given object.
It writes a text file containing the line information it found for that wavelength range.
    * The possible lines the code looks for are all those in the line_catalog.txt file in the spectrum folder
'''

############################################################################################################################################

''' OBJECTS OF THE SAMPLE '''

# name of the object
#                 0           1           2            3         4        5          6          7         8
objects_list =['arp252', 'iiizw107', 'iras08208', 'iras08339', 'mrk5', 'mrk960', 'mrk1087', 'mrk1199', 'ngc1741', 
               'pox4', 'sbs0218', 'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol9', 'tol1457']
#                9         10        11         12         13         14        15       16

# corresponding redshift
#             0        1        2        3         4         5        6        7         8
z_list = [0.032989, 0.01972, 0.04678, 0.19113, 0.002695, 0.021371, 0.02877, 0.013454, 0.01348, 
          0.01201, 0.05842, 0.046240, 0.013642, 0.002010, 0.0076, 0.01195, 0.01763]
#             9       10        11       12        13       14       15       16


############################################################################################################################################

# Choose parameters to run script

# 1) Select a number from objects_list, i = :
#       arp252 = 0,  iiizw107 = 1,  iras08208 = 2,  iras08339 = 3,  mrk5 = 4,  mrk960 = 5, mrk1087 = 6,  mrk1199 = 7,  ngc1741 = 8,  
#       pox4 =9,  sbs0218 = 10,  sbs0948 = 11, sbs0926 = 12,  sbs1054 = 13,  sbs1319 = 14,  tol9 =15,  tol1457 = 16
object_number = 9
object_name = objects_list[object_number]
z = z_list[object_number]

# 2) use all 3 files for NUV, optical, and NIR? Type which ones to use: nuv=0, opt=1, nir=2
specs = [1]

# 3) Do you want to use Vacuum wavelengths?
vacuum = False

# 4) Do you want to normalize the spectra to the continuum?
normalize = False

# 5) Choose how many sigmas to clip from the continuum array
sigmas_away = 3

# in case I want to use a specific order for the polynomial, else it will be determined by the algorithm
order = 3

# 6) What is the width of the window to use to find local continuum?
window = 600

# 7) Do you want to see the plots of the fitted continuum?
plot = True

# 8) write the text file with the line net fluxes and equivalent widths?
text_table = True

# Set width of Halpha in order to properly correct for reddening
Halpha_width = 43.


############################################################################################################################################

# Path of the text files of wavelengths and fluxes for the objects. 
# This works as long as the folder structure is always the same. This file is assumed to be in
#            /Users/home_direcotry/Documents/AptanaStudio3/src/
# so I want to go one back to find the results folder
results_path = "../results/"
# but just to make sure we are in the right place, get and work with full paths
full_results_path = os.path.abspath(results_path)
text_files_path = os.path.join(full_results_path, "1Dspecs/")
results4object_path = os.path.join(full_results_path, object_name)

add_str = "_selectedspecs"
data, full_file_list = spectrum.loadtxt_from_files(object_name, add_str, specs, text_files_path)
# alternative for when files have been corrected in splot
altern = '../results/sbs1319/sbs1319_opt_corr1.txt'
#f, w = numpy.loadtxt(altern, skiprows=5, usecols=(1,2), unpack=True)
#data = [numpy.array([w,f])]

# Terminations used for the lines text files
spectrum_region = ["_nuv", "_opt", "_nir"]

for d, s in zip(data, specs):
    print 'Working with:  %s' % full_file_list[s]
    # Correct spectra for redshift and calculate continuum with a polynomial of nth order
    object_spectra, fitted_continuum, err_cont_fit = spectrum.fit_continuum(object_name, d, z, order=order, sigmas_away=sigmas_away, window=window, plot=plot, z_correct=True, normalize=normalize)

    # Obtain the lines net fluxes and EWs
    new_file_name = object_name+"_lineinfo"+spectrum_region[s]+".txt"
    lineinfo_text_file = os.path.join(results4object_path, new_file_name)
    
    object_lines_info = spectrum.find_lines_info(object_spectra, fitted_continuum, err_cont_fit, lineinfo_text_file, Halpha_width=Halpha_width, text_table=text_table, vacuum=vacuum)
    print ''
print sigmas_away, 'sigmas_away'    
print 'Code finished!'
    
    
