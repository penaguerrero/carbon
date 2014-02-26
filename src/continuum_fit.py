#import numpy
import os
from science import spectrum


############################################################################################################################################

''' OBJECTS OF THE SAMPLE '''

# name of the object
#                 0           1           2            3         4        5          6        7         8
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                 9           10         11         12         13       14        15         16         17


# Choose parameters to run script

# 1) Select a number from objects_list, i = :
object_number = 12

# 2) use all 3 files for NUV, optical, and NIR? Type which ones to use: nuv=0, opt=1, nir=2
specs = [2]

# 3) Do you want to use Vacuum wavelengths?
vacuum = False

# 4) Do you want to normalize the spectra to the continuum?
normalize = False

# 5) Choose how many sigmas to clip from the continuum array
sigmas_away = 2

# in case I want to use a specific order for the polynomial, else it will be determined by the algorithm
order = 1

# 6) What is the width of the window to use to find local continuum?
window = 550

# 7) Do you want to see the plots of the fitted continuum?
plot = True

# 8) write the text file with the line wavelengths, fluxes, and fitted continuum?
text_table = False

# Set width of Halpha in order to properly correct for reddening
#Halpha_width = 40.


############################################################################################################################################

object_name = objects_list[object_number]

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
#altern = '../results/sbs1319/sbs1319_opt_corr1.txt'
#f, w = numpy.loadtxt(altern, skiprows=5, usecols=(1,2), unpack=True)
### To get altern files run: 1. correct_spec script, 2.rspectext, 3.splot, 4.correct with j and save with i, 5.wspectext
#altern = '../results/tol9/tol9_opt_corr1.txt'
#w, f = numpy.loadtxt(altern, unpack=True)
#data = [numpy.array([w,f])]

# Terminations used for the lines text files
spectrum_region = ["_nuv", "_opt", "_nir"]

for d, s in zip(data, specs):
    print 'Working with:  %s' % full_file_list[s]
    
    # Correct spectra for redshift and calculate continuum with a polynomial of nth order
    z = 0.0  # this is non-important because we are NOT correcting for redshift
    object_spectra, fitted_continuum, err_cont_fit = spectrum.fit_continuum(object_name, d, z, order=order, sigmas_away=sigmas_away, window=window, plot=plot, z_correct=False, normalize=normalize)
    wavs, fluxs = object_spectra
    _, cont_fluxs = fitted_continuum

    if text_table == True:
        # Write text file of wavelengths, fluxes, and fitted continuum
        new_file_name = object_name+spectrum_region[s]+".txt"
        file_name = os.path.join(results4object_path, new_file_name)
        txt = open(file_name, 'w+')
        print >> txt, '{:<10} {:>30} {:>35}'.format('Wavelength [A]', 'Flux [ergs/s/cm$^2$/$\AA$]', 'Continuum [ergs/s/cm$^2$/$\AA$]')
        for w, f, c in zip(wavs, fluxs, cont_fluxs):
            print >> txt, '{:<7.3f} {:>25.5e} {:>32.5e}'.format(w, f, c)
        txt.close()

    ### If needing to create a text file with only wavelengths and fluxes for splot change splot_text to true and change the corresponding part
    # of the spectra to appear in the title of the text file
    splot_text = False
    part_of_spec = 'opt'
    if splot_text == True:
        name_out_file = os.path.join(results4object_path, object_name+"_"+part_of_spec+"spec.txt")
        fout = open(name_out_file, 'w+')
        wavs, fluxs = object_spectra
        for w, f in zip(wavs, fluxs):
            fout.write('{:<4.5f} {:>20.10e}\n'.format(w, f))
        fout.close()

    
print 'Code finished!'

