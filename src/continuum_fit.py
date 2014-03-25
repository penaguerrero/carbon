import numpy
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
plot = False

# 8) write the text file with the line wavelengths, fluxes, and fitted continuum?
text_table = True

# Want to see the quiasi-final spectrum?  (i.e. correct for redshift and rebin)
correct_redshift = False
rebin = False

############################################################################################################################################

# Set width of Halpha in order to properly correct for reddening
#Halpha_width = 40.

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
#altern = '../results/tol9/tol9_opt21_fix.txt'
#w, f = numpy.loadtxt(altern, unpack=True)
#data = [numpy.array([w,f])]

# Terminations used for the lines text files
spectrum_region = ["_nuv", "_opt", "_nir"]

z_list = [0.01985, 0.019581, 0.02813, 0.01354, 0.002695, 0.02346, 0.013631, 0.01201, 0.05842,
          0.046240, 0.013642, 0.002010, 0.0076, 0.01732, 0.01199, 0.032989, 0.04678, 0.002031]
# original phase 2 z for iras08339 0.019113, mrk960 0.021371, ngc1741 0.013631, tol1457 0.01763, tol9 0.010641

'''The STIS data handbook gives a dispersion of 0.60, 1.58, 2.73, and 4.92 Angstroms per pixel for grating settings G140L, 
G230L, G430L, and G750L, respectively. The resolution element is ~1.5 pixels. '''
originals = [1.58, 2.73, 4.92]
desired_disp_list = [2.0, 4.0, 5.0]

for d, s in zip(data, specs):
    print 'Working with:  %s' % full_file_list[s]
    
    if rebin == True:
        # This mode is just to show how much the spectrum will be rebinned
        d = spectrum.rebin2AperPix(originals[s], desired_disp_list[s], d)
        text_table = False

    # Correct spectra for redshift and calculate continuum with a polynomial of nth order
    if correct_redshift == True:
        z = z_list[object_number]
        text_table = False
    else:
        z = 0.0  # this is non-important because we are NOT correcting for redshift
    object_spectra, fitted_continuum, err_cont_fit = spectrum.fit_continuum(object_name, d, z, order=order, sigmas_away=sigmas_away, window=window, plot=plot, z_correct=correct_redshift, normalize=normalize)
    wavs, fluxs = object_spectra
    _, cont_fluxs = fitted_continuum

    if text_table == True:
        # Write text file of wavelengths, fluxes, and fitted continuum
        if normalize == False:
            new_file_name = object_name+spectrum_region[s]+".txt"
        else:
            new_file_name = object_name+spectrum_region[s]+"_normalized.txt"
        file_name = os.path.join(results4object_path, new_file_name)
        txt = open(file_name, 'w+')
        print >> txt, 'Percentage Error in continuum fitting  =', err_cont_fit
        print >> txt, '{:<10} {:>30} {:>35}'.format('Wavelength [A]', 'Flux [ergs/s/cm$^2$/$\AA$]', 'Continuum [ergs/s/cm$^2$/$\AA$]')
        for w, f, c in zip(wavs, fluxs, cont_fluxs):
            print >> txt, '{:<7.3f} {:>25.5e} {:>32.5e}'.format(w, f, c)
        txt.close()

    ### If needing to create a text file with only wavelengths and fluxes for splot change splot_text to true and change the corresponding part
    # of the spectra to appear in the title of the text file
    splot_text = False
    part_of_spec = 'nir'
    if splot_text == True:
        name_out_file = os.path.join(results4object_path, object_name+"_"+part_of_spec+"spec.txt")
        fout = open(name_out_file, 'w+')
        wavs, fluxs = object_spectra
        for w, f in zip(wavs, fluxs):
            fout.write('{:<4.5f} {:>20.10e}\n'.format(w, f))
        fout.close()

    
c = 3.0e5 #km/s
z = z_list[object_number]
velocity = c * z
print 'v = c * z = 3e5 * %0.5f = %0.3f' % (z, velocity)
H0 = 75#67.8 #+-0.77 km/sec/Mpc (from Planck mission)
distance = velocity / H0
print 'Estimated distance assuming H0 = %f:   d[Mpc] = %0.3f' % (H0, distance)
    
    
print 'Code finished!'

