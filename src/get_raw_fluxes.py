import numpy
import os
import pyneb as pn 
from science import spectrum

'''
This code finds the net fluxes (or raw intensities) and equivalent widths from the redshift-corrected spectra of the given object.
It writes a text file containing the line information it found for that wavelength range.
    * The possible lines the code looks for are all those in the line_catalog.txt file in the spectrum folder
'''

############################################################################################################################################

''' OBJECTS OF THE SAMPLE '''

# name of the object
#                 0           1           2            3         4        5          6        7         8
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                 9           10         11         12         13       14        15         16         17

# corresponding redshift
# from observations
#             0        1        2        3         4         5        6        7         8
z_list = [0.01985, 0.019113, 0.02877, 0.01354, 0.002695, 0.021371, 0.01348, 0.01201, 0.05842,
          0.046240, 0.013642, 0.002010, 0.0076, 0.01763, 0.010641, 0.032989, 0.04678, 0.002031]
#             9       10        11       12        13       14       15       16       17
# tol9=14 is found to have a much different z than that of phase2.... previous value = 0.01195
# taken from phase 2 document

# A_B:  Background reddening correction, taken from NED
#              0     1      2       3      4     5       6    7       8
A_B_list = [0.259, 0.397, 0.272, 0.231, 0.364, 0.101, 0.221, 0.170, 0.155, 
            0.054, 0.135, 0.089, 0.061, 0.680, 0.281, 0.210, 0.137, 0.039]
#             9    10      11      12     13    14      15    16     17
# A_V:  Background reddening correction, taken from NED
#              0     1      2       3      4     5       6    7       8
A_V_list = [0.199, 0.305, 0.209, 0.178, 0.279, 0.077, 0.169, 0.130, 0.119, 
            0.042, 0.104, 0.068, 0.047, 0.522, 0.216, 0.161, 0.105, 0.030]
#             9    10      11      12     13    14      15    16     17

############################################################################################################################################

# Choose parameters to run script

# 1) Select a number from objects_list, i = :
object_number = 12
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
order = 1

# 6) What is the width of the window to use to find local continuum?
window = 550

# 7) Do you want to see the plots of the fitted continuum?
plot = True

# 8) write the text file with the line net fluxes and equivalent widths?
text_table = False

# Set width of Halpha in order to properly correct for reddening
Halpha_width = 35.


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
#altern = '../results/sbs1319/sbs1319_opt_corr1.txt'
#f, w = numpy.loadtxt(altern, skiprows=5, usecols=(1,2), unpack=True)
### To get altern files run: 1. correct_spec script, 2.rspectext, 3.splot, 4.correct with j and save with i, 5.wspectext
#altern = '../results/tol9/tol9_opt_corr1.txt'
#w, f = numpy.loadtxt(altern, unpack=True)
#data = [numpy.array([w,f])]

# Terminations used for the lines text files
spectrum_region = ["_nuv", "_opt", "_nir"]

'''
The STIS data handbook gives a dispersion of 0.60, 1.58, 2.73, and 4.92 Angstroms per pixel for grating settings G140L, 
G230L, G430L, and G750L, respectively. The resolution element is ~1.5 pixels.
'''
# Rebin the spectra to the following dispersions for grating settings G140L, G230L, G430L, and G750L, respectively
# Since I am not using G140L, I will remove it from desired_disp_list = [0.9, 2.4, 4.1, 7.4]
desired_disp_list = [2.4, 4.1, 7.4]

#######
''' FUNCTION TAKEN FROM PYNEB '''
def CCM89(wave, R_V):
    """
    Cardelli 1989
    """
    x = 1e4 / numpy.asarray([wave]) # inv microns
    a = numpy.zeros_like(x)
    b = numpy.zeros_like(x)
    
    tt = (x > 0.3) & (x <= 1.1)
    a[tt] = 0.574 * x[tt] ** 1.61 
    b[tt] = -0.527 * x[tt] ** 1.61

    tt = (x > 1.1) & (x <= 3.3)
    yg = x[tt] - 1.82
    a[tt] = (1. + 0.17699 * yg - 0.50447 * yg ** 2. - 0.02427 * yg ** 3. + 0.72085 * yg ** 4. + 
             0.01979 * yg ** 5. - 0.7753 * yg ** 6. + 0.32999 * yg ** 7.)
    b[tt] = (0. + 1.41338 * yg + 2.28305 * yg ** 2. + 1.07233 * yg ** 3. - 5.38434 * yg ** 4. - 
             0.622510 * yg ** 5. + 5.3026 * yg ** 6. - 2.09002 * yg ** 7.)
    
    tt = (x > 3.3) & (x <= 5.9)
    a[tt] = 1.752 - 0.316 * x[tt] - 0.104 / ((x[tt] - 4.67) ** 2. + 0.341)
    b[tt] = -3.090 + 1.825 * x[tt] + 1.206 / ((x[tt] - 4.62) ** 2 + 0.263)
    
    tt = (x > 5.9) & (x <= 8.0)
    a[tt] = (1.752 - 0.316 * x[tt] - 0.104 / ((x[tt] - 4.67) ** 2. + 0.341) - 
             0.04473 * (x[tt] - 5.9) ** 2. - 0.009779 * (x[tt] - 5.9) ** 3.)
    b[tt] = (-3.090 + 1.825 * x[tt] + 1.206 / ((x[tt] - 4.62) ** 2. + 0.263) + 
             0.2130 * (x[tt] - 5.9) ** 2. + 0.1207 * (x[tt] - 5.9) ** 3.)
    
    tt = (x > 8.0) & (x < 10.0)
    a[tt] = (-1.073 - 0.628 * (x[tt] - 8) + 0.137 * (x[tt] - 8) ** 2. - 
             0.070 * (x[tt] - 8) ** 3.)
    b[tt] = (13.670 + 4.257 * (x[tt] - 8) - 0.420 * (x[tt] - 8) ** 2. + 
             0.374 * (x[tt] - 8) ** 3.)
    
    Xx = R_V * a + b
    return numpy.squeeze(Xx)
#######


for d, s in zip(data, specs):
    print 'Working with:  %s' % full_file_list[s]
    
    # Rebin the spectra to the corresponding dispersion
    desired_dispersion = desired_disp_list[s]
    rebinned_arr = spectrum.rebin_spec2disp(desired_dispersion, d)
    
    # Determine the corresponding E(B-V) value for each object
    ebv = A_B_list[object_number] - A_V_list[object_number]
    
    # Correct for background reddening
    bacground_corr_flux = []
    wavenelgths, fluxes = d
    Rv = A_V_list[object_number] / ebv
    for wav, flx in zip(wavenelgths, fluxes):
        Xx = CCM89(wav, R_V=Rv)
        corr_flx = Xx * flx
        print wav, flx, corr_flx
        bacground_corr_flux.append(corr_flx)
    d = numpy.array([d[0], bacground_corr_flux])
    
    # Correct spectra for redshift and calculate continuum with a polynomial of nth order
    object_spectra, fitted_continuum, err_cont_fit = spectrum.fit_continuum(object_name, d, z, order=order, sigmas_away=sigmas_away, window=window, plot=plot, z_correct=True, normalize=normalize)

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

    # Obtain the lines net fluxes and EWs
    new_file_name = object_name+"_lineinfo"+spectrum_region[s]+".txt"
    lineinfo_text_file = os.path.join(results4object_path, new_file_name)
    # Now obtain the continuum and equivalent widths
    faintObj = True  # --> use this option if object is VERY faint and want to use thiner widths for emission lines
    object_lines_info = spectrum.find_lines_info(object_spectra, fitted_continuum, err_cont_fit, lineinfo_text_file, Halpha_width=Halpha_width, text_table=text_table, vacuum=vacuum, faintObj=faintObj)
    print ''
    
c = 3.0e5 #km/s
velocity = c * z
print 'v = c * z = 3e5 * %0.5f = %0.3f' % (z, velocity)
H0 = 75#67.8 #+-0.77 km/sec/Mpc (from Planck mission)
distance = velocity / H0
print 'Estimated distance assuming H0 = %f:   d[Mpc] = %0.3f' % (H0, distance)
print sigmas_away, 'sigmas_away'    

print 'Code finished!'
    
    
