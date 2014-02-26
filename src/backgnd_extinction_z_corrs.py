import numpy
import os
import pyneb as pn 
from science import spectrum

############################################################################################################################################


'''  Choose parameters to run script  '''

# name of the object
#                 0           1           2            3         4        5          6        7         8
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                 9           10         11         12         13       14        15         16         17

object_number = 12

# use this option if object is VERY faint and want to use thinner widths for emission lines
faintObj = False  

# Set width of Halpha in order to properly correct for reddening
Halpha_width = 40.

# Write the text file with line info?
text_table = True


############################################################################################################################################


# corresponding redshifts
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

# Desired Angstroms per pixel
desired_disp_list = [2.4, 8.0, 7.4]

# Object to be analyzed
object_name = objects_list[object_number]
z = z_list[object_number]

# Go into the directory of each object. Files assumed to be in: /Users/home_direcotry/Documents/AptanaStudio3/src/
results_path = "../results/"
# but just to make sure we are in the right place, get and work with full paths
full_results_path = os.path.abspath(results_path)
# Go into the object folder
results4object_path = os.path.join(full_results_path, object_name)

# Gather the nuv, opt, and nir files into a single list of wavelengths, fluxes, and continuum fluxes
specs = [0, 1, 2]
object_file = os.path.join(results4object_path, object_name)
text_file_list, _ = spectrum.get_obj_files2use(object_file, specs)
wavs_nuv = []
flxs_nuv = []
cont_nuv = []
wavs_opt = []
flxs_opt = []
cont_opt = []
wavs_nir = []
flxs_nir = []
cont_nir = []
continuum_flxs = []
for s, tf in zip(specs, text_file_list): 
    wavs, flxs, cont_flxs = numpy.loadtxt(tf, skiprows=1, unpack=True)
    for w, f, c in zip(wavs, flxs, cont_flxs):
        continuum_flxs.append(c)
        if s == 0:
            wavs_nuv.append(w)
            flxs_nuv.append(f)
            cont_nuv.append(c)
        elif s == 1:
            wavs_opt.append(w)
            flxs_opt.append(f)
            cont_opt.append(c)
        elif s == 2:
            wavs_nir.append(w)
            flxs_nir.append(f)
            cont_nir.append(c)

# Check if the continuum becomes negative and if it does add a constant to the fluxes so that it dosen't
'''
if min(continuum_flxs) < 0.0:
    constant = numpy.fabs(min(continuum_flxs)) + numpy.fabs(min(continuum_flxs))*0.1  
    print 'Adding a constant to the flux so that the continuum does not become negative. Constant = ', constant
    for f, c in zip(flxs_nuv, cont_nuv):
        new_f = f + constant
        new_c = c + constant
        idx = flxs_nuv.index(f)
        flxs_nuv[idx] = new_f
        cont_nuv[idx] = new_c
    for f, c in zip(flxs_opt, cont_opt):
        new_f = f + constant
        new_c = c + constant
        idx = flxs_opt.index(f)
        flxs_opt[idx] = new_f
        cont_opt[idx] = new_c
    for f, c in zip(flxs_nir, cont_nir):
        new_f = f + constant
        new_c = c + constant
        print c, new_c
        idx = flxs_nir.index(f)
        flxs_nir[idx] = new_f
        cont_nir[idx] = new_c
'''
# Gather the info into a list of numpy arrays
nuv = numpy.array([wavs_nuv, flxs_nuv])
opt = numpy.array([wavs_opt, flxs_opt])
nir = numpy.array([wavs_nir, flxs_nir])
data = []
data.append(nuv)
data.append(opt)
data.append(nir)
nuv_cont = numpy.array([wavs_nuv, cont_nuv])
opt_cont = numpy.array([wavs_opt, cont_opt])
nir_cont = numpy.array([wavs_nir, cont_nir])
cont_data = []
cont_data.append(nuv_cont)
cont_data.append(opt_cont)
cont_data.append(nir_cont)

# Terminations used for the lines text files
spectrum_region = ["_nuv", "_opt", "_nir"]

'''The STIS data handbook gives a dispersion of 0.60, 1.58, 2.73, and 4.92 Angstroms per pixel for grating settings G140L, 
G230L, G430L, and G750L, respectively. The resolution element is ~1.5 pixels. '''
originals = [1.58, 2.73, 4.92]

for d, cd, s in zip(data, cont_data, specs):
    # Rebin the spectra to the corresponding dispersion
    desired_dispersion = desired_disp_list[s]
    #rebinned_arr = spectrum.rebin_spec2disp(desired_dispersion, d)
    #rebinned_cont = spectrum.rebin_spec2disp(desired_dispersion, cd)
    rebinned_arr = spectrum.rebin2AperPix(originals[s], desired_disp_list[s], d)
    rebinned_cont = spectrum.rebin2AperPix(originals[s], desired_disp_list[s], cd)
    
    #print '    *** Wavelengths corrected for redshift.'
    w_corr = rebinned_arr[0] / (1+float(z))
    # rebinned and z-corrected data
    object_spectra = numpy.array([w_corr, rebinned_arr[1]]) 
    contum_spectra = numpy.array([w_corr, rebinned_cont[1]]) 
    
    # Obtain the lines net fluxes and EWs
    new_file_name = object_name+"_lineinfo"+spectrum_region[s]+".txt"
    lineinfo_text_file = os.path.join(results4object_path, new_file_name)
    # Now obtain the continuum and equivalent widths
    object_lines_info = spectrum.find_lines_info(object_spectra, contum_spectra, lineinfo_text_file, Halpha_width=Halpha_width, text_table=text_table, vacuum=False, faintObj=faintObj)
    print 'This are the lines in the ', spectrum_region[s]
    print ''
    #raw_input('    press enter to continue...')
    
    # Determine the corresponding E(B-V) value for each object
    ebv = A_B_list[object_number] - A_V_list[object_number]
    

print 'Code finished!'

