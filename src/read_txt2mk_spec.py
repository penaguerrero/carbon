import numpy as np
import os
import string
#import lineid_plot
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from science import spectrum


#######################################################################################################################

'''
This script reads a 2-column text file of Angstroms and Flux to create a plot of the specified region.
If the region is not specified.
Output is 2 plots: 
            - plot from 1600-2000A = enhancement with marked OIII] and CIII] lines 
            - plot from 1600-10,000A = whole spectrum
'''

# name of the object
#                  0            1           2          3         4        5         6         7         8
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                 9           10         11         12         13       14        15         16         17
object_number = 4

# Specify the region to enhance in the form: [4600, 5100]  or leave as None to enhance 1620-2000A.
enchance = None #[4600, 5100]

# Specify the type of image to be saved
save_figs = False
img_format = '.jpg'
 
#######################################################################################################################

#### FUNCTIONS

def find_overlapping_region(w_nuv, w_opt, w_nir, f_nuv, f_opt, f_nir, region):
    # disregard the first 50 Angstroms of the opt and nir because they are too noisy
    new_w_opt = []
    new_f_opt = []
    for wi, fi in zip(w_opt, f_opt):
        if wi >= (w_opt[0]+50.0):
            new_w_opt.append(wi)
            new_f_opt.append(fi)
    new_w_nir = []
    new_f_nir = []
    for wi, fi in zip(w_nir, f_nir):
        if wi >= (w_nir[0]+50.0):
            new_w_nir.append(wi)
            new_f_nir.append(fi)
    # find the overlapping regions
    if region == 'nuv2opt':
        _, overlap1_idx = spectrum.find_nearest(w_nuv, new_w_opt[0])
        _, overlap2_idx = spectrum.find_nearest(new_w_opt, w_nuv[-1])
        overlap_wavs1 = w_nuv[overlap1_idx:-1]
        overlap_flxs1 = f_nuv[overlap1_idx:-1]
        overlap_wavs2 = new_w_opt[0:overlap2_idx]
        overlap_flxs2 = new_f_opt[0:overlap2_idx]
    if region == 'opt2nir':
        _, overlap1_idx = spectrum.find_nearest(new_w_opt, new_w_nir[0])
        _, overlap2_idx = spectrum.find_nearest(new_w_nir, new_w_opt[-1])
        overlap_wavs1 = new_w_opt[overlap1_idx:-1]
        overlap_flxs1 = new_f_opt[overlap1_idx:-1]
        overlap_wavs2 = new_w_nir[0:overlap2_idx]
        overlap_flxs2 = new_f_nir[0:overlap2_idx]
    # find the average wavelength and flux of the two spectra
    overlap_avg_wavs = []
    overlap_avg_flxs = []
    for w1, w2, f1, f2 in zip(overlap_wavs1, overlap_wavs2, overlap_flxs1, overlap_flxs2):
        wavg = (w1 + w2) / 2.0
        overlap_avg_wavs.append(wavg)
        favg = (f1 + f2) / 2.0
        overlap_avg_flxs.append(favg)
    return (overlap_avg_wavs, overlap_avg_flxs)

def find_ticks(main_ticks_repeated, minor_ticks_repeated):
    majorLocator   = MultipleLocator(main_ticks_repeated)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator   = MultipleLocator(minor_ticks_repeated)
    return majorLocator, majorFormatter, minorLocator
    
    
#### CODE 

object_name = objects_list[object_number]

# Read the text files and record the wavelengths and fluxes
regions_list = ['NUV', 'OPT', 'NIR']
wavs = []
flxs = []
for reg in regions_list:
    file2use = object_name+'_'+reg+'_smooth3.txt'
    path2txt = 'carbon/results/'+object_name+'/'+file2use
    path2txt_list = string.split(os.getcwd(), sep='carbon')
    path2txt_path = os.path.join(path2txt_list[0], path2txt)
    w, f = np.loadtxt(path2txt_path, unpack=True)
    wavs.append(w)
    flxs.append(f)

# Create the whole spectrum array
w_nuv, w_opt, w_nir = wavs
f_nuv, f_opt, f_nir = flxs
nuv2opt_overlap_avg_wavs, nuv2opt_overlap_avg_flxs = find_overlapping_region(w_nuv, w_opt, w_nir, f_nuv, f_opt, f_nir, region='nuv2opt')
opt2nir_overlap_avg_wavs, opt2nir_overlap_avg_flxs = find_overlapping_region(w_nuv, w_opt, w_nir, f_nuv, f_opt, f_nir, region='opt2nir')
# get all the wavelengths into a single array
whole_spec_wavs = []
whole_spec_flxs = []
# add all the nuv wavelengths and fluxes before the first overlapping region
for wi, fi in zip(w_nuv, f_nuv):
    if (wi > w_nuv[0]+50.0) and (wi < nuv2opt_overlap_avg_wavs[0]):
        whole_spec_wavs.append(wi)
        whole_spec_flxs.append(fi)
# add the first overlapping region
for wi, fi in zip(nuv2opt_overlap_avg_wavs, nuv2opt_overlap_avg_flxs):
    whole_spec_wavs.append(wi)
    whole_spec_flxs.append(fi)
# add the optical region before the second overlapping region
for wi, fi in zip(w_opt, f_opt):
    if (wi > nuv2opt_overlap_avg_wavs[-1]) and (wi < opt2nir_overlap_avg_wavs[0]):
        whole_spec_wavs.append(wi)
        whole_spec_flxs.append(fi)
# add the second overlapping region
for wi, fi in zip(opt2nir_overlap_avg_wavs, opt2nir_overlap_avg_flxs): 
    whole_spec_wavs.append(wi)
    whole_spec_flxs.append(fi)
# add the nir wavelengths after the second overlaping region
for wi, fi in zip(w_nir, f_nir):
    if (wi > opt2nir_overlap_avg_wavs[-1]) and (wi <= w_nir[-1]-50.0):
        whole_spec_wavs.append(wi)
        whole_spec_flxs.append(fi)
whole_spec_arr = np.array([whole_spec_wavs, whole_spec_flxs])

# Make the whole spectrum plot
full_names_list = ['IIIZw 107', 'IRAS 08339+6517', 'MRK 1087', 'MRK 1199', 'MRK 5', 'MRK 960', 
                   'NGC 1741', 'POX 4', 'SBS 0218+003', 'SBS 0948+532', 'SBS 0926+606', 'SBS 1054+365', 
                   'SBS 1319+579', 'TOL 1457-262', 'TOL 9', 'ARP 252', 'IRAS 08208+2816', 'SBS 1415+437']
full_name = full_names_list[object_number]
plt.figure(1, figsize=(15, 8))
#plt.title(object_name)
plt.xlabel('Wavelength  [$\AA$]')
plt.ylabel('Flux  [ergs/s/cm$^2$/$\AA$]')
plt.xlim(1600.0, 9700.0)
majorLocator, majorFormatter, minorLocator = find_ticks(1000, 100)
ax = plt.gca() 
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)   # for the minor ticks, use no labels; default NullFormatter
ymax = max(whole_spec_arr[1]) + max(whole_spec_arr[1])*0.12
plt.ylim(-1.5e-15, ymax)
fmax = max(whole_spec_arr[1]) - max(whole_spec_arr[1])*0.2 
plt.text(8000, fmax, full_name, fontsize=18)
plt.plot(whole_spec_arr[0], whole_spec_arr[1], 'k')
plt_name = object_name+'_wholespec'+img_format
path2plt = 'carbon/cloudy_tests/pycloudy/'+plt_name
path2plt_list = string.split(os.getcwd(), sep='carbon')
destination = os.path.join(path2plt_list[0], path2plt)
if save_figs:
    plt.savefig(os.path.abspath(destination))
    print 'Plot ', os.path.abspath(destination), 'saved!'
plt.show()


# Make the enhancement plot
plt.figure(1, figsize=(15, 8))
plt.xlabel('Wavelength  [$\AA$]')
plt.ylabel('Flux  [ergs/s/cm$^2$/$\AA$]')
default_enhancement = False
if enchance == None:
    default_enhancement = True
    enchance = [1620.0, 2020.0]
else:
    enchance = [float(enchance[0]), float(enchance[1])]
w = whole_spec_arr[0][(whole_spec_arr[0] >= enchance[0]) & (whole_spec_arr[0] <= enchance[1])]
f = whole_spec_arr[1][(whole_spec_arr[0] >= enchance[0]) & (whole_spec_arr[0] <= enchance[1])]
plt.plot(w, f, 'k')
plt.xlim(enchance[0]+20.0, enchance[1]-20.0)
ax = plt.gca() 
majorLocator, majorFormatter, minorLocator = find_ticks(50, 10)
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)
ymax = max(f) + max(f)*0.12
ymin = min(f)-min(f)*0.1
plt.ylim(ymin, ymax)
fmax = max(f) - max(f)*0.2 
if default_enhancement:
    line_label = ['O III]', 'C III]']
    lines_wave = [1662.0, 1908.5]
    for l, lid in zip(lines_wave, line_label):
        wl, idx = spectrum.find_nearest(w, l)
        ftxt = ymax*0.8
        line_length_parts = ((ymax - ymin)/2.0)/10.0
        flmin = f[idx] + 3.0*line_length_parts
        flmax = flmin + 2.0*line_length_parts    
        ftxt = flmax + line_length_parts
        #print wl, f[idx], flmin, flmax, lid, flmax+flmax*0.1
        plt.vlines(wl, flmin, flmax, colors='r', linestyles='solid')   
        plt.text(l-5.0, ftxt, lid, fontsize=14)
        #plt.ylim(ymin, ymax+ftxt/2)
    #lineid_plot.plot_line_ids(w, f, lines_wave, line_label, color='k')
else:
    plt.plot(w, f, 'k')
plt_name = object_name+'_'+str(enchance[0])+'-'+str(enchance[1])+img_format
path2plt = 'carbon/cloudy_tests/pycloudy/'+plt_name
path2plt_list = string.split(os.getcwd(), sep='carbon')
destination = os.path.join(path2plt_list[0], path2plt)
if save_figs:
    plt.savefig(os.path.abspath(destination))
    print 'Plot ', os.path.abspath(destination), 'saved!'
plt.show()


print ' Code finished!'
