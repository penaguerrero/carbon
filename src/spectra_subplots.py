import numpy as np
import os
import matplotlib
from matplotlib import pyplot as plt
import subplots_dictionary as spd
from science import spectrum

'''
This scripts creates sub-plots for the following regions, identifying important lines:
    UV    -> 1500-2100: O III] and C III] lines
    blue  -> 3700-5100: [O II], [O III] (4363, 4959, 5007), Hd, Hg, Hb, [Ne III], [S II])
    red   -> 6250 - 6750: [O I], [S III], Ha+[N II], He I, [S II])
    IR    -> 9000 - > 9600: [SIII]
'''

#                  0            1           2          3         4        5         
objects_list = ['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 
#                  6          7         8         9          10         11          
                'ngc1741', 'pox4', 'sbs0218', 'sbs0948', 'sbs0926', 'sbs1054', 
#                  12         13       14        15         16         17
                'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']

# Set script parameters
object_number = 8       # object to create sub-plots for
save_fig = False         # Specify the type of image to be saved
img_format = '.jpg'     # format for resulting tile of sub-plots


#################################################################################################


# corresponding redshifts from observations
#             0        1        2        3         4     50.021371     6        7         8
z_list = [0.01985, 0.019581, 0.02813, 0.01354, 0.002695, 0.02346, 0.013631, 0.01201, 0.05842,
          0.046240, 0.013642, 0.002010, 0.0073, 0.01763, 0.01199, 0.03207, 0.04678, 0.002031]
#             9       10        11       12        13       14       15       16       17

#                       0      1     2     3     4      5     6     7     8
perform_rebin_list = [True, False, True, True, False, True, True, True, True, 
                      False, False, True, True, True, True, True, True, True]
#                      9      10    11    12     13    14    15     16     17
# Desired Angstroms per pixel
# Mrk 1087: 2.0, 8.0, 8.0
'''The STIS data handbook gives a dispersion of 0.60, 1.58, 2.73, and 4.92 Angstroms per pixel for grating settings G140L, 
G230L, G430L, and G750L, respectively. The resolution element is ~1.5 pixels. '''
originals = [1.58, 2.73, 4.92]
or1, or2, or3 = originals[0], originals[1], originals[2]
#                                    0                         1         2:or1*1.9, or2*1.85, or3*1.5           3                           4                5
desired_disp_listoflists = [[or1*2.0, or2*2.0, or3*2.0], [1.6, 6.5, 5.0], [or1*1.2, or2*1.51, or3*1.2], [or1*2.1, or2*2.1, or3*2.1], [2.0, 3.0, 5.0], [1.7, 5.6, 9.8], 
                            #        6              7                 8                9               10               11
                            [or1*2.0, or2*2.0, or3*2.0], [1.8, 4.0, 6.0], [or1, or2, or3], [1.6, or2*2., 5.0], [1.6, 3.5, 5.6], [1.6, 3.9, 5.0],
#                                   12             13                14               15               16               17                            
                            [1.6, 5.0, 5.0], [1.6, 4.5, 5.0], [1.6, 3.0, 5.0], [1.7, 5.6, or3], [1.6, 3.1, 5.1], [1.6, 8.3, 11.0]]


# Object to be analyzed
object_name = objects_list[object_number]
z = z_list[object_number]

# Go into the directory of each object. Files assumed to be in: /Users/home_direcotry/Documents/AptanaStudio3/src/
results_path = "../results/"
# but just to make sure we are in the right place, get and work with full paths
full_results_path = os.path.abspath(results_path)
# Go into the object folder
results4object_path = os.path.join(full_results_path, object_name)
perform_rebin = perform_rebin_list[object_number]
desired_disp_list = desired_disp_listoflists[object_number]

# Gather the nuv, opt, and nir files into a single list of wavelengths, fluxes, and continuum fluxes
specs = [0, 1, 2]
object_file = os.path.join(results4object_path, object_name)
text_file_list, _ = spectrum.get_obj_files2use(object_file, specs)
wavs_nuv, flxs_nuv, cont_nuv = [], [], []
wavs_opt, flxs_opt, cont_opt = [], [], []
wavs_nir, flxs_nir, cont_nir = [], [], []
continuum_flxs, all_err_cont_fit = [], []
for s, tf in zip(specs, text_file_list): 
    wavs, flxs, cont_flxs = np.loadtxt(tf, skiprows=2, unpack=True)
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
    # get error of the continuum
    err_cont_fit = spectrum.get_err_cont_fit(tf)
    all_err_cont_fit.append(err_cont_fit)    
# Gather the info into a list of np arrays
nuv = np.array([wavs_nuv, flxs_nuv])
opt = np.array([wavs_opt, flxs_opt])
nir = np.array([wavs_nir, flxs_nir])
data = []
data.append(nuv)
data.append(opt)
data.append(nir)
nuv_cont = np.array([wavs_nuv, cont_nuv])
opt_cont = np.array([wavs_opt, cont_opt])
nir_cont = np.array([wavs_nir, cont_nir])
cont_data = []
cont_data.append(nuv_cont)
cont_data.append(opt_cont)
cont_data.append(nir_cont)
object_wavelengths, object_fluxes, object_continuum = [], [], []
for d, cd, s in zip(data, cont_data, specs):
    # Rebin the spectra to the corresponding dispersion
    desired_dispersion = desired_disp_list[s]
    if perform_rebin:
        rebinned_arr = spectrum.rebin2AperPix(originals[s], desired_disp_list[s], d)
        rebinned_cont = spectrum.rebin2AperPix(originals[s], desired_disp_list[s], cd)
    else:
        rebinned_arr = d
        rebinned_cont = cd
    #  Wavelengths corrected for redshift
    w_corr = rebinned_arr[0] / (1+float(z))
    if s == 0:
        corrected_nuv_wav, corrected_nuv_flx = w_corr, rebinned_arr[1]
        corr_cont_nuv_flx = rebinned_cont[1]
    if s == 1:
        corrected_opt_wav, corrected_opt_flx = w_corr, rebinned_arr[1]
        corr_cont_opt_flx = rebinned_cont[1]
    if s == 2:
        corrected_nir_wav, corrected_nir_flx = w_corr, rebinned_arr[1]
        corr_cont_nir_flx = rebinned_cont[1]

full_names_list = ['IIIZw 107', 'IRAS 08339+6517', 'MRK 1087', 'MRK 1199', 'MRK 5', 'MRK 960', 
                'NGC 1741', 'POX 4', 'SBS 0218+003', 'SBS 0948+532', 'SBS 0926+606', 'SBS 1054+365',  
                'SBS 1319+579', 'TOL 1457-262', 'TOL 9', 'ARP 252', 'IRAS 08208+2816', 'SBS 1415+437']

def fill_list(list_to_check, reference_list):
    if len(list_to_check) == 1:
        for _ in reference_list:
            list_to_check.append(list_to_check[0])

def create_subplt(plt_data, lineymin=None):
    '''
    This function parses the specific conditions to create a sub-plot.
    '''
    # unfold data
    corrected_wav, corrected_flx, plt_number, xlims, lines, legends, fmin, subxcoord_list, subycoord_list, side_list, lineymax_list, idx_list = plt_data
    plt.subplot(2, 2, plt_number)
    xarr = corrected_wav[(corrected_wav >= xlims[0]) & (corrected_wav <= xlims[1])]
    yarr = corrected_flx[(corrected_wav >= xlims[0]) & (corrected_wav <= xlims[1])]
    plt.xlim(xlims[0], xlims[1])
    fmax = max(yarr)
    plt.ylim(fmin, fmax+fmax*0.35)
    plt.plot(xarr, yarr, 'k')#, corrected_nuv_wav, corr_cont_nuv_flx, 'r--')
    # make sure the list are the same length
    fill_list(subxcoord_list, lines)
    fill_list(subycoord_list, lines)
    fill_list(idx_list, lines)
    fill_list(side_list, lines)
    fill_list(lineymax_list, lines)   
    # Annotate the points 5 _points_ above and to the left of the vertex
    for line, leg, sl, sx, sy, lineymax, i in zip(lines, legends, side_list, subxcoord_list, 
                                                      subycoord_list, lineymax_list, idx_list):
        _, idx = spectrum.find_nearest(xarr, line)
        x, y = xarr[idx+i], yarr[idx+i]
        if lineymin is None:
            lineymin = 0.1
        subxcoord, subycoord, side = sx, sy, sl
        if line == 6563:
            leg, line = '[N II] + Ha', ''
        vlinex, vliney = [x, x], [y+y*lineymin, y+y*lineymin+lineymax]
        plt.plot(vlinex, vliney, 'r')
        plt.annotate('{} {}'.format(leg, line), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')

# Make the plots
font = {#'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}
matplotlib.rc('font', **font)
fig, axs = plt.subplots(2, 2, figsize=(12,10))
fig.suptitle(full_names_list[object_number])
fig.subplots_adjust(hspace=0.30)
xlab = 'Wavelength [$\AA$]'
ylab = 'Flux [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]'
# Set common labels
fig.text(0.5, 0.03, xlab, ha='center', va='center')
fig.text(0.05, 0.5, ylab, ha='center', va='center', rotation='vertical')
plt.minorticks_on()

regions2plots = ['NUV', 'Blue', 'Red', 'NIR']
for reg in regions2plots:
    # get data from dictionary
    plt_number = spd.galaxies_subplt_dict[object_name][reg]['plt_number']
    xlims = spd.galaxies_subplt_dict[object_name][reg]['xlims']
    lines = spd.galaxies_subplt_dict[object_name][reg]['lines']
    legends = spd.galaxies_subplt_dict[object_name][reg]['legends']
    fmin = spd.galaxies_subplt_dict[object_name][reg]['fmin']
    subxcoord_list = spd.galaxies_subplt_dict[object_name][reg]['subxcoord_list']
    subycoord_list = spd.galaxies_subplt_dict[object_name][reg]['subycoord_list']
    side_list = spd.galaxies_subplt_dict[object_name][reg]['side_list']
    lineymax_list = spd.galaxies_subplt_dict[object_name][reg]['lineymax_list']
    idx_list = spd.galaxies_subplt_dict[object_name][reg]['idx_list']
    # use appropriate wavelength and flux arrays
    if reg == 'NUV':
        corrected_wav, corrected_flx = corrected_nuv_wav, corrected_nuv_flx
    if reg == 'Blue':
        corrected_wav, corrected_flx = corrected_opt_wav, corrected_opt_flx
    if (reg == 'Red') or (reg == 'NIR'):
        corrected_wav, corrected_flx = corrected_nir_wav, corrected_nir_flx
    # compress data
    plt_data = [corrected_wav, corrected_flx, plt_number, xlims, lines, legends, fmin, 
                 subxcoord_list, subycoord_list, side_list, lineymax_list, idx_list]
    # create subplots
    create_subplt(plt_data, lineymin=0.05)

if save_fig:
    subplot_destination = os.path.join(results_path, 'plots/'+object_name+'_suplots.jpg')
    fig.savefig(subplot_destination)
    print 'Figure ', subplot_destination, ' saved!'

plt.show()

