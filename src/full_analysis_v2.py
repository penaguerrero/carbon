import numpy as np
import os
import metallicity
from science import spectrum
#import pyneb as pn 
#from glob import glob

"""
This script corrects for reddening assuming Case B, according to the ratio of
either 6563/4861=2.86 or 4340/4861=0.47
* NOTE: To run, this script requires PyNeb to be installed (Luridiana, V.,
        Morisset, C., & Shaw, R. A. 2015, A&A, 573, 42 
        http://www.iac.es/proyecto/PyNeb/)
"""

############################################################################################################################################


'''  Choose parameters to run script  '''

# name of the object
#                  0            1           2          3         4        5         6         7         8
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                 9           10         11         12         13       14        15         16         17
object_number = 2

# Is this a TEST run?
TEST_run = False

# If only wanting to perform the reddening and redshift correction set to True: first round of corrections
first_redcorr = False

# Do you want to use C_Hbeta to correct for extinction?   (if set to false the values of A_V and A_B will be used)
use_Chbeta = True

# Choose case
case = 'B'


############################################################################################################################################


### CODE

# Write the text file with line info?
create_txt_lineinfo = False
# Write the text file with temperatures, densities, and abundances?
create_txt_temdenabunds = False


# Skip the finding of the line info and go directly to gather the spectra?
# this line tells which set of line info did I take
#                            0     1     2     3      4       5      6      7     8
use_given_lineinfo_list = [False, True, False, False, False, False, False, False, False, 
                           True, False, False, False, False, True, False, False, False]
#                            9     10     11     12    13     14     15    16    17
use_given_lineinfo = use_given_lineinfo_list[object_number]

# In case of wanting to use a specific temperature and/or density (accepts values with errors as lists)
#                       0                1               2               3              4             5               6            7        
forceTeO3_list = [[10900.,11900.], [8400.,9600.], [7100.,1000.], [6850.,7950.], [12500.,13800.], [9600.,11000.], [9600.,1000.], [12500.,13500.],
                  [12600.,13600], None, [13700,14300.], None, [13200.,14200.], [14300., 15300.], [7600.,8600.], [8450.,9500.], [9700., 10500.], [15500.,16600.]]
#                       8          9          10         11         12              13               14             15               16           17
forceTeO3 = forceTeO3_list[object_number]

#                        0               1         2     3     4     5     6     7     8           
forceTeO2_list = [[10500.,11500.], [8000.,9900.], None, None, None, None, None, None, None,
                  None, None, None, None, None, [8300.,9000.], [8900.,9900.], None, [13850.,14250.]]
#                   9         10          11         12          13        14             15         16        17 
forceTeO2 = forceTeO2_list[object_number]

#                   0           1             2              3        4         5             6           7           8           
forceNe_list = [[200,900], [650., 1000.], [120., 3000.], [100, 2100], [375, 1000], [100., 3000.], [100,300], [600.,180.], [140., 80],
                [250.,80],  [100,200.], None, None, [400., 1000], [180.,60.], [300,100], 100., [100, 200.]]
#                   9          10        11    12           13       14         15      16         17 
forceNe = forceNe_list[object_number]

# use this option if object is VERY faint and want to use thinner widths for emission lines
#                  0      1     2     3      4     5     6      7      8
faintObj_list = [False, False, True, False, True, True, False, True, True, 
                 False, False, False, False, True, True, True, False, True]
#                  9      10    11    12     13    14     15     16     17
faintObj = faintObj_list[object_number]

# use this option if wanting to rebin or not
#                       0      1     2     3     4      5     6     7     8
perform_rebin_list = [True, False, True, True, False, True, True, True, True, 
                      False, False, True, True, True, True, True, True, True]
#                      9      10    11    12     13    14    15     16     17
perform_rebin = perform_rebin_list[object_number]

# Do deblend of lines?
deblend4363 = False
deblend3727 = False
deblend6563 = False

# Used values of Halpha_width in order to properly correct for reddening
#                     0    1    2    3    4    5    6    7   8     9   10   11   12   13   14   15   16   17
Halpha_width_list = [40., 28., 26., 25., 27., 28., 30., 40., 35., 27., 27., 30., 43., 28., 30., 40., 50., 35.]
Halpha_width = Halpha_width_list[object_number]

# Set theoretical ratios
I_theo_HaHb = 2.86   # Halpha/Hbeta
#I_theoHgHb = 0.47   # Hgamma/Hbeta

# Found values of EWabsHbeta and C_Hbeta in case the E(B-V) and Rv values are not known
#                   02.86            1           2            3        42.68          5*           
combos_list = [[2.0, 0.64], [2.0, 0.22], [2.0, 0.14], [2.5, 0.35], [1.0, 0.19], [0.5, 0.03], 
#                   6             7            8        90.5, 0.3      10*          11
               [2.0, 0.1], [2.0, 0.22], [2.0, 0.38], [0.5, 0.57], [0.07, 0.1], [0.8, 0.01], 
#                    12          13          140.63      15 0.54          16          17 
               [0.001, 0.8], [1.4, 0.57], [7.5, 0.58], [2.0, 0.35], [3.4, 0.13], [1.00, 0.0001]]
combo = combos_list[object_number]
# Set initial value of EWabsHbeta (this is a guessed value taken from HII regions)
# for HII region type objects typical values are 2.0-4.0 
EWabsHbeta = combo[0]
# Set value for extinction
# for HII region type objects there is no restriction to max but values MUST be positive
C_Hbeta = combo[1]

# Desired Angstroms per pixel
# Mrk 1087: 2.0, 8.0, 8.0
'''The STIS data handbook gives a dispersion of 0.60, 1.58, 2.73, and 4.92 Angstroms per pixel for grating settings G140L, 
G230L, G430L, and G750L, respectively. The resolution element is ~1.5 pixels. '''
originals = [1.58, 2.73, 4.92]
or1 = originals[0]
or2 = originals[1]
or3 = originals[2]
#                                    0                         1         2:or1*1.9, or2*1.85, or3*1.5           3                           4                5
desired_disp_listoflists = [[or1*2.0, or2*2.0, or3*2.0], [1.6, 6.5, 5.0], [or1*1.2, or2*1.51, or3*1.2], [or1*2.1, or2*2.1, or3*2.1], [2.0, 3.0, 5.0], [1.7, 5.6, 9.8], 
                            #        6              7                 8                9               10               11
                            [or1*2.0, or2*2.0, or3*2.0], [1.8, 4.0, 6.0], [or1, or2, or3], [1.6, or2*2., 5.0], [1.6, 3.5, 5.6], [1.6, 3.9, 5.0],
#                                   12             13                14               15               16               17                            
                            [1.6, 5.0, 5.0], [1.6, 4.5, 5.0], [1.6, 3.0, 5.0], [1.7, 5.6, or3], [1.6, 3.1, 5.1], [1.6, 8.3, 11.0]]
desired_disp_list = desired_disp_listoflists[object_number]

# corresponding redshifts from observations
#             0        1        2        3         4     50.021371     6        7         8
z_list = [0.01985, 0.019581, 0.02813, 0.01354, 0.002695, 0.02346, 0.013631, 0.01201, 0.05842,
          0.046240, 0.013642, 0.002010, 0.0073, 0.01763, 0.01199, 0.03207, 0.04678, 0.002031]
#             9       10        11       12        13       14       15       16       17
# original phase 2 z for iras08339 0.019113, mrk960 0.021371, ngc1741 0.013631, tol1457 0.01763, tol9 0.010641
# tol9=14 is found to have a much different z than that of phase2.... previous value = 0.01195
# taken from phase 2 document

# A_B:  Background reddening correction, taken from NED
#              0      1      2      3      4      5      6      7      8
A_B_list = [0.259, 0.397, 0.272, 0.231, 0.364, 0.101, 0.221, 0.170, 0.155, 
            0.054, 0.135, 0.089, 0.061, 0.680, 0.281, 0.210, 0.137, 0.039]
#              9     10     11     12     13     14     15     16     17
# A_V:  Background reddening correction, taken from NED
#              0      1      2      3      4      5      6      7      8
A_V_list = [0.199, 0.305, 0.209, 0.178, 0.279, 0.077, 0.169, 0.130, 0.119, 
            0.042, 0.104, 0.068, 0.047, 0.522, 0.216, 0.161, 0.105, 0.030]
#              9     10     11     12     13     14     15     16     17

'''
# Obtained Av values from PyNeb
>>> l=[0.59, 0.20, 0.13, 0.32, 0.17, 0.03, 0.09, 0.20, 0.35, 0.53, 0.09, 0.01, 0.74, 0.53, 0.54, 0.32, 0.13, 0.00]
>>> import numpy as np
>>> av_arr = np.array(l)
>>> av_avg = sum(l)/float(len(l))
>>> print av_avg
0.276111111111
'''

# Object to be analyzed
object_name = objects_list[object_number]
z = z_list[object_number]

# Go into the directory of each object. Files assumed to be in: /Users/home_direcotry/Documents/AptanaStudio3/src/
results_path = "../results/"
if TEST_run:
    results_path = "../results/TESTS/"
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
all_err_cont_fit =[]
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

# Terminations used for the lines text files
spectrum_region = ["_nuv", "_opt", "_nir"]

# This is the error in the flux obtained by STIS, according to the Data Handbook the typical accuracies for 
# spectroscopic l mode is about 2%, 5%, and 5% for NUV, opt, and NIR - tables 4.1 and 4.2
err_stis_list = [0.02, 0.05, 0.05]

def make_lineinfo_file(object_spectra, contum_spectra, Halpha_width, text_table, vacuum, faintObj, linesinfo_file_name, do_errs):
    object_lines_info = spectrum.find_lines_info(object_spectra, contum_spectra, Halpha_width, text_table, 
                                             vacuum, faintObj, linesinfo_file_name, do_errs)
    return object_lines_info

object_wavelengths = []
object_fluxes = []
object_continuum = []

if use_given_lineinfo == False:
    for d, cd, s in zip(data, cont_data, specs):
        # Rebin the spectra to the corresponding dispersion
        desired_dispersion = desired_disp_list[s]
        ##rebinned_arr = spectrum.rebin_spec2disp(desired_dispersion, d)
        ##rebinned_cont = spectrum.rebin_spec2disp(desired_dispersion, cd)
        if perform_rebin:
            rebinned_arr = spectrum.rebin2AperPix(originals[s], desired_disp_list[s], d)
            rebinned_cont = spectrum.rebin2AperPix(originals[s], desired_disp_list[s], cd)
        else:
            rebinned_arr = d
            rebinned_cont = cd
        
        #print '    *** Wavelengths corrected for redshift.'
        w_corr = rebinned_arr[0] / (1+float(z))
        object_spectra = np.array([w_corr, rebinned_arr[1]]) 
        contum_spectra = np.array([w_corr, rebinned_cont[1]])     # this is in case we want to rebin the continuum 
        # Save the redshift corrected whole spectra with the rebinned fluxes and continua
        for wa, fl, co in zip(w_corr, rebinned_arr[1], rebinned_cont[1]):
            object_wavelengths.append(wa)
            object_fluxes.append(fl)
            object_continuum.append(co)
    
        ''' # This commented section is used when we do not want the continuum to be rebinned, instead interpolate...
        cont = []
        for wrd in w_corr:
            fcd = np.interp(wrd, cd[0], cd[1])
            cont.append(fcd)
        contum_spectra = np.array([w_corr, cont]) 
        #print 'shapes of arrays after rebin:  object_spectra=', np.shape(object_spectra), '   contum_spectra=', np.shape(contum_spectra)
        '''    
        # Obtain the lines net fluxes and EWs
        new_file_name = object_name+"_lineinfo"+spectrum_region[s]+".txt"
        lineinfo_text_file = os.path.join(results4object_path, new_file_name)
        err_continuum = all_err_cont_fit[s] / 100.
        err_stis = err_stis_list[s]
        err_lists = [err_stis, err_continuum]
        # Now obtain the continuum and equivalent widths and write the _lineinfo files
        vacuum = False
        # in case manual changes need to be done to the line fluxes change use_mod_lineinfo_files to True
        if create_txt_lineinfo == False:
            print ' lineinfo file created already...  using previous files...'
            lines_info = spectrum.readlines_from_lineinfo(lineinfo_text_file)
            object_lines_info = [lines_info[0], lines_info[1], lines_info[6], lines_info[7], lines_info[8], lines_info[9]]
        else:    
            print ' created a lineinfo file for region', spectrum_region[s]
            object_lines_info = make_lineinfo_file(object_spectra, contum_spectra, Halpha_width, create_txt_lineinfo, vacuum, faintObj, lineinfo_text_file, err_lists)
        # line_info: 0=catalog_wavs_found, 1=central_wavelength_list, 2=width_list, 3=net_fluxes_list, 4=continuum_list, 5=EWs_list
    
        print 'There are ', len(object_lines_info[0]), ' lines in the ', spectrum_region[s]
        err_fluxes, err_continuum, err_ews = spectrum.get_lineinfo_uncertainties(object_spectra, contum_spectra, Halpha_width=Halpha_width, faintObj=faintObj, 
                                                                                 err_instrument=err_stis, err_continuum=err_continuum)
        '''    
        print '{:<12} {:>8} {:>15} {:>7} {:>15} {:>14} {:>5} {:>9} {:>8} {:>5}'.format('Wavelength', 'Flux', 'Flux err', '% err', 'Continuum', 'Continuum err', '% err', 'EW', 'EW err', '% err')
        for w, f, ef, c, ec, ew, eew in zip(object_lines_info[0], object_lines_info[3], object_lines_info[6], object_lines_info[4], err_continuum, object_lines_info[5], object_lines_info[7]):
            efp = (ef * 100.) / np.abs(f)
            ecp = (ec * 100.) / np.abs(c)
            eewp = (eew * 100.) / np.abs(ew)
            print '{:<10.2f} {:>14.5e} {:>12.5e} {:>6.1f} {:>16.5e} {:>12.5e} {:>6.1f} {:>10.2f} {:>6.2f} {:>6.1f}'.format(w, f, ef, efp, c, ec, ecp, ew, eew, eewp)
        '''
        make_text_file_errors = True
        if make_text_file_errors:
            err_file = os.path.join(results4object_path, object_name+"_lineerrs"+spectrum_region[s]+".txt")
            errf = open(err_file, 'w+')
            print >> errf, '{:<12} {:>8} {:>15} {:>7} {:>15} {:>14} {:>5} {:>9} {:>8} {:>5}'.format('Wavelength', 'Flux', 'Flux err', '% err', 'Continuum', 'Continuum err', '% err', 'EW', 'EW err', '% err')
            for w, f, ef, c, ec, ew, eew in zip(object_lines_info[0], object_lines_info[3], err_fluxes, object_lines_info[4], err_continuum, object_lines_info[5], err_ews):
                efp = (ef * 100.) / np.abs(f)
                ecp = (ec * 100.) / np.abs(c)
                eewp = (eew * 100.) / np.abs(ew)
                print >> errf, '{:<10.2f} {:>14.5e} {:>12.5e} {:>6.1f} {:>16.5e} {:>12.5e} {:>6.1f} {:>10.2f} {:>6.2f} {:>6.1f}'.format(w, f, ef, efp, c, ec, ecp, ew, eew, eewp)
            errf.close()
        print ''
        #raw_input('    press enter to continue...')
    
    # Gather all the *_lineinfo.txt files into a single file with the name defined below    
    add_str = "_lineinfo"
    text_file_list, _ = spectrum.get_obj_files2use(object_file, specs, add_str=add_str)
    # Read the observed lines from the table of lines_info.txt and normalize to Hbeta
    add_str = "_lineerrs"
    errs_files, _ = spectrum.get_obj_files2use(object_file, specs, add_str=add_str)
    
    # Define the file name for for all the lines
    name_out_file = os.path.join(results4object_path, object_name+"_linesNUV2NIR.txt")    
    cols_in_file, flxEW_errs, all_err_fit = spectrum.gather_specs(text_file_list, name_out_file, reject=50.0, start_w=None, create_txt=True, err_cont_fit=True, errs_files=errs_files)
else:
    # Load the text files with the measured lines and equivalent widths and write the corresponding file
    nameoutfile = os.path.join(results4object_path, object_name+"_measuredLI_linesNUV2NIR.txt") 
    divided_by_continuum = True   # Change to true if when fitted continuum used division to normalize  
    if object_number == 3:
        divided_by_continuum = False   
    cols_in_file, flxEW_errs = metallicity.use_measured_lineinfo_files(object_file, faintObj, Halpha_width, specs, data, cont_data, err_stis_list, 
                                                                       all_err_cont_fit, divided_by_continuum=divided_by_continuum,
                                                                       reject=1.0, start_w=None, create_txt=True, name_out_file=nameoutfile)

catalog_wavelength, observed_wavelength, element, ion, forbidden, how_forbidden, width, flux, continuum, EW = cols_in_file
err_flux, err_EW, cont_errs = flxEW_errs

# Determine the corresponding E(B-V) value for each object
av = A_V_list[object_number]
if use_Chbeta:
    ebv=0.0
else:
    ebv = A_B_list[object_number] - av
print 'This is the E(B-V) = ', ebv

# Do reddening correction 
# Choose the one we intend to use 
#redlaw= 'S79 H83'     ## Seaton (1979: MNRAS 187, 73)  and  Howarth (1983, MNRAS 203, 301) Galactic law
#redlaw = 'CCM89'      ## Cardelli, Clayton & Mathis 1989, ApJ 345, 245
#redlaw = 'CCM89 Bal07'## Cardelli, Clayton & Mathis 1989, ApJ 345, 245, modified by Balgrave et al 2007, 
#                         ApJ, 655, 299 for wavlength range 1250<lambda<3030.
redlaw = 'F99'        ## Fitzpatrick 1999 PASP, 11, 63
#redlaw = 'oD94'       ## O'Donnell 1994, ApJ, 422, 1580
#redlaw = 'LMCG03'     ## Gordon et al. (2003, ApJ, 594,279)
#redlaw = 'B07'        ## Blagrave et al 2007, ApJ, 655, 299
cHbeta = 0.434*C_Hbeta
# Names of the text files with the results with errors
if use_Chbeta:
    if use_given_lineinfo:
        name_extension = '_measuredLI_RedCor_CHbeta.txt'
        tfile2ndRedCor = os.path.join(results4object_path, object_name+"_Case"+case+"_measuredLI_2ndRedCor_CHbeta.txt")
    else:
        name_extension = '_RedCor_CHbeta.txt'
        tfile2ndRedCor = os.path.join(results4object_path, object_name+"_Case"+case+"_2ndRedCor_CHbeta.txt")
else:
    if use_given_lineinfo:
        name_extension = '_measuredLI_RedCor_Ebv.txt'
    else:
        name_extension = '_RedCor_Ebv.txt'
    tfile2ndRedCor = None
RedCor_file = os.path.join(results4object_path, object_name+name_extension)
print '     Reddening corrected intensities wrote on file:', RedCor_file

# Do reddening correction
opt_Hlines = True
kk = metallicity.BasicOps(results4object_path, redlaw, cols_in_file, I_theo_HaHb, EWabsHbeta, cHbeta, av, ebv, 
                          RedCor_file, do_errs=flxEW_errs, opt_Hlines=opt_Hlines)
normfluxes, Idered, I_dered_norCorUndAbs, errs_normfluxes, perc_errs_normfluxes, errs_Idered, perc_errs_I_dered, flambdas = kk.do_ops()
#flambdas = metallicity.find_flambdas(cHbeta, catalog_wavelength, I_dered_norCorUndAbs, normfluxes)
if use_Chbeta:
    idx4861 = catalog_wavelength.index(4861.33)
    idx6563 = catalog_wavelength.index(6562.820)
    new_ratio = Idered[idx6563]/Idered[idx4861]
    print '  new Halpha/Hbeta = %0.2f' % new_ratio
    #raw_input()

# Do collisional excitation correction
# If still fitting the C_Hbeta and EW_abs_Hbeta the code will stop here
if first_redcorr == True:
    exit()
# Write the first round of reddeding correction in pyneb readable format or correct for collisional excitation if
# asked to use CHbeta instead of E(B-V).
Hlines = None
# SBS1319       Halpha, H5,   H6,   H7,    H8,     H9,     H10,    H11,    H12
#theoCE_caseA = [2.80, 0.47, 0.265, 0.164, 0.109, 0.0760, 0.0553, 0.0415, 0.0320]
#theoCE_caseB = [2.85, 0.469, 0.260, 0.160, 0.105, 0.733, 0.0532, 0.0398, 0.0306]
# HOWERVER H8 and H9 are contaminated by HeI (H8 is also contaminated with [NeIII]), and H12 is too weak.
theoCE_caseA = [2.80, 0.47, 0.265, 0.164, 0.0553, 0.0415]
theoCE_caseB = [2.85, 0.469, 0.260, 0.160, 0.0532, 0.0398]
if case == 'A':
    theoCE = theoCE_caseA
elif case == 'B':
    theoCE = theoCE_caseB    
do_errs = None
writeouts=create_txt_temdenabunds
verbose = False
# This is only for pueposes of determining the metallicity Z of the object -- Non important but additional info:
# From Lopez-Sanchez % Esteban 2009, except for 1(Lopez-Sanchez&Esteban2006), 2(Lopez-SanchezEtal2004), and 6(Izotov&Yhuan1998)
#                    0         1         2         3         4        5         6         7         8
He_value_list = [0.087096, 0.081283, 0.0903325, 0.06166, 0.081283, 0.093325, 0.087, 0.075858, 0.087096,
                 0.075858, 0.087096, 0.075858, 0.087096, 0.089125, 0.085114, 0.093325, 0.123027, 0.058884]
#                    9        10        11        12         13       14        15         16         17
He_value = He_value_list[object_number]
#
advops = metallicity.AdvancedOps(results4object_path, redlaw, cols_in_file, I_theo_HaHb, EWabsHbeta, cHbeta, av, ebv, RedCor_file, do_errs,
                                 case, use_Chbeta, normfluxes, flambdas, Idered, I_dered_norCorUndAbs, theoCE, He_value, writeouts, verbose, tfile2ndRedCor)

lines_pyneb_matches = advops.perform_advanced_ops(forceTeH=forceTeO3, forceTeL=forceTeO2, forceNe=forceNe, theoCE=theoCE)


print '\n Code finished!'

