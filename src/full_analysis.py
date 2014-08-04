import numpy
import os
#import pyneb as pn 
import metallicity
from science import spectrum
#from uncertainties import unumpy

############################################################################################################################################


'''  Choose parameters to run script  '''

# name of the object
#                  0            1           2          3         4        5*        6         7         8
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                 9           10         11         12         13       14        15         16         17
object_number = 17

# Write the text file with line info?
create_txt_lineinfo = True

# If only wanting to perform the reddening and redshift correction set to True: first round of corrections
first_redcorr = False

# Do you want to use C_Hbeta to correct for extinction?   (if set to false the values of A_V and A_B will be used)
use_Chbeta = False

# Write the text file with temperatures, densities, and abundances?
create_txt_temdenabunds = True

# Choose case
case = 'B'

############################################################################################################################################

# Skip the finding of the line info and go directly to gather the spectra?
#                            0     1     2     3     4     5     6     7     8
use_given_lineinfo_list = [True, True, True, False, True, False, True, True, True, 
                           True, True, True, True, True, True, True, True, True]
#                            9     10    11    12    13    14    15    16    17
use_given_lineinfo = use_given_lineinfo_list[object_number]

# In case of wanting to use a specific temperature and/or density (accepts values with errors as lists)
#                      0                1                 2                3          4 (based on Na4)     5*    6     7     8           
forceTe_list = [[10900.,15000.],  [9000., 16000.], [10500.,12500.], [10100.,12100.], [15000.0, 17000.0], None, None, None, 11900.0,
                None,    16200., None, None, None, None, None, None, None]
#                   9      10*     11    12    13    14    15    16    17 
forceTe = forceTe_list[object_number]
#forceTe = None
#                   0      1             2              3        4     5    6     7     8           
forceNe_list = [None, [100., 650.], [100., 3000.], [800, 1000], None, None, None, None, [100., 655],
                None,    None, None, None, None, None, None, None, [100, 200.]]
#                 9      10*     11    12    13    14    15    16    17 
forceNe = forceNe_list[object_number]

# use this option if wanting to rebin or not
#                       0      1      2     3      4      5     6     7     8
perform_rebin_list = [False, False, True, False, False, True, True, True, False, 
                      False, False, True, True, True, True, True, True, True]
#                      9      10    11    12     13    14    15     16     17
perform_rebin = perform_rebin_list[object_number]

# use this option if object is VERY faint and want to use thinner widths for emission lines
#                  0      1     2     3      4     5     6      7      8
faintObj_list = [False, True, True, False, True, True, False, True, False, 
                 False, False, True, False, True, False, False, False, False]
#                  9      10    11    12     13    14     15     16     17
faintObj = faintObj_list[object_number]

# Do deblend of lines?
deblend4363 = False
deblend3727 = False
deblend6563 = False

# Used values of Halpha_width in order to properly correct for reddening
#                     0    1    2    3    4    5    6    7   8     9   10   11   12   13   14   15   16   17
Halpha_width_list = [40., 28., 28., 25., 33., 28., 30., 40., 35., 27., 27., 30., 43., 30., 30., 40., 50., 35.]
Halpha_width = Halpha_width_list[object_number]

# Set theoretical Halpha/Hbeta ratio
I_theo_HaHb = 2.86 

# Found values of EWabsHbeta and C_Hbeta in case the E(B-V) and Rv values are not known
#                   0            1           2            3            4          5*            6             7            8           
combos_list = [[2.0, 0.6], [2.0, 2.3], [2.0, 2.8], [0.7, 0.6], [2.0, 3.1], [0.1, 0.001], [0.8, 0.07], [2.6, 0.96], [2.5, 2.35], 
               [1.0, 1.15], [0.01, 0.01], [2.0, 2.], [0.8, 1.15], [2.5, 1.8], [2.5, 2.7], [2.7, 3.86], [1.6, 1.7], [0.001, 0.0001]]
#                   9            10*          11         12          13          14          15            16          17 
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
#                                    0              1                 2                 3              4                5
desired_disp_listoflists = [[1.6, 3.0, 5.0], [1.6, 6.5, 5.0], [1.6, 6.5, 5.0], [2.0, 5.0, 5.0], [2.0, 3.0, 5.0], [1.7, 5.6, 9.8], 
                            #        6              7                 8                9               10               11
                            [1.6, 5.6, 7.4], [1.8, 4.0, 6.0], [1.6, 3.0, 5.0], [1.6, 3.0, 5.0], [1.6, 3.5, 5.6], [1.6, 4.2, 5.0],
#                                   12             13                14               15               16               17                            
                            [1.6, 5.0, 5.0], [1.6, 4.5, 5.0], [1.6, 3.0, 5.0], [1.7, 5.6, 9.8], [1.6, 3.1, 5.1], [1.6, 8.3, 11.0]]
desired_disp_list = desired_disp_listoflists[object_number]

# corresponding redshifts from observations
#             0        1        2        3         4         5        6        7         8
z_list = [0.01985, 0.019581, 0.02813, 0.01354, 0.002695, 0.02346, 0.013631, 0.01201, 0.05842,
          0.046240, 0.013642, 0.002010, 0.0073, 0.01763, 0.01199, 0.032989, 0.04678, 0.002031]
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
all_err_cont_fit =[]
for s, tf in zip(specs, text_file_list): 
    wavs, flxs, cont_flxs = numpy.loadtxt(tf, skiprows=2, unpack=True)
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

# This is the error in the flux obtained by STIS, according to the Data Handbook the typical accuracies for 
# spectroscopic l mode is about 2%, 5%, and 5% for NUV, opt, and NIR - tables 4.1 and 4.2
err_stis_list = [0.015, 0.05, 0.05]

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
        object_spectra = numpy.array([w_corr, rebinned_arr[1]]) 
        contum_spectra = numpy.array([w_corr, rebinned_cont[1]])     # this is in case we want to rebin the continuum 
        # Save the redshift corrected whole spectra with the rebinned fluxes and continua
        for wa, fl, co in zip(w_corr, rebinned_arr[1], rebinned_cont[1]):
            object_wavelengths.append(wa)
            object_fluxes.append(fl)
            object_continuum.append(co)
    
        ''' # This commented section is used when we do not want the continuum to be rebinned, instead interpolate...
        cont = []
        for wrd in w_corr:
            fcd = numpy.interp(wrd, cd[0], cd[1])
            cont.append(fcd)
        contum_spectra = numpy.array([w_corr, cont]) 
        #print 'shapes of arrays after rebin:  object_spectra=', numpy.shape(object_spectra), '   contum_spectra=', numpy.shape(contum_spectra)
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
            print ' lineinfo file crated already...  using previous files...'
            lines_info = spectrum.readlines_from_lineinfo(lineinfo_text_file)
            object_lines_info = [lines_info[0], lines_info[1], lines_info[6], lines_info[7], lines_info[8], lines_info[9]]
        else:    
            print ' created a lineinfo file for region', spectrum_region[s]
            object_lines_info = make_lineinfo_file(object_spectra, contum_spectra, Halpha_width, create_txt_lineinfo, vacuum, faintObj, lineinfo_text_file, err_lists)
        # line_info: 0=catalog_wavs_found, 1=central_wavelength_list, 2=width_list, 3=net_fluxes_list, 4=continuum_list, 5=EWs_list
        if s == 1:
            if deblend3727:
                idx3726 = object_lines_info[0].index(3726.030)
                idx3729 = object_lines_info[0].index(3728.820)
                idx3727 = object_lines_info[0].index(3727.000)
                '''
                # Debled 3727 assuming that tot = 3729+3726. From ionic we know that 3726 = 0.7*3729
                o1 = object_lines_info[3][idx3726]
                tot = object_lines_info[3][idx3727]
                print 'I(3727) = %0.3e   ew = %0.3f' % (tot, tot/object_lines_info[4][idx3727])
                print 'I(3726) = %0.3e    I(3729) = %0.3e' % (o1, o1*0.7)
                print '   ew = %0.3f         ew = %0.3f ' % (o1/object_lines_info[4][idx3726], (o1*0.7)/object_lines_info[4][idx3729])
                '''
                ### Calculate fluxes of 3726 and 3739   -- this works as long as the sum of the individual lines add up to better than 85% of the total 3727 measurement
                # assume that real3726 + real3729 = real3727 and that measured3726 + measured3729 = measured3727
                # then we need a constant K that allows the fluxes of real3726 and real3729 to be comparable with measured3726 and measured3729
                measured3727 = object_lines_info[3][idx3727]
                measured3726 = object_lines_info[3][idx3726]
                measured3729 = object_lines_info[3][idx3729]
                K = measured3727 / (measured3726 + measured3729)
                real3726 = K * measured3726
                real3729 = K * measured3729
                # insert these fluxes and new equivalent wids in the corresponding places
                print 'PREVIOUS flux:  3726 = %0.3e    3729 = %0.3e    sum = %0.3e' % (object_lines_info[3][idx3726], object_lines_info[3][idx3729], object_lines_info[3][idx3726]+object_lines_info[3][idx3729])
                print '           ew:         %0.3f              %0.3f             %0.3f' % (object_lines_info[5][idx3726], object_lines_info[5][idx3729], object_lines_info[5][idx3726]+object_lines_info[5][idx3729])
                object_lines_info[3][idx3729] = real3729
                object_lines_info[5][idx3726] = real3726 / object_lines_info[4][idx3726]
                object_lines_info[5][idx3729] = real3729 / object_lines_info[4][idx3729]
                print '   NEW   flux:  3726 = %0.3e    3729 = %0.3e    sum = %0.3e' % (object_lines_info[3][idx3726], object_lines_info[3][idx3729], object_lines_info[3][idx3726]+object_lines_info[3][idx3729])
                print '           ew:         %0.3f              %0.3f             %0.3f' % (object_lines_info[5][idx3726], object_lines_info[5][idx3729], object_lines_info[5][idx3726]+object_lines_info[5][idx3729])
                raw_input()
        if s == 2:
            if deblend6563:
                # Deblend the N2 lines from 6563. This assumes that n1+H+n2 = tot, and we know that n2 is approx 3n1
                # hence, n1 = (tot - H) / 4
                idx6548 = object_lines_info[0].index(6548.030)
                idx6563 = object_lines_info[0].index(6562.820)
                idx6565 = object_lines_info[0].index(6565.00)
                idx6584 = object_lines_info[0].index(6583.410)
                H = object_lines_info[3][idx6563]
                tot = object_lines_info[3][idx6565]
                n1 = (tot - H) / 4.0
                print 'I(sum of 6548, 6563, 6584) = %0.3e' % tot
                print 'I(6548) = %0.3e    I(6563) = %0.3e    I(6584) = %0.3e' % (n1, H, n1*3)
                print '     ew = %0.3f           ew = %0.3f           ew = %0.3f' % (n1/object_lines_info[4][idx6548], H/object_lines_info[4][idx6563], (n1*3)/object_lines_info[4][idx6584])
                raw_input()
    
        print 'There are ', len(object_lines_info[0]), ' lines in the ', spectrum_region[s]
        err_fluxes, err_continuum, err_ews = spectrum.get_lineinfo_uncertainties(object_spectra, contum_spectra, Halpha_width=Halpha_width, faintObj=faintObj, 
                                                                                 err_instrument=err_stis, err_continuum=err_continuum)
        '''    
        print '{:<12} {:>8} {:>15} {:>7} {:>15} {:>14} {:>5} {:>9} {:>8} {:>5}'.format('Wavelength', 'Flux', 'Flux err', '% err', 'Continuum', 'Continuum err', '% err', 'EW', 'EW err', '% err')
        for w, f, ef, c, ec, ew, eew in zip(object_lines_info[0], object_lines_info[3], object_lines_info[6], object_lines_info[4], err_continuum, object_lines_info[5], object_lines_info[7]):
            efp = (ef * 100.) / numpy.abs(f)
            ecp = (ec * 100.) / numpy.abs(c)
            eewp = (eew * 100.) / numpy.abs(ew)
            print '{:<10.2f} {:>14.5e} {:>12.5e} {:>6.1f} {:>16.5e} {:>12.5e} {:>6.1f} {:>10.2f} {:>6.2f} {:>6.1f}'.format(w, f, ef, efp, c, ec, ecp, ew, eew, eewp)
        '''
        make_text_file_errors = True
        if make_text_file_errors:
            err_file = os.path.join(results4object_path, object_name+"_lineerrs"+spectrum_region[s]+".txt")
            errf = open(err_file, 'w+')
            print >> errf, '{:<12} {:>8} {:>15} {:>7} {:>15} {:>14} {:>5} {:>9} {:>8} {:>5}'.format('Wavelength', 'Flux', 'Flux err', '% err', 'Continuum', 'Continuum err', '% err', 'EW', 'EW err', '% err')
            for w, f, ef, c, ec, ew, eew in zip(object_lines_info[0], object_lines_info[3], err_fluxes, object_lines_info[4], err_continuum, object_lines_info[5], err_ews):
                efp = (ef * 100.) / numpy.abs(f)
                ecp = (ec * 100.) / numpy.abs(c)
                eewp = (eew * 100.) / numpy.abs(ew)
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
    cols_in_file, flxEW_errs = metallicity.use_measured_lineinfo_files(object_file, faintObj, Halpha_width, specs, data, cont_data, err_stis_list, 
                                                                       all_err_cont_fit, reject=50.0, start_w=None, create_txt=True, name_out_file=nameoutfile)

catalog_wavelength, observed_wavelength, element, ion, forbidden, how_forbidden, width, flux, continuum, EW = cols_in_file
err_flux, err_EW, cont_errs = flxEW_errs
'''
# Deblend lines
lines2deblend = [4359.340, 4363.210]
# the width is the same for all the components
idx4363 = catalog_wavelength.index(lines2deblend[1])
width_of_lines = width[idx4363]
tot_flx = flux[idx4363]
tot_cont = continuum[idx4363]
# Turn into arrays the whole object data
whole_object_data = numpy.array([object_wavelengths, object_fluxes])
whole_object_continua = numpy.array([object_wavelengths, object_continuum])
#spectrum.deblend_line(whole_object_data, whole_object_continua, catalog_wavelength, observed_wavelength, tot_flx, tot_cont, lines2deblend, width_of_lines, plot_fit=True)
#spectrum.deblend_with_gauss(whole_object_data, whole_object_continua, catalog_wavelength, observed_wavelength, tot_flx, tot_cont, lines2deblend, width_of_lines, plot_fit=True)
#exit()
'''
# Determine the corresponding E(B-V) value for each object
av = A_V_list[object_number]
if use_Chbeta:
    ebv=0.0
else:
    ebv = A_B_list[object_number] - av
    print 'This is the E(B-V) = ', ebv

# Do reddening correction 
# Choose the one we intend to use 
#redlaw= 'S 79 H 83'   ## Seaton (1979: MNRAS 187, 73)  and  Howarth (1983, MNRAS 203, 301) Galactic law
redlaw = 'CCM 89'     ## Cardelli, Clayton & Mathis 1989, ApJ 345, 245
#redlaw = 'oD 94'      ## O'Donnell 1994, ApJ, 422, 1580
#redlaw = 'LMC G 03'   ## Gordon et al. (2003, ApJ, 594,279)
#redlaw = 'B 07'       ## Blagrave et al 2007, ApJ, 655, 299
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
kk = metallicity.BasicOps(object_name, redlaw, cols_in_file, I_theo_HaHb, EWabsHbeta, cHbeta, av, ebv, RedCor_file, do_errs=flxEW_errs)
normfluxes, Idered, I_dered_norCorUndAbs, errs_normfluxes, perc_errs_normfluxes, errs_Idered, perc_errs_I_dered = kk.do_ops()
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
advops = metallicity.AdvancedOps(object_name, redlaw, cols_in_file, I_theo_HaHb, EWabsHbeta, cHbeta, av, ebv, RedCor_file, do_errs,
                                 case, use_Chbeta, theoCE, writeouts, verbose, tfile2ndRedCor)

lines_pyneb_matches = advops.perform_advanced_ops(forceTe=forceTe, forceNe=forceNe, theoCE=theoCE, )


print '\n Code finished!'

