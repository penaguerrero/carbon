from __future__ import division

import os
import pyfits
import numpy
import string
import copy
import math
import collections
import pyneb as pn
import PIL.Image as Image
from uncertainties import unumpy
from science import spectrum
from matplotlib import pyplot
from collections import OrderedDict
 

'''
This program contains various classes that together determine the object's metallicity.
'''

def round_wavs(catalog_wavelength, decimals=None):
    ### Round all catalog lines to make it easier to find lines
    rounded_catalog_wavelength = []
    for item in catalog_wavelength:
        if decimals != None:
            roundwav = numpy.round(item, decimals=decimals)
        else:
            roundwav = numpy.round(item)
        rounded_catalog_wavelength.append(roundwav)
    return rounded_catalog_wavelength

def use_measured_lineinfo_files(object_file, faintObj, Halpha_width, specs, data, cont_data, err_stis_list, all_err_cont_fit, reject=0.0, start_w=None, create_txt=False, name_out_file=None):
    all_wavs, all_flxs, all_cont, all_ews, all_ferrs, all_cerrs, all_ewerrs = read_measured_lineinfo_files(object_file, specs, data, cont_data, err_stis_list, all_err_cont_fit)
    data_arrays_list = [all_wavs, all_flxs, all_cont, all_ews, all_ferrs, all_cerrs, all_ewerrs]
    wavelengths, fluxes, continua, ews, errfluxes, errcontinua, errews = gather_measured_lineinfo(data_arrays_list, reject, start_w)
    data_arrays_list = [wavelengths, fluxes, continua, ews, errfluxes, errcontinua, errews]
    cols_in_file, flxEW_errs = find_matching_lines_with_catalog(data_arrays_list, faintObj, Halpha_width, all_err_cont_fit, create_txt, name_out_file)
    # cols_in_file contains: catalog_wavs_found, central_wavelength_list, found_element, found_ion, found_ion_forbidden, found_ion_how_forbidden, 
    #                        width_list, net_fluxes_list, continuum_list, EWs_list 
    # flxEW_errs contains: errs_net_fluxes, errs_ews, errs_conts
    return cols_in_file, flxEW_errs
  
def read_measured_lineinfo_files(object_file, specs, data, cont_data, err_stis_list, all_err_cont_fit):
    # object_file is the full path of the object in study
    # specs is a list of the spectral regions to be considered (i.e. specs=[0,1,2] )
    # cont_data is a list of numpy arrays of wavs and continuum fluxes for the specs in study
    # err_stis_list is the list with the 3 documented STIS percentage errors
    # all_err_cont_fit is the list of the 3 percentage errors from the continuum fit
    add_str = '_measured_lineinfo'
    text_file_list, _ = spectrum.get_obj_files2use(object_file, specs, add_str=add_str)
    all_wavs = []
    all_flxs = []
    all_cont = []
    all_ews = []
    all_ferrs = []
    all_cerrs = []
    all_ewerrs = []
    #print 'err_stis_list, all_err_cont_fit: ', err_stis_list, all_err_cont_fit
    for i in specs:
        # Load the text files with the measured lines and equivalent widths
        mwavs, normconts, normflxs, mews = numpy.loadtxt(text_file_list[i], unpack=True)
        mews = mews * -1   # this is to take the measured equivalent widths to the same formalism as the rest of the work
        realcont = cont_data[i] # this is the numpy array of wavs and flux continuum before normalization 
        # get the the flux and continuum values un-normalized and determine the errors in the fluxes and equivalent widths
        measuredflxs = numpy.array([])
        measuredconts = numpy.array([])
        measuredflxserrs = numpy.array([])
        measuredcontserrs = numpy.array([])
        measuredewserrs = numpy.array([])
        for mw, nf, nc, mew in zip(mwavs, normflxs, normconts, mews):
            c = numpy.interp(mw, realcont[0], realcont[1])
            mc = nc * c
            mf = nf * numpy.abs(c)
            measuredconts = numpy.append(measuredconts, mc)
            measuredflxs = numpy.append(measuredflxs, mf)
            # errors
            err_contfit = all_err_cont_fit[i]/100.0
            err_stis = err_stis_list[i]
            cabserrdue2me = numpy.abs(nc - 1.0)                           # human error
            cabserrdue2fit = err_contfit * numpy.abs(nc)                  # error due to the fit
            # total error in the continuum 
            cerr = numpy.sqrt(cabserrdue2me**2 + cabserrdue2fit**2) * numpy.abs(c)   
            fabserr_due2instrument = err_stis * numpy.abs(nf)             # intrinsic instrument error
            # I estimate that my error in determining where the line begins/ends is ~15%, hence I have to go little further in wavelength to determine
            # how much that affects the flux.
            object_wavsflx = data[i]
            leftf = numpy.interp(mw-3.0, object_wavsflx[0], object_wavsflx[1]) / numpy.abs(c)
            rightf = numpy.interp(mw+3.0, object_wavsflx[0], object_wavsflx[1]) / numpy.abs(c)
            #fabserrdue2me = numpy.sqrt((leftf-nf)**2 + (rightf-nf)**2) * 0.15 
            fabserrdue2me = (numpy.abs(leftf-nf) + numpy.abs(rightf-nf))/2 * 0.15   # this is assuming that the error is symmetric
            ferr = numpy.sqrt(fabserr_due2instrument**2 + fabserrdue2me**2) * numpy.abs(c)
            ewerr = mew * numpy.sqrt((ferr/mf)**2 + (cerr/mc)**2)
            #print 'normflux=',nf, '  fabserr_due2instrument=', fabserr_due2instrument, '  fabserrdue2me=', fabserrdue2me
            #print 'normcont=',nc,'  cabserrdue2me=', cabserrdue2me, '  cabserrdue2fit=', cabserrdue2fit
            #print mw, 'f=',mf,'+-',ferr, (ferr*100)/mf, '  c=',mc,'+-',cerr,(cerr*100.)/mc, '  ew=', mew,'+-', ewerr, (ewerr*100.)/mew
            measuredflxserrs = numpy.append(measuredflxserrs, ferr)
            measuredcontserrs = numpy.append(measuredcontserrs, cerr)
            measuredewserrs = numpy.append(measuredewserrs, ewerr)
            #raw_input()
        all_wavs.append(mwavs)
        all_flxs.append(measuredflxs)
        all_cont.append(measuredconts)
        all_ews.append(mews)
        all_ferrs.append(measuredflxserrs)
        all_cerrs.append(measuredcontserrs)
        all_ewerrs.append(measuredewserrs)
    return all_wavs, all_flxs, all_cont, all_ews, all_ferrs, all_cerrs, all_ewerrs   # these are all lists of numpy arrays

def gather_measured_lineinfo(data_arrays_list, reject=0.0, start_w=None):
    '''This function gathers the infro read from the measured lineinfo files and returns 7 single numpy arrays for all the
    observed spectral range.'''
    all_wavs, all_flxs, all_cont, all_ews, all_ferrs, all_cerrs, all_ewerrs = data_arrays_list
    wavelengths = numpy.array([])
    fluxes = numpy.array([])
    continua = numpy.array([])
    ews = numpy.array([])
    errfluxes = numpy.array([])
    errcontinua = numpy.array([])
    errews = numpy.array([])
    # The input are expected to be lists of numpy arrays
    for i in range(len(all_wavs)):  # loop through the lists
        for _ in zip(all_wavs[i], all_flxs[i], all_cont[i], all_ews[i]):
            if start_w == None:
                ini_wav = all_wavs[i][0]+reject
            else:
                # Make sure to start the right array and leave the others with the regular reject
                if start_w <= 2000.0:
                    if i == 0:  # this is for the NUV spec
                        ini_wav = start_w
                    else:
                        ini_wav = all_wavs[i][0]+reject                
                if (start_w > 2000.0) and (start_w < 5000.0):
                    if i == 1:  # this is for the OPT spec
                        ini_wav = start_w
                    else:
                        ini_wav = all_wavs[i][0]+reject                
                if start_w >= 5000.0:  # this is for the NIR spec
                    if i == 2:
                        ini_wav = start_w
                    else:
                        ini_wav = all_wavs[i][0]+reject
        initial_wav, initial_idx = spectrum.find_nearest(all_wavs[i], ini_wav)
        _, ending_idx = spectrum.find_nearest(all_wavs[i], all_wavs[i][-1]+reject)
        wavelengths = numpy.append(wavelengths, initial_wav)
        fluxes = numpy.append(fluxes, all_flxs[i][initial_idx])
        continua = numpy.append(continua, all_cont[i][initial_idx])
        ews = numpy.append(ews, all_ews[i][initial_idx])
        errfluxes = numpy.append(errfluxes, all_ferrs[i][initial_idx])
        errcontinua = numpy.append(errcontinua, all_cerrs[i][initial_idx])
        errews = numpy.append(errews, all_ewerrs[i][initial_idx])
        # In order to have a single array, every single item must be added to the array 
        for j in range(len(all_wavs[i][initial_idx : ending_idx])):
            idx = initial_idx+1 + j 
            wavelengths = numpy.append(wavelengths, all_wavs[i][idx])
            fluxes = numpy.append(fluxes, all_flxs[i][idx])
            continua = numpy.append(continua, all_cont[i][idx])
            ews = numpy.append(ews, all_ews[i][idx])
            errfluxes = numpy.append(errfluxes, all_ferrs[i][idx])
            errcontinua = numpy.append(errcontinua, all_cerrs[i][idx])
            errews = numpy.append(errews, all_ewerrs[i][idx])
    # Choose the right intensities if lines are repeated
    # Because the sensibility of the detector is better with the g430l than with the g230l, we want to keep OPT the lines.
    # Because the sensibility of the detector is better with the g750l than with the g430l, we want to keep NIR the lines.
    rounded_obs_wavs = round_wavs(wavelengths, decimals=0)
    repeated_lines = collections.Counter(rounded_obs_wavs)
    rl_list = [i for i in repeated_lines if repeated_lines[i]>1]
    list_of_indeces_to_remove =[]
    for i in range(len(rounded_obs_wavs)):
        if rounded_obs_wavs[i] in rl_list:
            if (rounded_obs_wavs[i] > 3015.0) and (rounded_obs_wavs[i] < 3150.0) or (rounded_obs_wavs[i] > 5000.0) and (rounded_obs_wavs[i] < 5600.0):
                list_of_indeces_to_remove.append(i)
    # Remove the first occurence
    wavelengths = numpy.delete(wavelengths, list_of_indeces_to_remove)
    fluxes = numpy.delete(fluxes, list_of_indeces_to_remove)
    continua = numpy.delete(continua, list_of_indeces_to_remove)
    ews = numpy.delete(ews, list_of_indeces_to_remove)   
    errfluxes = numpy.delete(errfluxes, list_of_indeces_to_remove)
    errcontinua = numpy.delete(errcontinua, list_of_indeces_to_remove)
    errews = numpy.delete(errews, list_of_indeces_to_remove)   
    return wavelengths, fluxes, continua, ews, errfluxes, errcontinua, errews
            
def find_matching_lines_with_catalog(measured_data, faintObj, Halpha_width, all_err_cont_fit, create_txt=False, name_out_file=None):
    ''' Find the lines that match the measured lines
    # err_stis_list is the documented error for that spectral region
    # all_err_cont_fit is the list of errors of the fitted continua
    # faintObj is a true or false
    # Halpha_width is the width of Halpha
    # name_out_file is the name of the output file
    # create_txt is expected to be a True or False argument
    '''
    # Read the line_catalog file, assuming that the path is the same:
    # '/Users/name_of_home_directory/Documents/AptanaStudio3/science/science/spectrum/lines_catalog.txt'
    line_catalog_path = os.path.abspath('../../science/science/spectrum/lines_catalog.txt')
    measured_wavs, measured_flxs, measured_conts, measured_ews, errfluxes, errcontinua, errews = measured_data 
    # Define the columns of the catalog file
    catalog_wavelength = []
    element = []
    ion =[]
    forbidden = []
    how_forbidden = []
    transition = []
    strong_line = []
    cols_in_file = [catalog_wavelength, element, ion, forbidden, how_forbidden, transition, strong_line]
    # Define the list of the files to be read
    text_file_list = [line_catalog_path]
    # Read the catalog files
    data, widths_faintObj, widths_strongObj  = spectrum.readlines_from_lineinfo(text_file_list, cols_in_file)
    catalog_wavelength, element, ion, forbidden, how_forbidden, transition, strong_line = data
    # If the wavelength is grater than 2000 correct the theoretical air wavelengths to vacuum using the IAU
    # standard for conversion from air to vacuum wavelengths is given in Morton (1991, ApJS, 77, 119). To
    # correct find the refraction index for that wavelength and then use:
    #       wav_vac / wav_air -1 = n - 1
    # (To start I usded NIST, n=0.999271)
    wavs_air = []
    wavs_vacuum = []
    for w in catalog_wavelength:
        # separate air and vacuum wavelengths into 2 lists
        if w < 2000.0:
            # For keeping all vacuum wavelengths
            wavs_vacuum.append(w)
            # For converting vaccuum to air
            wav_refraction_index = spectrum.n4airvac_conversion(w)
            #print 'Refraction index  n = %f' % (wav_refraction_index)
            wair = w / (2 - wav_refraction_index)
            wavs_air.append(wair)
        elif w >= 2000.0:
            # For converting to vacuum wavelengths
            wav_refraction_index = spectrum.n4airvac_conversion(w)
            wvac = w * (2 - wav_refraction_index)
            wavs_vacuum.append(wvac)
            wavs_air.append(w)
    # Determine the strength of the lines: 
    width = []
    for sline in strong_line:
        if faintObj == True: 
            if sline == "nw":
                s = widths_faintObj[0]
            elif sline == "no":
                s = widths_faintObj[1]
            elif sline == "weak":
                s = widths_faintObj[2]
            elif sline == "medium":
                s = widths_faintObj[3]
            elif sline == "yes":
                s = widths_faintObj[4]
            elif sline == "super":
                s = widths_faintObj[5]
            elif sline == "Halpha":
                s = Halpha_width
            width.append(s)
        else:
            if sline == "nw":
                s = widths_strongObj[0]
            elif sline == "no":
                s = widths_strongObj[1]
            elif sline == "weak":
                s = widths_strongObj[2]
            elif sline == "medium":
                s = widths_strongObj[3]
            elif sline == "yes":
                s = widths_strongObj[4]
            elif sline == "super":
                s = widths_strongObj[5]
            elif sline == "Halpha":
                s = Halpha_width
            width.append(s)
    # Search in the object given for the lines in the lines_catalog
    lines_catalog = (wavs_air, wavs_vacuum, element, ion, forbidden, how_forbidden, transition, width)
    net_fluxes_list = []
    EWs_list = []
    central_wavelength_list =[]
    catalog_wavs_found = []
    continuum_list =[]
    width_list = []
    found_element = []
    found_ion = []
    found_ion_forbidden = []
    found_ion_how_forbidden = []
    errs_net_fluxes = []
    errs_conts = []
    errs_ews = [] 
    # but choose the right wavelength column
    vacuum = False  # wavelengths are expected to be in AIR
    if vacuum == True:
        use_wavs = 1
        use_wavs_text = '# Used VACUUM wavelengths to find lines from line_catalog.txt'
    if vacuum == False:
        use_wavs = 0
        use_wavs_text = '# Used AIR wavelengths to find lines from line_catalog.txt'
    #print 'vacuum was set to %s, %s' % (vacuum, use_wavs_text)
    for i in range(len(lines_catalog[0])):
        # find the line in the catalog that is closest to a 
        line_looked_for = lines_catalog[use_wavs][i]
        nearest2line = spectrum.find_nearest_within(measured_wavs, line_looked_for, 10.0)
        nearest2line_idx = numpy.where(measured_wavs == nearest2line)
        idx = nearest2line_idx[0]
        if nearest2line > 0.0:  
            catalog_wavs_found.append(line_looked_for)
            central_wavelength = float(measured_wavs[idx])
            line_width = lines_catalog[7][i]
            round_line_looked_for = numpy.round(line_looked_for, decimals=0)
            if (round_line_looked_for == 1907.0) or (round_line_looked_for == 1909.0):
                if faintObj:
                    line_width = 3.0
                else:
                    line_width = 7.0
            if (line_looked_for ==  4267.15) or (line_looked_for == 4640.0) or (line_looked_for == 4650.0):
                line_width = 5.0 
            f = float(measured_flxs[idx])
            c = float(measured_conts[idx])
            e = float(measured_ews[idx])
            ferr = float(errfluxes[idx])
            errs_net_fluxes.append(ferr)
            errs_conts.append(float(errcontinua[idx]))
            errs_ews.append(float(errews[idx]))
            print '\n Looking for:',  line_looked_for
            print nearest2line, '   flux =', f,'+-',ferr, '   cont =', c, '   ew =', e
            #if line_looked_for ==  2144.0:
            #    raw_input()
            #if line_looked_for ==  4650.0:
            #    raw_input()
            width_list.append(line_width)
            central_wavelength_list.append(central_wavelength)
            continuum_list.append(c)
            net_fluxes_list.append(f)
            EWs_list.append(e) 
            found_element.append(lines_catalog[2][i])
            found_ion.append(lines_catalog[3][i])
            found_ion_forbidden.append(lines_catalog[4][i])
            found_ion_how_forbidden.append(lines_catalog[5][i])
    # Create the table of Net Fluxes and EQWs
    if name_out_file != None:
        if create_txt:
            #linesinfo_file_name = raw_input('Please type name of the .txt file containing the line info. Use the full path.')
            txt_file = open(name_out_file, 'w+')
            print >> txt_file,  use_wavs_text
            print >> txt_file,   '# Positive EW = emission        Negative EW = absorption' 
            print >> txt_file,   '# Percentage Errors of Continuum Fits: NUV, Opt, NIR = %0.2f, %0.2f, %0.2f' % (all_err_cont_fit[0], all_err_cont_fit[1], all_err_cont_fit[2])
            print >> txt_file,   '#    NUV: wav <= 2000,   Opt: 2000 > wav < 5000,   NIR: wav >= 5000'
            print >> txt_file,  ('{:<12} {:<12} {:>10} {:<4} {:<9} {:>8} {:<9} {:<12} {:<10} {:<7} {:<12} {:<7} {:<6} {:<6}'.format('# Catalog WL', 'Observed WL', 'Element', 
                                                                                                                                    'Ion', 'Forbidden', 'HowForb', 'Width[A]', 
                                                                                                                                    'Flux [cgs]', 'FluxErr', '%Err', 
                                                                                                                                    'Continuum [cgs]', 'EW [A]', 'EWErr', '%Err'))
            for cw, w, e, i, fd, h, s, F, Fe, C, ew, ewe in zip(catalog_wavs_found, central_wavelength_list, found_element, found_ion, found_ion_forbidden, 
                                                       found_ion_how_forbidden, width_list, net_fluxes_list, errs_net_fluxes, continuum_list, EWs_list, errs_ews):
                Fep = (Fe * 100.) / numpy.abs(F)
                ewep = (ewe * 100.) / numpy.abs(ew)
                print >> txt_file,  ('{:<12.3f} {:<12.3f} {:>10} {:<6} {:<8} {:<8} {:<6} {:>12.3e} {:>10.3e} {:>6.1f} {:>13.3e} {:>10.3f} {:>6.3f} {:>6.1f}'.format(cw, w, e, i, fd, h, s, F, Fe, Fep, C, ew, ewe, ewep))
            txt_file.close()
            print 'File   %s   writen!' % name_out_file
    elif create_txt == False:
        print '# Positive EW = emission        Negative EW = absorption' 
        print  ('{:<12} {:<12} {:>10} {:<4} {:<9} {:>8} {:<9} {:<12} {:<10} {:<7} {:<12} {:<7} {:<6} {:<6}'.format('# Catalog WL', 'Observed WL', 'Element', 
                                                                                                                                'Ion', 'Forbidden', 'HowForb', 'Width[A]', 
                                                                                                                                'Flux [cgs]', 'FluxErr', '%Err', 
                                                                                                                                'Continuum [cgs]', 'EW [A]', 'EWErr', '%Err'))
        for cw, w, e, i, fd, h, s, F, Fe, C, ew, ewe in zip(catalog_wavs_found, central_wavelength_list, found_element, found_ion, found_ion_forbidden, 
                                                   found_ion_how_forbidden, width_list, net_fluxes_list, errs_net_fluxes, continuum_list, EWs_list, errs_ews):
            Fep = (Fe * 100.) / numpy.abs(F)
            ewep = (ewe * 100.) / numpy.abs(ew)
            print ('{:<12.3f} {:<12.3f} {:>10} {:<6} {:<8} {:<8} {:<6} {:>12.3e} {:>10.3e} {:>6.1f} {:>13.3e} {:>10.3f} {:>6.3f} {:>6.1f}'.format(cw, w, e, i, fd, h, s, F, Fe, Fep, C, ew, ewe, ewep))
    cols_in_file = [catalog_wavs_found, central_wavelength_list, found_element, found_ion, found_ion_forbidden, found_ion_how_forbidden, width_list, net_fluxes_list, continuum_list, EWs_list] 
    flxEW_errs = [errs_net_fluxes, errs_ews, errs_conts]
    return cols_in_file, flxEW_errs

def normalize_lines(rounded_catalog_wavelength, element, ion, forbidden,
                        how_forbidden, observed_wavelength, flux, intensities, EW, continuum):
    Hb_idx = rounded_catalog_wavelength.index(4861.0)
    catalog_lines = []
    element_lines = []
    ion_lines = []
    forbidden_lines = []
    how_forbidden_lines = []
    wavs_lines = []
    normfluxes = []
    norm_intensities = []
    EW_lines = []
    calc_cont = []
    for i in range(0, len(flux)):
        catalog_lines.append(rounded_catalog_wavelength[i])
        element_lines.append(element[i])
        ion_lines.append(ion[i])
        forbidden_lines.append(forbidden[i])
        how_forbidden_lines.append(how_forbidden[i])
        wavs_lines.append(observed_wavelength[i])
        norm_flux = flux[i] / flux[Hb_idx] * 100.
        normfluxes.append(norm_flux)
        normI = intensities[i] / intensities[Hb_idx] * 100.
        norm_intensities.append(normI)
        EW_lines.append(EW[i])
        calc_cont.append(continuum[i])
    return (catalog_lines, wavs_lines, element_lines, ion_lines, forbidden_lines, 
            how_forbidden_lines, normfluxes, calc_cont, norm_intensities, EW_lines)    

def find_emission_lines(rounded_catalog_wavelength, element, ion, forbidden,
                        how_forbidden, observed_wavelength, flux, intensities, EW, continuum):
    catalog_emission_lines = []
    element_emission_lines = []
    ion_emission_lines = []
    forbidden_emission_lines = []
    how_forbidden_emission_lines = []
    wavs_emission_lines = []
    positive_normfluxes = []
    positive_norm_intensities = []
    EW_emission_lines = []
    pos_calc_cont = []
    for i in range(len(flux)):
        if flux[i] > 0.00000:
            catalog_emission_lines.append(rounded_catalog_wavelength[i])
            element_emission_lines.append(element[i])
            ion_emission_lines.append(ion[i])
            forbidden_emission_lines.append(forbidden[i])
            how_forbidden_emission_lines.append(how_forbidden[i])
            wavs_emission_lines.append(observed_wavelength[i])
            positive_normfluxes.append(flux[i])
            positive_norm_intensities.append(intensities[i])
            EW_emission_lines.append(EW[i])
            pos_calc_cont.append(continuum[i])
    #print '  ***  There are ', len(positive_norm_intensities), ' emission lines in this object!'
    return (catalog_emission_lines, wavs_emission_lines, element_emission_lines, ion_emission_lines, forbidden_emission_lines, 
            how_forbidden_emission_lines, positive_normfluxes, pos_calc_cont, positive_norm_intensities, EW_emission_lines)

def find_flambdas(cHbeta, catalog_wavelength, I_dered_noCorUndAbs, normfluxes):
    # Finding the f_lambda values
    all_flambdas = []
    for Icor, Iobs, w in zip(I_dered_noCorUndAbs, normfluxes, catalog_wavelength):
        if Iobs < 0.0:   # this is to avoid nan values
            all_flambdas.append(0.0)
        else:
            f12 = (numpy.log10(Icor) - numpy.log10(Iobs)) / cHbeta
            all_flambdas.append(f12)
    # remove the zeros in order to interpolate
    wavs = []
    pos_flambdas = []
    for w, f in zip(catalog_wavelength, all_flambdas):
        if f != 0.0:
            wavs.append(w)
            pos_flambdas.append(f)
    # now replace the 0.0 for interpolated values
    flambdas = []
    for w, f in zip(catalog_wavelength, all_flambdas):
        if f is 0.0:
            newf = numpy.interp(w, wavs, pos_flambdas)
            flambdas.append(newf)
        else:
            flambdas.append(f)
    return flambdas

def write_RedCorfile(object_name, path_and_name_outfile, cols_in_file):
    # Columns to be written in the text file: 
    catalog_wavelength, flambdas, element, ion, forbidden, how_forbidden, normfluxes, errs_normfluxes, perc_errs_normfluxes, Idered, errs_Idered, perc_errs_I_dered, EW, err_EW, perrEW = cols_in_file
    RedCor_file = open(path_and_name_outfile, 'w+')
    print >> RedCor_file, '#{:<12} {:>8} {:>8} {:<5} {:>6} {:>6} {:>10} {:>9} {:>6} {:>14} {:>8} {:>5} {:>11} {:>8} {:>5}'.format('Wavelength', 'flambda', 'Element', 'Ion', 'Forbidden', 'How much', 'Flux', 'FluxErr', '%err', 'Intensity', 'IntyErr', '%err', 'EW', 'EWerr', '%err')
    for w, fl, ele, ion, forb, hforb, f, ef, epf, i, ei, epi, ew, eew, eewp in zip(catalog_wavelength, flambdas, element, ion, forbidden, how_forbidden, normfluxes, errs_normfluxes, perc_errs_normfluxes, Idered, errs_Idered, perc_errs_I_dered, EW, err_EW, perrEW):
        print >> RedCor_file, '{:<10.2f} {:>9.3f} {:>8} {:<6} {:>6} {:>8} {:>14.3f} {:>8.3f} {:>6.1f} {:>13.3f} {:>8.3f} {:>6.1f} {:>12.2f} {:>6.2f} {:>6.1f}'.format(w, fl, ele, ion, forb, hforb, f, ef, epf, i, ei, epi, ew, eew, eewp)
    RedCor_file.close()
    print 'Text file with corrected intensities written in: ', path_and_name_outfile

def read_IonTotAbs(file_name):
    ''' This function reads the ionic and total abundances text files and returns a dictionary with:
    element, abundance, abund error, %err, logAbund, logErr, X/O, X/Oerr '''
    # the file columns are for the following elements in this specific order:
    default_elements = ['O', 'N', 'Ne', 'S', 'Cl', 'Ar', 'Fe', 'C']
    # these are the columns of interest:
    elements = []
    abund = []
    abunderr = []
    perr = []
    logabund = [] 
    loerr = []
    XO = []
    XOerr = []
    cols = [elements, abund, abunderr, perr, logabund, loerr, XO, XOerr]
    f = open(file_name, 'r')
    for row in f:
        line_data = row.split()
        # the common thing between the lines of interest is the number of columns 
        if (len(line_data) == 8) and ('#' not in row):
            # Split each row into columns
            line_data = row.split()
            # append the element into each column in the cols_in_file
            for item, col in zip(line_data, cols):
                if item not in default_elements:
                    item = float(item)
                col.append(item)
    f.close()
    # create the dictionary for this object
    elem_abun = OrderedDict()
    for el, a, ae, pe, la, le, xo, xoe in zip(elements, abund, abunderr, perr, logabund, loerr, XO, XOerr):
        elem_abun[el] = [a, ae, pe, la, le, xo, xoe]
    return elem_abun


class OneDspecs:
    '''
    This uses the 1d fits extraction files from all three filters (g230l, g430l, and g750l), 
    in order to let me choose which one has the best S/N. 
    '''
    def __init__(self, oneDspecs, selectedspecs_file, plot_name, img_format, full_name):
        self.oneDspecs = oneDspecs 
        self.wavs = []
        self.flxs = []
        self.colors = ['b', 'r', 'g', 'c', 'k', 'm', 'y']
        self.specs = []
        self.counter = 0
        self.xlolim = 1500
        self.xuplim = 10500
        self.ylolim = -1.5e-13
        self.yuplim = 1.5e-13
        self.indx_selected_wavs_and_flxs = []
        self.wf_used_specs = []
        self.get_specs()
        self.plot_chose_from(plot_name, img_format)
        self.selecting_specs(selectedspecs_file, plot_name, img_format)
        self.smooth_spec(plot_name, img_format, full_name)
        
    def get_specs(self):
        print 'Found %i files in the oneDspecs list.' % len(self.oneDspecs)
        for spec in self.oneDspecs:
            s = pyfits.open(spec)
            s.info()
            print('This will be spectrum number %i' % self.counter)            
            # read the data and store the wavelengths and fluxes
            t = s['sci'].data
            self.specs.append(repr(self.counter)+'='+spec+'_'+repr('sci'))
            w = t.WAVELENGTH.flatten()
            f = t.FLUX.flatten()
            self.wavs.append(w)
            self.flxs.append(f)
            # Count the spectra you find with the sci extension
            self.counter = self.counter + 1   

    def do_plot(self, plot_name, img_format, specs_list, x_axs, y_axs, xlolim, xuplim, ylolim, yuplim, save=True):
        fig1 = pyplot.figure(1, figsize=(15, 8))
        # Making the plot
        j = 0
        k = 0
        curves = []
        pyplot.xlim(xlolim, xuplim)
        pyplot.ylim(ylolim, yuplim)
        for i in range(len(x_axs)):
            #print(len(x_axs[i]), len(y_axs[i]), len(specs_list[i]))
            if i > 6:
                j = k
                c, = pyplot.plot(x_axs[i], y_axs[i], self.colors[j])
                k = k + 1
                if k > 6:
                    k = 0
            else:
                c, = pyplot.plot(x_axs[i], y_axs[i], self.colors[j])
            curves.append(c)
            #print j, k
            j = j + 1
        pyplot.figure(1).legend(curves, specs_list, labelspacing=0.2, bbox_to_anchor=(0.9, 0.9))
        if save == True:
            fig1.savefig(plot_name+img_format)
        pyplot.show(fig1)

    def save_plot(self, plot_name, img_format, specs_list, x_axs, y_axs, xlolim, xuplim, ylolim, yuplim):
        pyplot.ioff()
        save_plt = raw_input('Do you want to save this plot? (n)  y  ')
        if save_plt == 'n' or save_plt == '':
            os.remove(plot_name+img_format)
            print('Plot not saved.')
        else:
            axis_scale = raw_input('Do you want to modify axis scale? (n) y  ')
            if axis_scale == 'y':
                x_axis = raw_input('Want to modify x-axis? (n) y  ')
                if x_axis == 'y':
                    xlo = raw_input('Enter lower x-limit: ')
                    xlolim = float(xlo)
                    xup = raw_input('Enter upper x-limit: ')
                    xuplim = float(xup)
                else:
                    print('x-lolim: %i,  x-uplim: %i' % (xlolim, xuplim))
                y_axis = raw_input('Want to modify y-axis? (n) y  ')
                if y_axis == 'y':
                    ylo = raw_input('Enter lower y-limit: ')
                    ylolim = float(ylo)
                    yup = raw_input('Enter upper y-limit: ')
                    yuplim = float(yup)
                else:
                    print('y-lolim: %e,  y-uplim: %e' % (ylolim, yuplim))
                os.remove(plot_name+img_format)
                self.do_plot(plot_name, img_format, specs_list, x_axs, y_axs, xlolim, xuplim, ylolim, yuplim)
            else:
                print('x-lolim: %i,  x-uplim: %i, y-lolim: %e,  y-uplim: %e' % (xlolim, xuplim, ylolim, yuplim))
            #plot_name = raw_input('Enter plot name (including directory)')
            print('Plot %s was saved!' % plot_name+img_format)
      
    def plot_chose_from(self, plot_name, img_format):
        xlolim = self.xlolim
        xuplim = self.xuplim
        ylolim = self.ylolim
        yuplim = self.yuplim
        keep_prevFig = raw_input('Is there is a previous figure you want to keep?  (n)  y   ')
        if keep_prevFig == 'y':
            self.do_plot(plot_name, img_format, self.specs, self.wavs, self.flxs, xlolim, xuplim, ylolim, yuplim, save=False)
        else:
            self.do_plot(plot_name, img_format, self.specs, self.wavs, self.flxs, xlolim, xuplim, ylolim, yuplim, save=True)
            pyplot.ioff()
            self.save_plot(plot_name, img_format, self.specs, self.wavs, self.flxs, xlolim, xuplim, ylolim, yuplim)
    
    def selecting_specs(self, selectedspecs_file, plot_name, img_format):
        # Now choose the spectra to use for building the whole spectrum
        print('List which spectra do you want to keep for NUV, optical, and NIR (i.e. 0,2,5)')
        indxs = raw_input()
        indx_list = string.split(indxs, sep=',')
        for i in indx_list:
            self.indx_selected_wavs_and_flxs.append(int(i))
        # Do the plot of the selected files, save the text file, and save the eps
        xlolim = self.xlolim
        xuplim = self.xuplim
        ylolim = self.ylolim
        yuplim = self.yuplim
        # this is assuming there are ONLY 3 spectra that were selected, 1 for nuv, 1 for opt, and 1 for nir.
        w_nuv = copy.deepcopy(self.wavs[self.indx_selected_wavs_and_flxs[0]])
        w_opt = copy.deepcopy(self.wavs[self.indx_selected_wavs_and_flxs[1]])
        w_nir = copy.deepcopy(self.wavs[self.indx_selected_wavs_and_flxs[2]])
        f_nuv = copy.deepcopy(self.flxs[self.indx_selected_wavs_and_flxs[0]])
        f_opt = copy.deepcopy(self.flxs[self.indx_selected_wavs_and_flxs[1]])
        f_nir = copy.deepcopy(self.flxs[self.indx_selected_wavs_and_flxs[2]])
        # store the wavelengths and fluxes of the spectra selected in wf_used_specs            
        wf_nuv = numpy.array([w_nuv, f_nuv])
        wf_opt = numpy.array([w_opt, f_opt])
        wf_nir = numpy.array([w_nir, f_nir])
        # wf_used_specs is a tuple of arrays composed by wavelengths and fluxes
        wf = (wf_nuv, wf_opt, wf_nir)
        wf_arr = numpy.array(wf)
        self.wf_used_specs.append(wf_arr) 
        # make a list of the spectra used to make things prettier in the plot
        used_specs = []    
        for i in self.indx_selected_wavs_and_flxs:
            used_specs.append(self.specs[i])              
        # saving or not saving the text files
        write_file_wf = raw_input('Is there a wavelength-flux set of texts files you want to keep? (n)  y  ')
        if write_file_wf == 'n' or write_file_wf == '':
            #txtout = open(selectedspecs_file, 'w+')
            txtout0 = open(selectedspecs_file+'_nuv.txt', 'w+')
            txtout1 = open(selectedspecs_file+'_opt.txt', 'w+')
            txtout2 = open(selectedspecs_file+'_nir.txt', 'w+')
            # Writting the text file of wavelengths and fluxes of the selected spectra
            for i in range(len(w_opt)):
                txtout0.write("%.2f    %.3e \n" % (w_nuv[i], f_nuv[i]))
                txtout1.write("%.2f    %.3e \n" % (w_opt[i], f_opt[i]))
                txtout2.write("%.2f    %.3e \n" % (w_nir[i], f_nir[i]))
            txtout0.close()
            txtout1.close()
            txtout2.close()            
        # now do the plotting
        plot_selecSpecs = plot_name+"_selecSpectra"
        keep_prevFig = raw_input('Is there is a previous figure you want to keep?  (n)  y   ')
        # plot the smooth spectra
        if keep_prevFig == 'y':
            self.do_plot(plot_selecSpecs, img_format, used_specs, wf_arr[:,0], wf_arr[:,1], xlolim, xuplim, ylolim, yuplim, save=False)
        else:
            self.do_plot(plot_selecSpecs, img_format, used_specs, wf_arr[:,0], wf_arr[:,1], xlolim, xuplim, ylolim, yuplim, save=True)
            self.save_plot(plot_selecSpecs, img_format, used_specs, wf_arr[:,0], wf_arr[:,1], xlolim, xuplim, ylolim, yuplim)
        # used_specs is the list of the used spectra that will appear in the plot
        # self.wf_used_specs is an array that contains the arrays of wavs and fluxes of the used spectra
        return(self.wf_used_specs, used_specs)
        
    def find_overlapping_region(self, w_nuv, w_opt, w_nir, f_nuv, f_opt, f_nir, region):
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
    
    def plot_rebinned_spec(self, whole_spec_arr, full_name):
        pyplot.figure(1, figsize=(15, 8))
        pyplot.xlabel('Wavelength  [$\AA$]')
        pyplot.ylabel('Flux  [ergs/s/cm$^2$/$\AA$]')
        pyplot.xlim(1640.0, 10100.0)
        ymax = max(whole_spec_arr[1]) + max(whole_spec_arr[1])*0.12
        pyplot.ylim(-1.5e-15, ymax)
        fmax = max(whole_spec_arr[1]) - max(whole_spec_arr[1])*0.2 
        pyplot.text(8000, fmax, full_name, fontsize=18)
        pyplot.plot(whole_spec_arr[0], whole_spec_arr[1], 'k')
    
    def smooth_spec(self, plot_name, img_format, full_name):
        w_nuv = copy.deepcopy(self.wf_used_specs[0][0][0])
        f_nuv = copy.deepcopy(self.wf_used_specs[0][0][1])
        w_opt = copy.deepcopy(self.wf_used_specs[0][1][0])
        f_opt = copy.deepcopy(self.wf_used_specs[0][1][1])
        w_nir = copy.deepcopy(self.wf_used_specs[0][2][0])
        f_nir = copy.deepcopy(self.wf_used_specs[0][2][1])
        # show the entire spectra with smoothened overlapping regions
        # first smoothen the overlapping regions of the spectra and make a single nice spectrum to plot 
        nuv2opt_overlap_avg_wavs, nuv2opt_overlap_avg_flxs = self.find_overlapping_region(w_nuv, w_opt, w_nir, f_nuv, f_opt, f_nir, region='nuv2opt')
        opt2nir_overlap_avg_wavs, opt2nir_overlap_avg_flxs = self.find_overlapping_region(w_nuv, w_opt, w_nir, f_nuv, f_opt, f_nir, region='opt2nir')
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
        original_whole_spec_arr = numpy.array([whole_spec_wavs, whole_spec_flxs])
        # rebin the arrays to make a prettier spectra
        print 'Do you want to rebin?  If no press enter, if yes type rebinning factor:  '
        desired_rebin_factor = raw_input('(value < 1  will reduce the resolution,  0.4 generally gives a good "clean" spectrum)  ')
        if (desired_rebin_factor == '') or (desired_rebin_factor == 'n'):
            whole_spec_arr = original_whole_spec_arr
            self.plot_rebinned_spec(whole_spec_arr, full_name)
            pyplot.show()
        else:
            rebin_factor = (1, float(desired_rebin_factor))
            whole_spec_arr = spectrum.rebin(original_whole_spec_arr, rebin_factor)
            self.plot_rebinned_spec(whole_spec_arr, full_name)
            pyplot.show()
            rebin_again = raw_input('    Do you want to rebin again from the original spectrum?  (n)  y   ')
            end_while = False
            if (rebin_again == 'n') or (rebin_again == ''):
                end_while = True
            while end_while == False:
                print '    Type rebinning factor:  '
                desired_rebin_factor = raw_input('    (value < 1  will reduce the resolution,  0.4 generally gives a good "clean" spectrum) ')
                rebin_factor = (1, float(desired_rebin_factor))
                whole_spec_arr = spectrum.rebin(original_whole_spec_arr, rebin_factor)
                self.plot_rebinned_spec(whole_spec_arr, full_name)
                pyplot.show()
                rebin_again = raw_input('    Do you want to rebin again from the original spectrum?  (n)  y   ')                 
                if (rebin_again == 'n') or (rebin_again == ''):
                    end_while = True
        self.plot_rebinned_spec(whole_spec_arr, full_name)
        save_plt = raw_input('Do you want to save the smoothened whole spectrum image?  (n)  y   ')
        if save_plt == 'y':
            destination = plot_name+'_wholespec'+img_format
            pyplot.savefig(destination)
            print('Plot %s was saved!' % destination)
        elif save_plt == 'n':
            print('Plot not saved.')

        
        '''
        # Do the plot of the selected files, save the text file, and save the eps
        xlolim = self.xlolim
        xuplim = self.xuplim
        ylolim = self.ylolim
        yuplim = self.yuplim
        w = []
        f = []
        used_specs = []    
        txtout = open(selectedspecs_file, 'w+')
        write_file_wf = raw_input('Is there a wavelength-flux text file you want to keep? (n)  y  ')
        for i in self.indx_selected_wavs_and_flxs:
            wi = self.wavs[i]
            fi = self.flxs[i]
            w.append(wi)
            f.append(fi)
            used_specs.append(self.specs[i])
            # Writting the text file of wavelengths and fluxes of the selected spectra
            if write_file_wf == 'n' or write_file_wf == '':
                for i in range(0, len(wi)):
                    txtout.write("%.2f    %.3e \n" % (wi[i], fi[i]))
            # store the wavelengths and fluxes of the spectra selected in wf_used_specs
            wf = numpy.array([wi, fi])
            # wf_used_specs is a list of pairs of arrays of wavelengths and fluxes
            self.wf_used_specs.append(wf)       

        plot_name = plot_name+img_format
        plot_selecSpecs = plot_name.replace(img_format, "_selecSpectra"+img_format)
        keep_prevFig = raw_input('Is there is a previous figure you want to keep?  (n)  y   ')
        if keep_prevFig == 'y':
            self.do_plot(plot_selecSpecs, used_specs, w, f, xlolim, xuplim, ylolim, yuplim, save=False)
        else:
            self.do_plot(plot_selecSpecs, used_specs, w, f, xlolim, xuplim, ylolim, yuplim, save=True)
            self.save_plot(plot_selecSpecs, used_specs, w, f, xlolim, xuplim, ylolim, yuplim)        
        txtout.close()
        return(self.wf_used_specs)
        '''
        
    def make_splot_readable_files(self):
        pass


class BasicOps:
    '''
    This class gathers all the basic operations after we have the 1D spectra. These operations are:
    - reddening correction
    - underlying stellar absorption correction
    - line intensity measurement, equivalent widths, and FWHM
    '''
    def __init__(self, object_name, redlaw, cols_in_file, I_theo_HaHb, EWabsHbeta, cHbeta, av, ebv, 
                 path_and_name_RedCoroutfile, do_errs=None):
        self.object_name = object_name
        self.path_and_name_RedCoroutfile = path_and_name_RedCoroutfile
        if not isinstance(redlaw, str):
            print 'redlaw should be a string, got ', type(redlaw)
        if not isinstance(cols_in_file, list):
            print 'cols_in_file should be a list, got ', type(cols_in_file)
        if not isinstance(I_theo_HaHb, float):
            print 'I_theo_HaHb should be a float, got ', type(I_theo_HaHb)
        if not isinstance(EWabsHbeta, float):
            print 'EWabsHbeta should be a float, got ', type(EWabsHbeta)
        if not isinstance(EWabsHbeta, float):
            print 'EWabsHbeta should be a float, got ', type(EWabsHbeta)
        if not isinstance(cHbeta, float):
            print 'cHbeta should be a float, got ', type(cHbeta)
        if not isinstance(av, float):
            print 'av should be a float, got ', type(av)
        if not isinstance(ebv, float):
            print 'ebv should be a float, got ', type(ebv)
        self.redlaw = redlaw
        # Variables in cols_in_file: catalog_wavelength, observed_wavelength, element, ion, forbidden, how_forbidden, width, flux, continuum, EW
        self.catalog_wavelength = cols_in_file[0]
        self.observed_wavelength = cols_in_file[1]
        self.element = cols_in_file[2]
        self.ion = cols_in_file[3]
        self.forbidden = cols_in_file[4]
        self.how_forbidden = cols_in_file[5]
        self.width = cols_in_file[6]
        self.flux = cols_in_file[7]
        self.continuum = cols_in_file[8]
        self.EW = cols_in_file[9]
        self.I_theo_HaHb = I_theo_HaHb
        self.EWabsHbeta = EWabsHbeta
        self.cHbeta = cHbeta
        self.corr_undelyingAbs_EWs = []
        self.intensities_corr_undelyingAbs = []
        self.normfluxes = []
        self.Idered = []
        self.I_dered_norCorUndAbs = []
        self.av = av
        self.ebv = ebv
        # If dealing with uncertainties is desired, do_errs will be different than None
        self.errs_list = None
        if do_errs != None:
            self.errs_list = do_errs
        # Adopt predefined atomic data set
        pn.atomicData.getPredefinedDataFileDict()   # same as  pn.atomicData.getDataFile()
        #pn.atomicData.setDataFileDict('IRAF_09')   # Adopt IRAF atomic data set
        #pn.atomicData.resetDataFileDict            # Reset the atomic data to the predefined one
        #pn.atomicData.printDirForAllFiles()        # This prints all the paths for all the atomic data files
        #pn.atomicData.printAllSources()            # Returns the bibliographic reference of the source paper
        #raw_input()
        
    def underlyingAbsCorr(self, return_Is=False, EWabsHbeta=None):
        catalog_wavelength = self.catalog_wavelength
        continuum = self.continuum
        flux = self.flux
        if EWabsHbeta == None:
            EWabsHbeta = self.EWabsHbeta
        ### Step 1 of first iteration of reddening correction: Assume that there is NO collisional excitation
        ### get the EW_abs of the H and He lines with respect to EW_abs(Hbeta)
        Hline_and_EWs, Heline_and_EWs = spectrum.readlines_EWabsRelHbeta()
        ### Both H and He lines are going to be used to correct for underlying absorption
        undabs_wav = []
        undabs_EW = []
        for w, ew in zip(Hline_and_EWs[0], Hline_and_EWs[1]):
            undabs_wav.append(w)
            undabs_EW.append(ew)
        for w, ew in zip(Heline_and_EWs[0], Heline_and_EWs[1]):
            undabs_wav.append(w)
            undabs_EW.append(ew)
        #lines_undabs_and_EW = [undabs_wav, undabs_EW]
        ### Now add a 0.000 to all other observed lines
        corr_undelyingAbs_EWs = []
        for w in catalog_wavelength:
            if int(w) in undabs_wav:
                w_idx_undabs = undabs_wav.index(int(w))
                e = undabs_EW[w_idx_undabs]
            else:
                e = 0.000
            corr_undelyingAbs_EWs.append(e)
        ### Remove UNDERLYING ABSORPTION for optical lines to get Intensities
        corr_intensities = []
        for EWabsLine,cont,flx in zip(corr_undelyingAbs_EWs, continuum, flux):
            I = EWabsLine * EWabsHbeta * cont + flx
            corr_intensities.append(I)
        self.intensities_corr_undelyingAbs = corr_intensities
        if return_Is:
            return corr_intensities
    
    def Halpha2Hbeta_dered(self, av=None, ebv=None, cHbeta=None, fluxes=None, intensities=None):
        ''' Function to dered and obtain the Halpha/Hbeta ratio nicely printed along with the redden values (for comparison). '''
        redlaw = self.redlaw
        I_theo_HaHb = self.I_theo_HaHb
        catalog_wavelength = self.catalog_wavelength
        rounded_catalog_wavelength = round_wavs(catalog_wavelength)
        observed_wavelength = self.observed_wavelength
        element = self.element
        ion = self.ion
        forbidden = self.forbidden 
        how_forbidden = self.how_forbidden
        continuum = self.continuum
        EW = self.EW
        # These if statements are so that in the next reddening round different values than the initials are to be used.
        if av == None:
            ebv = 0.0
        if cHbeta == None:
            cHbeta = self.cHbeta
        if fluxes == None:
            flux = self.flux
        else:
            norm_fluxes = fluxes
        if intensities == None:
            intensities = self.intensities_corr_undelyingAbs
            norm_data = normalize_lines(rounded_catalog_wavelength, element, ion, forbidden,
                                    how_forbidden, observed_wavelength, flux, intensities, EW, continuum)
            catalog_lines, _, _, _, _, _, norm_fluxes, _, norm_intensities, _ = norm_data
        else:
            catalog_lines = rounded_catalog_wavelength
            norm_intensities = intensities
        # Normalized fluxes are the raw flux normalized to Hbeta. 
        # Intensities are already corrected for underlying absorption and are also normalized to Hbeta.
        Idered = []
        I_dered_noCorUndAbs = []
        # Print all available extinction laws
        #pn.RedCorr().printLaws()
        # or define a new one
        #RC.UserFunction = my_X
        #RC.law = 'user'
        # Plot the available laws
        #RC.plot(laws='all')
        #pyplot.show()        
        if ebv != 0.0:
            rv = av / ebv
            # Define a reddening correction object
            RC = pn.RedCorr(E_BV=ebv, R_V=rv, law=redlaw)
            # The relation between cHbeta and E(B-V) is:   E(B-V) = 0.61 * cHbeta + 0.024 * cHbeta ** 2.
            # Once one of the two parameters is defined, the other is directly obtained by:
            cHbeta = RC.cHbeta
            print 'THIS IS THE PYNEB OBTAINED C_Hbeta =', cHbeta/0.434
        else:
            # Obtain the reddening corrected intensities based on the given law and c(Hb)
            RC = pn.RedCorr(law=redlaw, cHbeta=cHbeta)
        for w, nF, nI in zip(catalog_lines, norm_fluxes, norm_intensities):   
            I_dered = nI * RC.getCorrHb(w)
            Idered.append(I_dered)
            # Obtain the reddening corrected intensities WITHOUT the correction due to 
            # underlying absorption based on the given law and c(Hb)
            IdnUA = nF * RC.getCorrHb(w)
            I_dered_noCorUndAbs.append(IdnUA)
            #print w, 'norm_Flux=', nI, '  dered_factor=',RC.getCorrHb(w) , ' dered_I=', I_dered
            #print '   FluxCorUndAbs', nF, '  I_dered_norCoUndAbs=', IdnUA
        # Find observed Halpha/Hbeta ratio
        Halpha_idx = catalog_lines.index(6563)
        Hbeta_idx = catalog_lines.index(4861)
        Halpha = norm_fluxes[Halpha_idx]
        Hbeta = norm_fluxes[Hbeta_idx]
        raw_ratio = Halpha/Hbeta
        I_Halpha = Idered[Halpha_idx]
        I_Hbeta = Idered[Hbeta_idx]
        I_obs_HaHb = I_Halpha/I_Hbeta
        print ''
        print 'C_Hbeta = %0.3f  -->  c_Hbeta = %0.3f' % (cHbeta/0.434, cHbeta) 
        print ' Intensities corrected for reddening and underlying absorption.'
        print '            Using', redlaw, '                   Normalized fluxes before extinction correction'
        print catalog_lines[Halpha_idx], '    ', I_Halpha, '                  ', Halpha
        print catalog_lines[Hbeta_idx], '    ', I_Hbeta, '                          ', Hbeta
        print 'theoretical ratio Ha/Hb = %0.2f' % (I_theo_HaHb)
        #print '      observed Ha/Hb = %0.3f           raw Ha/Hb = %0.3f' % (I_obs_HaHb, raw_ratio)
        print '      observed Ha/Hb =', numpy.round(I_obs_HaHb, decimals=2), '           raw Ha/Hb =', numpy.round(raw_ratio, decimals=2)
        self.Idered = Idered
        self.I_dered_norCorUndAbs = I_dered_noCorUndAbs
        self.normfluxes = norm_fluxes
        return self.normfluxes, self.Idered, self.I_dered_norCorUndAbs
    
    def get_uncertainties(self, data2process=None):
        ''' This is the function that actually runs the extintion law to deredd the lines. It also obtains the uncertainties.'''
        if data2process == None:
            wavs = self.catalog_wavelength
            flux = self.flux
            norm_flx = self.normfluxes
            norm_Idered = self.Idered
            errs_Flx, errs_EW, _ = self.errs_list   #errs_Flx, errs_EW, cont_errs = self.errs_list  
        else:
            wavs, flux, errs_Flx, norm_flx, norm_Idered, errs_EW = data2process
        errs_Idered = []
        perc_errs_I_dered = []
        errs_normfluxes = []
        perc_errs_normfluxes = []
        # Find observed Hbeta 
        Hbeta_idx = self.catalog_wavelength.index(4861.330)
        Hbeta = self.flux[Hbeta_idx]
        err_Hbeta = errs_Flx[Hbeta_idx]
        for w, F, eF, nF, nI in zip(wavs, flux, errs_Flx, norm_flx, norm_Idered):
            # Assuming that the errors are gaussian (which should be a good enough approximation for STIS)
            ''' the error of Flux/Hbeta is calculated using: 
            F(5007) = X +- deltaX
            F(Hbeta) = Y +- deltaY
            R = F(5007) / F(Hbeta) = X/Y
            deltaR**2 = (partial derivative of R with respect to X * deltaY)**2 + (partial derivative of R with respect to Y * deltaX)**2
            partial derivative of R with respect to X = 1/Y = R * 1/X
            partial derivative of R with respect to Y = X/Y * 1/Y = R * 1/Y
            deltaR**2 = R**2 ( (deltaY/X)**2 + (deltaX/Y)**2 ) 
            '''
            tot_err_nF = numpy.abs(nF) * numpy.sqrt((err_Hbeta/Hbeta)*(err_Hbeta/Hbeta) + (eF/F)*(eF/F) )
            errs_normfluxes.append(tot_err_nF)
            perc_tot_err_nF = (tot_err_nF * 100.) / numpy.abs(nF)
            perc_errs_normfluxes.append(perc_tot_err_nF)
            RC = pn.RedCorr()
            # Obtain the reddening corrected intensities based on the given law and c(Hb)
            if self.ebv != 0.0:
                rv = self.av / self.ebv
                RC = pn.RedCorr(E_BV=self.ebv, R_V=rv, law=self.redlaw, cHbeta=self.cHbeta)
            else:
                RC = pn.RedCorr(law=self.redlaw, cHbeta=self.cHbeta)
            errIdered = tot_err_nF * numpy.abs(RC.getCorrHb(w))
            errs_Idered.append(errIdered)
            perc_errIdered = (errIdered * 100.) / numpy.abs(nI)
            perc_errs_I_dered.append(perc_errIdered)
            #print w, ' NormFlux=', nF, ' err=',tot_err_nF, ' err%=', perc_tot_err_nF, '   I=', nI, ' err=', errIdered, ' err%=', perc_errIdered 
        #print 'min(perc_errs_I_dered)', min(perc_errs_I_dered)
        return errs_normfluxes, perc_errs_normfluxes, errs_Idered, perc_errs_I_dered, errs_EW
    
    def do_ops(self):
        self.underlyingAbsCorr()
        normfluxes, Idered, I_dered_norCorUndAbs = self.Halpha2Hbeta_dered(self.av, self.ebv)
        if self.errs_list != None:
            errs_normfluxes, perc_errs_normfluxes, errs_Idered, perc_errs_I_dered, errs_EW = self.get_uncertainties()
            flambdas = find_flambdas(self.cHbeta, self.catalog_wavelength, I_dered_norCorUndAbs, normfluxes)
            perrEW = []
            for ew, ewe in zip(self.EW, errs_EW):
                erp = ewe*100.0 / ew
                perrEW.append(erp)  
            cols_2write_in_file = [self.catalog_wavelength, flambdas, self.element, self.ion, self.forbidden, self.how_forbidden, normfluxes, errs_normfluxes, perc_errs_normfluxes, Idered, errs_Idered, perc_errs_I_dered, self.EW, errs_EW, perrEW]
            write_RedCorfile(self.object_name, self.path_and_name_RedCoroutfile, cols_2write_in_file)
            return normfluxes, Idered, I_dered_norCorUndAbs, errs_normfluxes, perc_errs_normfluxes, errs_Idered, perc_errs_I_dered
        else:
            return normfluxes, Idered, I_dered_norCorUndAbs

    
class AdvancedOps(BasicOps):    
    #def __init__(self, object_name, cHbeta, case, use_Chbeta, theoCE, writeouts=False, verbose=False):
    def __init__(self, object_name, redlaw, cols_in_file, I_theo_HaHb, EWabsHbeta, cHbeta, av, ebv, RedCor_file1, do_errs, #variables from parent class
                 case, use_Chbeta, theoCE, writeouts=False, verbose=False, tfile2ndRedCor=None):  #variables from child class
        # Initialize the inherited class
        BasicOps.__init__(self, object_name, redlaw, cols_in_file, I_theo_HaHb, EWabsHbeta, cHbeta, av, ebv, RedCor_file1, do_errs)
        # Inputs:
        self.cHbeta = cHbeta
        self.case = case                        # this is the Case used through out the class
        self.verbose = verbose                  # if True print midpoints in order to know what is the code working on
        self.writeouts = writeouts              # write or not text files with outputs (temps, densities, and abundances)
        # Variables defined in the class
        self.lines_pyneb_matches = []
        # Define the string added to the name of the text files to be written according to the reddening process used
        self.use_Chbeta = use_Chbeta
        if self.use_Chbeta:
            self.RedCorType = 'CHbeta'
        else:
            self.RedCorType = 'Ebv'
        self.used_measured_info = False
        if 'measuredLI' in RedCor_file1:
            self.used_measured_info = True
        # If performing the collisional excitation correction, tfile2ndRedCor is the name of the text file with the corrected intensities.
        self.tfile2ndRedCor = tfile2ndRedCor


    def writeRedCorrFile(self, verbose=None):
        ''' This function writes a text file in a pyneb readable format... Necessary to find the temperatures and densities. '''
        object_name = self.object_name
        if verbose == None:
            verbose = self.verbose
        input_file = object_name+'_RedCor_'+self.RedCorType+'.txt'
        if self.used_measured_info:
            input_file = object_name+'_measuredLI_RedCor_'+self.RedCorType+'.txt'
        path_object = '../results/'+object_name
        path_inputfile = os.path.join(path_object, input_file)
        # define the columns to be filled when reading the file
        self.wavelength = []
        self.flambda = []
        self.element = []
        self.ion = []
        self.forbidden = []
        self.howforb = []
        self.flux = []
        self.errflux = []
        self.percerrflx = []
        self.intensity = []
        self.errinten = []
        self.percerrinten = []
        self.ew = []
        self.errew = []
        self.percerrew = []
        AdvOpscols_in_file = [self.wavelength, self.flambda, self.element, self.ion, self.forbidden, self.howforb, self.flux, self.errflux, self.percerrflx, self.intensity, self.errinten, self.percerrinten, self.ew, self.errew, self.percerrew]
        # Read the info from the text file that contains IDs, dered info, and errors (this is the object_redCor.txt)
        self.AdvOpscols_in_file = spectrum.readlines_from_lineinfo(path_inputfile, AdvOpscols_in_file)
        if self.used_measured_info:
            self.pynebIDstxt = os.path.join(path_object, object_name+"_measuredLI_pynebIDs.txt")
        else:
            self.pynebIDstxt = os.path.join(path_object, object_name+"_pynebIDs.txt")
        tf = open(self.pynebIDstxt, 'w+')
        if verbose == True:
            print '{:<15} {:>15} {:>15}'.format('Ion_line', 'Intensity', 'Abs Error')
            print 'For  cHbeta = %0.3f' % self.cHbeta
        print >> tf, 'cHbeta = %0.3f' % self.cHbeta
        for w, el, io, Ic, er in zip(self.wavelength, self.element, self.ion, self.intensity, self.errinten):
            intw = int(numpy.round(w, decimals=0))
            cw = str(intw)
            el = str(el)
            io = str(io)
            wavid = cw+'A'
            pynebid = el+io
            if verbose == True:
                print w, el, io, Ic, er
                print 'pynebid =', pynebid, '    wavid =', wavid
            lineID = pynebid + '_' + wavid
            matching_line = [cw, Ic, er]
            if pynebid in pn.LINE_LABEL_LIST:
                if wavid in pn.LINE_LABEL_LIST[pynebid]:
                    print >> tf, '{:<15} {:>15.3f} {:>15.3f}'.format(lineID, Ic, er)
                    self.lines_pyneb_matches.append(matching_line)
                else:
                    pynebid = el+io+'_'+wavid+'+'
                    if pynebid in pn.BLEND_LIST:
                        print >> tf, '{:<15} {:>15.3f} {:>15.3f}'.format(pynebid, Ic, er)
                        self.lines_pyneb_matches.append(matching_line)
                    else:
                        continue
        tf.close()
        if verbose:
            print 'Pyneb IDs found!  Lines written in file  %s' % self.pynebIDstxt
        return self.lines_pyneb_matches 

    def get_tempsdens(self):
        ''' 
        This is a VERY big function that unfortunately cannot be broken... It determines ALL posible temperatures and densities 
        that are avaiable to veryfy with IRAF. The abundance determination can be turned off by setting iontotabs to False.
        '''
        RedCorType = self.RedCorType
        if self.used_measured_info:
            out_file = self.object_name+'_measuredLI_TempDens_'+RedCorType+'.txt'
        else:
            out_file = self.object_name+'_TempDens_'+RedCorType+'.txt'
        path_object = '../results/'+self.object_name
        fullpath_outfile = os.path.join(path_object, out_file)
        if self.writeouts:
            outf = open(fullpath_outfile, 'w+')
        # Determine first approximation of temperatures and densities
        # Define an Observation object and assign it to name 'obs'
        self.obs = pn.Observation()
        # read data from file created specifically for pyneb reading
        self.obs.readData(self.pynebIDstxt, fileFormat='lines_in_rows', corrected=True, errIsRelative=False)
        ###obs.readData(self.pynebIDstxt, fileFormat='lines_in_rows_err_cols', corrected=True)    # this option is not working
        # Intensities
        # O
        self.I_3727 = [0.0, 0.0]
        self.I_4959 = [0.0, 0.0]
        self.I_5007 = [0.0, 0.0]
        # He 1
        self.I_5876 = 0.0
        self.I_7065= 0.0                               
        self.I_4686 = 0.0 
        # C 3
        self.I_1661 = [0.0, 0.0]        
        self.I_1666 = [0.0, 0.0] 
        # N2               
        self.I_2144 = [0.0, 0.0] 
        for line in self.obs.lines:
            if self.verbose == True:            
                print 'line.wave', line.wave, '     line.corrIntens', line.corrIntens
            # Carbon
            if line.wave == 1907:           # C 3
                I_1907 = line.corrIntens
            elif line.wave == 1909:
                self.I_1909 = line.corrIntens 
            elif line.wave == 2326:         # C 2
                I_2326 = line.corrIntens
            elif line.wave == 2328:
                I_2328 = line.corrIntens 
            # Nitrogen
            elif line.wave == 1749:         # N 3
                self.I_1749 = line.corrIntens                      
            elif line.wave == 1752:
                self.I_1752 = line.corrIntens            
            elif line.wave == 2144:         # N 2
                self.I_2144 = line.corrIntens                      
            elif line.wave == 5755:
                I_5755 = line.corrIntens            
            elif line.wave == 6548:
                I_6548 = line.corrIntens            
            #elif line.wave == 6584:   # not using it because most of the time is blended with Halpha
            #    I_6584 = line.corrIntens            
            elif line.wave == 5198:         # N 1
                I_5198 = line.corrIntens                      
            elif line.wave == 5200:
                I_5200 = line.corrIntens                      
            # Oxygen 2 and 3
            elif line.wave == 1661:
                self.I_1661 = line.corrIntens
            elif line.wave == 1666:
                self.I_1666 = line.corrIntens
            elif line.wave == 4363:
                I_4363 = line.corrIntens
                print '4363 has an intensity of', I_4363
            elif line.wave == 4959:
                self.I_4959 = line.corrIntens
                print '4959 has an intensity of', self.I_4959
            elif line.wave == 5007:
                self.I_5007 = line.corrIntens
                print '5007 has an intensity of', self.I_5007
            elif line.wave == 3726:
                I_3726 = line.corrIntens
                print '3726 has an intensity of', I_3726
            elif line.wave == 3729:
                I_3729 = line.corrIntens
                print '3729 has an intensity of', I_3729
            elif line.wave == 3727:
                self.I_3727 = line.corrIntens
                print '3727 has an intensity of', self.I_3727
            #elif line.wave == 7320:
            #    I_7320 = line.corrIntens
            #    print '7320 has an intensity of', I_7320
            elif line.wave == 7330:
                I_7330 = line.corrIntens
                print '7330 has an intensity of', I_7330
            elif line.label == 'O2_7325+':
                I_7325 = line.corrIntens
                print '3725 has an intensity of', I_7325
            # Neon
            elif line.wave == 2975:         # Ne 5
                I_2975 = line.corrIntens            
            elif line.wave == 3346:
                I_3346 = line.corrIntens            
            #elif line.wave == 3426:
            #    I_3426 = line.corrIntens            
            elif line.wave == 2423:         # Ne 4
                I_2423 = line.corrIntens            
            elif line.wave == 2425:
                I_2425 = line.corrIntens            
            elif line.wave == 3342:         # Ne 3
                I_3342 = line.corrIntens            
            elif line.wave == 3869:
                I_3869 = line.corrIntens            
            #elif line.wave == 3969:
            #    I_3969 = line.corrIntens   
            # Aluminum
            elif line.wave == 2661:         # Al 2
                I_2661 = line.corrIntens            
            elif line.wave == 2670:
                I_2670 = line.corrIntens            
            # Silicon
            elif line.wave == 1883:         # Si 3
                I_1883 = line.corrIntens            
            elif line.wave == 1892:
                I_1892 = line.corrIntens            
            elif line.wave == 2335:         # Si 2
                I_2335 = line.corrIntens            
            elif line.wave == 2345:
                I_2345 = line.corrIntens            
            # Sodium
            elif line.wave == 2569:         # Na 6
                I_2569 = line.corrIntens                        
            elif line.wave == 2871:
                I_2871 = line.corrIntens                      
            elif line.wave == 2970:
                I_2970 = line.corrIntens                   
            elif line.wave == 2805:         # Na 4
                I_2805 = line.corrIntens                        
            elif line.wave == 3242:
                I_3242 = line.corrIntens                      
            elif line.wave == 3362:
                I_3362 = line.corrIntens  
            # Magnesium
            elif line.wave == 2418:         # Mg 5
                I_2418 = line.corrIntens                               
            elif line.wave == 2783:
                I_2783 = line.corrIntens  
            elif line.wave == 2928:
                I_2928 = line.corrIntens  
            # Sulphur
            #elif line.wave == 4068:         # S 2
            #    I_4068 = line.corrIntens    # not using it because it is the weakest of the auroral lines
            elif line.wave == 4076:
                I_4076 = line.corrIntens  
            elif line.wave == 6716:
                I_6716 = line.corrIntens  
            elif line.wave == 6731:
                I_6731 = line.corrIntens  
            elif line.wave == 6312:         # S 3
                I_6312 = line.corrIntens
                print '6312 has an intensity of', I_6312
            elif line.wave == 9069:
                I_9069 = line.corrIntens
                print '9069 has an intensity of', I_9069
            elif line.wave == 9531:
                I_9531 = line.corrIntens
                print '9531 has an intensity of', I_9531
            # Chlorine
            elif line.wave == 5518:         # Cl 3
                I_5518 = line.corrIntens
                print '5518 has an intensity of', I_5518
            elif line.wave == 5538:
                I_5538 = line.corrIntens
                print '5538 has an intensity of', I_5538
            # Argon
            elif line.wave == 4626:         # Ar 5
                I_4626= line.corrIntens                               
            elif line.wave == 6435:
                I_6435 = line.corrIntens  
            elif line.wave == 7006:
                I_7006 = line.corrIntens  
            #elif line.wave == 2854:         # Ar 4
            #    I_2854= line.corrIntens                               
            elif line.wave == 2868:
                I_2868 = line.corrIntens  
            #elif line.wave == 4711:
            #    I_4711 = line.corrIntens  
            elif line.wave == 4740:
                I_4740 = line.corrIntens  
            elif line.wave == 5192:         # Ar 3
                I_5192= line.corrIntens                               
            elif line.wave == 7136:
                I_7136 = line.corrIntens  
            elif line.wave == 7751:
                I_7751 = line.corrIntens  
            # Potasium
            elif line.wave == 2515:         # K 5
                I_2515= line.corrIntens                               
            elif line.wave == 2495:
                I_2495 = line.corrIntens  
            elif line.wave == 4123:
                I_4123 = line.corrIntens  
            elif line.wave == 4163:
                I_4163 = line.corrIntens  
            elif line.wave == 6223:
                I_6223 = line.corrIntens
            elif line.wave == 6349:
                I_6349 = line.corrIntens
            elif line.wave == 4511:         # K 4
                I_4511= line.corrIntens                               
            elif line.wave == 6102:
                I_6102 = line.corrIntens  
            elif line.wave == 6796:
                I_6796 = line.corrIntens  
            # Helium
            elif line.wave == 5876:         # He 1
                self.I_5876 = line.corrIntens
            elif line.wave == 7065:         # He 1
                self.I_7065= line.corrIntens                               
            elif line.wave == 4686:         # He 2
                self.I_4686 = line.corrIntens  
        # simultaneously compute temperature and density from pairs of line ratios
        # First of all, a Diagnostics object must be created and initialized with the relevant diagnostics.
        diags = pn.Diagnostics()   # Instantiate the Diagnostics class
        diags.getAllDiags()  # see what Diagnostics exist
        # temperature determination from an intensity ratio
        # explore some specific atom in the atoms collection
        # TEMPERATURES
        print ' \n *** FIRST estimations of TEMPERATURES:'
        print '   {:<8} {:<25} {:<10} {:<15} {:<15}'.format('Ion', 'Line Ratio', 'Ioniz Zone', 'Temp [K]', 'Temp+err')
        if self.writeouts:
            print >> outf, '#{:<8} {:<24} {:<14} {:<11} {:<11} {:<11}'.format('Ion', 'Line Ratio', 'Ioniz Zone', 'Temp/Dens', 'Temp/Dens+err', 'Error')
            print >> outf, '# FIRST estimations of TEMPERATURES:'
        self.TO3 = [0.0, 0.0]
        try:
            self.O3 = pn.Atom("O", "3")
            #self.tem_diag_O3 = '(L(4959)+L(5007)) / L(4363)'
            #print 'ratio of O3 = ', (I_4959+I_5007)/I_4363
            self.tem_diag_O3 = 'L(5007) / L(4363)'
            self.strongO3 = self.I_5007
            self.TO3 = self.O3.getTemDen(self.strongO3/I_4363, den=100., to_eval=self.tem_diag_O3)
            print '   {:<8} {:<25} {:<10} {:<15} {:<15}'.format('[O 3]','5007/4363', 'Medium', self.TO3[0], self.TO3[1])
            if self.writeouts:
                Terr = numpy.abs(self.TO3[0] - self.TO3[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[O 3]','5007/4363', 'Medium', self.TO3[0], self.TO3[1], Terr)
            #self.O3.plotGrotrian(tem=1e4, den=1e2, thresh_int=1e-3, unit = 'eV')
        except Exception as e:
            (NameError,),e
        try:
            self.N2 = pn.Atom("N", "2")
            tem_diag_N2 = 'L(6548) / L(5755)'    #'(L(6548)+L(6584)) / L(5755)'
            self.strongN2 = I_6548             #(I_6548 + I_6584)
            self.temN2 = self.N2.getTemDen(self.strongN2/I_5755, den=100.0, to_eval=tem_diag_N2) 
            print '   {:<8} {:<25} {:<10} {:<15} {:<15}'.format('[N 2]','6548/5755', 'Low', self.temN2[0], self.temN2[1])
            if self.writeouts:
                Terr = numpy.abs(self.temN2[0] - self.temN2[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[N 2]','6548/5755', 'Low', self.temN2[0], self.temN2[1], Terr)
        except Exception as e:
            (NameError,),e
        try:
            self.temO2 = [0.0, 0.0]
            self.O2 = pn.Atom("O", "2")
            tem_diag_O2 = '(L(3726)+L(3729)) / (L(7329) + L(7330))'
            self.temO2 = self.O2.getTemDen(self.I_3727/I_7330, den=100.0, to_eval=tem_diag_O2) 
            print '   {:<8} {:<25} {:<10} {:<15} {:<15}'.format('[O 2]','3727/7325', 'Low', self.temO2[0], self.temO2[1])
            if self.writeouts:
                Terr = numpy.abs(self.temO2[0] - self.temO2[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[O 2]','3727/7325', 'Low', self.temO2[0], self.temO2[1], Terr)
        except Exception as e:
            (NameError,),e
        try:
            self.Ne5 = pn.Atom("Ne", "5")
            tem_diag_Ne5 = 'L(3346) / L(2975)'  #'(L(3426)+L(3346)) / L(2975) '
            self.strongNe5 = I_3346             #(I_3426 + I_3346)
            self.temNe5 = self.Ne5.getTemDen(self.strongNe5/I_2975, den=100.0, to_eval=tem_diag_Ne5) 
            print '   {:<8} {:<25} {:<10} {:<15} {:<15}'.format('[Ne 5]','3346/2975', 'High', self.temNe5[0], self.temNe5[1])
            if self.writeouts:
                Te = self.temNe5
                Terr = numpy.abs(Te[0] - Te[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[Ne 5]','(3426+3346)/2975', 'High', Te[0], Te[1], Terr)
        except Exception as e:
            (NameError,),e
        try:
            self.Ne3 = pn.Atom("Ne", "3")
            tem_diag_Ne3 = 'L(3869)'    #'(L(3869)+L(3969)) / L(3342) '
            self.strongNe3 = I_3869     #(I_3869 + I_3969)
            self.temNe3 = self.Ne3.getTemDen(self.strongNe3/I_3342, den=100.0, to_eval=tem_diag_Ne3) 
            print '   {:<8} {:<25} {:<10} {:<15} {:<15}'.format('[Ne 3]','3869/3342', 'High', self.temNe3[0], self.temNe3[1])
            if self.writeouts:
                Te = self.temNe3
                Terr = numpy.abs(Te[0] - Te[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[Ne 3]','(3869+3969)/3342', 'High', Te[0], Te[1], Terr)
        except Exception as e:
            (NameError,),e
        try:
            self.Na6 = pn.Atom("Na", "6")
            tem_diag_Na6 = '(L(2871)+L(2970)) / L(2569) '
            self.strongNa6 = (I_2871 + I_2970)
            self.temNa6 = self.Na6.getTemDen(self.strongNa6/I_2569, den=100.0, to_eval=tem_diag_Na6) 
            print '   {:<8} {:<25} {:<10} {:<15} {:<15}'.format('[Na 6]','(2871+2970)/2569', 'High', self.temNa6[0], self.temNa6[1])
            if self.writeouts:
                Te = self.temNa6
                Terr = numpy.abs(Te[0] - Te[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[Na 6]','(2871+2970)/2569', 'High', Te[0], Te[1], Terr)
        except Exception as e:
            (NameError,),e
        try:
            self.Na4 = pn.Atom("Na", "4")
            tem_diag_Na4 = '(L(3242)+L(3362)) / L(2805) '
            self.strongNa4 = (I_3242 + I_3362)
            self.temNa4 = self.Na4.getTemDen(self.strongNa4/I_2805, den=100.0, to_eval=tem_diag_Na4) 
            print '   {:<8} {:<25} {:<10} {:<15} {:<15}'.format('[Na 4]','(3242+3362)/2805', 'Medium', self.temNa4[0], self.temNa4[1])
            if self.writeouts:
                Te = self.temNa4
                Terr = numpy.abs(Te[0] - Te[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[Na 4]','(3242+3362)/2805', 'Medium', Te[0], Te[1], Terr)
        except Exception as e:
            (NameError,),e
        try:
            self.Mg5 = pn.Atom("Mg", "5")
            tem_diag_Mg5 = '(L(2783)+L(2928)) / L(2418)'
            self.strongMg5 = I_2783     # (I_2783 + I_2928)
            self.temMg5 = self.Mg5.getTemDen((self.strongMg5+I_2928)/I_2418, den=100.0, to_eval=tem_diag_Mg5) 
            print '   {:<8} {:<25} {:<10} {:<15} {:<15}'.format('[Mg 5]','(2783+2928)/2418', 'High', self.temMg5[0], self.temMg5[1])
            if self.writeouts:
                Te = self.temMg5
                Terr = numpy.abs(Te[0] - Te[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[Mg 5]','(2783+2928)/2418', 'High', Te[0], Te[1], Terr)
        except Exception as e:
            (NameError,),e
        self.TS3 = [0.0, 0.0]
        try:
            self.S3 = pn.Atom("S", "3")
            S3ratio = I_6312 / I_9531 
            self.strongS3 = I_9531
            #print 'ratio of S3 = ', S3ratio
            self.TS3 = self.S3.getTemDen(S3ratio, den=100., wave1=6312, wave2=9531)
            print '   {:<8} {:<25} {:<10} {:<15} {:<15}'.format('[S 3]','6312/9531', 'High', self.TS3[0], self.TS3[1])
            if self.writeouts:
                Te = self.TS3
                Terr = numpy.abs(Te[0] - Te[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[S 3]','6312/9531', 'High', Te[0], Te[1], Terr)
        except Exception as e:
            (NameError,),e
        try:
            self.S2 = pn.Atom("S", "2")
            tem_diag_S2 = '(L(6716)+L(6731)) / L(4076)'    # Not using 4068 because is too weak
            self.strongS2 = (I_6716 + I_6731)
            self.temS2 = self.S2.getTemDen(self.strongS2/I_4076, den=100.0, to_eval=tem_diag_S2) 
            print '   {:<8} {:<25} {:<10} {:<15} {:<15}'.format('[S 2]','(6716+6731)/4076', 'Low', self.temS2[0], self.temS2[1])
            if self.writeouts:
                Te = self.temS2
                Terr = numpy.abs(Te[0] - Te[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[S 2]','(6716+6731)/4076', 'Low', Te[0], Te[1], Terr)
        except Exception as e:
            (NameError,),e
        try:
            self.Ar5 = pn.Atom("Ar", "5")
            tem_diag_Ar5 = '(L(6435) + L(7006)) / L(4626) '
            self.strongAr5 = (I_6435)
            self.temAr5 = self.Ar5.getTemDen((self.strongAr5 + I_7006)/I_4626, den=100.0, to_eval=tem_diag_Ar5) 
            print '   {:<8} {:<25} {:<10} {:<15} {:<15}'.format('[Ar 5]','(6435+7006)/4626', 'High', self.temAr5[0], self.temAr5[1])
            if self.writeouts:
                Te = self.temAr5
                Terr = numpy.abs(Te[0] - Te[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[Ar 5]','(6435+7006)/4626', 'High', Te[0], Te[1], Terr)
        except Exception as e:
            (NameError,),e
        try:
            self.Ar4 = pn.Atom("Ar", "4")    # not using 4711 and 2854 because they are VERY likely to be blended
            tem_diag_Ar4 = 'L(4740) / L(2868) '
            self.strongAr4 = I_4740
            self.temAr4 = self.Ar4.getTemDen(self.strongAr4/I_2868, den=100.0, to_eval=tem_diag_Ar4) 
            print '   {:<8} {:<25} {:<10} {:<15} {:<15}'.format('[Ar 4]','4740/2868', 'High', self.temAr4[0], self.temAr4[1])
            if self.writeouts:
                Te = self.temAr4
                Terr = numpy.abs(Te[0] - Te[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[Ar 4]','4740/2868', 'High', Te[0], Te[1], Terr)
        except Exception as e:
            (NameError,),e
        try:
            self.Ar3 = pn.Atom("Ar", "3")
            tem_diag_Ar3 = '(L(7136) + L(7751)) / L(5192) '
            self.strongAr3 = (I_7136 + I_7751)
            self.temAr3 = self.Ar3.getTemDen(self.strongAr3/I_5192, den=100.0, to_eval=tem_diag_Ar3) 
            print '   {:<8} {:<25} {:<10} {:<15} {:<15}'.format('[Ar 3]','(7136+7751)/5192', 'Medium', self.temAr3[0], self.temAr3[1])
            if self.writeouts:
                Te = self.temAr3
                Terr = numpy.abs(Te[0] - Te[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[Ar 3]','(7136+7751)/5192', 'Medium', Te[0], Te[1], Terr)
        except Exception as e:
            (NameError,),e
        try:
            self.K5 = pn.Atom("K", "5")
            tem_diag_K5 = '(L(4123)+L(4163)) / (L(2515) + L(2495))'
            self.strongK5 = (I_4123 + I_4163)
            self.temK5 = self.K5.getTemDen(self.strongK5/(I_2515 + I_2495), den=100.0, to_eval=tem_diag_K5) 
            print '   {:<8} {:<25} {:<10} {:<15} {:<15}'.format('[K 5]','(4123+4163)/(2515+2495)', 'High', self.temK5[0], self.temK5[1])
            if self.writeouts:
                Te = self.temK5
                Terr = numpy.abs(Te[0] - Te[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[K 5]','(4123+4163)/(2515+2495)', 'High', Te[0], Te[1], Terr)
        except Exception as e:
            (NameError,),e
        try:
            self.K4 = pn.Atom("K", "4")
            tem_diag_K4 = '(L(6102)+L(6796)) / L(4511)'
            self.strongK4 = (I_6102 + I_6796)
            self.temK4 = self.K4.getTemDen(self.strongK4/I_4511, den=100.0, to_eval=tem_diag_K4) 
            print '   {:<8} {:<25} {:<10} {:<15} {:<15}'.format('[K 4]','(6102+6796)/4511', 'High', self.temK4[0], self.temK4[1])
            if self.writeouts:
                Te = self.temK4
                Terr = numpy.abs(Te[0] - Te[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[K 4]','(6102+6796)/4511', 'High', Te[0], Te[1], Terr)
        except Exception as e:
            (NameError,),e
            
        # DENSITIES
        print ' \n *** FIRST estimations of DENSITIES:'
        print '   {:<8} {:<12} {:<10} {:<15} {:<15}'.format('Ion', 'Line Ratio', 'Ioniz Zone', 'Density [cm^-3]', 'Dens+err')
        if self.writeouts:
            print >> outf, '#'
            print >> outf, '# FIRST estimations of DENSITIES:'
            #print >> outf, '#{:<8} {:<24} {:<14} {:<11} {:<11}'.format('Ion', 'Line Ratio', 'Ioniz Zone', 'Den [cm^-3]', 'Den Err [cm^-3]')
        try:
            self.C3 = pn.Atom("C", "3")
            den_diag_C3 = 'L(1907) / L(1909)'
            self.strongC3 = I_1907
            self.denC3 = self.C3.getTemDen(I_1907/self.I_1909, tem=10000.0, to_eval=den_diag_C3) 
            print '   {:<8} {:<12} {:<10} {:<15} {:<15}'.format('C 3]','1907/1909', 'Medium', self.denC3[0], self.denC3[1])
            if self.writeouts:
                den = self.denC3
                Derr = numpy.abs(den[0] - den[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('C 3]','1907/1909', 'Medium', den[0], den[1], Derr)
        except Exception as e:
            (NameError,),e        
        try:
            self.C2 = pn.Atom("C", "2")
            den_diag_C2 = 'L(2326) / L(2328)'
            self.strongC2 = I_2328
            self.denC2 = self.C2.getTemDen(I_2326/I_2328, tem=10000.0, to_eval=den_diag_C2) 
            print '   {:<8} {:<12} {:<10} {:<15} {:<15}'.format('C 2]','2326/2328', 'Medium', self.denC2[0], self.denC2[1])
            if self.writeouts:
                den = self.denC2
                Derr = numpy.abs(den[0] - den[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('C 2]','2326/2328', 'Medium', den[0], den[1], Derr)
        except Exception as e:
            (NameError,),e        
        try:
            self.N3 = pn.Atom("N", "3")
            den_diag_N3 = 'L(1749) / L(1752)'
            self.strongN3 = self.I_1752
            self.denN3 = self.N3.getTemDen(self.I_1749/self.I_1752, tem=10000.0, to_eval=den_diag_N3) 
            print '   {:<8} {:<12} {:<10} {:<15} {:<15}'.format('N 3]','1749/1752', 'Medium', self.denN3[0], self.denN3[1])
            if self.writeouts:
                den = self.denN3
                Derr = numpy.abs(den[0] - den[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('N 3]','1749/1752', 'Medium', den[0], den[1], Derr)
        except Exception as e:
            (NameError,),e        
        try:
            self.N1 = pn.Atom("N", "1")
            den_diag_N1 = 'L(5198) / L(5200)'
            self.strongN1 = I_5200
            self.denN1 = self.N1.getTemDen(I_5198/I_5200, tem=10000.0, to_eval=den_diag_N1) 
            print '   {:<8} {:<12} {:<10} {:<15} {:<15}'.format('[N 1]','5198/5200', 'Low', self.denN1[0], self.denN1[1])
            if self.writeouts:
                den = self.denN1
                Derr = numpy.abs(den[0] - den[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[N 1]','5198/5200', 'Low', den[0], den[1], Derr)
        except Exception as e:
            (NameError,),e        
        try:
            O2ratio = 'L(3729) / L(3726)'
            self.strongO2 = I_3729
            self.denO2 = self.O2.getTemDen(I_3729 / I_3726, tem=10000.0, to_eval=O2ratio) 
            print '   {:<8} {:<12} {:<10} {:<15}'.format('[O 2]','3729/3726', 'Medium', self.denO2[0], self.denO2[1])
            if self.writeouts:
                den = self.denO2
                Derr = numpy.abs(den[0] - den[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[O 2]','3729/3726', 'Medium', den[0], den[1], Derr)
        except Exception as e:
            (NameError,),e
        try:
            self.Ne4 = pn.Atom("Ne", "4")
            den_diag_Ne4 = 'L(2423) / L(2425)'
            self.strongNe4 = I_2425
            self.denNe4 = self.Ne4.getTemDen(I_2423/I_2425, tem=10000.0, to_eval=den_diag_Ne4) 
            print '   {:<8} {:<12} {:<10} {:<15} {:<15}'.format('[Ne 4]','2423/2425', 'High', self.denNe4[0], self.denNe4[1])
            if self.writeouts:
                den = self.denNe4
                Derr = numpy.abs(den[0] - den[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[Ne 4]','2423/2425', 'High', den[0], den[1], Derr)
        except Exception as e:
            (NameError,),e        
        try:
            self.Al2 = pn.Atom("Al", "2")
            den_diag_Al2 = 'L(2661) / L(2670)'
            self.strongAl2 = I_2670
            self.denAl2 = self.Al2.getTemDen(I_2661/I_2670, tem=10000.0, to_eval=den_diag_Al2) 
            print '   {:<8} {:<12} {:<10} {:<15} {:<15}'.format('[Al 2]','2661/2670', 'Low', self.denAl2[0], self.denAl2[1])
            if self.writeouts:
                den = self.denAl2
                Derr = numpy.abs(den[0] - den[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[Al 2]','2661/2670', 'Low', den[0], den[1], Derr)
        except Exception as e:
            (NameError,),e        
        try:
            self.Si3 = pn.Atom("Si", "3")
            den_diag_Si3 = 'L(1883) / L(1892)'
            self.strongSi3 = I_1883
            self.denSi3 = self.Si3.getTemDen(I_1883/I_1892, tem=10000.0, to_eval=den_diag_Si3) 
            print '   {:<8} {:<12} {:<10} {:<15} {:<15}'.format('Si 3]','1883/1892', 'Low', self.denSi3[0], self.denSi3[1])
            if self.writeouts:
                den = self.denSi3
                Derr = numpy.abs(den[0] - den[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('Si 3]','1883/1892', 'Low', den[0], den[1], Derr)
        except Exception as e:
            (NameError,),e        
        try:
            self.Si2 = pn.Atom("Si", "2")
            den_diag_Si2 = 'L(2335) / L(2345)'
            self.strongSi2 = I_2335
            self.denSi2 = self.Si2.getTemDen(I_2335/I_2345, tem=10000.0, to_eval=den_diag_Si2) 
            print '   {:<8} {:<12} {:<10} {:<15} {:<15}'.format('[Si 2]','2335/2345', 'Low', self.denSi2[0], self.denSi2[1])
            if self.writeouts:
                den = self.denSi2
                Derr = numpy.abs(den[0] - den[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('Si 2]','2335/2345', 'Low', den[0], den[1], Derr)
        except Exception as e:
            (NameError,),e        
        try:
            self.denS2 = [0.0, 0.0]
            den_diag_S2 = 'L(6716) / L(6731)'
            self.denS2 = self.S2.getTemDen(I_6716/I_6731, tem=10000.0, to_eval=den_diag_S2) 
            print '   {:<8} {:<12} {:<10} {:<15} {:<15}'.format('[S 2]','6716/6731', 'Low', self.denS2[0], self.denS2[1])
            if self.writeouts:
                den = self.denS2
                Derr = numpy.abs(den[0] - den[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[S 2]','6716/6731', 'Low', den[0], den[1], Derr)
        except Exception as e:
            (NameError,),e        
        try:
            self.Cl3 = pn.Atom("Cl", "3")
            Cl3ratio = I_5538 / I_5518 
            self.strongCl3 = I_5518
            self.dCl3 = self.Cl3.getTemDen(Cl3ratio, tem=10000.0, wave1=5538, wave2=5518)
            print '   {:<8} {:<12} {:<10} {:<15} {:<15}'.format('[Cl 3]','5538/5518', 'Medium', self.dCl3[0], self.dCl3[1])
            if self.writeouts:
                den = self.dCl3
                Derr = numpy.abs(den[0] - den[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[Cl 3]','5538/5518', 'Medium', den[0], den[1], Derr)
        except Exception as e:
            (NameError,),e
        try:
            den_diag_K5 = 'L(6223) / L(6349)'
            self.denK5 = self.K5.getTemDen(I_6223/I_6349, tem=10000.0, to_eval=den_diag_K5) 
            print '   {:<8} {:<12} {:<10} {:<15} {:<15}'.format('[K 5]','6223/6349', 'High', self.denK5[0], self.denK5[1])
            if self.writeouts:
                den = self.denK5
                Derr = numpy.abs(den[0] - den[1])
                print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('[K 5]','6223/6349', 'High', den[0], den[1], Derr)
        except Exception as e:
            (NameError,),e
            
        ### Density measurement from [Fe III] lines -- taken from Peimbert, Pena-Guerrero, Peimbert (2012, ApJ, 753, 39)
        if self.TO3[0] != 0.0:
            if not math.isnan(self.TO3[0]):
                # Iron
                I4986 = 0.0
                I4987 = 0.0
                I4658 = 0.0
                if line.wave == 4986:         # Fe 3
                    I4986 = line.corrIntens                               
                elif line.wave == 4987:
                    I4987 = line.corrIntens  
                elif line.wave == 4658:
                    I4658 = line.corrIntens  
                if self.writeouts:
                    print >> outf, '# '
                    print >> outf, '# Theoretically obtained values'
                if (I4986 != 0.0) and (I4987 != 0) and (I4658 !=0):
                    log_Fe3den = 2 - ( (numpy.log10((I4986+I4987)/I4658) - 0.05 - 0.25*(numpy.log10(self.TO3[0]-4))) / (0.66 - 0.18*(numpy.log10(self.TO3[0]-4))) )
                    self.Fe3den = 10**(log_Fe3den)
                    print '   Density measured from [Fe3] lines:', self.Fe3den
                    if self.writeouts:
                        den = self.Fe3den
                        Derr = numpy.abs(den[0] - den[1])
                        print >> outf,'{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('ne[Fe 3]','Peimbert et al 2012', 'High', den[0], den[1], Derr)
                    self.densities.append(den)
                else:
                    print '   {:<8} {:<12} {:<10} {:<15} {:<15}'.format('[Fe 3]','Peimbert et al 2012', 'High','nan', 'nan', 'nan')
                    if self.writeouts:
                        print >> outf, '{:<8} {:<25} {:<14} {:<11} {:<11} {:<11}'.format('ne[Fe 3]','Peimbert et al 2012', 'High', 'nan', 'nan', 'nan')
    
                ### Following analysis presented in Pena-Guerrero, Peimbert, Peimbert, Ruiz (2012, ApJ, 746, 115) and
                ### Peimbert, Pena-Guerrero, Peimbert (2012, ApJ, 753, 39).
                # If the [O II] temperature was not obtained directly from observations, get an estimate temperature of OII from OII:
                # 1) using equation of Peimbert, Peimbert, & Luridiana (2002, ApJ, 565, 668) - Based on data of Stasinska's models.
                self.TO2pei = 2430. + self.TO3 * (1.031 - self.TO3/54350.)
                TO2pei_err = numpy.abs(self.TO2pei[0] - self.TO2pei[1])
                print 'This is the theoretically obtained temperature of O2 from Peimbert et al. 2002 = ', self.TO2pei
                if self.writeouts:
                    print >> outf, '{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('Te[O 2]','Peimbert et al 2002', 'Low', self.TO2pei[0], self.TO2pei[1], TO2pei_err)
                print ' * this theoretical relation works fine if Te[OIII] > 12,000'
                print ' * for comparison, from observations Te[O III] = ', self.TO3
                # 2) using equation of Garnett, D. R. 1992, AJ, 103, 1330
                self.TO2gar = 0.7 * self.TO3 + 3000.
                TO2gar_err = numpy.abs(self.TO2gar[0] - self.TO2gar[1])
                print 'Theoretically obtained temperature of O2 from Garnet 1992 = ', self.TO2gar
                if self.writeouts:
                    print >> outf, '{:<8} {:<25} {:<14} {:<11.2f} {:<11.2f} {:<10.2f}'.format('Te[O 2]','Garnett 1992', 'Low', self.TO2gar[0], self.TO2gar[1], TO2gar_err)
            else:
                self.TO2gar = [0.0, 0.0]
        # Make sure that the temperatures and densities file closes properly
        if self.writeouts:
            outf.close()
            
    def define_TeNe_HighLow_andVLow(self, forceTeH=None, forceTeL=None, forceNe=None):
        print '\n These are the temperatures and density that will be used to determine abundances...'
        if (forceTeH != None):
            if type(forceTeH) is not list:
                if (forceTeL != None) and type(forceTeL) is not list:
                    teL = forceTeL
                else:
                    teHgarnett = forceTeH*0.7 + 3000.0     # from Garnett(1992) 
                    teL = teHgarnett  
                tLfromTS3 = forceTeH*0.83 + 1700.0     # from Garnett(1992), if TS3 is measured 
                # error is not defined, use a default value of 10% of the given temperature
                perr = 0.1
                teH = forceTeH
                teHerr = teH + (teH * perr)
                teVL = tLfromTS3  
            else:
                teH = forceTeH[0]
                if (forceTeL != None) and type(forceTeL) is list:
                    teL = forceTeL[0]
                else:                
                    teHgarnett = teH*0.7 + 3000.0     # from Garnett(1992) 
                    teL = teHgarnett
                tLfromTS3 = teH*0.83 + 1700.0     # from Garnett(1992), if TS3 is measured 
                # error is defined, determine the percentage to make it absolute error
                if forceTeH[1] < forceTeH[0]:
                    abserr = forceTeH[1]
                else:
                    abserr = forceTeH[1] - teH
                perr = (abserr * 1.0) / teH
                teHerr = teH + abserr
                #teL = forceTeH[0] * 0.85         # 85% of the high temperature
                #teVL = forceTeH[0] * 0.75        # 75% of the high temperature
                teVL = tLfromTS3  
            te_high = [teH, teHerr]
            teLerr = teL + (teL*perr)
            te_low = [teL, teLerr]
            teVLerr = teVL + (teVL * perr)
            te_verylow = [teVL, teVLerr]
            print '   High ionization degree temperature forced to be:            ', te_high[0], '+-', te_high[1]-te_high[0]
            print '   Low ionization degree temperature = t_high*0.7 + 3000.0:    ', te_low[0], '+-', te_low[1]-te_low[0]
            print '   Very Low ionization degree temperature = t_high*0.85+ 1700.0', te_verylow[0], '+-', te_verylow[1]-te_verylow[0]
            print '   * Error is calculated from percentage error in Te_high:    %0.2f' % (perr*100.0), '%'
        else:
            te_high = [10000.0, 11000.0]
            #te_low = [9000.000, 10000.0]
            teHgarnett = te_high[0]*0.7 + 3000.0     # from Garnett(1992) 
            te_low = [teHgarnett, teHgarnett+1000.0]
            #te_verylow = [8500.0, 9500.0]
            tLfromTS3 = te_high[0]*0.83 + 1700.0     # from Garnett(1992), if TS3 is measured 
            te_verylow = [tLfromTS3, tLfromTS3+1000.0]
            if math.isnan(self.temO2[0]):
                print '   Te[O 2] not available, using default value from Garnett (1992): t_high*0.7 + 3000.0'
            elif self.temO2[0] != 0.0:
                te_low = self.temO2
                print '   Te_low (O2) =', te_low
            else:
                te_low = [teHgarnett, teHgarnett+1000.0]
            if math.isnan(self.TO3[0]):
                if math.isnan(self.temO2[0]):
                    print '   Te[O 3] not available, using default value:   te_high = 10,000 +- 1000.0'
            else:
                if (self.TO3[0] > 20000.0):
                    te_hi = te_low[0] + (te_low[0] * 0.1)
                    err = te_hi*0.3
                    te_high = [te_hi, err]
                    print '   Te[O 3] is not available but te_low=', te_low,', using:   te_high = te_low+(0.1*te_low) +- 0.2*errte_low'
                else:
                    te_high = self.TO3
                    print '   Te_high (O3) =', te_high
                # make sure the low ionization temperature has a reasonable value compared to the high one
                diff_temps = te_low[0] - te_high[0]
                if (diff_temps > 3000.0) or math.isnan(self.temO2[0]):
                    print '   Te_low (O2-Garnet92) =', self.TO2gar
                    te_low = self.TO2gar
            if math.isnan(self.TS3[0]):
                print '   Te[S 3] not available, using default value from Garnett (1992): t_high*0.83 + 1700.0'
            elif self.TS3[0] <= self.TO3[0]:
                    te_verylow = self.TS3
                    print '   Te_Verylow (S3) = ', te_verylow
            elif self.TS3[0] > self.TO3[0]:
                #teVlow = self.TO3[0]*0.85
                teVlow = tLfromTS3
                te_verylow = [teVlow, teVlow+2000.0]
                print '   Te[S 3] > Te[O 3], using default value:   te_Verylow = Te[O 3]*0.85 +-  % err of Te[S 3] =', te_verylow
        print 'Densitiy being used for calculations:'
        if forceNe != None:
            if type(forceNe) is not list:
                # error is not defined, use a default value of 50% of the given temperature
                perr = 0.50
                ne = forceNe
            else:
                ne = forceNe[0]
                if forceNe[1] < forceNe[0]:
                    abserr = forceNe[1]
                else:
                    abserr = forceNe[1] - forceNe[0]
                perr = (abserr * 1.0) / ne
            dens_err = ne + (ne * perr)
            dens = [ne, dens_err]
            print '   Density forced to be: ', dens[0], '+-', dens[1]-dens[0], '\n'
        else:        
            dens = [150.0, 200.0]
            if math.isnan(self.denS2[0]):
                print 'ne[S 2] not available, using default value: 150.0 +- 50.0 \n'
            elif self.denS2[0] != 0.0:
                #ne = self.denO2
                ne = self.denS2
                if numpy.abs(ne[1] - ne[0]) < 20.0:
                    ne[1] = ne[0] + ne[0]*0.1
                dens = ne
                print 'ne =', dens, '\n'
        print ''
        # Define the temperatures and density for the class. 
        self.te_high = te_high        
        self.te_low = te_low        
        self.te_verylow = te_verylow 
        self.dens = dens         
        
    def get_iontotabs(self, data2use=None):
        ''' With the available temperatures determine ionic and total abundances
        data2use = is only different to None when correcting for collisional excitation. It is a list of
                    [lines_info, dereddening_info, uncertainties_info]  where
                    lines_info = [catalog_lines, wavs_lines, element_lines, ion_lines, forbidden_lines,  
                                  how_forbidden_lines, norm_fluxes, norm_intensities, EW_lines]
                    dereddening_info = [EWabsHbeta, C_Hbeta, norm_Idered, I_dered_norCorUndAbs, flambdas]
                    uncertainties_info = [percent_Iuncert, absolute_Iuncert, S2N]
        '''
        # Define the string added to the name of the text files to be written according to the reddening process used
        RedCorType = self.RedCorType
        # Define the high and lo ionization zones temperatures and densities
        print 'Temperatures being used for estimation of abundances:'
        te_high = self.te_high        
        te_low = self.te_low        
        te_verylow = self.te_verylow       
        dens = self.dens         
        print 'te_high =', te_high[0], '+-', te_high[1] - te_high[0]
        print 'te_low =', te_low[0], '+-',te_low[1] - te_low[0]
        print 'te_verylow = ', te_verylow[0], '+-', te_verylow[1] - te_verylow[0]
        print 'density =', dens[0], '+-', dens[1] - dens[0]
        #raw_input(' ***  press enter to continue')
                
        print '\n  Calculating abundances.... \n'
            
        # Define all atoms to make calculations
        all_atoms = pn.getAtomDict()
        #print 'all_atoms', all_atoms
        # Determine all available abundances -- THIS IS FOR COMPARISON PURPOSES ONLY
        ab_high_list = []
        ab_low_list = []
        ab_very_low = []
        line_label_list = []
        line_I_list = []
        try:
            for line in self.obs.lines:
                if line.atom in all_atoms:
                    abH = all_atoms[line.atom].getIonAbundance(line.corrIntens, te_high, dens, to_eval=line.to_eval)
                    abL = all_atoms[line.atom].getIonAbundance(line.corrIntens, te_low, dens, to_eval=line.to_eval)
                    abVL = all_atoms[line.atom].getIonAbundance(line.corrIntens, te_verylow, dens, to_eval=line.to_eval)
                    ab_high_list.append(abH)
                    ab_low_list.append(abL)
                    ab_very_low.append(abVL)
                    line_label_list.append(line.label)
                    line_I_list.append(line.corrIntens)
                    #print line.label, abH, abL
                else:
                    pn.log_.warn('line from %s not used because ion not found' % line.atom, calling='full_analysis.py')
            pn.log_.timer('Ending full_analysis.py', calling='full_analysis.py')
        except (RuntimeError, TypeError, NameError):
            pass
         
        # recombination atoms
        He1 = pn.RecAtom('He', 1)
        He2 = pn.RecAtom('He', 2)
        pn.atomicData.getDataFile(data_type='rec')
        print 'HELIUM ABUNDANCES...'
        try:
            abhe1 = He1.getIonAbundance(int_ratio=self.I_5876, tem=te_high, den=dens, wave=5876)
            print 'He1 abundance =', abhe1
        except (RuntimeError, TypeError, NameError):
            print 'Could not find He1 abundance.'
            pass
        try:
            abhe2 = He1.getIonAbundance(int_ratio=self.I_7065, tem=te_high, den=dens, wave=7065)
            print 'He1 abundance =', abhe2
        except (RuntimeError, TypeError, NameError):
            print 'Could not find He1 abundance.'
            pass
        try:
            abhe3 = He2.getIonAbundance(int_ratio=self.I_4686, tem=te_high, den=dens, wave=4686)
            print 'He2 abundance =', abhe3
        except (RuntimeError, TypeError, NameError):
            print 'Could not find He2 abundance.'
            pass
        #raw_input(' ***  press enter to continue')
        
        # ions of zones of high and medium ionization degree are combined
        ab_high = ['Ar3', 'Ar4', 'Ar5', 'C2', 'C3', 'Ca5', 'Cl2', 'Cl3', 'Cl4', 'Fe3', 'K4', 'K5', 'Mg5', 
                   'N3', 'Ne3', 'Ne4', 'Ne5', 'Na4', 'Na6', 'Ne3', 'Ne5', 'O3', 'S3', 'He2']
        ab_low = ['Al2', 'N1', 'N2', 'O1', 'O2', 'S2', 'Si2', 'Si3', 'He1']
        
        # create two lists: 1) all atoms list, 2) empty list with all ions to be filled with the total ionic abundances
        atoms_list = []
        totabs_ions_list = []
        for ah in ab_high:
            atoms_list.append(ah)
            empty_list_for_ionab_and_error = []
            totabs_ions_list.append(empty_list_for_ionab_and_error)
        for al in ab_low:
            atoms_list.append(al)
            empty_list_for_ionab_and_error = []
            totabs_ions_list.append(empty_list_for_ionab_and_error)
        # sort those atoms alphabetically and from low to high ionization degree
        sorted_atoms = sorted(set(atoms_list))
        
        # Get abundances only from strong (or more easily measurable) lines and appropriate temperature
        #                     Al2     Ar3     Ar4     Ar5     C2      C3      Ca5     Cl2     Cl3     Cl4     Fe3 
        strong_lines_list = ['2670', '7751', '4740', '6435', '2328', '1907', '6087', '9124', '5518', '7531', '4987',
                             # K4     K5      Mg5     N1      N2      N3      Na4     Na6     Ne3     Ne4     Ne5    
                             '6796', '4163', '2783', '5200', '6548', '1752', '3362', '2970', '3869', '2425', '3346',
                             # O1     O2     O3      S2      
                             '6300', '3727', '5007', '6731',#'6716', 
                             # S3      Si2     Si3
                             '9531', '2345', '1892',
                             # He1     He2
                             '7065', '4686'] 
        for label in line_label_list:
            for line in strong_lines_list:
                if line in label:
                    i = line_label_list.index(label)
                    kk = label.split("_")
                    ion = kk[0]
                    if ion in ab_high:
                        ab = ab_high_list[i]
                    elif ion in ab_low:
                        ab = ab_low_list[i]
                    idx = sorted_atoms.index(ion)
                    totabs_ions_list[idx] = ab
                    #print ion, line, ab, sorted_atoms[idx], idx
        # add a 0.0 where no abundances could be determined
        for ion, ab in zip(sorted_atoms, totabs_ions_list):
            if len(ab)<2:
                idx = sorted_atoms.index(ion)
                totabs_ions_list[idx] = [0.0, 0.0]
                logab = 0.0
            else:
                logab = 12+numpy.log10(ab[0])
            print ion, ab, '    12+log(X+i/H+) =', logab
        
        # Define the dictionary with the elements that I want to have total abundances for
        #elements = ['Ar', 'Cl', 'N', 'Ne', 'O', 'S',   # using pyneb ICF with traditional way
        #            'C', 'He',                         # using different ICF not in pyneb 
        #            'Al', 'Ca', 'K', 'Mg', 'Na','Si']  # pyneb does not have ICF for HII regions for these elements
        elem_abun = OrderedDict()
        
        ### TOTAL abundances with ICFs in pyneb
        icf = pn.ICF()
        icf.getAvailableICFs() # get all available ICFs
        #icf.printAllICFs(type_=['HII']) # print out the ICFs only for HII regions
        #r = icf.getReference('Ial06_22a')
        self.atom_abun = OrderedDict()
        #atom_abun = {}
        for ion, ionabund in zip(sorted_atoms, totabs_ions_list):
            self.atom_abun[ion] = ionabund
        #print self.atom_abun
        # use a specific recipy to determine abundances
        #elem_abunTPP85 = icf.getElemAbundance(self.atom_abun, icf_list=['TPP85']) 
        #elem_abunIal06_22a = icf.getElemAbundance(self.atom_abun, icf_list=['Ial06_22a'])
        #icf.getExpression('Ial06_22a') # returns the analytical expression of the icf identified by the label TPP85
        #icf.getReference('Ial06_22a') # returns the bibliographic reference of the source paper
        #icf.getURL('Ial06_22a') # returns the ADS URL of the source paper
        #icf.printInfo('Ial06_20a')
        #print 'elem_abunIal06_22a', elem_abunIal06_22a#['Ial06_22a']
        
        # Oxygen
        print '\n OXYGEN'
        Otot = self.atom_abun['O2'][0] + self.atom_abun['O3'][0]
        O23sq = self.atom_abun['O2'][1]**2 + self.atom_abun['O3'][1]**2
        #O23sq = (Otot * 0.35)**2
        Ototerr = numpy.sqrt(O23sq)
        O_errp = (Ototerr*100) / Otot
        #print '    absOtot = ', Otot, '+-', Ototerr, ', which is', O_errp,'% of error'
        print '    absOtot = %0.3e +- %0.3e (~%0.3f percent)' % (Otot, Ototerr, O_errp)
        print '    Assuming that O+++ contributes less than 1% to Otot, hence Otot = O+ + O++  '
        # For such assumption see Lopez-Sanchez & Esteban (2009) and Izotov et al. (2006) 
        logOtot = 12+numpy.log10(Otot)
        #logOtoterr = numpy.log10(numpy.abs(100+O_errp) / numpy.abs(100-O_errp))/2.0
        logOtoterr = Ototerr / (2.303 * Otot)        # equivalent method to above
        R = Otot / Otot
        Ratio = numpy.log10(R)
        Ratioerr = numpy.sqrt((Ototerr/Otot)**2 + (Ototerr/Otot)**2) / (2.303 )
        elem_abun['O'] = [Otot, Ototerr, O_errp, logOtot, logOtoterr, Ratio, Ratioerr]        
        print '    O_tot = %0.2f +- %0.2f ' % (logOtot, logOtoterr)
        logOtot_sun = 8.66 #+-0.05    taken from Asplund et al. 2005
        logOtot_sun_err = 0.05
        OH = logOtot - logOtot_sun
        OHerr = numpy.sqrt(logOtoterr**2 + logOtot_sun_err**2)
        print '    [O/H]  =  log(O/H) - log(O/H)_sun  =  %0.2f +- %0.2f' % (OH, OHerr)
        print '    * considering the existence of O+3 due to detection of HeII 4686, [O/H] should be higher'
        print '      by about 0.01 to 0.02 dex (Lopez-Sanchez & Esteban 2009).'
        
        # Nitrogen
        print '\n NITROGEN'
        print '    Assuming ICF(N) from Peimbert+Costero 69 = (Otot / O+) * N+'
        N_tot = (Otot / self.atom_abun['O2'][0]) * self.atom_abun['N2'][0]
        N_toterr = numpy.sqrt(O23sq**2 * (self.atom_abun['O3'][0]/Otot)**2 * (self.atom_abun['N2'][0]/self.atom_abun['O2'][0])**2 +
                             self.atom_abun['N2'][1]**2)
        Nicf = Otot / self.atom_abun['O2'][0]
        N_errp = (N_toterr/N_tot)*100
        logele = 12+numpy.log10(N_tot)
        #logeleerr = numpy.log10((100+N_errp) / (100-N_errp))/2.0
        logeleerr = N_toterr / (2.303 * N_tot)        # equivalent method to above
        #print 'Nicf=%0.2f ,   Ntot = %0.2f +- %0.2f' % (Nicf, 12+numpy.log10(Ntot), numpy.log10((Ntot+Ntoterr) / (Ntot-Ntoterr))/2.0)
        print '    ICF(N)=%0.2f ,   N_tot = %0.2f +- %0.2f' % (Nicf, logele, logeleerr)
        R = N_tot / Otot
        Ratio = numpy.log10(R)
        Ratioerr = numpy.sqrt((N_toterr/N_tot)**2 + (Ototerr/Otot)**2) / (2.303)
        elem_abun['N'] = [N_tot, N_toterr, N_errp, logele, logeleerr, Ratio, Ratioerr]
        print '    N/O = %0.2f +- %0.2f' % (Ratio, Ratioerr)
        if self.I_2144[0] != 0.0:
            abN2 = self.N2.getIonAbundance(int_ratio=self.I_2144, tem=te_low, den=dens, wave=2144)
            print 'N2 abundance from 2144A line =', abN2
            N_tot = (Otot / self.atom_abun['O2'][0]) * abN2[0]
            N_toterr = numpy.sqrt(O23sq**2 * (self.atom_abun['O3'][0]/Otot)**2 * (abN2[0]/self.atom_abun['O2'][0])**2 + abN2[1]**2)
            Nicf = Otot / self.atom_abun['O2'][0]
            N_errp = (N_toterr/N_tot)*100
            logele = 12+numpy.log10(N_tot)
            logeleerr = numpy.log10((100+N_errp) / (100-N_errp))/2.0
            #print 'Nicf=%0.2f ,   Ntot = %0.2f +- %0.2f' % (Nicf, 12+numpy.log10(Ntot), numpy.log10((Ntot+Ntoterr) / (Ntot-Ntoterr))/2.0)
            print '    ICF(N)=%0.2f ,   N_tot = %0.2f +- %0.2f' % (Nicf, logele, logeleerr)
            R = N_tot / Otot
            Ratio = numpy.log10(R)
            Ratioerr = numpy.sqrt((N_toterr/N_tot)**2 + (Ototerr/Otot)**2) / (2.303)
            #elem_abun['N'] = [N_tot, N_toterr, N_errp, logele, logeleerr, Ratio, Ratioerr]
            print '    N/O = %0.2f +- %0.2f' % (Ratio, Ratioerr)        
        else:
            print 'Could not find N2 abundance from 2144A line.'
            pass
        #raw_input(' ***  press enter to continue')
        
        # Neon
        print '\n NEON'
        print '    Assuming ICF(Ne) from Peimbert+Costero 69 =  (Otot / O++) * Ne++'
        Ne_icf = Otot / self.atom_abun['O3'][0]
        Ne_tot = self.atom_abun['Ne3'][0] * Ne_icf
        Ne_toterr = numpy.sqrt((O23sq * (self.atom_abun['O2'][0]/Otot)**2) * (self.atom_abun['Ne3'][0]/self.atom_abun['O3'][0])**2 +
                               self.atom_abun['Ne3'][1]**2)
        Ne_errp = (Ne_toterr/Ne_tot)*100
        logele = 12+numpy.log10(Ne_tot)
        logeleerr = numpy.log10((100+Ne_errp) / (100-Ne_errp))/2.0
        print '    ICF(Ne)=%0.2f ,   Ne_tot = %0.2f +- %0.2f' % (Ne_icf, logele, logeleerr)
        R = Ne_tot / Otot
        Ratio = numpy.log10(R)
        Ratioerr = numpy.sqrt((Ne_toterr/Ne_tot)**2 + (Ototerr/Otot)**2) / (2.303 )
        elem_abun['Ne'] = [Ne_tot, Ne_toterr, Ne_errp, logele, logeleerr, Ratio, Ratioerr]
        print '    Ne/O = %0.2f +- %0.2f' % (Ratio, Ratioerr)
        
        # Sulphur
        print '\n SULPHUR'
        if (self.atom_abun['S2'][0] > 0.0) and (self.atom_abun['S3'][0] > 0.0): 
            print '    Assuming ICF(S) from Garnett 89:  (S+ + S++)/Stot = [1 - (1 - O+/Otot)^alpha] ^ 1/alpha'
            ofrac = self.atom_abun['O2'][0]/Otot
            print '    *  O+ / Otot =:', ofrac
            #OpOtot = float(raw_input('Enter O+/Otot: '))
            if ofrac <= 0.15:
                correspondingSvalue2OpOtot = -0.25
            elif ofrac > 0.15 and ofrac < 0.3:
                correspondingSvalue2OpOtot = -0.1
            elif ofrac >= 0.3:
                correspondingSvalue2OpOtot = -0.05
            S_icf = 1.0 / 10**(correspondingSvalue2OpOtot)
            S_tot = (self.atom_abun['S2'][0] + self.atom_abun['S3'][0]) * S_icf
            S_icf_perr = 0.2 # this is the percentage error of the ICF taken from the min scale in Fig 7 of Garnett 89 = 0.05
            # the measured value in the x-asis is -0.25 +- 0.05, thus 20% of the measured value
            S_toterr = numpy.sqrt( ((S_icf_perr*S_icf)**2) * (self.atom_abun['S2'][0] + self.atom_abun['S3'][0])**2 +
                                   (self.atom_abun['S2'][1]**2 + self.atom_abun['S3'][1]**2) * S_icf**2 )
            S_errp = (S_toterr/S_tot)*100
            logele = 12+numpy.log10(S_tot)
            logeleerr = numpy.log10((100+S_errp) / (100-S_errp))/2.0
            print '    S_toterr = ', S_toterr, '=', S_errp, '%'
            print '    ICF(S)=%0.2f ,   S_tot = %0.2f +- %0.2f' % (S_icf, logele, logeleerr)
            # For comparison also clculate abundances with Izotov et al (2006)
            if logOtot <= 7.2:
                rule = 'Ial06_20a'
            elif (logOtot > 7.2) and (logOtot < 8.2):
                rule = 'Ial06_20b'
            elif logOtot >= 8.2:
                rule = 'Ial06_20c'
            S_Ial06 = icf.getElemAbundance(self.atom_abun, icf_list=[rule])
            abS = S_Ial06[rule][0]
            erS = S_Ial06[rule][1]
            print '    Abundances with:  Garnett 89        =  %0.3e +- %0.3e' % (S_tot, S_toterr) 
            print '                      Izotov et al. 06  =  %0.3e +- %0.3e' % (abS, erS)
            print '        they compare as: Izotov/Garnett = ', S_tot/abS
            print '                          Garnett + err = ', S_tot + S_toterr
            print '                          Izotov - err  = ', abS - erS
            R = S_tot / Otot
            Ratio = numpy.log10(R)
            Rerr = R * numpy.sqrt((S_toterr/S_tot)**2 + (Ototerr/Otot)**2)
            Ratioerr = Rerr / (2.303 * R)
            elem_abun['S'] = [S_tot, S_toterr, S_errp, logele, logeleerr, Ratio, Ratioerr]   
            print '    S/O = %0.2f +- %0.2f' % (Ratio, Ratioerr)
        else:
            elem_abun['S'] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]   
            print 'Not possible to determine the S abundance: S2=', self.atom_abun['S2'], '  S3=', self.atom_abun['S3']
        
        # Chlorine
        print '\n CHLORINE'
        if self.atom_abun['Cl3'][0] > 0.0:
            # Use Izotov et al (2006) just in case conditions for Peimbert et al (2005) are not met
            if logOtot <= 7.2:
                rule = 'Ial06_21a'
            elif (logOtot > 7.2) and (logOtot < 8.2):
                rule = 'Ial06_21b'
            elif logOtot >= 8.2:
                rule = 'Ial06_21c'
            ab_Ial06 = icf.getElemAbundance(self.atom_abun, icf_list=[rule])
            ab = ab_Ial06[rule][0]
            er = ab_Ial06[rule][1]
            # There are 2 conditions to use the ICF from Peimbert et al (2005):
            # 1. Diff between S3 and S2 is between 2 - 10 %
            # 2. Both Ar3 and Ar4 must be positive
            diffS = numpy.abs(numpy.abs(self.atom_abun['S3'][0]) - numpy.abs(self.atom_abun['S2'][0]))
            diffSper =  numpy.abs(self.atom_abun['S3'][0])/diffS
            if (diffSper >= 0.02) and (diffSper <= 0.1):
                print '    Percentage difference between S3 and S2 =', diffSper, '%'
                if (self.atom_abun['Ar4'][0] > 0.0) and (self.atom_abun['Ar3'][0] > 0.0):
                    print '    Assuming ICF(Cl) from Peimbert, Peimbert Ruiz (2005):  ((S+/S++)*Cl++ + Cl++ + (Ar+++/Ar++)*Cl++) / Cl++'
                    cl3 = self.atom_abun['Cl3'][0]
                    S23 = self.atom_abun['S2'][0] / self.atom_abun['S3'][0]
                    Ar43 = self.atom_abun['Ar4'][0]/self.atom_abun['Ar3'][0]
                    Cl_icf = (S23*cl3 + cl3 + Ar43*cl3) / cl3 
                    Cl_tot = cl3 * Cl_icf
                    S23err = S23 * numpy.sqrt( (self.atom_abun['S2'][1]/self.atom_abun['S2'][0])**2 + (self.atom_abun['S3'][1]/self.atom_abun['S3'][0])**2 )
                    Ar43err =  Ar43 * numpy.sqrt( (self.atom_abun['Ar4'][1]/self.atom_abun['Ar4'][0])**2 + (self.atom_abun['Ar3'][1]/self.atom_abun['Ar3'][0])**2 )
                    er1 = (S23*cl3) * numpy.sqrt((S23err/S23)**2 + (self.atom_abun['Cl3'][1]/self.atom_abun['Cl3'][0])**2)
                    er2 = (Ar43*cl3) * numpy.sqrt((Ar43err/Ar43)**2 + (self.atom_abun['Cl3'][1]/self.atom_abun['Cl3'][0])**2)
                    er3 = numpy.sqrt(er1**2 + er2**2 + self.atom_abun['Cl3'][1]**2)
                    clicferr = Cl_icf * numpy.sqrt((er3/(S23*cl3 + cl3 + Ar43*cl3))**2 + (self.atom_abun['Cl3'][1]/cl3)**2)
                    Cl_toterr = Cl_tot * numpy.sqrt((self.atom_abun['Cl3'][1]/cl3)**2 + (clicferr/Cl_icf)**2)
                    Cl_errp = (Cl_toterr/Cl_tot)*100
                    logele = 12+numpy.log10(Cl_tot)
                    logeleerr = numpy.log10((100+Cl_errp) / (100-Cl_errp))/2.0
                    print '    Cl_toterr = ', Cl_toterr, '=', Cl_errp, '%'
                    print '    ICF(Cl)=%0.2f ,   Cl_tot = %0.2f +- %0.2f' % (Cl_icf, logele, logeleerr)
                    # For comparison also clculate abundances with Izotov et al (2006)
                    print '    Abundances with: Peimbert et al. 05 =  %0.3e +- %0.3e' % (Cl_tot, Cl_toterr) 
                    print '                     Izotov et al. 06   =  %0.3e +- %0.3e' % (ab, er)
                    print '       they compare as: Izotov/Peimbert = ', Cl_tot/ab
                    print '                         Peimbert + err = ', Cl_tot + Cl_toterr
                    print '                          Izotov - err  = ', ab - er
                    R = Cl_tot / Otot
                    Ratio = numpy.log10(R)
                    Rerr = R * numpy.sqrt((Cl_toterr/Cl_tot)**2 + (Ototerr/Otot)**2)
                    elem_abun['Cl'] = [ab, er, Cl_errp, logele, logeleerr, Ratio, Ratioerr]
                    Ratioerr = Rerr / (2.303 * R)
            else:
                print '    Assuming ICF(Cl) from Izotov et al. (2006)'
                perr = er/ab *100
                print '    Cl_toterr = ', er, '=', perr, '%'
                logele = 12+numpy.log10(ab)
                logeleerr = numpy.log10((100+perr) / (100-perr))/2.0
                print '    Cl_tot = %0.2f +- %0.2f' % (logele, logeleerr)
                R = ab /Otot
                Ratio = numpy.log10(R)
                Rerr = R * numpy.sqrt((er/ab)**2 + (Ototerr/Otot)**2)
                Ratioerr = Rerr / (2.303 * R)
                elem_abun['Cl'] = [ab, er, perr, logele, logeleerr, Ratio, Ratioerr]
            print '    Cl/O = %0.2f +- %0.2f' % (Ratio, Ratioerr)
        else:
            elem_abun['Cl'] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            print '    No Cl abundance available.'            
        
        # Argon
        print '\n ARGON'
        if logOtot <= 7.2:
            rule = 'Ial06_22a'
        elif (logOtot > 7.2) and (logOtot < 8.2):
            rule = 'Ial06_22b'
        elif logOtot >= 8.2:
            rule = 'Ial06_22c'
        ab_Ial06 = icf.getElemAbundance(self.atom_abun, icf_list=[rule])
        ab = ab_Ial06[rule][0]
        er = ab_Ial06[rule][1]
        # If possible, use Peimbert et al (2005)
        if (self.atom_abun['Ar4'][0] > 0.0) and (self.atom_abun['Ar3'][0] > 0.0) and (self.atom_abun['S2'][0] > 0.0) and (self.atom_abun['S3'][0] > 0.0):
            print '    Assuming ICF(Ar) from Peimbert, Peimbert Ruiz (2005):  (S+/S++)*Ar++ * (Ar++ + Ar+++)/(Ar++ + Ar+++)'
            Ar34 = self.atom_abun['Ar3'][0] + self.atom_abun['Ar4'][0]
            S23 = self.atom_abun['S2'][0] / self.atom_abun['S3'][0]
            Ar_icf = (S23 * self.atom_abun['Ar3'][0] + Ar34)/Ar34 
            Ar_tot = Ar34 * Ar_icf
            S23err = S23 * numpy.sqrt( (self.atom_abun['S2'][1]/self.atom_abun['S2'][0])**2 + (self.atom_abun['S3'][1]/self.atom_abun['S3'][0])**2 )
            Ar34err =  numpy.sqrt( self.atom_abun['Ar3'][1]**2 + self.atom_abun['Ar4'][1]**2 )
            er1 = (S23 * self.atom_abun['Ar3'][0]) * numpy.sqrt((S23err/S23)**2 + (self.atom_abun['Ar3'][1]/self.atom_abun['Ar3'][0])**2)
            er2 = numpy.sqrt(er1**2 + self.atom_abun['Ar3'][1]**2 + self.atom_abun['Ar4'][1]**2)
            Ar_icferr = Ar_icf * numpy.sqrt((er2/(S23 * self.atom_abun['Ar3'][0] + Ar34))**2 + (Ar34err/Ar34)**2)
            Ar_toterr = Ar_tot * numpy.sqrt((Ar34err/Ar34)**2 + (Ar_icferr/Ar_icf)**2)
            Ar_errp = (Ar_toterr/Ar_tot)*100
            print '    Ar_toterr = ', Ar_toterr, '=', Ar_errp, '%'
            logele = 12+numpy.log10(Ar_tot)
            logeleerr = numpy.log10((100+Ar_errp) / (100-Ar_errp))/2.0
            print '    ICF(Ar)=%0.2f ,   Ar_tot = %0.2f +- %0.2f' % (Ar_icf, logele, logeleerr)
            # For comparison also clculate abundances with Izotov et al (2006) when we have both
            if logOtot <= 7.2:
                rule = 'Ial06_23a'
            elif (logOtot > 7.2) and (logOtot < 8.2):
                rule = 'Ial06_23b'
            elif logOtot >= 8.2:
                rule = 'Ial06_23c'
            ab_Ial06 = icf.getElemAbundance(self.atom_abun, icf_list=[rule])
            ab = ab_Ial06[rule][0]
            er = ab_Ial06[rule][1]
            #logele = 12+numpy.log10(ab)
            #perr = er/ab *100
            #logeleerr = numpy.log10((100+perr) / (100-perr))/2.0
            #elem_abun['Ar'] = [ab, er, logele, logeleerr]   
            print '    Abundances with: Peimbert et al. 05 =  %0.3e +- %0.3e' % (Ar_tot, Ar_toterr) 
            print '                     Izotov et al. 06   =  %0.3e +- %0.3e' % (ab, er)
            print '       they compare as: Izotov/Peimbert = ', Ar_tot/ab
            print '                         Peimbert + err = ', Ar_tot + Ar_toterr
            print '                          Izotov - err  = ', ab - er
            R = Ar_tot / Otot
            Ratio = numpy.log10(R)
            Rerr = R * numpy.sqrt((Ar_toterr/Ar_tot)**2 + (Ototerr/Otot)**2)
            Ratioerr = Rerr / (2.303 * R)
            elem_abun['Ar'] = [ab, er, Ar_errp, logele, logeleerr, Ratio, Ratioerr]
        else:
            print '    Assuming ICF(Ar) from Izotov et al. (2006)'
            perr = er/ab * 100.0
            print '    Ar_toterr = ', er, '=', perr, '%'
            logele = 12+numpy.log10(ab)
            logeleerr = numpy.log10((100+perr) / (100-perr))/2.0
            print '    Ar_tot = %0.2f +- %0.2f' % (logele, logeleerr)
            R = ab / Otot
            Ratio = numpy.log10(R)
            Rerr = R * numpy.sqrt((er/ab)**2 + (Ototerr/Otot)**2)
            Ratioerr = Rerr / (2.303 * R)
            elem_abun['Ar'] = [ab, er, perr, logele, logeleerr, Ratio, Ratioerr]
        print '    Ar/O = %0.2f +- %0.2f' % (Ratio, Ratioerr)
        
        # Iron
        print '\n IRON'
        if self.atom_abun['Fe3'][0] > 0.0:
            # Use Izotov et al (2006)
            if logOtot <= 7.2:
                rule = 'Ial06_41a'
            elif (logOtot > 7.2) and (logOtot < 8.2):
                rule = 'Ial06_24b'
            elif logOtot >= 8.2:
                rule = 'Ial06_24c'
            ab_Ial06 = icf.getElemAbundance(self.atom_abun, icf_list=[rule])
            ab = ab_Ial06[rule][0]
            er = ab_Ial06[rule][1]
            print '    Assuming ICF(Fe) from Izotov et al. (2006)'
            perr = er/ab * 100.0
            print '    Fe_toterr = ', er, '=', perr, '%'
            logele = 12+numpy.log10(ab)
            logeleerr = numpy.log10((100+perr) / (100-perr))/2.0
            print '    Fe_tot = %0.2f +- %0.2f' % (logele, logeleerr)
            R = ab / Otot
            Ratio = numpy.log10(R)
            Rerr = R * numpy.sqrt((er/ab)**2 + (Ototerr/Otot)**2)
            Ratioerr = Rerr / (2.303 * R)
            elem_abun['Fe'] = [ab, er, perr, logele, logeleerr, Ratio, Ratioerr]
            print '    Fe/O = %0.2f +- %0.2f' % (Ratio, Ratioerr)
        else:
            elem_abun['Fe'] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            print '    No Fe abundance available.'  
                      
        '''
        # Aluminum STILL PENDING BECAUSE NO SOURCE OF ICF FOR HII REGIONS.
        print '\n ALUMINUM'
        '''
        
        # METHODS NOT IN PYNEB
        
        # Corrected Auroral Line Method (Pena-Guerrero et al. 2012b):
        Ocalm = 1.0825 * elem_abun['O'][3] - 0.375
        errOcalm = 1.0825 * elem_abun['O'][4]
        print '\nO_direct = %0.2f +- %0.2f    O_CALM = %0.2f +- %0.2f' % (elem_abun['O'][3], elem_abun['O'][4], Ocalm, errOcalm)
        
        # Recalibrated R23 Method O abundances (Pena-Guerrero et al. 2012b):
        #if (self.I_3727[0] != 0.0) and (self.I_4959[0] != 0.0) and (self.I_5007[0] != 0.0):
        if (self.I_3727[0] != 0.0) and (self.I_5007[0] != 0.0):
            self.I_4959[0] = self.I_5007[0]/3.0
            R23 = (self.I_3727[0] + self.I_4959[0] + self.I_5007[0])/ 100.0
            errR23 = numpy.sqrt(self.I_3727[1]**2 + (self.I_5007[1]/3.)**2 + self.I_5007[1]**2)/ 100.0
            P = (self.I_4959[0] + self.I_5007[0]) / (self.I_3727[0] + self.I_4959[0] + self.I_5007[0])
            e1 = numpy.sqrt((self.I_5007[1]/3.)**2 + self.I_5007[1]**2)
            e2 = numpy.sqrt(self.I_3727[1]**2 + (self.I_5007[1]/3.)**2 + self.I_5007[1]**2)
            errP = P * numpy.sqrt((e1/(self.I_4959[0] + self.I_5007[0]))**2 + (e2/(self.I_3727[0] + self.I_4959[0] + self.I_5007[0]))**2)
            errP2 = 2*P*errP
            print 'R_23 = %0.2f +- %0.2f        P = %0.2f +- %0.2f' % (R23, errR23, P, errP)
            if Ocalm >= 8.55:#elem_abun['O'][3] >= 8.25:
                Orrm = (R23 + 1837.0 + 2146.0*P + 850.0*P*P) / (209.5 + 201.7*P + 43.98*P*P + 1.793*R23)
                e1 = numpy.sqrt(errR23**2 + 2146.0*errP**2 + 850.0*errP2**2)
                e2 = numpy.sqrt(201.7*errP**2 + 43.98*errP2**2 + 1.793*errR23**2)  
                errOrrm = Orrm * numpy.sqrt((e1/(R23 + 1837.0 + 2146.0*P + 850.0*P*P))**2 + (e2/(209.5 + 201.7*P + 43.98*P*P + 1.793*R23))**2)
                print 'Upper branch object:       O_RRM = %0.2f +- %0.2f' % (Orrm, errOrrm)
            elif Ocalm <= 8.29:#elem_abun['O'][3] <= 8.00:
                Orrm = (R23 + 90.73 + 94.58*P - 5.26*P*P) / (14.81 + 5.52*P + 5.81*P*P - 0.252*R23)
                e1 = numpy.sqrt(errR23**2 + 94.58*errP**2 + 5.26*errP2**2)
                e2 = numpy.sqrt(5.52*errP**2 + 5.81*errP2**2 + 0.252*errR23**2)  
                errOrrm = Orrm * numpy.sqrt((e1/(R23 + 90.73 + 94.58*P - 5.26*P*P))**2 + (e2/(14.81 + 5.52*P + 5.81*P*P - 0.252*R23))**2)
                print 'Lower branch object:       O_RRM = %0.2f +- %0.2f' % (Orrm, errOrrm)
            else:
                print 'O abundance in the degeneracy zone... Unable to derermine O_RRM.'
            raw_input()
        
        # Now calculate the ionic abundances of C^{++}/O^{++}, N^{++}, and C/O according to Garnett et al. (1995)
        # Equation 2 for C^{++}/O^{++}
        te_used = te_high
        tc = te_used[0]/10000.0         # central temp
        tpluserr = te_used[1]/10000.0   # central temp plus error
        if self.I_1666[0] <= 0.0:
            I_1663 = self.I_1661
        else:
            I_1663 = self.I_1666
        I_1909 = self.I_1909
        C2toO2 = 0.089 * numpy.exp(-1.09/tc) * (I_1909[0]/I_1663[0])
        # error calculation
        uplimC2toO2_temp = 0.089 * numpy.exp(-1.09/tpluserr) * (I_1909[0]/I_1663[0])
        lolimC2toO2_temp = 0.089 * numpy.exp(-1.09/(tc-numpy.abs(tpluserr-tc))) * (I_1909[0]/I_1663[0])
        uplimC2toO2_lines = 0.089 * numpy.exp(-1.09/tc) * ((I_1663[0] + I_1663[1])/(I_1909[0] + I_1909[1]))
        lolimC2toO2_lines = 0.089 * numpy.exp(-1.09/tc) * ((I_1663[0] - I_1663[1])/(I_1909[0] - I_1909[1]))
        err_C2toO2_temp = numpy.sqrt((uplimC2toO2_temp - C2toO2)**2 + (C2toO2 - lolimC2toO2_temp)**2)
        err_C2toO2_lines = numpy.sqrt((uplimC2toO2_lines - C2toO2)**2 + (C2toO2 - lolimC2toO2_lines)**2)
        #C2toO2_err = numpy.sqrt(err_C2toO2_temp**2 + err_C2toO2_lines**2)
        C2toO2_err = (err_C2toO2_temp + err_C2toO2_lines) / 2.0
        #if C2toO2 < 0.0:
        #    C2toO2 = 0.0
        # Equation 3 for C^{++}/O^{++}
        #I_1759 = self.I_1749
        # Find the fraction of O++/Otot
        op = sorted_atoms.index('O2')
        opp = sorted_atoms.index('O3')
        Ofrac = totabs_ions_list[opp][0] / Otot
        print '\n X(O++) = %0.3f' % Ofrac
        print ' C++/O++ = %0.3f +- %0.03f' % (C2toO2, C2toO2_err)
        # now determine the metallicity of the object: Z
        # We know that 1/X = (X + Y + Z)/X = 1 + Y/X + Z/X
        # the ammount of helium is obtained from the HELIO10 code
        #***
        Tavg, delt = self.get_avgTandDelta4helio10(Otot, totabs_ions_list[op][0], totabs_ions_list[opp][0], te_low[0], te_high[0])
        print ' This is the weighted temperature between O+ and O++: ', Tavg
        print '    and this is the value of detla: ', delt
        self.get_logIs4_helio10()
        wavsHe, logInts, logErrs = self.get_logIs4_helio10()
        for w, lI, le in zip(wavsHe, logInts, logErrs):
            print 'For', w, 'log Intensity =', lI, '+-', le
        #***
        helio = 0.278 #0.298
        YoverX = 4 * (helio)
        # Z/X = Mass(C, N, O, ...)/Mass(H) ~ 2M(O)/M(H) = 2m(O^16)/m(H^1)*n(O)/n(H) = 32 * Otot
        ZoverX = 32 * Otot
        oneoverX = 1 + YoverX + ZoverX
        X = 1 / oneoverX
        # now solve for Z from Z/X 
        Z = ZoverX * X
        print ' The metallicity Z =', Z
        #raw_input(' ***  press enter to continue')
        # use this metallicity in figure 2 of Garnett ('95) to find X(C++)/X(O++)
        XCXO = 0.9#0.86
        icfC = 1/XCXO #1 / (C2toO2 * Ofrac)
        print ' ICF(C) using Garnett (1995) = %0.3f' % (icfC)
        cpp = sorted_atoms.index('C3')
        totC_garnett = totabs_ions_list[cpp] * icfC     
        logele_gar = 12+numpy.log10(totC_garnett[0])
        logeleerr_gar = totC_garnett[1] / (2.303 * totC_garnett[0])
        percent_err = totC_garnett[1]*100 / totC_garnett[0]
        #logeleerr_gar = numpy.log10((100+percent_err)/(100-percent_err)) / 2.0
        print '  C++ = ', totabs_ions_list[cpp], '%err=', percent_err
        print ' total C = ', totC_garnett
        print ' 12+log(C) = %0.2f +- %0.2f' % (logele_gar, logeleerr_gar)
        # now determine carbon abundance with my thesis method
        if data2use != None:
            if self.use_Chbeta:
                # data2use = CHbeta, catalog_wavelength, flambdas, element, ion, norm_fluxes, errflux, percerrflx, norm_intenties, errinten
                CHbeta, catalog_wavelength, flambdas, _, _, norm_fluxes, _, _, norm_intenties, errinten = data2use
        else:
            #AdvOpscols_in_file = [self.wavelength, self.flambda, self.element, self.ion, self.forbidden, self.howforb, self.flux, self.errflux, self.percerrflx, self.intensity, self.errinten, self.percerrinten, self.ew, self.errew, self.percerrew]
            CHbeta = self.cHbeta / 0.434
            catalog_wavelength = self.AdvOpscols_in_file[0]
            flambdas = self.AdvOpscols_in_file[1]
            ion = self.AdvOpscols_in_file[3]
            norm_fluxes = self.AdvOpscols_in_file[7]
            norm_intenties = self.AdvOpscols_in_file[10]
            errinten = self.AdvOpscols_in_file[11]        
        # We are interested in the flambda values in order to obtain the C(lambda) to correct the flux of that wavelength range
        rounded_catwavs = round_wavs(catalog_wavelength)
        idx1661 = rounded_catwavs.index(1661.0)
        idx1907 = rounded_catwavs.index(1907.0)
        deltaflambda = numpy.abs(flambdas[idx1907] - flambdas[idx1661])
        avgC = 10**(CHbeta*deltaflambda)
        IC3IO3_ratio = norm_fluxes[idx1907] / norm_fluxes[idx1661] * avgC
        # Since we do not trust the 1666 measurement enough, we will use the optical part to back it up
        if 4959.0 in rounded_catwavs:
            idxO3 = rounded_catwavs.index(4959.0)
            # With the corresponding temperature and density, we need a theoretical ratio of the intensity of 1661/4959
            ionic_ratio = 6.771e-24 / 1.222e-21     # from ionic with Te=10,000 and ne=100
        else:
            idxO3 = rounded_catwavs.index(5007.0)
            # With the corresponding temperature and density, we need a theoretical ratio of the intensity of 1661/5007
            ionic_ratio = 6.771e-24 / 3.531e-21     # from ionic with Te=10,000 and ne=100
        I1661 = norm_intenties[idxO3] * ionic_ratio
        err_I1661 = errinten[idxO3] * ionic_ratio
        # now from the IC3IO3_ratio equation, solve for the I1661 intensity and use the optical 4959 to find C3] 1907
        corrI1907 = IC3IO3_ratio * I1661
        err_corrI1907 = IC3IO3_ratio * err_I1661 
        I_1907 = numpy.array([corrI1907, err_corrI1907])
        C3_thesis = self.C3.getIonAbundance(I_1907, te_used, dens, to_eval='L(1907)')
        print '         This is the ionic abundance of C++ with my thesis method', C3_thesis 
        # now, using the correction factor given by Garnett '95
        Ctot_thesis = icfC*C3_thesis
        logele_thes = 12+numpy.log10(Ctot_thesis[0])
        #logeleerr = numpy.log10((100+C_errp) / (100-C_errp))/2.0
        logeleerr_thes = Ctot_thesis[1] / (2.303 * Ctot_thesis[0])
        print '         total abundance with thesis method, Ctot =', Ctot_thesis
        print '         12+log(Ctot_thesis) = %0.2f +- %0.2f' % (logele_thes, logeleerr_thes)
        
        # But only print the correct value: if CHbeta take the one from my thesis, else use the ebv one
        if data2use != None:
            Ctot = Ctot_thesis
            logele = logele_thes
            logeleerr = logeleerr_thes
        else:
            Ctot = totC_garnett
            logele = logele_gar
            logeleerr = logeleerr_gar
        C_errp = Ctot[1]*100/Ctot[0]
        R = Ctot[0] / Otot
        print Ctot[0], '/', Otot, '=', R
        Ratio = numpy.log10(R)
        Ratioerr = numpy.sqrt((Ctot[1]/Ctot[0])**2 + (Ototerr/Otot)**2) / (2.303)            
        elem_abun['C'] = [Ctot[0], Ctot[1], C_errp, logele, logeleerr, Ratio, Ratioerr]
        print ' C/O = %0.2f +- %0.2f\n' % (Ratio, Ratioerr)
        #raw_input(' ***  press enter to continue')
        
        te_used = te_low#te_high
        tc = te_used[0]/10000.0         # central temp
        tpluserr = te_used[1]/10000.0   # central temp plus error
        I_1752 = self.I_1752
        N2toO2 = 0.212 * numpy.exp(-0.43/tc) * (I_1752[0]/I_1663[0])
        # error calculation
        uplimN2toO2_temp =  0.212 * numpy.exp(-0.43/tpluserr) * (I_1752[0]/I_1663[0])
        lolimN2toO2_temp =  0.212 * numpy.exp(-0.43/(tc-numpy.abs(tpluserr-tc))) * (I_1752[0]/I_1663[0])
        uplimN2toO2_lines =  0.212 * numpy.exp(-0.43/tc) * ((I_1752[0] + I_1752[1])/(I_1663[0] + I_1663[1]))
        lolimN2toO2_lines =  0.212 * numpy.exp(-0.43/tc) * ((I_1752[0] - I_1752[1])/(I_1663[0] - I_1663[1]))
        err_N2toO2_temp = numpy.sqrt((uplimN2toO2_temp - N2toO2)**2 + (N2toO2 - lolimN2toO2_temp)**2)
        err_N2toO2_lines = numpy.sqrt((uplimN2toO2_lines - N2toO2)**2 + (N2toO2 - lolimN2toO2_lines)**2)
        #N2toO2_err = numpy.sqrt(err_N2toO2_temp**2 + err_N2toO2_lines**2)
        N2toO2_err = (err_N2toO2_temp + err_N2toO2_lines) / 2.0
        print ' N++/O++ = %0.3f +- %0.03f' % (N2toO2, N2toO2_err)
        #raw_input(' ***  press enter to continue')
        
        self.t2_RLs(dens, totabs_ions_list[opp][0], totabs_ions_list[cpp][0])
            
        # Write results in text file
        if self.writeouts:
            if self.used_measured_info:
                out_file = self.object_name+'_measuredLI_IonicTotAbundances_'+RedCorType+'.txt'
            else:
                out_file = self.object_name+'_IonicTotAbundances_'+RedCorType+'.txt'
            path_object = '../results/'+self.object_name
            fullpath_outfile = os.path.join(path_object, out_file)
            outf = open(fullpath_outfile, 'w+')
            print >> outf, ('{:<45} {:>10} {:>6}'.format('# Temperatures used [K]', 'High =', int(te_high[0]))+'{:>30} {:>6}'.format('Low =', int(te_low[0]))+
                            '{:>35} {:>6}'.format('Very_low =', int(te_verylow[0])))
            print >> outf, ('{:<9} {:>13} {:>13} {:>10} {:>17} {:>18} {:>17} {:>18} {:>17}'.format('# Line_label', 'Intensity', 'percent_err', 'abH', 'abH_err', 
                                                                                            'abL', 'abL_err', 'abVL', 'abVL_err'))
            for l, I, ab1, ab2, ab3 in zip(line_label_list, line_I_list, ab_high_list, ab_low_list, ab_very_low):
                percent_I = (I[1] * 100.)/I[0]
                print >> outf, ('{:<9} {:>15.3f} {:>8.3f} {:>20.5e} {:>15.5e} {:>20.5e} {:>15.5e} {:>20.5e} {:>15.5e}'.format(l, I[0], percent_I, ab1[0], ab1[1], 
                                                                                                                              ab2[0], ab2[1], ab3[0], ab3[1]))
            print >> outf, '#####'
            print >> outf, '# IONIC ABUNDANCES'
            print >> outf, ('{:<6} {:>15} {:>12} {:>10} {:>7}'.format('# Ion', 'abundance', 'abs error', 'LOGabund', 'LOGerr'))
            for a, abund in zip(sorted_atoms, totabs_ions_list):
                if abund[0] > 0.0:
                    ionab = abund[0]
                    ionerr = abund[1]
                    logionab = 12 + numpy.log10(ionab)
                    erlogionerr = (numpy.log10(ionab+ionerr) - numpy.log10(ionab-ionerr)) /2
                else:
                    ionab = 0.0
                    ionerr = 0.0
                    logionab = 0.0
                    erlogionerr = 0.0
                print >> outf, ('{:<6} {:>15.3e} {:>12.3e} {:>10.2f} {:>7.2f}'.format(a, ionab, ionerr, logionab, erlogionerr))
            print >> outf, '#####'
            print >> outf, '# RATIOS -- using equations from Garnett et al. (1995)'
            print >> outf, ('{:<6} {:>11} {:>11} {:>10}'.format('# Ion ratio', 'Line ratio', 'Value', 'Abs err'))
            print >> outf, ('{:<6} {:>15} {:>12.3f} {:>10.3f}'.format('C++/O++', '1909/1666', C2toO2, C2toO2_err))
            print >> outf, ('{:<6} {:>15} {:>12.3f} {:>10.3f}'.format('N++/O++', '1752/1666', N2toO2, N2toO2_err))
            '''
            npp = sorted_atoms.index('N2')
            opp = sorted_atoms.index('O3')
            n2too2 = totabs_ions_list[npp][0]/totabs_ions_list[opp][0]
            n2too2_err = n2too2 * numpy.sqrt((totabs_ions_list[npp][1]/totabs_ions_list[npp][0])**2 + (totabs_ions_list[opp][1]/totabs_ions_list[opp][0])**2)
            print >> outf, ('{:<6} {:>15} {:>12.3f} {:>10.3f}'.format('N+/O++', '6548/5007', n2too2, n2too2_err))
            '''
            # print total abundances in file
            print >> outf, '#####'
            print >> outf, '# TOTAL ABUNDANCES '
            print >> outf, ('{:<5} {:>12} {:>12} {:>7} {:>10} {:>7} {:>10} {:>7}'.format('# Element', 'abundance', 'abs error', '% err', 
                                                                                   'LOGabund', 'LOGerr', 'X/O', 'err'))
            for ele, vals in elem_abun.items():
                # ratio and error with absolute values
                #XO = vals[0] / Otot
                #XOerr = XO * numpy.sqrt((Ototerr/Otot)**2 + (vals[1]/vals[0])**2)
                print >> outf, ('{:<6} {:>15.3e} {:>12.3e} {:>7.0f} {:>10.2f} {:>7.2f} {:>10.2f} {:>7.2f}'.format(ele, vals[0], vals[1], 
                                                                                                                  vals[2], vals[3], vals[4], 
                                                                                                                  vals[5], vals[6]))
        # Make sure that the temperatures and densities file closes properly
        if self.writeouts:
            outf.close()
        
    def corr_ColExcit(self, Idered_norm, Hlines=None, verbose=False):
        '''
        This function is the second iteration of reddening correction. It fits the observed H lines given to the theoretical
        ones found with INTRAT by Storey & Hummer (1995).
        # Hlines = list of the hydrogen wavelengths to look for (typically from Halpha to H12).
        '''
        # Hydrogen lines to be considered for correction
        #        Halpha,  H5,   H6,   H7,   H8,   H9,  H10,   H11,  H12
        #Hlines = [6563, 4340, 4101, 3967, 3889, 3835, 3798, 3771, 3750]
        # HOWERVER H8 and H9 are contaminated by HeI (H8 is also contaminated with [NeIII], and H12 is too weak.
        if Hlines == None:
            Hlines = [6563.0, 4340.0, 4102.0, 3967.0, 3798.0, 3771.0]
        ### For collisional excitation, interpolate from Table 1 of Peimbert, Luridiana, Peimbert (2007, ApJ, 666, 636)
        # Table 1
        # Objects = NGC346, NGC2363, Haro29, SBS0335-052, IZw18
        TeOII = [12600.0, 13800.0, 14000.0, 15600.0, 15400.0]
        xalpha = [0.011, 0.037, 0.033, 0.086, 0.070]
        #xbeta = [0.007, 0.027, 0.021, 0.066, 0.053]
        # x_lambda = I_col/I_tot
        # now do the interpolation according to the [OII] temperature
        xL = []
        if self.TO2gar[0] == 0.0:
            TeLow = self.te_low[0]
        else:
            TeLow = self.TO2gar[0]
        xalpha_interp = numpy.interp(TeLow, TeOII, xalpha)
        xL.append(xalpha_interp)
        xbeta = xalpha_interp * 0.67
        xL.append(xbeta)
        # Halpha, H5,  H6,  H7,  H10,  H11
        L = [3.0, 5.0, 6.0, 7.0, 10.0, 11.0]
        for l in L:
            xl = xalpha_interp / ( 2**((l-2.0)/3.0) )
            xL.append(xl)
        # For the following hydrogen lines recalculate the intensity correcting for collisional exitation
        norm_IcorrCE = []   # intensities corrected for collisional excitation
        obs_ratios = []     # observed ratios normalized to Hbeta --> these are the ones to be compared with the INTRAT ones
        found_Hlines = []   # the Hlines that were found in the observations
        rounded_wavelengths = round_wavs(self.wavelength)
        for w, el, I in zip(rounded_wavelengths, self.element, Idered_norm):
            #print 'w, el, I', w, el, I
            for h, l in zip(Hlines, xL):
                #print 'h, l', h, l
                if (w == h) and (el == 'H'):
                    #print 'FOUND ONE!', h
                    found_Hlines.append(h)
                    newI = I * (1-l)
                    normI = newI/100.
                    if verbose == True:
                        print w, 'before', I, '  after collisionally excited corrected and norm to Hbeta:', newI, '  ratio2Hbeta', normI
                    norm_IcorrCE.append(newI)
                    obs_ratios.append(normI)
                    if w == 4102:
                        norm_H6theo = normI
        # Return the corrected intensities, the observed normalized intensities, the Hlines that were looked for, 
        # the found Hlines, and the specific ratio of H6 to Hbeta.
        return norm_IcorrCE, obs_ratios, Hlines, found_Hlines, norm_H6theo
    
    def find_Chi_of_CE(self, Idered_norm, theoCE, Hlines=None, verbose=False):
        ''' This function finds the best combination of CHbeta and EW_of_Hbeta.
        # Idered_norm = normalized and dereddend intensities
        # theoCE = theoretical hydrogen intensities corrected for collisional excitation found with INTRAT (see below) '''
        # Correct for collisional excitation
        IcorrCE, obs_ratios, Hlines, found_Hlines, norm_H6theo = self.corr_ColExcit(Idered_norm, Hlines, verbose)
        # Recalculate the intensities of the most prominent hydrogen lines (Halpha through H12) to match them 
        # with the theoretical ratios given by INTRAT (Storey & Hummer, 1995, MNRAS, 272, 41).
        uncert = []
        rounded_wavelengths = round_wavs(self.wavelength)
        for H, Ic, obsr in zip(found_Hlines, IcorrCE, obs_ratios):
            idx = Hlines.index(H)
            if H in rounded_wavelengths:
                H_index = rounded_wavelengths.index(H)
                u = (1 - (obsr / theoCE[idx])) * 100.0 / self.percerrinten[H_index]
                I = self.intensity[H_index]
                if verbose:
                    print H, 'theo_ratio =', theoCE[idx], 'obs_ratio', obsr, '   IcorrCE =', Ic, '   Idered = ', I, '   percent_Iuncert[H_index]=', self.percerrinten[H_index]
                    print '   % err_collisional_excitation_corr = ', u 
                uncert.append(u)
        # In order to find the best combination of C_Hbeta and EWabsHbeta determine chi squared
        sqs = []
        for u in uncert:
            us = u * u
            sqs.append(us)
        Chi_sq = sum(sqs)
        if verbose:
            print 'Chi_sq = ', Chi_sq
        return Chi_sq, norm_H6theo
        
    def redcor2(self, theoCE, verbose=False, em_lines=False):
        '''This function is only used if C(Hbeta) was used for reddening correction. If E(B-v) was used instead, this
        function will be skiped.'''
        # Variables obtained from parent class
        I_theo_HaHb = self.I_theo_HaHb
        C_Hbeta = self.cHbeta / 0.434
        EWabsHbeta = self.EWabsHbeta
        continuum = self.continuum
        #AdvOpscols_in_file = [self.wavelength, self.flambda, self.element, self.ion, self.forbidden, self.howforb, self.flux, self.errflux, self.percerrflx, self.intensity, self.errinten, self.percerrinten, self.ew, self.errew, self.percerrew]
        flux = self.AdvOpscols_in_file[6]
        observed_wavelength = self.AdvOpscols_in_file[0]
        intensities = self.AdvOpscols_in_file[9]
        # Variables used and defined for this function
        number_iterations = 14 #this number must be even
        EWabsHbeta_increase = 0.1
        C_Hbeta_increase = 0.01
        Chi_sq_models = []
        EWabsHbeta_values = []
        EWabsHbeta_values.append(EWabsHbeta)
        C_Hbeta_values = []
        C_Hbeta_values.append(C_Hbeta)
        rounded_wavelengths = round_wavs(self.wavelength)
        Halpha_idx = rounded_wavelengths.index(6563.)
        Hbeta_idx = rounded_wavelengths.index(4861.)
        H6_idx = rounded_wavelengths.index(4102)
        dif_TheoObs_H6Hb_values = []        
        Hlines=None
        # Find the the Chi squared of the first round dereddened intensities 
        I_obs_H6Hb = rounded_wavelengths[H6_idx] / rounded_wavelengths[Hbeta_idx]
        Chi_sq, I_theo_H6Hb = self.find_Chi_of_CE(intensities, theoCE, Hlines, verbose)
        Chi_sq_models.append(Chi_sq)
        dif_TheoObs_H6Hb = numpy.fabs(I_theo_H6Hb - I_obs_H6Hb) 
        dif_TheoObs_H6Hb_values.append(dif_TheoObs_H6Hb)
        I_obs_HaHb = self.intensity[Halpha_idx] / self.intensity[Hbeta_idx]
        print ' ***    I_theo_HaHb =', I_theo_HaHb, '   I_obs_HaHb =', I_obs_HaHb
        diff_HaHb_values = []
        diff_HaHb = numpy.fabs(I_theo_HaHb - I_obs_HaHb)
        diff_HaHb_values.append(diff_HaHb)
        # First, variate EWabsHbeta with C_Hbeta fixed
        for EWabsHbeta_iterations in range(0, number_iterations):
            print 'cHbeta = %0.3f        EWabsHbeta = % 0.3f' % (self.cHbeta, EWabsHbeta)
            print 'EWabsHbeta_iterations', EWabsHbeta_iterations
            if I_theo_HaHb < I_obs_HaHb:
                EWabsHbeta = EWabsHbeta + EWabsHbeta_increase
            elif I_theo_HaHb > I_obs_HaHb:
                EWabsHbeta = EWabsHbeta - EWabsHbeta_increase
                if EWabsHbeta < 0.0:
                    EWabsHbeta = 0.00001
            EWabsHbeta_values.append(EWabsHbeta)
            # We need to correct these fluxes for underlying absorption correction with the new value of EWabsHbeta
            intensities = self.underlyingAbsCorr(return_Is=True, EWabsHbeta=EWabsHbeta)
            # If asked to, keep only the positive fluxes = emission lines
            if em_lines:
                catalog_lines, _, element_lines, _, _, _, norm_fluxes, _, norm_intensities, _ = find_emission_lines(rounded_wavelengths, self.element, self.ion, self.forbidden, self.howforb, observed_wavelength, flux, intensities, self.EW, continuum)
            else:
                catalog_lines = rounded_wavelengths
                element_lines = self.AdvOpscols_in_file[2]
                norm_fluxes = flux
                norm_intensities = intensities
            # Dered again and find the Chi_squared of that model
            cHbeta = 0.434*C_Hbeta
            _, Idered_norm, _ = self.Halpha2Hbeta_dered(cHbeta=cHbeta, fluxes=norm_fluxes, intensities=norm_intensities)
            I_obs_H6Hb = rounded_wavelengths[H6_idx] / rounded_wavelengths[Hbeta_idx]
            Chi_sq, norm_H6theo = self.find_Chi_of_CE(Idered_norm, theoCE, Hlines, verbose)
            Chi_sq_models.append(Chi_sq)
            dif_TheoObs_HaHb = numpy.fabs(norm_H6theo - I_obs_H6Hb) 
            dif_TheoObs_H6Hb_values.append(dif_TheoObs_HaHb)
            I_obs_HaHb = Idered_norm[Halpha_idx] / Idered_norm[Hbeta_idx]
            print ' ***    I_theo_HaHb =', I_theo_HaHb, '   I_obs_HaHb =', I_obs_HaHb
            diff_HaHb = numpy.fabs(I_theo_HaHb - I_obs_HaHb)
            diff_HaHb_values.append(diff_HaHb)
            EWabsHbeta_iterations = EWabsHbeta_iterations + 1
        # Second, variate C_Hbeta with EWabsHbeta fixed
        EWabsHbeta = EWabsHbeta_values[0]
        for C_Hbeta_iterations in range(0, number_iterations):
            print 'cHbeta = %0.3f        EWabsHbeta = % 0.3f' % (cHbeta, EWabsHbeta)
            print 'C_Hbeta_iterations =', C_Hbeta_iterations
            if I_theo_HaHb < I_obs_HaHb:
                C_Hbeta = C_Hbeta + C_Hbeta_increase
            elif I_theo_HaHb > I_obs_HaHb:
                C_Hbeta = C_Hbeta - C_Hbeta_increase
                if C_Hbeta < 0.0:
                    C_Hbeta = 0.00001            
            C_Hbeta_values.append(C_Hbeta)
            cHbeta = 0.434*C_Hbeta
            intensities = self.underlyingAbsCorr(return_Is=True)
            # If asked to, keep only the epositive fluses = emission lines
            if em_lines:
                catalog_lines, _, element_lines, _, _, _, norm_fluxes, _, norm_intensities, _ = find_emission_lines(rounded_wavelengths, self.element, self.ion, self.forbidden, self.howforb, _, flux, intensities, _, continuum)
            else:
                catalog_lines = rounded_wavelengths
                element_lines = self.element
                norm_fluxes = flux
                norm_intensities = intensities
            # Dered again and find the Chi_squared of that model
            _, Idered_norm, _ = self.Halpha2Hbeta_dered(cHbeta=cHbeta, fluxes=norm_fluxes, intensities=norm_intensities)
            I_obs_H6Hb = rounded_wavelengths[H6_idx] / rounded_wavelengths[Hbeta_idx]
            Chi_sq, norm_H6theo = self.find_Chi_of_CE(Idered_norm, theoCE, Hlines, verbose)
            Chi_sq_models.append(Chi_sq)
            dif_TheoObs_HaHb = numpy.fabs(norm_H6theo - I_obs_H6Hb) 
            dif_TheoObs_H6Hb_values.append(dif_TheoObs_HaHb)
            I_obs_HaHb = Idered_norm[Halpha_idx] / Idered_norm[Hbeta_idx]
            print ' ***    I_theo_HaHb =',I_theo_HaHb, '   I_obs_HaHb =', I_obs_HaHb
            diff_HaHb = numpy.fabs(I_theo_HaHb - I_obs_HaHb)
            diff_HaHb_values.append(diff_HaHb)
            C_Hbeta_iterations = C_Hbeta_iterations + 1
        # With all 41 models find the one that has the smallest Chi_sq
        #print 'Chi_sq_models:', Chi_sq_models
        minChi = min(Chi_sq_models)
        minChi_idx = Chi_sq_models.index(minChi)
        #print 'minChi', minChi, 'minChi_idx', minChi_idx
        # Now find the model that has the closest observed Hdelta/Hbeta ratio to the theoretical one
        min_dif_TheoObs_H6Hb = min(dif_TheoObs_H6Hb_values)
        min_dif_TheoObs_H6Hb_idx = dif_TheoObs_H6Hb_values.index(min_dif_TheoObs_H6Hb)
        #print ' min_dif_TheoObs_H6Hb = ', min_dif_TheoObs_H6Hb, '   min_dif_TheoObs_HaHb_idx', min_dif_TheoObs_H6Hb_idx
        # Calculate the final dereddend values but keep in mind that model 0 is the first reddening iteration. 
        # Note that if there were number_iterations = 10,
        # model 5 through 9 is where EWabsHbeta varied and C_Hbeta was fixed at model 0, and
        # model 10 though 14 is where C_Hbeta varied and EWabsHbeta was fixed at model 0.
        #print 'LENGTHS OF CHbeta, EWabsHbeta, and diff_HaHb_values lists: ', len(C_Hbeta_values), len(EWabsHbeta_values), len(diff_HaHb_values)
        if minChi_idx == min_dif_TheoObs_H6Hb_idx:
            print 'min indeces are the same!'
        tolerance = 0.005
        if (I_obs_HaHb < I_theo_HaHb+tolerance) and (I_obs_HaHb > I_theo_HaHb-tolerance):
            print ' VALUE WITHIN TOLERANCE!'
            EWabsHbeta = EWabsHbeta_values[0]
            C_Hbeta = C_Hbeta_values[minChi_idx - number_iterations]
        elif (I_obs_HaHb > I_theo_HaHb+tolerance) or (I_obs_HaHb < I_theo_HaHb-tolerance):
            print ' VALUE STILL FAR... LOOKING FOR ALTERNATIVE...'
            min_diff_HaHb = min(diff_HaHb_values)
            min_diff_HaHb_idx = diff_HaHb_values.index(min_diff_HaHb)
            if min_diff_HaHb_idx >= number_iterations:
                idx = min_diff_HaHb_idx - number_iterations
                EWabsHbeta = EWabsHbeta_values[0]
                C_Hbeta = C_Hbeta_values[idx]
            else:
                EWabsHbeta = EWabsHbeta_values[min_diff_HaHb_idx]
                C_Hbeta = C_Hbeta_values[0]
        #print 'Chi_sq_models', Chi_sq_models
        #print 'dif_TheoObs_H6Hb_values', dif_TheoObs_H6Hb_values
        cHbeta = 0.434*C_Hbeta
        intensities = self.underlyingAbsCorr(return_Is=True, EWabsHbeta=EWabsHbeta)
        #AdvOpscols_in_file = [self.wavelength, self.flambda, self.element, self.ion, self.forbidden, self.howforb, self.flux, self.errflux, self.percerrflx, self.intensity, self.errinten, self.percerrinten, self.ew, self.errew, self.percerrew]
        catalog_lines = rounded_wavelengths
        #wavs_lines = observed_wavelength
        #element_lines = self.AdvOpscols_in_file[2]
        #ion_lines = self.AdvOpscols_in_file[3]
        #forbidden_lines = self.AdvOpscols_in_file[4]
        #how_forbidden_lines = self.AdvOpscols_in_file[5]
        norm_fluxes = flux
        norm_intensities = intensities
        #EW_lines = self.AdvOpscols_in_file[12]
        #lines_info = [catalog_lines, wavs_lines, element_lines, ion_lines, forbidden_lines, how_forbidden_lines, norm_fluxes, norm_intensities, EW_lines]
        # Dered again and find the Chi_squared of that model
        _, norm_Idered, I_dered_norCorUndAbs = self.Halpha2Hbeta_dered(cHbeta=cHbeta, fluxes=norm_fluxes, intensities=norm_intensities)
        flambdas = find_flambdas(cHbeta, catalog_lines, I_dered_norCorUndAbs, norm_fluxes)
        print 'Flambda values per wavelength'
        for w, l in zip(catalog_lines, flambdas):
            print w, l
        #dereddening_info = [EWabsHbeta, C_Hbeta, norm_Idered, I_dered_norCorUndAbs, flambdas]
        # Determine uncertainties    
        percent_Iuncert = self.percerrinten
        absolute_Iuncert = []
        S2N = []
        for I, perIu in zip(norm_Idered, percent_Iuncert):
            errI = I * (perIu/100.0)
            absolute_Iuncert.append(errI)
            sn = I / errI
            S2N.append(sn)
        #uncertainties_info = [percent_Iuncert, absolute_Iuncert, S2N]
        I_obs_HaHb = norm_Idered[Halpha_idx] / norm_Idered[Hbeta_idx]
        print ' ***    I_theo_HaHb =', I_theo_HaHb, '   I_obs_HaHb =', I_obs_HaHb
        print''
        print 'First iteration of reddening correction:   EWabsHbeta = %0.3f  C_Hbeta = %0.3f' % (EWabsHbeta_values[0], C_Hbeta_values[0])
        print '    The best combination was:              EWabsHbeta = %0.3f  C_Hbeta = %0.3f' % (EWabsHbeta, C_Hbeta)
        print '                                                     this means cHbeta = %0.3f' % (cHbeta)
        # Now define the columns to be written in the text file
        cols_2write_in_file = [catalog_lines, flambdas, self.element, self.ion, self.forbidden, self.how_forbidden, norm_fluxes, self.errflux, self.percerrflx, norm_Idered, absolute_Iuncert, percent_Iuncert, self.ew, self.errew, self.percerrew]
        return (C_Hbeta, cols_2write_in_file)
        
    def get_avgTandDelta4helio10(self, Otot, Op, Opp, TOp, TOpp):
        '''This function gets the temperature to use in the HELIO10 program.
        # Otot =  sum of O ionic abundances
        # Op = O+ = O2 abundance
        # Opp = O++ = O3 abundance
        # TOp = Temperature of O2
        # TOpp = Temperature of O3
        '''        
        # Assuming that N(O++) + N(O+) = 100%
        perOp = Op*1.0 / Otot
        perOpp = Opp*1.0 / Otot
        # The temperature is T = T(O++) * %  +  T(O+) * %
        avgT = (TOp*perOp + TOpp*perOpp) / 10000.0
        print ' fraction of O+ = ', perOp, '     fraction of O++ = ', perOpp
        # now get the value of delta using the equation given by Antonio Peimbert
        delt = 0.65 * (Op/Otot)
        return avgT, delt
        
    def get_logIs4_helio10(self):
        ''' This function obtains the log(Intensities) +- err for the helium lines needed in HELIO10. '''
        object_name = self.object_name
        use_Chbeta = self.use_Chbeta
        # Go into the directory of each object. Files assumed to be in: /Users/home_direcotry/Documents/AptanaStudio3/src/
        results_path = "../results/"
        # but just to make sure we are in the right place, get and work with full paths
        full_results_path = os.path.abspath(results_path)
        # Go into the object folder and get the file
        results4object_path = os.path.join(full_results_path, object_name)
        if use_Chbeta:
            if self.used_measured_info:
                RedCor_file = os.path.join(results4object_path, object_name+"_measuredLI_RedCor_CHbeta.txt")
            else:
                RedCor_file = os.path.join(results4object_path, object_name+"_RedCor_CHbeta.txt")
        else:
            if self.used_measured_info:
                RedCor_file = os.path.join(results4object_path, object_name+"_measuredLI_RedCor_Ebv.txt")
            else:
                RedCor_file = os.path.join(results4object_path, object_name+"_RedCor_Ebv.txt")
        # Load the data of interest
        wavs, intensities, errs = numpy.loadtxt(RedCor_file, skiprows=(1), usecols=(0,9,10), unpack=True)
        w4He = [2945.0, 3188.0, 3614.0, 3819.0, 3889.0, 3965.0, 4026.0, 4121.0, 4388.0, 4438.0, 
                4471.0, 4713.0, 4922.0, 5016.0, 5048.0, 5876.0, 6678.0, 7065.0, 7281.0, 9464.0]
        rounded_wavs = round_wavs(wavs)
        logInts = []
        logErrs = []
        wavsHe = []
        for w, I, e in zip(rounded_wavs, intensities, errs):
            if w in w4He:
                if I < 0.0:
                    logI = -1.000
                    loge = 99.000
                else:
                    logI = numpy.log10(I/100.0)
                    loge = 0.5 * numpy.log10((I+e)/(I-e))
                    if w == 3889.0:
                        logI = numpy.log10((I-10.6)/100.0)
                        loge = 0.5 * numpy.log10((I-10.6+e)/(I-10.6-e))
                wavsHe.append(w)
                logInts.append(logI)
                logErrs.append(loge)
        return wavsHe, logInts, logErrs 
                
        
    def t2_RLs(self, dens, Opp, Cpp):
        '''This function determines the value of t^2 using recombination lines.'''
        print '\n***   Recombination Lines...' 
        # Oxygen recombination lines (ORLs)
        # The O++ abundance can be determined from the equation of recombination lines intensities:
        # dens(O++)    alpha_effHbeta * h nu_ORLs * dens(O++) * n_e * V_eff * integral(4*pi*radius^2)
        # --------- = --------------------------------------------------------------------------------
        # dens(H+)     alpha_effORLs * h nu_Hbeta * dens(H+)  * n_e * V_eff * integral(4*pi*radius^2)
        # the electron densities, volumes, and the integrals cancel out because the observed region is the 
        # same, and since nu is inversely proportional to lambda, hence:
        #              alpha_effHbeta * lambda_Hbeta * dens(O++)
        #           = ------------------------------------------- ,   now solving for dens(O++)/dens(H+)
        #              alpha_effHbeta * lambda_ORLs * dens(O++)
        # dens(O++)/dens(H+) = I(ORLs)/I(Hbeta) * alpha_eff(Hbeta)/alpha_eff(ORLs) * lambda_Hbeta/lambda_ORLs
        # The value of alpha_eff(Hbeta) can be obtained from Osterbrock & Ferland (2006, pp22, Table 2.1)
        # and Storey & Hummer (1995, MNRAS, 272, 41):  2.59e-13 cm^3 s^-1 for T=10,000
        # but our temperature differs from 10,000 K, hence to find alpha_effHbeta=aHbeta*t_4^b:
        t_4 = self.te_high[0] / 10000.0 # temperatures are in Kelvin
        err_t4 = self.te_high[1] / 10000.0
        aHbeta = 3.03e-14  # This alpha_Hbeta(10,000)
        bHbeta = -0.9033 #= log[alpha_Hbeta(20,000)/alpha_Hbeta(10,000)] / log(20,000/10,000) = log[1.62e-14/3.03e-14]/log(20000/10000)
        alpha_effHbeta =  aHbeta * t_4**bHbeta # cm^3 s^-1
        # The value of alpha_eff(ORLs) can be obtained from Storey (1994, A&A, 282, 999), equation 7:
        # alpha_effORLs = 10^-14 * a * t_4^b (1 + c * (1 - t_4) + d * (1-t_4^2))
        # The values for constants a, b, c, and d can be obtained from table 3 for 4652 Angstroms, Case B (Case A is commented): 
        a = 36.2 #34.9
        b = -0.736 #-0.749
        c = 0.033 #0.023
        d = 0.077 #0.074
        alpha_effORLs = 1e-14 * a * t_4**b * (1 + c * (1 - t_4) + d * (1 - t_4**2))
        err_alpha_effORLs = 1e-14 * a * err_t4**b * (1 + c * (1 - err_t4) + d * (1 - err_t4**2))
        # Now find the sum of the observed lines of the multiplet 1:
        # lambda_ORLs / lambda_Hbeta = ((4639+4642+4649+4651)/4) / 4861 = 4645.25/4861 = 0.956
        lambdaORLs_over_lambdaHbeta = 0.956
        # assuming that each line is right within 3.5 angstroms,
        err_obslines = 1/4.0 * numpy.sqrt(4.0 * 3.5**2)
        err_lambdaORLs_over_lambdaHbeta = lambdaORLs_over_lambdaHbeta * numpy.sqrt((err_obslines/4645.25)**2 + (3.5/4861.0)**2)
        # Now obtain the sum of the intensities of the ORLs observed
        # Using equations 8 through 11 from Peimbert, Peimbert, & Ruiz (2005, ApJ, 634, 1056)
        n_crit = 2800.0 #+-500  for HII regions,  for PN it is 1325+-300
        dens_ratio = 1.0 + dens[0]/n_crit
        I4651p74_over_Isum = (0.101 + 0.144 / dens_ratio) * 0.844
        I4639p62p96_over_Isum = (0.201 + 0.205 / dens_ratio) * 0.455       
        I4642p76_over_Isum = (0.301 - 0.057 / dens_ratio) * 0.742 
        I4649_over_Isum = 0.397 - 0.292 / dens_ratio 
        # These equations give out a percentage of the total, IsumORLs, which is assumed to be 100%. Hence, the 
        # observed doublets have to be summed up:
        # blending of 4638.86 and 4641.81, known as 4639+42
        percent_4639p42 = I4639p62p96_over_Isum + I4642p76_over_Isum
        # blending of 4649.13 and 4659.84, known as 4649+51
        percent_4649p51 = I4649_over_Isum + I4651p74_over_Isum
        # Now, the sum of these two gives the percentage observed of the total multiplet, percent_ORLs_obs
        percent_ORLs_obs = percent_4639p42 + percent_4649p51
        # Now get the actual values for the intensities from the measured data:
        # The multiplet 1 has eight lines (4639, 42, 49, 51, 62, 74,76, and 4696) of which only two doublets are usually observed:
        # the blending of 4638.86 and 4641.81, known as 4639+42 and the blending of 4649.13 and 4659.84, known as 4649+51
        #AdvOpscols_in_file = [self.wavelength, self.flambda, self.element, self.ion, self.forbidden, self.howforb, self.flux, self.errflux, self.percerrflx, self.intensity, self.errinten, self.percerrinten, self.ew, self.errew, self.percerrew]
        I_4639p42 = [0.0, 0.0]
        I_4649p51 = [0.0, 0.0]
        for w, F, err in zip(self.AdvOpscols_in_file[0], self.AdvOpscols_in_file[6], self.AdvOpscols_in_file[7]):
            if w == 4640.0:
                if F > 0.0:
                    I_4639p42 = [F, err]
            elif w == 4650.0:
                if F > 0.0:        
                    I_4649p51 = [F, err]
            elif w == 4861.33:
                Hbeta = [F, err]
        I_ORLs_obs = I_4639p42[0] + I_4649p51[0]        
        if I_ORLs_obs == 0.0:
            densOpp_over_densHp = 0.0
            err_densOpp_over_densHp = 0.0
            log_densOpp_over_densHp = 0.0
            logeleerr = 0.0
        else:
            err_I_ORLs_obs = numpy.sqrt(I_4639p42[1]**2 + I_4649p51[1]**2)
            # So now we have that I_ORLs_obs +- err_I_ORLs_obs represents percent_ORLs_obs of the I_total of the multiplet, therefore
            # we have to multiply it for 1/percent_ORLs_obs in order to correct for the ~30% we do not see:
            # lets call I_ORLs_obs * (1/percent_ORLs_obs) = ORLs
            ORLs = I_ORLs_obs * (1.0/percent_ORLs_obs)
            IORLs_over_Hbeta =  ORLs / Hbeta[0]
            # errors in the intensity of the multiplet 1
            err_ORLs = err_I_ORLs_obs * 1.0/percent_ORLs_obs
            err_IORLs_over_Hbeta = IORLs_over_Hbeta * numpy.sqrt((err_ORLs/ORLs)**2 + (Hbeta[1]/Hbeta[0])**2) 
            # Now we can proceed with the O++ abundance determination:
            densOpp_over_densHp = IORLs_over_Hbeta * alpha_effHbeta/alpha_effORLs * lambdaORLs_over_lambdaHbeta
            # error in ionic abundance
            err_densOpp_over_densHp = densOpp_over_densHp * lambdaORLs_over_lambdaHbeta * numpy.sqrt((err_IORLs_over_Hbeta/IORLs_over_Hbeta)**2 + ((alpha_effHbeta/err_alpha_effORLs)/(alpha_effHbeta/alpha_effORLs))**2)
            log_densOpp_over_densHp = 12 + numpy.log10(densOpp_over_densHp)
            logeleerr = err_densOpp_over_densHp / (2.303 * densOpp_over_densHp)
            print 't_4 = ', t_4        
            #print 'percent_4639p42 =', I4639p62p96_over_Isum, '+', I4642p76_over_Isum
            #print 'percent_4649p51 =', I4649_over_Isum, '+', I4651p74_over_Isum
            print 'percent_ORLs_obs = ', percent_4639p42, '+', percent_4649p51, '=', percent_ORLs_obs
            print 'I_ORLs_obs * (1.0/percent_ORLs_obs) = ', I_ORLs_obs, '*', '1/',percent_ORLs_obs, '=', ORLs
            print 'IORLs_over_Hbeta =', IORLs_over_Hbeta
            print 'alpha_effHbeta/alpha_effORLs = ', alpha_effHbeta, '/', alpha_effORLs, '=', alpha_effHbeta/alpha_effORLs
            #print 'lambdaORLs_over_lambdaHbeta = ', lambdaORLs_over_lambdaHbeta
        print '\n The O++ abundance determined from RLs is: %0.3e +- %0.3e' % (densOpp_over_densHp, err_densOpp_over_densHp)
        print '                           12+log(O++/H+) = %0.2f +- %0.2f' % (log_densOpp_over_densHp, logeleerr)
        # To find approximate value of ADF, we divide the abundances found with CELs/RLs:
        approxADF = densOpp_over_densHp/Opp
        print ' The approximate value of ADF is: '
        print '            O++_RLs / O++_CELs =', densOpp_over_densHp, '/', Opp, ' =', approxADF
 
        '''
        # An alternative method (not so accurate) is using a proportionallity relation with the well known 30 Dor:
        # [ I(ORLmultiplet)/I(Hbeta) / N(O++)/N(H+) ]_object = [I(ORLmultiplet)/I(Hbeta) / N(O++)/N(H+) ]_30Dor * [T^0.9/T^0.8],
        # ORLabund_obj = ORLabund_30Dor * TRL
        # where the last term, TRL, is the temperature dependence of the oxygen lines to hydrogen lines, 
        # TRL = T^0.9 / T^0.8 = T^0.1, since it is very small, we will disregard it.  
        ORLsumoverHbeta_30Dor = 0.003385   # sum of the 8 oxygen recombination lines in 30 Dor
        OppoverHp_30Dor = 0.0002888        # abundance of O++ with respect to H+ in 30 Dor
        ORLabund_30Dor = ORLsumoverHbeta_30Dor/OppoverHp_30Dor   # this turns into a constant for the equation, hence:
        ORLabund_object = ORLs * (1.0/ORLabund_30Dor)
        print '       --> other method O++ abundance:', ORLabund_object
        '''
        
        # Carbon recombination lines (CRLs)
        # In the optical, there is only one line that is measurable:
        #CRLs = 4267.0, 7231.0
        # From A. R. Davey, P. J. Storey, & R. Kisielius (2000, A&ASS, 142,85), ecuation 3:
        # alpha_recC = 10e-4 * a * t_4**f * (1 + b*(1-t_4) + c*(1-t_4)**2 + d*(1-t_4)**3)
        # The values of constants a, b, c, d, and f were taken from table 5 for T_e = 5000 - 20,000K, Case B
        a = 27.586
        b = -0.055
        c = -0.039
        d = -0.208
        f = -1.1416
        err = 0.15 #%
        alpha_CRL = 1e-14 * a * t_4**f * (1 + b*(1-t_4) + c*(1-t_4)**2 + d*(1-t_4)**3)
        err_alpha_CRL =  1e-14 * a * err_t4**f * (1 + b*(1-err_t4) + c*(1-err_t4)**2 + d*(1-err_t4)**3)
        # Following the same recipy as for oxygen,
        lambdaCRL_over_lambdaHbeta = 4267.15/4861.33
        F_4267 = [0.0, 0.0]
        densCpp_over_densHp = 0.0
        err_densCpp_over_densHp = 0.0
        log_densCpp_over_densHp = 0.0
        logeleerr = 0.0
        for w, F, err in zip(self.AdvOpscols_in_file[0], self.AdvOpscols_in_file[6], self.AdvOpscols_in_file[7]):
            if w == 4267.15:
                if F > 0.0:
                    F_4267 = [F, err]
                else:
                    print  'This object had no 4267 CII line...'
        percent_CRL_obs = 1.
        CRLs = F_4267[0] * (1.0/percent_CRL_obs)
        CRLs_over_Hbeta =  CRLs / Hbeta[0]
        if F_4267[0] != 0.0:
            err_CRLs_over_Hbeta = numpy.sqrt((F_4267[1]/F_4267[0])**2 + (Hbeta[1]/Hbeta[0])**2)
            densCpp_over_densHp = CRLs_over_Hbeta * alpha_effHbeta/alpha_CRL * lambdaCRL_over_lambdaHbeta
            err_densCpp_over_densHp = densCpp_over_densHp * lambdaCRL_over_lambdaHbeta * numpy.sqrt((err_CRLs_over_Hbeta/CRLs_over_Hbeta)**2 + ((alpha_effHbeta/err_alpha_CRL)/(alpha_effHbeta/alpha_CRL))**2) 
            log_densCpp_over_densHp = 12 + numpy.log10(densCpp_over_densHp)
            logeleerr = err_densCpp_over_densHp / (2.303 * densCpp_over_densHp)
        print '\n'
        print 'CRLs =', F_4267[0], '* 1.0/', percent_CRL_obs, '=', F_4267[0] * (1.0/percent_CRL_obs)
        print 'CRLs_over_Hbeta =', CRLs_over_Hbeta
        print 'alpha_effHbeta/alpha_CRL = ', alpha_effHbeta, '/', alpha_CRL, '=', alpha_effHbeta/alpha_CRL
        print 'lambdaCRL_over_lambdaHbeta = ', lambdaCRL_over_lambdaHbeta
        print '\n The C++ abundance determined from RLs is: %0.3e  +- %0.3e' % (densCpp_over_densHp, err_densCpp_over_densHp)
        print '                           12+log(C++/H+) = %0.2f +- %0.2f' % (log_densCpp_over_densHp, logeleerr)
        approxADF_C = densCpp_over_densHp / Cpp
        print ' The approximate value of ADF is: '
        print '            C++_RLs / C++_CELs =', densCpp_over_densHp, '/', Cpp, ' =', approxADF_C
    
        
    def perform_advanced_ops(self, forceTeH, forceTeL, forceNe, theoCE):
        lines_pyneb_matches = self.writeRedCorrFile()
        self.get_tempsdens()
        self.define_TeNe_HighLow_andVLow(forceTeH, forceTeL, forceNe)
        if self.use_Chbeta:
            print '    Performing second iteration of extinction correction ... \n'
            CHbeta, cols_2write_in_file = self.redcor2(theoCE, verbose=False, em_lines=False)
            write_RedCorfile(self.object_name, self.tfile2ndRedCor, cols_2write_in_file)
            # cols_2write_in_file contains the following columns:
            # catalog_wavelength, flambdas, element, ion, forbidden, how_forbidden, norm_fluxes, errflux, percerrflx, 
            # norm_Idered, percent_Iuncert, absolute_Iuncert, EW, err_EW, perrEW
            catalog_wavelength, flambdas, element, ion, _, _, norm_fluxes, errflux, percerrflx, norm_intenties, errinten, _, _, _, _ = cols_2write_in_file
            data2use = [CHbeta, catalog_wavelength, flambdas, element, ion, norm_fluxes, errflux, percerrflx, norm_intenties, errinten]
        else:
            data2use = None
        self.get_iontotabs(data2use)
        return lines_pyneb_matches
