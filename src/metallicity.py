from __future__ import division

import os
import pyfits
import numpy
import string
import copy
import math
import pyneb as pn
import PIL.Image as Image
from uncertainties import unumpy
from science import spectrum
from matplotlib import pyplot
from collections import OrderedDict
 

'''
This program contains various classes that together determine the object's metallicity.
'''

def round_cat_wavs(catalog_wavelength):
    ### Round all catalog lines to make it easier to find lines
    rounded_catalog_wavelength = []
    for item in catalog_wavelength:
        rw = numpy.round(item)
        rounded_catalog_wavelength.append(rw)
    return rounded_catalog_wavelength

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
    Hb_idx = rounded_catalog_wavelength.index(4861.0)
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
            norm_flux = flux[i] / flux[Hb_idx] * 100.
            positive_normfluxes.append(norm_flux)
            normI = intensities[i] / intensities[Hb_idx] * 100.
            positive_norm_intensities.append(normI)
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
    def __init__(self, redlaw, cols_in_file, I_theo_HaHb, EWabsHbeta, cHbeta, av, ebv, do_errs=None):
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
        
    def underlyingAbsCorr(self):
        catalog_wavelength = self.catalog_wavelength
        continuum = self.continuum
        flux = self.flux
        EWabsHbeta = self.EWabsHbeta
        ### Step 1 of first iteration of reddening correction: Assume that there is no collisional excitation
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
    
    def Halpha2Hbeta_dered(self, av, ebv):
        ### Function to dered and obtain the Halpha/Hbeta ratio nicely printed along with the unreddend values to compare.
        redlaw = self.redlaw
        cHbeta = self.cHbeta
        I_theo_HaHb = self.I_theo_HaHb
        catalog_wavelength = self.catalog_wavelength
        rounded_catalog_wavelength = round_cat_wavs(catalog_wavelength)
        observed_wavelength = self.observed_wavelength
        element = self.element
        ion = self.ion
        forbidden = self.forbidden 
        how_forbidden = self.how_forbidden
        flux = self.flux
        continuum = self.continuum
        EW = self.EW
        intensities = self.intensities_corr_undelyingAbs
        norm_data = normalize_lines(rounded_catalog_wavelength, element, ion, forbidden,
                                how_forbidden, observed_wavelength, flux, intensities, EW, continuum)
        catalog_lines, _, _, _, _, _, normfluxes, _, norm_intensities, _ = norm_data
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
        for w, nF, nI in zip(catalog_lines, normfluxes, norm_intensities):   
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
        Halpha = normfluxes[Halpha_idx]
        Hbeta = normfluxes[Hbeta_idx]
        raw_ratio = Halpha/Hbeta
        I_Halpha = Idered[Halpha_idx]
        I_Hbeta = Idered[Hbeta_idx]
        I_obs_HaHb = I_Halpha/I_Hbeta
        print ''
        print 'cHbeta = %0.5f' % cHbeta
        print ' Intensities corrected for reddening and underlying absorption.'
        print '            Using', redlaw, '                   Normalized fluxes before extinction correction'
        print catalog_lines[Halpha_idx], '    ', I_Halpha, '                  ', Halpha
        print catalog_lines[Hbeta_idx], '    ', I_Hbeta, '                          ', Hbeta
        print 'theoretical ratio Ha/Hb = %0.2f' % (I_theo_HaHb)
        #print '      observed Ha/Hb = %0.3f           raw Ha/Hb = %0.3f' % (I_obs_HaHb, raw_ratio)
        print '      observed Ha/Hb =', numpy.round(I_obs_HaHb, decimals=2), '           raw Ha/Hb =', numpy.round(raw_ratio, decimals=2)
        self.Idered = Idered
        self.I_dered_norCorUndAbs = I_dered_noCorUndAbs
        self.normfluxes = normfluxes
        return self.normfluxes, self.Idered, self.I_dered_norCorUndAbs
    
    def get_uncertainties(self):
        errs_Flx, errs_EW, cont_errs = self.errs_list
        errs_Idered = []
        perc_errs_I_dered = []
        errs_normfluxes = []
        perc_errs_normfluxes = []
        # Find observed Hbeta 
        Hbeta_idx = self.catalog_wavelength.index(4861.330)
        Hbeta = self.flux[Hbeta_idx]
        err_Hbeta = errs_Flx[Hbeta_idx]
        for w, F, I, nF, eF, eC in zip(self.catalog_wavelength, self.flux, self.Idered, self.normfluxes, errs_Flx, cont_errs):
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
            tot_err_nF = numpy.abs(nF) * numpy.sqrt((err_Hbeta/Hbeta)*(err_Hbeta/Hbeta) + (eF/F)*(eF/F) )#+ (eC/C)*(eC/C) )
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
            perc_errIdered = (errIdered * 100.) / numpy.abs(I)
            perc_errs_I_dered.append(perc_errIdered)
            #print w, ' NormFlux=', nF, ' err=',tot_err_nF, ' err%=', perc_tot_err_nF, '   I=', I, ' err=', errIdered, ' err%=', perc_errIdered
        #print 'min(perc_errs_I_dered)', min(perc_errs_I_dered)
        return errs_normfluxes, perc_errs_normfluxes, errs_Idered, perc_errs_I_dered
    
    def do_ops(self):
        self.underlyingAbsCorr()
        normfluxes, Idered, I_dered_norCorUndAbs = self.Halpha2Hbeta_dered(self.av, self.ebv)
        if self.errs_list != None:
            errs_normfluxes, perc_errs_normfluxes, errs_Idered, perc_errs_I_dered = self.get_uncertainties()
            return normfluxes, Idered, I_dered_norCorUndAbs, errs_normfluxes, perc_errs_normfluxes, errs_Idered, perc_errs_I_dered
        else:
            return normfluxes, Idered, I_dered_norCorUndAbs

    
class AdvancedOps:
    def __init__(self, object_name, cHbeta, case, writeouts=False, verbose=False):
        # Inputs:
        self.object_name = object_name
        self.cHbeta = cHbeta
        self.case = case                        # this is the Case used through out the class
        self.verbose = verbose                  # if True print midpoints in order to know what is the code working on
        self.writeouts = writeouts              # write or not text files with outputs (temps, densities, and abundances)
        # Variables defined in the class
        self.lines_pyneb_matches = []
    
    def writeRedCorrFile(self):
        ''' This function writes a text file in a pyneb readable format... Necessary to find the temperatures and densities. '''
        object_name = self.object_name
        verbose = self.verbose
        input_file = object_name+'_RedCor.txt'
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
        cols_in_file = [self.wavelength, self.flambda, self.element, self.ion, self.forbidden, self.howforb, self.flux, self.errflux, self.percerrflx, self.intensity, self.errinten, self.percerrinten, self.ew, self.errew, self.percerrew]
        # Read the info from the text file that contains IDs, dered info, and errors (this is the object_redCor.txt)
        cols_in_file = spectrum.readlines_from_lineinfo(path_inputfile, cols_in_file)
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
        print 'Pyneb IDs found!  Lines written in file  %s' % self.pynebIDstxt
        return self.lines_pyneb_matches 

    def get_tempsdens(self):
        ''' 
        This is a VERY big function that unfortunately cannot be broken... It determines ALL posible temperatures and densities 
        that are avaiable to veryfy with IRAF. The abundance determination can be turned off by setting iontotabs to False.
        '''
        out_file = self.object_name+'_TempDens.txt'
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
            elif line.wave == 5755:         # N 2
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
                I_4959 = line.corrIntens
                print '4959 has an intensity of', I_4959
            elif line.wave == 5007:
                I_5007 = line.corrIntens
                print '5007 has an intensity of', I_5007
            elif line.wave == 3726:
                I_3726 = line.corrIntens
                print '3726 has an intensity of', I_3726
            elif line.wave == 3729:
                I_3729 = line.corrIntens
                print '3729 has an intensity of', I_3729
            elif line.wave == 3727:
                I_3727 = line.corrIntens
                print '3727 has an intensity of', I_3727
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
        try:
            self.O3 = pn.Atom("O", "3")
            #self.tem_diag_O3 = '(L(4959)+L(5007)) / L(4363)'
            #print 'ratio of O3 = ', (I_4959+I_5007)/I_4363
            self.tem_diag_O3 = 'L(5007) / L(4363)'
            self.strongO3 = I_5007
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
            self.O2 = pn.Atom("O", "2")
            tem_diag_O2 = '(L(3726)+L(3729)) / (L(7329) + L(7330))'
            self.temO2 = self.O2.getTemDen(I_3727/I_7330, den=100.0, to_eval=tem_diag_O2) 
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
        # Make sure that the temperatures and densities file closes properly
        if self.writeouts:
            outf.close()
        print ''
        
    def get_iontotabs(self, forceTeH=None, forceNe=None):
        ''' With the available temperatures determine ionic and total abundances
        forceTeH = specific high ionization region temperature to be used to determine abundances.
                    It can be either one value or a list of value and error.
        forceNe = specific density to be used 
                    It can be either one value or a list of value and error.
        '''
        # Define the high and lo ionization zones temperatures and densities
        print 'Temperatures being used for estimation of abundances:'
        if (forceTeH != None):
            if type(forceTeH) is not list:
                # error is not defined, use a default value of 25% of the given temperature
                perr = 0.15
                teH = forceTeH
                teHerr = teH + (teH * perr)
                teL = forceTeH * 0.9         # 90% of the high temperature
                teVL = forceTeH * 0.8        # 80% of the high temperature
            else:
                teH = forceTeH[0]
                # error is defined, determine the percentage to make it absolute error
                if forceTeH[1] < forceTeH[0]:
                    abserr = forceTeH[1]
                else:
                    abserr = forceTeH[1] - teH
                perr = (abserr * 1.0) / teH
                teHerr = teH + abserr
                teL = forceTeH[0] * 0.9         # 90% of the high temperature
                teVL = forceTeH[0] * 0.8        # 80% of the high temperature
            te_high = [teH, teHerr]
            teLerr = teL + (teL*perr)
            te_low = [teL, teLerr]
            teVLerr = teVL + (teVL * perr)
            te_verylow = [teVL, teVLerr]
            print '   High ionization degree temperature forced to be:      ', te_high[0], '+-', te_high[1]-te_high[0]
            print '   Low ionization degree temperature    =   Te_high*0.9: ', te_low[0], '+-', te_low[1]-te_low[0]
            print '   Very Low ionization degree temperature = Te_high*0.8: ', te_verylow[0], '+-', te_verylow[1]-te_verylow[0]
            print '   * Error is calculated from percentage error in Te_high:    %0.2f' % (perr*100.0), '%'
        else:
            if math.isnan(self.temO2[0]):
                te_low = [9000.000, 9500.0]
            else:
                te_low = self.temO2
                print '   Te_low (O2) =', te_low
                
            if math.isnan(self.TO3[0]):
                te_high = [10000.0, 10500.0]
                print '   Te[O 3] not available, using default value:   te_high = 10,000 +- 500.0'
            else:
                te_high = self.TO3
                print '   Te_high (O3) =', te_high
                # make sure the loz ionization temperature has a lowe temperature than the high one
                if (te_low[0] > te_high[0]) or math.isnan(self.temO2[0]):
                    print '   Te_low (O2-Garnet92) =', self.TO2gar
                    te_low = self.TO2gar
            if te_low[0] == 9000.00:
                print '   Te[O 2] not available, using default value:   te_low = 9,000 +- 500.0'
                
            if math.isnan(self.TS3[0]):
                te_verylow = [8000.0, 8500.0]
            else:
                te_verylow = self.TS3
            
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
            if math.isnan(self.denS2[0]):
                dens = [100.0, 150.0]
                print 'ne[S 2] not available, using default value: 100.0 +- 50.0 \n'
            else:
                #ne = self.denO2
                ne = self.denS2
                if numpy.abs(ne[1] - ne[0]) < 20.0:
                    ne[1] = ne[0] + ne[0]*0.1
                dens = ne
                print 'ne =', dens, '\n'
                
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
         
        # ions of zones of high and medium ionization degree are combined
        ab_high = ['Ar3', 'Ar4', 'Ar5', 'C2', 'C3', 'Ca5', 'Cl2', 'Cl3', 'Cl4', 'Fe3', 'K4', 'K5', 'Mg5', 
                   'N3', 'Ne3', 'Ne4', 'Ne5', 'Na4', 'Na6', 'Ne3', 'Ne5', 'O3', 'S3']
        ab_low = ['Al2', 'N1', 'N2', 'O1', 'O2', 'S2', 'Si2', 'Si3']
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
        strong_lines_list = ['2670', '7751', '4740', '6435', '2328', '1907', '6087', '9124', '5538', '7531', '4987',
                             # K4     K5      Mg5     N1      N2      N3      Na4     Na6     Ne3     Ne4     Ne5     
                             '6796', '4163', '2783', '5200', '6548', '1752', '3362', '2970', '3869', '2425', '3346',
                             # O1     O2     O3      S2      S3      Si2     Si3
                             '6300', '3727', '5007', '6731', '9531', '2345', '1892']
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
            print ion, ab, logab
        
        #Oab = self.O3.getIonAbundance(4.15, te_high, dens, to_eval='L(1661)')
        #print Oab        

        # Now calculate the ionic abundances of C^{++}/O^{++}, N^{++}, and C/O according to Garnett et al. (1995)
        # Equation 2 for C^{++}/O^{++}
        tc = te_high[0]/10000.0         # central temp
        tpluserr = te_high[1]/10000.0   # central temp plus error
        if self.I_1666[0] < 0.0:
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
        #if N2toO2 < 0.0: 
        #    N2toO2 = 0.0
        print '\nC2/O2 = ', C2toO2, '+-', C2toO2_err
        print 'N2/O2 = ', N2toO2, '+-', N2toO2_err,'\n'

        # Write results in text file
        if self.writeouts:
            out_file = self.object_name+'_IonicTotAbundances.txt'
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
            print >> outf, '# RATIOS -- using equations from Garnet et al. (1995)'
            print >> outf, ('{:<6} {:>15} {:>12.3f} {:>10.3f}'.format('C2/O2', '1909/1666', C2toO2, C2toO2_err))
            print >> outf, ('{:<6} {:>15} {:>12.3f} {:>10.3f}'.format('N2/O2', '1752/1666', N2toO2, N2toO2_err))
            n2 = sorted_atoms.index('N2')
            o2 = sorted_atoms.index('O2')
            n2too2 = totabs_ions_list[n2][0]/totabs_ions_list[o2][0]
            n2too2_err = n2too2 * numpy.sqrt((totabs_ions_list[n2][1]/totabs_ions_list[n2][0])**2 + (totabs_ions_list[o2][1]/totabs_ions_list[o2][0])**2)
            print >> outf, ('{:<6} {:>15} {:>12.3f} {:>10.3f}'.format('N2/O2', '6548/5007', n2too2, n2too2_err))
        
        ### TOTAL abundances
        icf = pn.ICF()
        icf.getAvailableICFs() # get all available ICFs
        icf.printAllICFs(type_=['HII']) # print out the ICFs only for HII regions
        #r = icf.getReference('Ial06_22a')
        self.atom_abun = OrderedDict()
        #atom_abun = {}
        for ion, ionabund in zip(sorted_atoms, totabs_ions_list):
            self.atom_abun[ion] = ionabund
        print self.atom_abun
        # use a specific recipy to determine abundances
        #elem_abunTPP85 = icf.getElemAbundance(self.atom_abun, icf_list=['TPP85']) 
        elem_abunIal06_22a = icf.getElemAbundance(self.atom_abun, icf_list=['Ial06_22a'])
        #icf.getExpression('Ial06_22a') # returns the analytical expression of the icf identified by the label TPP85
        #icf.getReference('Ial06_22a') # returns the bibliographic reference of the source paper
        #icf.getURL('Ial06_22a') # returns the ADS URL of the source paper
        print elem_abunIal06_22a['Ial06_22a']

        #elements = ['Al', 'Ar', 'C', 'Ca', 'Cl', 'K', 'Mg', 'N', 'Na', 'Ne', 'O', 'S', 'Si']
        elem_abun = OrderedDict()
        # Oxygen
        #Otot_test = self.atom_abun['O2'] + self.atom_abun['O3']
        Otot = self.atom_abun['O2'][0] + self.atom_abun['O3'][0]
        O23sq = self.atom_abun['O2'][1]**2 + self.atom_abun['O3'][1]**2
        Ototerr = numpy.sqrt(O23sq)
        O_errp = (Ototerr/Otot)*100
        elem_abun['O'] = [Otot, Ototerr]
        print ' Assuming that  Otot = O+ + O++:  ', Otot, '+-', Ototerr, '( which is ', O_errp, ' %)'
        #print 'Otot_test[1]', Otot_test[1], '    Ototerr', Ototerr
        print 'O_tot = %0.2f +- %0.2f' % (12+numpy.log10(Otot), numpy.log10((100+O_errp) / (100-O_errp))/2.0)
        
        # Nitrogen
        N_tot = (Otot * self.atom_abun['N2'][0]) / self.atom_abun['O2'][0]
        print ' Assuming ICF(N) from Peimbert+Costero 69 = (Otot / O+) * N+'
        N_toterr = numpy.sqrt(O23sq**2 * (self.atom_abun['O3'][0]/Otot)**2 * (self.atom_abun['N2'][0]/self.atom_abun['O2'][0])**2 +
                             self.atom_abun['N2'][1]**2)
        elem_abun['N'] = [N_tot, N_toterr]
        Nicf = Otot / self.atom_abun['O2'][0]
        N_errp = (N_toterr/N_tot)*100
        #print 'Nicf=%0.2f ,   Ntot = %0.2f +- %0.2f' % (Nicf, 12+numpy.log10(Ntot), numpy.log10((Ntot+Ntoterr) / (Ntot-Ntoterr))/2.0)
        print 'ICF(N)=%0.2f ,   N_tot = %0.2f +- %0.2f' % (Nicf, 12+numpy.log10(N_tot), numpy.log10((100+N_errp) / (100-N_errp))/2.0)
        
        # Neon
        print ' Assuming ICF(Ne) from Peimbert+Costero 69 =  (Otot / O++) * Ne++'
        Ne_icf = Otot / self.atom_abun['O3'][0]
        Ne_tot = self.atom_abun['Ne3'][0] * Ne_icf
        Ne_toterr = numpy.sqrt((O23sq * (self.atom_abun['O2'][0]/Otot)**2) * (self.atom_abun['Ne3'][0]/self.atom_abun['O3'][0])**2 +
                               self.atom_abun['Ne3'][1]**2)
        elem_abun['Ne'] = [Ne_tot, Ne_toterr]
        Ne_errp = (Ne_toterr/Ne_tot)*100
        print 'ICF(Ne)=%0.2f ,   Ne_tot = %0.2f +- %0.2f' % (Ne_icf, 12+numpy.log10(Ne_tot), numpy.log10((100+Ne_errp) / (100-Ne_errp))/2.0)
        
        # Sulphur
        print ' Assuming ICF(S) from Garnett 89:  (S+ + S++)/Stot = [1 - (1 - O+/Otot)^alpha] ^ 1/alpha'
        print 'O+ / Otot =', self.atom_abun['O2'][0]/Otot
        #OpOtot = float(raw_input('Enter O+/Otot: '))
        OpOtot = -0.25
        S_icf = 1.0 / 10**(OpOtot)
        S_tot = (self.atom_abun['S2'][0] + self.atom_abun['S3'][0]) * S_icf
        S_icf_perr = 0.2 # this is the percentage error of the ICF taken from the min scale in Fig 7 of Garnett 89 = 0.05
        # the measured value in the x-asis is -0.25 +- 0.05, thus 20% of the measured value
        S_toterr = numpy.sqrt( ((S_icf_perr*S_icf)**2) * (self.atom_abun['S2'][0] + self.atom_abun['S3'][0])**2 +
                               (self.atom_abun['S2'][1]**2 + self.atom_abun['S3'][1]**2) * S_icf**2 )
        elem_abun['S'] = [S_tot, S_toterr]
        S_errp = (S_toterr/S_tot)*100
        print ' S_toterr = ', S_toterr, '=', S_errp, '%'
        print 'ICF(S)=%0.2f ,   S_tot = %0.2f +- %0.2f' % (S_icf, 12+numpy.log10(S_tot), numpy.log10((100+S_errp) / (100-S_errp))/2.0)
        
        # Argon
        print ' Assuming ICF(Ar) from Peimbert, Peimbert Ruiz (2005):  (S+/S++)*Ar++ * (Ar++ + Ar+++)/(Ar++ + Ar+++)'
        Ar34 = self.atom_abun['Ar3'][0] + self.atom_abun['Ar4'][0]
        Ar_icf = (self.atom_abun['S2'][0] / self.atom_abun['S3'][0]) * self.atom_abun['Ar3'][0]
        Ar_tot = Ar34 * Ar_icf
        S23err = (self.atom_abun['S2'][0] / self.atom_abun['S3'][0]) * numpy.sqrt( (self.atom_abun['S2'][1]/self.atom_abun['S2'][0])**2 + (self.atom_abun['S3'][1]/self.atom_abun['S3'][0])**2 )
        Ar34err =  numpy.sqrt( self.atom_abun['Ar3'][1]**2 + self.atom_abun['Ar4'][1]**2 )
        Ar_toterr = numpy.sqrt( Ar_icf**2 * (((self.atom_abun['S2'][0]/self.atom_abun['S3'][0]) / S23err)**2 + (Ar34/Ar34err)**2 )**2 +
                               (self.atom_abun['Ar3'][1]**2 + self.atom_abun['Ar4'][1]**2) * Ar_icf**2 )
        elem_abun['S'] = [Ar_tot, Ar_toterr]
        Ar_errp = (Ar_toterr/Ar_tot)*100
        print ' Ar_toterr = ', Ar_toterr, '=', Ar_errp, '%'
        print 'ICF(Ar)=%0.2f ,   Ar_tot = %0.2f +- %0.2f' % (Ar_icf, 12+numpy.log10(Ar_tot), numpy.log10((100+Ar_errp) / (100-Ar_errp))/2.0)
        
        # Make sure that the temperatures and densities file closes properly
        if self.writeouts:
            outf.close()
        
    def corr_ColExcit(self, TO2gar, catalog_lines, element_lines, Idered, Hlines, verbose=False):
        '''
        This function is the second iteration of reddening correction. It fits the observed H lines given to the theoretical
        ones found with INTRAT by Storey & Hummer (1995).
        # Idered = first iteration of reddening correction
        # Hlines = list of the hydrogen wavelengths to look for (typically from Halpha to H12).
        # theoCE = theoretical hydrogen intensities corrected for collisional excitation.
        '''
        ### For collisional excitation, interpolate from Table 1 of Peimbert, Luridiana, Peimbert (2007, ApJ, 666, 636)
        # Table 1
        # Objects = NGC346, NGC2363, Haro29, SBS0335-052, IZw18
        TeOII = [12600.0, 13800.0, 14000.0, 15600.0, 15400]
        xalpha = [0.011, 0.037, 0.033, 0.086, 0.070]
        #xbeta = [0.007, 0.027, 0.021, 0.066, 0.053]
        # x_lambda = I_col/I_tot
        # now do the interpolation according to the [OII] temperature
        xL = []
        xalpha_interp = numpy.interp(TO2gar, TeOII, xalpha)
        xL.append(xalpha_interp)
        xbeta = xalpha_interp * 0.67
        xL.append(xbeta)
        L = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]
        for l in L:
            xl = xalpha_interp / ( 2**((l-2.0)/3.0) )
            xL.append(xl)
        # For the following hydrogen lines recalculate the intensity correcting for collisional exitation
        norm_IcorrCE = [] #intensities corrected for collisional excitation
        obs_ratios = []
        found_Hlines = []
        for w, el, I in zip(catalog_lines, element_lines, Idered):
            for h, l in zip(Hlines, xL):
                if (w == h) and (el == 'H'):
                    found_Hlines.append(h)
                    newI = I * (1-l)
                    normI = newI/100.
                    if verbose == True:
                        print w, 'before', I, '  after collisionally excited corrected', newI, '  ratio2Hbeta', normI
                    norm_IcorrCE.append(newI)
                    obs_ratios.append(normI)
                    if w == 4102:
                        norm_H6theo = normI
        return norm_IcorrCE, obs_ratios, found_Hlines, norm_H6theo
    
    def find_Chi_of_CE(self, TO2gar, catalog_lines, element_lines, Idered, Hlines, theoCE, percent_Iuncert, verbose=False):
        # Correct for collisional excitation
        IcorrCE, obs_ratios, found_Hlines, norm_H6theo = self.corr_ColExcit(TO2gar, catalog_lines, element_lines, Idered, Hlines)
        # Recalculate the intensities of the most prominent hydrogen lines (Halpha through H12) to match them 
        # with the theoretical ratios given by INTRAT (Storey & Hummer, 1995, MNRAS, 272, 41).
        uncert = []
        for H, Ic, obsr in zip(found_Hlines, IcorrCE, obs_ratios):
            idx = Hlines.index(H)
            if H in catalog_lines:
                H_index = catalog_lines.index(H)
                u = (1 - (obsr / theoCE[idx])) * 100.0 / percent_Iuncert[H_index]
                I = Idered[H_index]
                if verbose == True:
                    print H, 'theo_ratio =', theoCE[idx], 'obs_ratio', obsr, '   Icorr =', Ic, '   Idered = ', I#, '   percent_Iuncert[H_index]=', percent_Iuncert[H_index]
                    print '   error de excitacion colisional = ', u 
                uncert.append(u)
        # In order to find the best combination of C_Hbeta and EWabsHbeta determine chi squared
        sqs = []
        for u in uncert:
            nu = u * u
            sqs.append(nu)
        Chi_sq = sum(sqs)
        return Chi_sq, norm_H6theo
        
    def redcor2(self, I_theo_HaHb, theoCE, Hlines, TO2gar, Idered, C_Hbeta, EWabsHbeta, catalog_lines, corr_undelyingAbs_EWs, 
                rounded_catalog_wavelength, element, ion, forbidden, how_forbidden, observed_wavelength, flux, 
                intensities, EW, continuum, all_err_cont_fit, em_lines=False):
        '''This function is only used if C(Hbeta) was used for reddening correction. If E(B-v) was used instead, this
        function will be skiped.'''
        number_iterations = 14 #this number must be even
        EWabsHbeta_increase = 0.1
        C_Hbeta_increase = 0.01
        Chi_sq_models = []
        EWabsHbeta_values = []
        EWabsHbeta_values.append(EWabsHbeta)
        C_Hbeta_values = []
        C_Hbeta_values.append(C_Hbeta)
        Halpha_idx = catalog_lines.index(6563.)
        Hbeta_idx = catalog_lines.index(4861.)
        H6_idx = catalog_lines.index(4102)
        dif_TheoObs_H6Hb_values = []
        # Find the faintest detected emission line: get rid of negative fluxes
        if em_lines:
            catalog_lines, _, element_lines, _, _, _, normfluxes, _, norm_intensities, _ = find_emission_lines(rounded_catalog_wavelength, element, ion, forbidden, how_forbidden, observed_wavelength, flux, intensities, EW, continuum)
        else:
            catalog_lines, _, element_lines, _, _, _, normfluxes, _, norm_intensities, _ = normalize_lines(rounded_catalog_wavelength, element, ion, forbidden, how_forbidden, observed_wavelength, flux, intensities, EW, continuum)        
        I_obs_H6Hb = catalog_lines[H6_idx] / catalog_lines[Hbeta_idx]
        # Determine uncertainties
        percent_Iuncert, _, _ = BasicOps.get_uncertainties(catalog_lines, normfluxes, all_err_cont_fit)
        Chi_sq, I_theo_H6Hb = self.find_Chi_of_CE(TO2gar, catalog_lines, element_lines, Idered, Hlines, theoCE, percent_Iuncert)
        Chi_sq_models.append(Chi_sq)
        dif_TheoObs_H6Hb = numpy.fabs(I_theo_H6Hb - I_obs_H6Hb) 
        dif_TheoObs_H6Hb_values.append(dif_TheoObs_H6Hb)
        I_obs_HaHb = Idered[Halpha_idx] / Idered[Hbeta_idx]
        print ' ***    I_theo_HaHb =', I_theo_HaHb, '   I_obs_HaHb =', I_obs_HaHb
        diff_HaHb_values = []
        diff_HaHb = numpy.fabs(I_theo_HaHb - I_obs_HaHb)
        diff_HaHb_values.append(diff_HaHb)
        # First, variate EWabsHbeta with C_Hbeta fixed
        for EWabsHbeta_iterations in range(0, number_iterations):
            print 'EWabsHbeta_iterations', EWabsHbeta_iterations
            if I_theo_HaHb < I_obs_HaHb:
                EWabsHbeta = EWabsHbeta + EWabsHbeta_increase
            elif I_theo_HaHb > I_obs_HaHb:
                EWabsHbeta = EWabsHbeta - EWabsHbeta_increase
                if EWabsHbeta < 0.0:
                    EWabsHbeta = 0.00001
            EWabsHbeta_values.append(EWabsHbeta)
            intensities = BasicOps.underlyingAbsCorr(EWabsHbeta, corr_undelyingAbs_EWs, continuum, flux)
            # Find the faintest detected emission line: get rid of negative fluxes
            if em_lines:
                catalog_lines, _, element_lines, _, _, _, normfluxes, _, norm_intensities, _ = find_emission_lines(rounded_catalog_wavelength, element, ion, forbidden, how_forbidden, observed_wavelength, flux, intensities, EW, continuum)
            else:
                catalog_lines, _, element_lines, _, _, _, normfluxes, _, norm_intensities, _ = normalize_lines(rounded_catalog_wavelength, element, ion, forbidden, how_forbidden, observed_wavelength, flux, intensities, EW, continuum)        
            # Determine uncertainties
            percent_Iuncert, _, _ = BasicOps.get_uncertainties(catalog_lines, normfluxes, all_err_cont_fit)
            # Dered again and find the Chi_squared of that model
            cHbeta = 0.434*C_Hbeta
            Idered, _ = self.Halpha2Hbeta_dered(I_theo_HaHb, cHbeta, catalog_lines, normfluxes, norm_intensities)
            I_obs_H6Hb = catalog_lines[H6_idx] / catalog_lines[Hbeta_idx]
            Chi_sq, norm_H6theo = self.find_Chi_of_CE(TO2gar, catalog_lines, element_lines, Idered, Hlines, theoCE, percent_Iuncert)
            Chi_sq_models.append(Chi_sq)
            dif_TheoObs_HaHb = numpy.fabs(norm_H6theo - I_obs_H6Hb) 
            dif_TheoObs_H6Hb_values.append(dif_TheoObs_HaHb)
            I_obs_HaHb = Idered[Halpha_idx] / Idered[Hbeta_idx]
            print ' ***    I_theo_HaHb =', I_theo_HaHb, '   I_obs_HaHb =', I_obs_HaHb
            diff_HaHb = numpy.fabs(I_theo_HaHb - I_obs_HaHb)
            diff_HaHb_values.append(diff_HaHb)
            EWabsHbeta_iterations = EWabsHbeta_iterations + 1
        # Second, variate C_Hbeta with EWabsHbeta fixed
        EWabsHbeta = EWabsHbeta_values[0]
        for C_Hbeta_iterations in range(0, number_iterations):
            print 'C_Hbeta_iterations =', C_Hbeta_iterations
            if I_theo_HaHb < I_obs_HaHb:
                C_Hbeta = C_Hbeta + C_Hbeta_increase
            elif I_theo_HaHb > I_obs_HaHb:
                C_Hbeta = C_Hbeta - C_Hbeta_increase
                if C_Hbeta < 0.0:
                    C_Hbeta = 0.00001            
            C_Hbeta_values.append(C_Hbeta)
            cHbeta = 0.434*C_Hbeta
            intensities = BasicOps.underlyingAbsCorr(EWabsHbeta_values[0], corr_undelyingAbs_EWs, continuum, flux)
            # Find the faintest detected emission line: get rid of negative fluxes
            if em_lines:
                catalog_lines, _, element_lines, _, _, _, normfluxes, _, norm_intensities, _ = find_emission_lines(rounded_catalog_wavelength, element, ion, forbidden, how_forbidden, observed_wavelength, flux, intensities, EW, continuum)
            else:
                catalog_lines, _, element_lines, _, _, _, normfluxes, _, norm_intensities, _ = normalize_lines(rounded_catalog_wavelength, element, ion, forbidden, how_forbidden, observed_wavelength, flux, intensities, EW, continuum)        
            percent_Iuncert, _, _ = BasicOps.get_uncertainties(catalog_lines, normfluxes, all_err_cont_fit)
            # Dered again and find the Chi_squared of that model
            Idered, _ = self.Halpha2Hbeta_dered(I_theo_HaHb, cHbeta, catalog_lines, normfluxes, norm_intensities)
            I_obs_H6Hb = catalog_lines[H6_idx] / catalog_lines[Hbeta_idx]
            Chi_sq, norm_H6theo = self.find_Chi_of_CE(TO2gar, catalog_lines, element_lines, Idered, Hlines, theoCE, percent_Iuncert)
            Chi_sq_models.append(Chi_sq)
            dif_TheoObs_HaHb = numpy.fabs(norm_H6theo - I_obs_H6Hb) 
            dif_TheoObs_H6Hb_values.append(dif_TheoObs_HaHb)
            I_obs_HaHb = Idered[Halpha_idx] / Idered[Hbeta_idx]
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
        # Calculate the final dereddend values but keep in mind that model 0 is the first reddening iteration, 
        # if there were number_iterations = 10
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
        intensities = BasicOps.underlyingAbsCorr(EWabsHbeta, corr_undelyingAbs_EWs, continuum, flux)
        # Find the faintest detected emission line: get rid of negative fluxes
        if em_lines:
            catalog_lines, wavs_lines, element_lines, ion_lines, forbidden_lines, how_forbidden_lines, normfluxes, calc_cont, norm_intensities, EW_lines = find_emission_lines(rounded_catalog_wavelength, element, ion, forbidden, how_forbidden, observed_wavelength, flux, intensities, EW, continuum)
        else:
            catalog_lines, wavs_lines, element_lines, ion_lines, forbidden_lines, how_forbidden_lines, normfluxes, calc_cont, norm_intensities, EW_lines = normalize_lines(rounded_catalog_wavelength, element, ion, forbidden, how_forbidden, observed_wavelength, flux, intensities, EW, continuum)        
        lines_info = [catalog_lines, wavs_lines, element_lines, ion_lines, forbidden_lines, how_forbidden_lines, normfluxes, calc_cont, norm_intensities, EW_lines]
        # Dered again and find the Chi_squared of that model
        norm_Idered, I_dered_norCorUndAbs = self.Halpha2Hbeta_dered(I_theo_HaHb, cHbeta, catalog_lines, normfluxes, norm_intensities)
        flambdas = find_flambdas(cHbeta, I_dered_norCorUndAbs, normfluxes)
        dereddening_info = [EWabsHbeta, C_Hbeta, norm_Idered, I_dered_norCorUndAbs, flambdas]
        # Determine uncertainties    
        percent_Iuncert, absolute_Iuncert, S2N = BasicOps.get_uncertainties(catalog_lines, normfluxes, all_err_cont_fit)
        uncertainties_info = [percent_Iuncert, absolute_Iuncert, S2N]
        I_obs_HaHb = norm_Idered[Halpha_idx] / norm_Idered[Hbeta_idx]
        print ' ***    I_theo_HaHb =', I_theo_HaHb, '   I_obs_HaHb =', I_obs_HaHb
        print''
        print 'First iteration of reddening correction:   EWabsHbeta = %0.3f  C_Hbeta = %0.3f' % (EWabsHbeta_values[0], C_Hbeta_values[0])
        print '    The best combination was:              EWabsHbeta = %0.3f  C_Hbeta = %0.3f' % (EWabsHbeta, C_Hbeta)
        print '                                                     this means cHbeta = %0.3f' % (cHbeta)
        return (lines_info, dereddening_info, uncertainties_info)
        
        
    def perform_advanced_ops(self, forceTe, forceNe):
        lines_pyneb_matches = self.writeRedCorrFile()
        self.get_tempsdens()
        self.get_iontotabs(forceTe, forceNe)
        return lines_pyneb_matches


class TemperatureStruct():
    def t2_helio(self):
        pass
    def t2_RLs(self):
        pass

