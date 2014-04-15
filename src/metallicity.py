import os
import pyfits
import numpy
import string
import copy
import pyneb as pn
import PIL.Image as Image
from uncertainties import unumpy
from science import spectrum
from matplotlib import pyplot
 

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

def find_flambdas(cHbeta, catalog_wavelength, I_dered_norCorUndAbs, normfluxes):
    # Finding the f_lambda values
    all_flambdas = []
    for Icor, Iobs, w in zip(I_dered_norCorUndAbs, normfluxes, catalog_wavelength):
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
        I_dered_norCorUndAbs = []
        # Obtain the reddening corrected intensities based on the given law and c(Hb)
        if ebv != 0.0:
            rv = av / ebv
            # Define a reddening correction object
            RC = pn.RedCorr(E_BV=ebv, R_V=rv, law=redlaw, cHbeta=cHbeta)
        else:
            # Define a reddening correction object
            RC = pn.RedCorr(law=redlaw, cHbeta=cHbeta)
        for w, nF, nI in zip(catalog_lines, normfluxes, norm_intensities):   
            I_dered = nI * RC.getCorrHb(w)
            print w, nI, RC.getCorrHb(w) , I_dered                                            # ********
            Idered.append(I_dered)
            # Obtain the reddening corrected intensities WITHOUT the correction due to 
            # underlying absorption based on the given law and c(Hb)
            IdnUA = nF * RC.getCorrHb(w)
            I_dered_norCorUndAbs.append(IdnUA)
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
        self.I_dered_norCorUndAbs = I_dered_norCorUndAbs
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

    
class CollisionalExcitationCorr:
    def __init__(self, object_name, cHbeta, case, verbose=False):
        # Inputs:
        self.object_name = object_name
        self.cHbeta = cHbeta
        self.case = case                        # this is the Case used through out the class
        self.verbose = verbose                  # if True print midpoints in order to know what is the code working on
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

    def get_temps(self):
        # Determine first approximation of temperatures and densities
        # Define an Observation object and assign it to name 'obs'
        obs = pn.Observation()
        # read data from file created specifically for pyneb reading
        obs.readData(self.pynebIDstxt, fileFormat='lines_in_rows', corrected=True, errIsRelative=False)
        # Intensities
        for line in obs.lines:
            if self.verbose == True:            
                print 'line.wave', line.wave, '     line.corrIntens', line.corrIntens
            if line.wave == 4363:
                I1 = line.corrIntens
                print '4363 has an intensity of', I1[0]
            elif line.wave == 5007:
                I2 = line.corrIntens
                print '5007 has an intensity of', I2[0]
            elif line.wave == 3726:
                IO21 = line.corrIntens
                print '3726 has an intensity of', IO21[0]
            elif line.wave == 3729:
                IO22 = line.corrIntens
                print '3729 has an intensity of', IO22[0]
            elif line.wave == 6312:
                IS31 = line.corrIntens
                print '6312 has an intensity of', IS31[0]
            elif line.wave == 9069:
                IS32 = line.corrIntens
                print '9069 has an intensity of', IS32[0]
            elif line.wave == 9531:
                IS33 = line.corrIntens
                print '9531 has an intensity of', IS33[0]
            elif line.wave == 5518:
                ICl31 = line.corrIntens
                print '5518 has an intensity of', ICl31[0]
            elif line.wave == 5538:
                ICl32 = line.corrIntens
                print '5538 has an intensity of', ICl32[0]
        # simultaneously compute temperature and density from pairs of line ratios
        # First of all, a Diagnostics object must be created and initialized with the relevant diagnostics.
        diags = pn.Diagnostics()   # Instantiate the Diagnostics class
        diags.getAllDiags()  # see what Diagnostics exist
        # temperature determination from an intensity ratio
        # explore some specific atom in the atoms collection
        self.TO3 = 0.0
        try:
            O3 = pn.Atom("O", "3")
            O3ratio = I1[0] / I2[0]
            print 'ratio of O3 = ', O3ratio
            self.TO3 = O3.getTemDen(O3ratio, den=100., wave1=4363, wave2=5007)
            print '  First estimation of temperature of O3 = ', self.TO3 
            O2 = pn.Atom("O", "2")
            O2ratio = IO22[0] / IO21[0]
            print 'ratio of O2 = ', O2ratio
            denO2 = O2.getTemDen(O2ratio, tem=self.TO3, wave1=3729, wave2=3726) 
            print '   First estimation of density of O2 = ', denO2
            S3 = pn.Atom("S", "3")
            S3ratio = IS31[0] / (IS32[0]) 
            print 'ratio of S3 = ', S3ratio
            TS3 = S3.getTemDen(S3ratio, den=100., wave1=6312, wave2=9532)
            print '   First estimation of temperature of S3 = ', TS3 
            Cl3 = pn.Atom("Cl", "3")
            Cl3ratio = ICl32[0] / (ICl31[0]) 
            dCl3 = Cl3.getTemDen(S3ratio, temp=self.TO3, wave1=5538, wave2=5518)
            print '   First estimation of density of Cl3 = ', dCl3
        except Exception as e:
            (NameError,),e
        ### Density measurement from [Fe III] lines -- taken from Peimbert, Pena-Guerrero, Peimbert (2012, ApJ, 753, 39)
        if self.TO3 != 0.0:
            I4986 = 0.0
            I4987 = 0.0
            I4658 = 0.0
            for w, i in zip(self.wavelength, self.intensity):
                if int(w) == 4986:
                    I4986 = i
                elif int(w) == 4987:
                    I4987 = i
                elif int(w) == 4658:
                    I4658 = i
            if (I4986 != 0.0) and (I4987 != 0) and (I4658 !=0):
                log_Fe3den = 2 - ( (numpy.log10((I4986+I4987)/I4658) - 0.05 - 0.25*(numpy.log10(self.TO3-4))) / (0.66 - 0.18*(numpy.log10(self.TO3-4))) )
                Fe3den = 10**(log_Fe3den)
                print 'Density measured from [Fe3] lines:', Fe3den
            else:
                print 'No [Fe3] density available.'

        # Simultaneously determine temps and densities
        try:
            tem_N2, den_tmp = diags.getCrossTemDen('[NII] 5755/6548', '[SII] 6731/6716', obs=obs)
            tem_Ar3, den_S2 = diags.getCrossTemDen('[ArIII] 5192/7300+', '[SII] 6731/6716', obs=obs)
            tem_O3, den_Cl3 = diags.getCrossTemDen('[OIII] 4363/5007', '[ClIII] 5538/5518', obs=obs)
            tem_O3, den_Ar4 = diags.getCrossTemDen('[OIII] 4363/5007', '[ArIV] 4740/4711', obs=obs)
            tem_O3, den_O2 = diags.getCrossTemDen('[OIII] 4363/5007', '[OII] 3926/3929', obs=obs)
            # Printout of physical conditions
            print 'den_O2: ', den_O2
            print 'den_S2: ', den_S2
            print 'tem_O3: ', tem_O3
            print 'den_Ar4: ', den_Ar4
        except Exception as e:
            (NameError,),e
        '''
        # Define all atoms to make calculations
        all_atoms = pn.getAtomDict()
        # Alternate way of computing T(OIII)
        if tem_O3 == 'NA':
            tem_O3 = all_atoms['O3'].getTemDen(i5007/i4363, den=100., wave1=5007, wave2=4363)
        # Include in diags the relevant line ratios
        diags.addDiag([
                        ## temperatures
                        #'[NII] 5755/6548',
                        '[OII] 7320/3737+',
                        '[OIII] 4363/5007',
                        '[ArIII] 5192/7136',
                        '[ArIII] 5192/7300+',
                        '[ArIV] 7230+/4720+',
                        #'[SIII] 6312/9532',
                        ## densities
                        #'[NI] 5198/5200',
                        #'[OII] 3729/3736',
                        '[ArIV] 4740/4711',
                        #'[SII] 4072+/6725+',
                        '[SII] 6731/6716',
                        ])
        '''
        ### Second iteration of extinction correction: Collisional Excitation
        ### Following analysis presented in Pena-Guerrero, Peimbert, Peimbert, Ruiz (2012, ApJ, 746, 115) & Peimbert, Pena-Guerrero, Peimbert (2012, ApJ, 753, 39).
        # If the [O II] temperature was not obtained directly from observations, get an estimate
        ### Get temperature of OII from OIII. 
        # Using equation of Peimbert, Peimbert, & Luridiana (2002, ApJ, 565, 668) - Based on data of Stasinska's models.
        self.TO2pei = 2430. + self.TO3 * (1.031 - self.TO3/54350.)
        print 'This is the theoretically obtained temperature of O2 from Peimbert etal 2002 = ', self.TO2pei
        print ' * for comparison, Temperature of [O III] = ', self.TO3
        # Using equation of Garnett, D. R. 1992, AJ, 103, 1330
        self.TO2gar = 0.7 * self.TO3 + 3000.
        print 'Theoretically obtained temperature of O2 from Garnet 1992 = ', self.TO2gar
        
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
        
    def perform_colexcit_corr(self):
        lines_pyneb_matches = self.writeRedCorrFile()
        self.get_temps()
        return lines_pyneb_matches

class UseTemdenAbund():
    def do_temden(self):
        pass
    def do_abund(self):
        pass

class TemperatureStruct():
    def t2_helio(self):
        pass
    def t2_RLs(self):
        pass

