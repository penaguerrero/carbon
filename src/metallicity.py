import os
import pyfits
import numpy
import string
import copy
import pyneb as pn
import PIL.Image as Image
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
    for i in range(len(flux)):
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

def find_flambdas(cHbeta, I_dered_norCorUndAbs, normfluxes):
    # Finding the f_lambda values
    flambdas = []
    for Icor, Iobs in zip(I_dered_norCorUndAbs, normfluxes):
        f12 = (numpy.log10(Icor) - numpy.log10(Iobs)) / cHbeta
        flambdas.append(f12)
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
            print('This will be spectruum number %i' % self.counter)            
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
    def __init__(self, law, cols_in_file, I_theo_HaHb, EWabsHbeta, cHbeta, av, z, ebv=None):
        self.law = law
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
        self.underlyingAbsCorr()
        self.Halpha2Hbeta_dered(av, ebv)


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
        ''' Function to dered and obtain the Halpha/Hbeta ratio nicely printed along with the unreddend values to compare. '''
        law = self.law
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
        for w, nF, nI in zip(catalog_lines, normfluxes, norm_intensities):    
            # Obtain the reddening corrected intensities based on the given law and c(Hb)
            if ebv != None:
                rv = av / ebv
                RC = pn.RedCorr(E_BV=ebv, R_V=rv, law=law, cHbeta=cHbeta)
            else:
                RC = pn.RedCorr(law=law, cHbeta=cHbeta)
            I_dered = nI * RC.getCorrHb(w)
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
        print '            Using', law, '                   Normalized fluxes before extinction correction'
        print catalog_lines[Halpha_idx], '    ', I_Halpha, '                  ', Halpha
        print catalog_lines[Hbeta_idx], '    ', I_Hbeta, '                          ', Hbeta
        print 'theoretical ratio Ha/Hb = %0.3f' % (I_theo_HaHb)
        print '      observed Ha/Hb = %0.3f           raw Ha/Hb = %0.3f' % (I_obs_HaHb, raw_ratio)
        return normfluxes, Idered, I_dered_norCorUndAbs


    
class CollisionalExcitationCorr():
    def collisional_excit_corr(self):
        pass

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

