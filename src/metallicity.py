import os
import pyfits
import numpy
import string
import copy
import PIL.Image as Image
from science import spectrum
from matplotlib import pyplot


'''
This program contains various classes that together determine the object's metallicity.
'''

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
        whole_spec_arr = numpy.array([whole_spec_wavs, whole_spec_flxs])
        # rebin the arrays to make a prettier spectra
        print 'Type rebinning factor:  '
        desired_rebin_factor = raw_input('(i.e. value<1 will reduce the resolution)  ')
        rebin_factor = (1, float(desired_rebin_factor))
        whole_spec_arr = spectrum.rebin(whole_spec_arr, rebin_factor)
        pyplot.figure(1, figsize=(15, 8))
        save_plt = raw_input('Do you want to save the smoothened whole spectrum image?  (n)  y   ')
        pyplot.xlabel('Wavelength  [$\AA$]')
        pyplot.ylabel('Flux  [ergs/s/cm$^2$/$\AA$]')
        fmax = max(whole_spec_arr[1]) - max(whole_spec_arr[1])*0.2 
        pyplot.text(8500, fmax, full_name, fontsize=18)
        pyplot.plot(whole_spec_arr[0], whole_spec_arr[1], 'k')
        if save_plt == 'y':
            destination = plot_name+'_wholespec'+img_format
            pyplot.savefig(destination)
            print('Plot %s was saved!' % destination)
        elif save_plt == 'n':
            print('Plot not saved.')
        pyplot.show()

        
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


class BasicOps(object):
    '''
    This class gathers all the basic operations after we have the 1D spectra. These operations are:
    - reddening correction
    - finding local continuum
    - line intensity measurement, equivalent widths, and FWHM
    '''
    def __init__(self, wf_used_specs):
        self.wf_used_specs = wf_used_specs
        self.n = 0.999721
        self.w_air = []
        self.w_vac = []
        self.w_theo = []
        self.w_obs = []
        self.w_findz = [3727, 4861.33, 4958.91, 5006.84]
        #self.wf_usedspecs = []  # list of wavelengths and fluxes arrays of used spectra
        self.z = 0.0
        self.convert_wair2vacuum()
        self.wavs_flxs_usedspecs(wf_used_specs)
        self.find_z()
        self.correct_for_z()
        self.find_lines()
        self.determine_line_info()
           
    '''
    def convert_wair2vacuum(self, w_air):
        print('By defaulf the value of the refraction index is n=0.999721') 
        print('    n was taken from NIST, using 1-n=lambda(vac)/lambda(air)-1')
        print('Do you want to use another value of n? If yes type it, else enter.')
        q = raw_input()
        if q != 0.0:
            n = float(q)
        else:
            n = self.n
        for wa in w_air:
            if wa > 2900.:
                wv = wa * (2 - n)
        return(wv)
    '''
    def wavs_flxs_usedspecs(self, wf_used_specs):
        '''
        This function does 2 things:
        1) punts in ONE 2D array the arrays of wvelength and flux for each spectrum
        2) creates the list that contains all those 2D arrays for the 3 used spectra (i.e. 3 pairs of wavs and flxs),
            this is the list that is returned.
        '''
        #for pair in wf_used_specs:
        #    w = pair[0]
        #    f = pair[1]
            
            
    def find_z(self, wf_used_specs):
        '''
        This function corrects the observed wavelength for redshift. It uses the result from function wavs_flxs_usedspecs.
        '''
        print('List wavelengths to use for reddening correction.')
        print('To use default (3727, 4861.33, 4958.91, 5006.84) press enter')
        list_findz = raw_input()
        w_findz = []
        if list_findz != 0:
            for item in list_findz:
                w = float(item)
                w_findz.append(w)
        else:
            w_findz = self.w_findz
        # Find these wavelengths in the spectrum
        # Plot the optical part
        fig = pyplot.figure(1, figsize=(10, 10))
        for wf in wf_used_specs:
            # wf is the pair of wavelength and flux arrays for each one of the selected spectra
            print 'Indicate where are the followig lines:'
            print(w_findz)
            x = len(w_findz)
            w_findz_in_obs = fig.ginput(x, timeout=-1)
            fig.plot(wf[0], wf[1])
            pyplot.draw()
            print(type(w_findz_in_obs))
            
        
        '''
        wz = [] # this is the list of the wavelengths in the spec that are closest to 
        for w in w_findz:
                
            nearest_w = spectrum.find_nearest(, w)
                
            
        
        
        zs = []
        wf_OBS = self.wavs_flxs_usedspecs
        for wobs in w_fidzOBS:
            # Now find z
            for wtheo in w_findz:
                z = (wobs/wtheo) - 1.0
                zs.append(z)
        self.z = sum(zs)/len(zs)
        return(self.z)
        '''
    def correct_for_z(self):
        pass
    
    '''
    def reddening_correction(self, w_air):
        for wobs in w_HST:
            if wobs < 2900.0:
                wobs = self.convert_wair2vacuum(wobs)
                self.w_vac.append(wobs)
            else:
                pass
        for wtheo in self.w_vac:
            self.z = (wobs/wtheo) - 1
        return(self.z)
    '''
    
    def find_lines(self):
        pass
    def determine_line_info(self):
        '''
        This method determines the line intensity, its EW, and its FWHM.
        '''
        pass
    
class Corrections2optical:
    def underlying_abs_corr(self):
        pass
    def Ierr_determination(self):
        '''this method determines the error in the line intensities. '''
        pass
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

