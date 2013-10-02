import os
import pyfits
import numpy
import string
#import spectrum
from matplotlib import pyplot


'''
This program contains various classes that together determine the object's metallicity.
'''

class OneDspecs:
    '''
    This uses the 1d fits extraction files from all three filters (g230l, g430l, and g750l), 
    in order to let me choose which one has the best S/N. 
    '''
    def __init__(self, oneDspecs, selectedspecs_file, plot_name):
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
        self.plot_chose_from(plot_name)
        self.selecting_specs(selectedspecs_file, plot_name)
        
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

    def do_plot(self, plot_name, specs_list, x_axs, y_axs, xlolim, xuplim, ylolim, yuplim, save=True):
        fig1 = pyplot.figure(1, figsize=(12, 8))
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
            fig1.savefig(plot_name)
        pyplot.show(fig1)

    def save_plot(self, plot_name, specs_list, x_axs, y_axs, xlolim, xuplim, ylolim, yuplim):
        pyplot.ioff()
        save_plt = raw_input('Do you want to save this plot? (n)  y  ')
        if save_plt == 'n' or save_plt == '':
            os.remove(plot_name)
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
                os.remove(plot_name)
                self.do_plot(plot_name, specs_list, x_axs, y_axs, xlolim, xuplim, ylolim, yuplim)
            else:
                print('x-lolim: %i,  x-uplim: %i, y-lolim: %e,  y-uplim: %e' % (xlolim, xuplim, ylolim, yuplim))
            #plot_name = raw_input('Enter plot name (including directory)')
            print('Plot %s was saved!' % plot_name)
      
    def plot_chose_from(self, plot_name):
        xlolim = self.xlolim
        xuplim = self.xuplim
        ylolim = self.ylolim
        yuplim = self.yuplim
        keep_prevFig = raw_input('Is there is a previous figure you want to keep?  (n)  y   ')
        if keep_prevFig == 'y':
            self.do_plot(plot_name, self.specs, self.wavs, self.flxs, xlolim, xuplim, ylolim, yuplim, save=False)
        else:
            self.do_plot(plot_name, self.specs, self.wavs, self.flxs, xlolim, xuplim, ylolim, yuplim, save=True)
            pyplot.ioff()
            self.save_plot(plot_name, self.specs, self.wavs, self.flxs, xlolim, xuplim, ylolim, yuplim)
    
    def selecting_specs(self, selectedspecs_file, plot_name):
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
        w_nuv = self.wavs[self.indx_selected_wavs_and_flxs[0]]
        w_opt = self.wavs[self.indx_selected_wavs_and_flxs[1]]
        w_nir = self.wavs[self.indx_selected_wavs_and_flxs[2]]
        f_nuv = self.flxs[self.indx_selected_wavs_and_flxs[0]]
        f_opt = self.flxs[self.indx_selected_wavs_and_flxs[1]]
        f_nir = self.flxs[self.indx_selected_wavs_and_flxs[2]]
        
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
        plot_selecSpecs = plot_name.replace(".eps", "_selecSpectra.eps")
        keep_prevFig = raw_input('Is there is a previous figure you want to keep?  (n)  y   ')
        if keep_prevFig == 'y':
            self.do_plot(plot_selecSpecs, used_specs, wf_arr[:,0], wf_arr[:,1], xlolim, xuplim, ylolim, yuplim, save=False)
        else:
            self.do_plot(plot_selecSpecs, used_specs, wf_arr[:,0], wf_arr[:,1], xlolim, xuplim, ylolim, yuplim, save=True)
            self.save_plot(plot_selecSpecs, used_specs, wf_arr[:,0], wf_arr[:,1], xlolim, xuplim, ylolim, yuplim)        
        return(self.wf_used_specs)
        
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

        plot_selecSpecs = plot_name.replace(".eps", "_selecSpectra.eps")
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

