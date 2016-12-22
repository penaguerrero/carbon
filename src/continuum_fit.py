import numpy
import os
from science import spectrum
import matplotlib
from matplotlib import pyplot as plt

############################################################################################################################################

''' OBJECTS OF THE SAMPLE '''

# name of the object
#                 0           1           2            3         4        5          6        7         8
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                 9           10         11         12         13       14        15         16         17

# Choose parameters to run script
# Select a number from objects_list, i = :
object_number = 1
# use all 3 files for NUV, optical, and NIR? Type which ones to use: nuv=0, opt=1, nir=2
specs = [2]
# set commented parameters, choose options 0 through 3
choose_conditions4textfiles = 3

############################################################################################################################################
'''
# write the text file with the line wavelengths, fluxes, and fitted continuum?
text_table = False
# Do you want to normalize the spectra to the continuum?
normalize = False
# Want to see the quasi-final spectrum?  (i.e. correct for redshift and rebin)
correct_redshift = True
rebin = True
### If needing to create a text file with only wavelengths and fluxes for splot change splot_text to true 
splot_text = True

# Do you want to see the plots of the fitted continuum?
plot = True
'''
divide_by_continuum = False   # Change to true if you want to divide by the continuum instead of subtract
vacuum = False               # Do you want to use Vacuum wavelengths?
nullfirst150 = True          # Do you want to consider the first and last 150 A? If yes, set to False
lastangstroms2omit = None    # How many Angstroms to omit at the end? Number or None
firstangstroms2omit = None#300.0  # How many A to omit at the beginning? Number or None
wins2omit = None#[[4567.0, 5161.0], [6338., 6840.]] # these must be z-corrected
if specs == [2]:   # make sure it is always true for the NIR spectrum
    lastangstroms2omit = 500.0

# lists to make faster the text files for splot:
#          text_table       normalize correct_redshift rebin   splot_text          plot
create_text_files = False
reg     = [create_text_files, False,     False,      False,   create_text_files,    True]
norm    = [create_text_files, True,      True,       False,   create_text_files,    True]
normreb = [create_text_files, True,      True,        True,   create_text_files,    True]
reb     = [create_text_files, False,     True,        True,   create_text_files,    True]
if choose_conditions4textfiles==0:
    conditions4textfiles = reg
elif choose_conditions4textfiles==1:
    conditions4textfiles = norm
elif choose_conditions4textfiles==2:
    conditions4textfiles = normreb
elif choose_conditions4textfiles==3:
    conditions4textfiles = reb
if conditions4textfiles:
    text_table = conditions4textfiles[0]
    normalize = conditions4textfiles[1]
    correct_redshift = conditions4textfiles[2]
    rebin = conditions4textfiles[3]
    splot_text = conditions4textfiles[4]
    plot = conditions4textfiles[5]
    
for s in specs:
    # Choose how many sigmas to clip from the continuum array
    sigmas_lsit = [[3,3,3], [2,3,1.5], [3,3,3], [3,3,3], [3,3,3], [3,3,3], [3,3,3], [3,3,3], 
                   [3,3,3],[2,3,3], [3,3,3], [2,3,3], [2,2,2], [3,3,3], [3,3,3], [3,3,3], [3,3,3], [3,3,3]]
    sigmas_away = sigmas_lsit[object_number][s]
    # in case I want to use a specific order for the polynomial, else it will be determined by the algorithm
    order_list = [[9,5,9], [8,3,5], [2,7,3], [7,2,3], [3,3,1], [5,5,5], 
                  [7,3,3], [7,8,1], [3,2,9], [7,2,1], [7,1,1], [7,1,1], 
                  [5,1,1], [3,1,1], [9,3,1], [9,9,7], [7,7,3], [7,5,1]]
    order = order_list[object_number][s]
    # What is the width of the window to use to find local continuum?
    window_list = [[150,90,80], [70,400,800], [100,80,500], [80,130,280], [150,250,100], [100,300,200], 
                   [70,350,350], [80,60,550], [100,250,650], [80,95,150], [50,350,350], [80,250,300], 
                   [70,550,550], [300,300,300], [50,250,250], [100,380,200], [50,250,500], [40,550,550]]
    window = window_list[object_number][s]

object_name = objects_list[object_number]

# Path of the text files of wavelengths and fluxes for the objects. 
# This works as long as the folder structure is always the same. This file is assumed to be in
#            /Users/home_direcotry/Documents/AptanaStudio3/src/
# so I want to go one back to find the results folder
results_path = "../results/"
# but just to make sure we are in the right place, get and work with full paths
full_results_path = os.path.abspath(results_path)
text_files_path = os.path.join(full_results_path, "1Dspecs/")
results4object_path = os.path.join(full_results_path, object_name)

add_str = "_selectedspecs"
data, full_file_list = spectrum.loadtxt_from_files(object_name, add_str, specs, text_files_path)
# alternative for when files have been corrected in splot
altern_file = False
### To get altern files run: 1. correct_spec script, 2.rspectext, 3.splot, 4.correct with j and save with i, 5.wspectext
##altern = '../results/sbs1319/sbs1319_opt_corr1.txt'
##f, w = numpy.loadtxt(altern, skiprows=5, usecols=(1,2), unpack=True)  ## OLD FILE
for s in specs:
    if (object_name=='iras08339') and (s==2):
        altern = '../results/iras08339/iras08339_nir_corr.txt'
        altern_file = True
    if (object_name=='mrk1087') and (s==1):
        altern = '../results/mrk1087/mrk1087_opt_corr.txt'
        altern_file = True
    if (object_name=='mrk1087') and (s==2):
        altern = '../results/mrk1087/mrk1087_nir_corr.txt'
        altern_file = True
    if (object_name=='mrk1199') and (s==1):
        altern = '../results/mrk1199/mrk1199_optspec_corr.txt'
        altern_file = True
    if (object_name=='mrk5') and (s==1):
        altern = '../results/mrk5/mrk5_optspec.txt'
        altern_file = True
    #if (object_name=='mrk960') and (s==1):
    #    altern = '../results/mrk960/mrk960_optspec_corr.txt'
    #    altern_file = True
    if (object_name=='sbs0218') and (s==1):
        #altern = '../results/sbs0218/sbs0218_opt_corr.txt'
        altern = '../results/sbs0218/sbs0218_optspec_corr.txt'
        altern_file = True
    if (object_name=='sbs0948') and (s==1):
        altern = '../results/sbs0948/sbs0948_optspec_corr.txt'
        altern_file = True
    if (object_name=='sbs0926') and (s==1):
        altern = '../results/sbs0926/sbs0926_optspec_corr.txt'
        altern_file = True
    if (object_name=='sbs0926') and (s==2):
        altern = '../results/sbs0926/sbs0926_nirspec_corr.txt'
        altern_file = True
    if (object_name=='pox4') and (s==1):
        altern = '../results/pox4/pox4_opt_corr.txt'
        altern_file = True
    if (object_name=='sbs1054') and (s==1):
        altern = '../results/sbs1054/sbs1054_optspec_corr.txt'
        altern_file = True
    if (object_name=='sbs1319') and (s==1):
        altern = '../results/sbs1319/sbs1319_optspec_corr.txt'
        altern_file = True
    if (object_name=='tol9') and (s==1):
        altern = '../results/tol9/tol9_opt21_fix.txt'
        altern_file = True
    if (object_name=='arp252') and (s==1):
        altern = '../results/arp252/arp252_optspec_corr.txt'
        altern_file = True
    if (object_name=='sbs1415') and (s==1):
        altern = '../results/sbs1415/sbs1415_optspec_cor.txt'
        altern_file = True
if altern_file:
    w, f = numpy.loadtxt(altern, unpack=True)
    data = [numpy.array([w,f])]

# Terminations used for the lines text files
spectrum_region = ["_nuv", "_opt", "_nir"]

z_list = [0.01985, 0.019581, 0.02813, 0.01354, 0.002695, 0.021371, 0.013631, 0.01201, 0.05842,
          0.046240, 0.013642, 0.002010, 0.0073, 0.01763, 0.01199, 0.032989, 0.04678, 0.002031]
# original phase 2 z for iras08339 0.019113, mrk960 0.021371, ngc1741 0.013631, tol1457 0.01763, tol9 0.010641
if wins2omit is not None:
    zshift_wins2omit = []
    if choose_conditions4textfiles == 0:
        for item in wins2omit:
            w0 = item[0] * (1.0+z_list[object_number])
            w1 = item[1] * (1.0+z_list[object_number])
            w_pair = [w0, w1]
            zshift_wins2omit.append(w_pair)
    else:
        zshift_wins2omit = wins2omit
else:
    zshift_wins2omit = None
        
'''The STIS data handbook gives a dispersion of 0.60, 1.58, 2.73, and 4.92 Angstroms per pixel for grating settings G140L, 
G230L, G430L, and G750L, respectively. The resolution element is ~1.5 pixels. '''
originals = [1.58, 2.73, 4.92]
or1 = originals[0]
or2 = originals[1]
or3 = originals[2]
#                                    0                          1                  2                     3                4                5
desired_disp_listoflists = [[or1*2.0, or2*2.0, or3*2.0], [1.6, 6.5, 5.0], [or1*1.2, or2*1.51, or3*1.2], [2.0, 5.0, 5.0], [2.0, 3.0, 5.0], [1.7, 5.6, 9.8], 
                            #        6                          7                 8                9                     10               11
                            [or1*2.0, or2*2.0, or3*2.0], [1.8, 4.0, 6.0], [1.6, 3.0, 5.0], [or1*2., or2*2., or3*2.], [1.6, or2*2., 5.6], [1.6, 3.9, 5.0],
#                                   12             13                14               15               16               17                            
                            [1.6, 5.0, 5.0], [1.6, 4.5, 5.0], [1.6, 3.0, 5.0], [1.7, 5.6, 9.8], [1.6, 3.1, 5.1], [1.6, 8.3, 11.0]]
desired_disp_list = desired_disp_listoflists[object_number]


for d, s in zip(data, specs):
    print 'Working with:  %s' % full_file_list[s]
    
    if rebin == True:
        # This mode is just to show how much the spectrum will be rebinned
        d = spectrum.rebin2AperPix(originals[s], desired_disp_list[s], d)
        text_table = False

    # Correct spectra for redshift and calculate continuum with a polynomial of nth order
    if correct_redshift == True:
        z = z_list[object_number]
        #text_table = False
        spectrum_region[s] = spectrum_region[s]+'_zcorr'
    else:
        z = 0.0  # this is non-important because we are NOT correcting for redshift
      
    object_spectra, fitted_continuum, err_cont_fit = spectrum.fit_continuum(object_name, d, z, order=order, sigmas_away=sigmas_away, 
                                                                            window=window, plot=plot, z_correct=correct_redshift, 
                                                                            normalize=normalize, nullfirst150=nullfirst150,
                                                                            divide_by_continuum=divide_by_continuum,
                                                                            lastangstroms2omit=lastangstroms2omit,
                                                                            firstangstroms2omit=firstangstroms2omit, 
                                                                            wins2omit=zshift_wins2omit)
    wavs, fluxs = object_spectra
    _, cont_fluxs = fitted_continuum
    '''
    # Create the sub-plots
    def mk_subplt(wavs, fluxs, xlims, lines, legends):
        font = {#'family' : 'normal',
                'weight' : 'normal',
                'size'   : 16}
        matplotlib.rc('font', **font)
        fig = plt.figure(1, figsize=(12, 10))
        xlab = 'Wavelength [$\AA$]'
        ylab = 'Flux [ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]'
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.minorticks_on()
        yarr = fluxs[(wavs >= xlims[0]) & (wavs <= xlims[1])]
        plt.xlim(xlims[0], xlims[1])
        plt.ylim(min(yarr)+min(yarr)*0.25, max(yarr)+max(yarr)*0.35)
        plt.plot(wavs, fluxs, 'k')#, wavs, cont_fluxs, 'r--')
        # Annotate the points 5 _points_ above and to the left of the vertex
        subxcoord, subycoord, side = 0, 85, 'center'
        for line, leg in zip(lines, legends):
            x, idx = spectrum.find_nearest(wavs, line)
            y = fluxs[idx]
            #plt.axvline(x=x, ymin=0.25, ymax=0.3, color='r')
            vlinex, vliney = [x, x], [y+y*0.10, y+y*0.10+0.3e-15]
            plt.plot(vlinex, vliney, 'r')
            plt.annotate('{} {}'.format(leg, line), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, 
                         textcoords='offset points')
        plt.show()    
    if s == 0:
        xlims = [1500.0, 2100.0]
        lines = [1661, 1666, 1909]
        legends = ['O III]', 'O III]', 'C III]']
        mk_subplt(wavs, fluxs, xlims, lines, legends)
    if s == 1:
        xlims = [3700.0, 5100.0]
        lines = [3727, 3869, 3967, 4069, 4076, 4102, 4340, 4363, 4861, 4959, 5007]
        legends = ['[O II]', '[Ne III]', '[Ne III]', '[S II]', '[S II]', 'Hg', 'Hd', 
                   '[O III]', 'Hb', '[O III]', '[O III]']
        mk_subplt(wavs, fluxs, xlims, lines, legends)
    if s == 2:
        xlims = [6250.0, 6750.0]
        lines = [6300, 6312, 6563, 6678, 6716]
        legends = ['[O I]', '[S III]', 'Ha', 'He I', '[S II]']
        mk_subplt(wavs, fluxs, xlims, lines, legends)
        xlims = [9000.0, 9600.0]
        lines = [9069, 9531]
        legends = ['[S III]', '[S III]']
        mk_subplt(wavs, fluxs, xlims, lines, legends)
        '''

    if text_table == True:
        # Write text file of wavelengths, fluxes, and fitted continuum
        if normalize == False:
            new_file_name = object_name+spectrum_region[s]+".txt"
        else:
            new_file_name = object_name+spectrum_region[s]+"_normalized.txt"
        file_name = os.path.join(results4object_path, new_file_name)
        txt = open(file_name, 'w+')
        print >> txt, 'Percentage Error in continuum fitting  =', err_cont_fit
        print >> txt, '{:<10} {:>30} {:>35}'.format('Wavelength [A]', 'Flux [ergs/s/cm$^2$/$\AA$]', 'Continuum [ergs/s/cm$^2$/$\AA$]')
        for w, f, c in zip(wavs, fluxs, cont_fluxs):
            print >> txt, '{:<7.3f} {:>25.5e} {:>32.5e}'.format(w, f, c)
        txt.close()
        print 'I recorded the wavelengths, fluxes, and error in the file: ', txt

    ### The following lines are just in case of wanting to create a text file to convert into readable splot fits
    if splot_text:
        if s == 0:
            part_of_spec = 'nuv'
        elif s == 1:
            part_of_spec = 'opt'
        elif s == 2:
            part_of_spec = 'nir'
        name_out_file = os.path.join(results4object_path, object_name+"_"+part_of_spec+"spec.txt")
        if normalize:
            name_out_file = os.path.join(results4object_path, object_name+"_"+part_of_spec+"spec_NORMALIZED.txt")
        if rebin:
            if normalize:
                name_out_file = os.path.join(results4object_path, object_name+"_"+part_of_spec+"spec_NORMandREBIN.txt")
            else:
                name_out_file = os.path.join(results4object_path, object_name+"_"+part_of_spec+"spec_REBINNED.txt")
        fout = open(name_out_file, 'w+')
        wavs, fluxs = object_spectra
        for w, f in zip(wavs, fluxs):
            fout.write('{:<4.5f} {:>20.10e}\n'.format(w, f))
        fout.close()

    
c = 3.0e5 #km/s
z = z_list[object_number]
velocity = c * z
print 'v = c * z = 3e5 * %0.5f = %0.3f' % (z, velocity)
H0 = 75#67.8 #+-0.77 km/sec/Mpc (from Planck mission)
distance = velocity / H0
print 'Estimated distance assuming H0 = %f:   d[Mpc] = %0.3f' % (H0, distance)
    
    
print 'Code finished!'

