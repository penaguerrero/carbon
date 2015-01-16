import numpy as np
import os
import metallicity
import string
from science import spectrum
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator
#from matplotlib.ticker import NullFormatter
#import matplotlib.gridspec as gridspec


#####################################################################################################################

'''
This script crates the following plots:
            - S/N vs difference of benchmark abundances - cloudy abundances.
  
Choose parameters to run script  
'''

# name of the object
#                  0            1           2          3         4        5         6         7         8
objects_list = ['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
                'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                  9           10         11         12         13       14        15         16         17
# Type number for specific analysis or 'all' for only intensities and S/N of the given lines 
object_number = 'all'

# List the lines of interest to find their S/N
lines2look4 = [1661, 1666, 1907, 1908, 1909, 3727, 4363, 4686, 4861, 5007, 9068]

# Specify the type of image to be saved
save_figs = False
img_format = '.jpg'


#####################################################################################################################


##### FUNCTIONS

def get_S2N(use_given_lineinfo, object_name):
    # Go into the directory of each object. Files assumed to be in: /Users/home_direcotry/Documents/AptanaStudio3/src/
    if use_given_lineinfo:
        file2use = object_name+'_RedCor_Ebv.txt'
    else:
        file2use = object_name+'_measuredLI_RedCor_Ebv.txt'
    path2txt = 'carbon/results/'+object_name+'/'+file2use
    path2txt_list = string.split(os.getcwd(), sep='carbon')
    path2txt_path = os.path.join(path2txt_list[0], path2txt)
    # read the numbers as numpy arrays
    float_data = np.loadtxt(path2txt_path, skiprows=1, usecols=(0,1,6,7,8,9,10,11,12,13,14), unpack=True)
    # float data = wave, flambda, flux, ferr, fperr, intensity, ierr, iperr, ew, ewerr, ewperr
    #now read the element, ion, forbidden, and how forbidden columns
    element, ion, forb, howforb = [], [], [], []
    tf = open(path2txt_path, "r")
    for line in tf.readlines():
        if '#' not in line:
            line_list = line.split()
            element.append(line_list[2])
            ion.append(line_list[3])
            forb.append(line_list[4])
            howforb.append(line_list[5])
    tf.close()
    str_data = [element, ion, forb, howforb]
    # Calculate S/N
    s2n = float_data[5]/float_data[7] * 100.0
    return str_data, float_data, s2n

def get_abunds(object_number, path_from_carbon_dir):
    abunds_path_list = string.split(os.getcwd(), sep='carbon')
    fname = os.path.join(abunds_path_list[0], path_from_carbon_dir)
    #print '-->  looking to read ', fname, '...'
    he, o, c, n, ne, s, to3, to2 = np.loadtxt(fname, skiprows=1, usecols=(2,3,4,5,6,7,8,9), unpack=True)
    abunds = [he[object_number], o[object_number], c[object_number], n[object_number], ne[object_number], s[object_number]]
    temps = [to3[object_number], to2[object_number]]
    # set the ratios with respect to O, theta
    abunds2H = []
    abunds2H.append(abunds[0])
    abunds2H.append(abunds[1])
    for i in range(2, len(abunds)):
        abunds2H.append(abunds[i]-abunds[1])
    abunds2H = np.array([abunds2H])
    temps = np.array([temps])
    return abunds2H, temps

def get_diff2benchmark(object_number):
    # Benchmark abundances
    Final_meas_abunds = 'carbon/results/Final_meas_abunds.txt'
    bench_abunds, bench_temps = get_abunds(object_number, Final_meas_abunds)
    # Chosen Cloudy abundances
    Cldy_best_abunds = 'carbon/results/Cloudy_best_abunds.txt'
    Cldy_abunds, Cldy_temps = get_abunds(object_number, Cldy_best_abunds)
    # Find the differences
    diffs_abunds = bench_abunds - Cldy_abunds
    diffs_temps = bench_temps - Cldy_temps
    return diffs_abunds, diffs_temps

def mk_plot(destination, save_figs, x, y, xlabel, ylabel):
    plt.figure(1, figsize=(12, 10))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.plot(x, y, 'k.')
    if save_figs:
        plt.savefig(destination)
    plt.show()

def get_list2plot(thelist, lastcol, midcol=0):
    ''' This function creates the lists to plot from either of the following lists:
         diffs_abunds_list, diffs_temps_list, lines_info_list 
         These have the form of [[[0.0, 0.0,...]]] '''
    list2plot = []
    for i in range(18):   #iterate of the number of objects in the sample
        list2plot.append(thelist[i][midcol][lastcol])
    return list2plot

def find_ticks(main_ticks_repeated, minor_ticks_repeated):
    majorLocator   = MultipleLocator(main_ticks_repeated)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator   = MultipleLocator(minor_ticks_repeated)
    return majorLocator, majorFormatter, minorLocator

def pretty_print_info(objects_list, float_data_all, s2n_all, sub_list):
    diffs_abunds_list, diffs_temps_list, lines_info_list = [], [], []
    for obj in sub_list:
        object_number = objects_list.index(obj)
        #wave, flambda, flux, ferr, fperr, intensity, ierr, iperr, ew, ewerr, ewperr = float_data_all[object_number]
        wave, _, _, _, _, intensity, _, _, _, _, _ = float_data_all[object_number]
        s2n = s2n_all[object_number]
        lines_info = []
        # get ONLY the intensity and S/N to print for all objects
        for l in lines2look4:
            if l in rounded_wavs:
                _, idx = spectrum.find_nearest(wave, l)
                linfo = [ intensity[idx], s2n[idx] ]
                lines_info.append(linfo)
        lines_info_list.append(lines_info)
        lines2print = ''
        for li in lines_info:
            lines2print = lines2print+'{:>8.2f} {:<8.2f}'.format(li[0], li[1])
        # Calculate the differences with the benchmark values
        diffs_abunds, diffs_temps = get_diff2benchmark(object_number)
        diffs_abunds_list.append(diffs_abunds)
        diffs_temps_list.append(diffs_temps)
        diffs_abundss2print = ''
        for d in diffs_abunds:
            diffs_abundss2print = diffs_abundss2print+'{:>6.2f} {:>5.2f} {:>5.2f} {:>5.2f} {:>5.2f} {:>5.2f}'.format(
                                                                                  d[0], d[1], d[2], d[3], d[4], d[5] )        
        # Print as a formated string the whole line per object
        diffs_temps2print = ''
        for dt in diffs_temps:
            diffs_temps2print = diffs_temps2print+'{:>6} {:>6}'.format(int(dt[0]), int(dt[1]))        
        print '{:<10}'.format(objects_list[object_number]), lines2print, diffs_abundss2print, diffs_temps2print
    return diffs_abunds_list, diffs_temps_list, lines_info_list


#### CODE

# Write the text file with line info?
create_txt_lineinfo = True

# Which set of line info did I use:  True = I measured the lines,   False = the code measured the lines
#                            0     1     2     3     4       5      6      7     8
use_given_lineinfo_list = [False, True, True, True, True, False, False, False, False, 
                           False, False, True, False, False, True, False, True, False]
#                            9     10     11     12    13     14     15    16    17
    
# Divide the sample into 2 lists: 
flatCO = ['mrk5', 'ngc1741', 'sbs0926', 'tol1457', 'tol9', 'iras08208']
slopeCO = ['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk960', 'pox4', 'sbs0218',
           'sbs0948', 'sbs1054', 'sbs1319', 'arp252', 'sbs1415']

if object_number == 'all':
    str_data_all, float_data_all, s2n_all = [], [], []
    for object_number in range(len(objects_list)):
        use_given_lineinfo = use_given_lineinfo_list[object_number]
        object_name = objects_list[object_number]
        str_data, float_data, s2n = get_S2N(use_given_lineinfo, object_name)
        str_data_all.append(str_data)
        float_data_all.append(float_data)
        s2n_all.append(s2n)
    rounded_wavs = metallicity.round_wavs(float_data[0])
    float_lines2look4 = []
    for l in lines2look4:
        if l in rounded_wavs:
            float_lines2look4.append(l)
    # Headers
    colhead0 = '{:<50}'.format('# ')
    colhead0 = colhead0+' {:<130} {:>40}'.format('Intensity/Hbeta  and  S/N  for the following lines', 'Differences with benchmark values')
    colhead = '{:8}'.format('#Galaxy')
    for line_of_interest in float_lines2look4:
        colhead = colhead+' {:>16}'.format(line_of_interest)
    colhead = colhead+' {:>8} {:>5} {:>5} {:>5} {:>5} {:>5} {:>7} {:>6}'.format('He', 'O', 'C', 'N', 'Ne', 'S', 'TeO3', 'TeO2')
    print colhead0
    print colhead
    # Create the lines to print 
    diffs_abunds_list, diffs_temps_list, lines_info_list = [], [], []
    print '\nFlat C/O vs O/H behavior'
    diffs_abunds_list1, diffs_temps_list1, lines_info_list1 = pretty_print_info(objects_list, float_data_all, s2n_all, flatCO)
    print '\nSlope C/O vs O/H behavior'
    diffs_abunds_list2, diffs_temps_list2, lines_info_list2 = pretty_print_info(objects_list, float_data_all, s2n_all, slopeCO)
    # combine information in the right order
    for obj in objects_list:
        if obj in flatCO:
            i = flatCO.index(obj)
            diffs_abunds_list.append(diffs_abunds_list1[i])
            diffs_temps_list.append(diffs_temps_list1[i])
            lines_info_list.append(lines_info_list1[i])
        elif obj in slopeCO:
            i = slopeCO.index(obj)
            diffs_abunds_list.append(diffs_abunds_list2[i])
            diffs_temps_list.append(diffs_temps_list2[i])
            lines_info_list.append(lines_info_list2[i])
    
    '''
    for object_number in range(len(objects_list)):
        wave, flambda, flux, ferr, fperr, intensity, ierr, iperr, ew, ewerr, ewperr = float_data_all[object_number]
        s2n = s2n_all[object_number]
        lines_info = []
        # get ONLY the intensity and S/N to print for all objects
        for l in lines2look4:
            if l in rounded_wavs:
                _, idx = spectrum.find_nearest(wave, l)
                linfo = [ intensity[idx], s2n[idx] ]
                lines_info.append(linfo)
        lines_info_list.append(lines_info)
        lines2print = ''
        for li in lines_info:
            lines2print = lines2print+'{:>8.2f} {:<8.2f}'.format(li[0], li[1])
        # Calculate the differences with the benchmark values
        diffs_abunds, diffs_temps = get_diff2benchmark(object_number)
        diffs_abunds_list.append(diffs_abunds)
        diffs_temps_list.append(diffs_temps)
        diffs_abundss2print = ''
        for d in diffs_abunds:
            diffs_abundss2print = diffs_abundss2print+'{:>6.2f} {:>5.2f} {:>5.2f} {:>5.2f} {:>5.2f} {:>5.2f}'.format(
                                                                                  d[0], d[1], d[2], d[3], d[4], d[5] )        
        # Print as a formated string the whole line per object
        diffs_temps2print = ''
        for dt in diffs_temps:
            diffs_temps2print = diffs_temps2print+'{:>6} {:>6}'.format(int(dt[0]), int(dt[1]))        
        print '{:<10}'.format(objects_list[object_number]), lines2print, diffs_abundss2print, diffs_temps2print
    '''
      
    # PLOTS
    path_list = string.split(os.getcwd(), sep='carbon')
    carbon_dir = 'carbon/cloudy_tests/pycloudy/'
    path_from_carbon_dir = carbon_dir+'S2Ntile'+img_format
    destination = os.path.join(path_list[0], path_from_carbon_dir)
    font = {#'family' : 'Vera Sans',
            'weight' : 'regular',
            'size'   : 14}
    matplotlib.rc('font', **font)
    '''
    # Note that, unlike matplotlib's subplot, the index starts from 0 in gridspec
    fig = plt.figure(1, figsize=(20, 16))
    # Set common labels
    fig.text(0.5, 0.06, 'S/N', ha='center', va='center')
    fig.text(0.035, 0.5, 'Difference with respect to benchmark values', ha='center', va='center', rotation='vertical')
    # Set conditions for each plot
    f1 = plt.subplot2grid( (2,3), (0, 0) )
    f1.set_xlim(-25.0, 120.0)
    f2 = plt.subplot2grid( (2,3), (0, 1) )
    f3 = plt.subplot2grid( (2,3), (0, 2) )
    f4 = plt.subplot2grid( (2,3), (1, 0) )
    f5 = plt.subplot2grid( (2,3), (1, 1) )
    f6 = plt.subplot2grid( (2,3), (1, 2) )
    # create a subplot that spans multiple cells
    #f7 = plt.subplot2grid( (3,3), (2, 0), colspan=3 )
    # Create the tile
    sn1909 = get_list2plot(lines_info_list, 1, 4)   # S/N of 1909
    diffc = get_list2plot(diffs_abunds_list, 2)   # Difference to C benchmark
    for i, obj in enumerate(objects_list):
        subxcoord = 6
        subycoord = -6
        side = 'left'
        f1.plot(sn1909, diffc, 'k.')
        if (obj == 'tol9') or (obj == 'iras08339') or (obj == 'arp252'):
            side = 'right'
            subxcoord = -4
        f1.annotate('{}'.format(obj), xy=(sn1909[i], diffc[i]), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
    # remove the first tick so that they do not overlap
    ax = plt.gca() 
    ax.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    plt.show()
    '''
    figs = plt.figure(1, figsize=(20, 12))
    figs.subplots_adjust(wspace=0.25, hspace=0.2)
    # Set common labels
    figs.text(0.5, 0.06, 'S/N', ha='center', va='center')
    figs.text(0.035, 0.5, 'Difference with respect to benchmark values [dex]', ha='center', va='center', rotation='vertical')
    # Set conditions for each plot
    f1 = figs.add_subplot(231)
    f1.set_title('1666 O III]')
    f1.set_xlim(-35.0, 120.0)
    f1.set_ylim(-0.6, 0.4)
    f1.set_ylabel('Difference to O')
    #f1.text(6.5, -30, 'stellar component')   
    # adding the enhancement of the stellar component  
    sn1666 = get_list2plot(lines_info_list, 1, 1)
    sn1666[3] = lines_info_list[3][1][0]
    diffo = get_list2plot(diffs_abunds_list, 1)   
    f1.plot(sn1666, diffo, 'k.')
    # remove the first tick so that they do not overlap
    #f1.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    
    f2 = figs.add_subplot(232)
    f2.set_title('4363 [O III]')
    sn4363 = get_list2plot(lines_info_list, 1, 6)
    f2.plot(sn4363, diffo, 'k.')
    f2.set_xlim(-70.0, 330.0)
    f2.set_ylim(-0.6, 0.4)
    #f2.xaxis.set_major_formatter( NullFormatter() )
    #f2.yaxis.set_major_formatter( NullFormatter() )
    f2.set_ylabel('Difference to O')
    #f2.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    
    f3 = figs.add_subplot(233)
    f3.set_title('5007 [O III]')
    sn5007 = get_list2plot(lines_info_list, 1, 9)
    f3.plot(sn5007, diffo, 'k.')
    f3.set_xlim(-1950.0, 11500.0)
    f3.set_ylim(-0.6, 0.4)
    f3.set_ylabel('Difference to O')
    #f3.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    
    f4 = figs.add_subplot(234)
    f4.set_title('1909 C III]')
    sn1909 = get_list2plot(lines_info_list, 1, 4)
    diffc = get_list2plot(diffs_abunds_list, 2)   
    f4.plot(sn1909, diffc, 'k.')
    f4.set_xlim(-35.0, 120.0)
    #f4.set_ylim(-0.75, 0.4)
    f4.set_ylabel('Difference to C')
    #f4.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    
    f5 = figs.add_subplot(235)
    f5.set_title('4686 He II')
    sn4686 = get_list2plot(lines_info_list, 1, 7)
    diffhe = get_list2plot(diffs_abunds_list, 0)   
    f5.plot(sn4686, diffhe, 'k.')
    f5.set_xlim(-83.0, 300.0)
    f5.set_ylabel('Difference to He')
    #f5.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    
    f6 = figs.add_subplot(236)
    f6.set_title('4686 He II')
    f6.plot(sn4686, diffo, 'k.')
    f6.set_xlim(-83.0, 300.0)
    f6.set_ylim(-0.6, 0.4)
    f6.set_ylabel('Difference to O')
    #f6.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    
    # Create the tile
    for i, obj in enumerate(objects_list):
        subxcoord = 5
        subycoord = -4
        side = 'left'
        # plot of S/N of 1666  VS  Difference to O benchmark
        if (obj == 'mrk1199') or (obj == 'sbs1319'):            
            subycoord = 4            
        if (obj == 'mrk5') or (obj == 'sbs1415'):
            subycoord = -10            
        if (obj == 'iras08339') or (obj == 'arp252') or (obj == 'sbs0218') or (obj == 'mrk960'):
            side = 'right'
            subxcoord = -4
        if (obj == 'iras08208'):
            side = 'right'
            subxcoord = -4
            subycoord = -10
        if (obj == 'sbs0948'):
            side = 'right'
            subxcoord = -2
            subycoord = -13            
        f1.annotate('{}'.format(obj), xy=(sn1666[i], diffo[i]), xytext=(subxcoord, subycoord), ha=side, 
                    textcoords='offset points', fontsize=13)
        # plot of S/N of 4363  VS  Difference to O benchmark
        subxcoord = 6
        subycoord = -7
        side = 'left'
        if (obj == 'mrk1199'):            
            subycoord = 4            
        if (obj == 'iras08339') or (obj == 'arp252') or (obj == 'iiizw107') or (obj == 'sbs0218') or (obj == 'iras08208'):
            side = 'right'
            subxcoord = -4
            subycoord = 0
        if (obj == 'ngc1741') or (obj == 'sbs1415') or (obj == 'pox4'):
            subycoord = -10
        if (obj == 'iras08339') or (obj == 'sbs1054') or (obj == 'tol1457') or (obj == 'sbs1319'):
            subycoord = -1
        if (obj == 'iras08208') or (obj == 'mrk960') or (obj == 'sbs0948'):
            subycoord = -5
        f2.annotate('{}'.format(obj), xy=(sn4363[i], diffo[i]), xytext=(subxcoord, subycoord), ha=side, 
                    textcoords='offset points', fontsize=13)
        # plot of S/N of 5007  VS  Difference to O benchmark
        subxcoord = 6
        subycoord = -4
        side = 'left'
        if (obj == 'sbs1415') or (obj == 'arp252') or (obj == 'iiizw107') or (obj == 'mrk960') or (obj == 'sbs0218') or (obj == 'sbs1054'):
            side = 'right'
            subxcoord = -4
        if (obj == 'ngc1741') or (obj == 'iras08208') or (obj == 'sbs0948'):
            side = 'right'
            subxcoord = -4
            subycoord = -8
        if (obj == 'sbs1054') or (obj == 'mrk1199'):
            subycoord = 0            
        f3.annotate('{}'.format(obj), xy=(sn5007[i], diffo[i]), xytext=(subxcoord, subycoord), ha=side, 
                    textcoords='offset points', fontsize=13)
        # plot of S/N of 1909  VS  Difference to C benchmark
        subxcoord = 6
        subycoord = -6
        side = 'left'
        if (obj == 'mrk1087'):            
            subycoord = -10            
        if (obj == 'iras08339') or (obj == 'iras08208'):
            side = 'right'
            subxcoord = -4
            subycoord = -4
        if (obj == 'tol1457') or (obj == 'ngc1741'):
            side = 'right'
            subxcoord = -4
            subycoord = 0
        if (obj == 'arp252'):
            side = 'right'
            subxcoord = -4
            subycoord = -8
        f4.annotate('{}'.format(obj), xy=(sn1909[i], diffc[i]), xytext=(subxcoord, subycoord), ha=side, 
                    textcoords='offset points', fontsize=13)
        # plot of S/N of 4686  VS  Difference to He benchmark
        subxcoord = -5
        subycoord = -2
        side = 'right'
        if (obj == 'sbs0948') or (obj == 'tol1457') or (obj == 'iras08208') or (obj == 'tol9'):
            side = 'left'
            subxcoord = 6
        if (obj == 'pox4'):
            side = 'left'
            subxcoord = 4
            subycoord = -10
        f5.annotate('{}'.format(obj), xy=(sn4686[i], diffhe[i]), xytext=(subxcoord, subycoord), ha=side, 
                    textcoords='offset points', fontsize=13)
        # plot of S/N of 4686  VS  Difference to O benchmark
        subxcoord = -4
        subycoord = -4
        side = 'right'
        if (obj == 'sbs0218') or (obj == 'sbs0926'):            
            subycoord = -7          
        if (obj == 'sbs1319'):            
            subycoord = 1            
        if (obj == 'sbs1415') or (obj == 'iras08208') or (obj == 'mrk5') or (obj == 'pox4'):
            side = 'left'
            subxcoord = 4
            subycoord = -8
        if (obj == 'tol1457') or (obj == 'mrk1199'):
            side = 'left'
            subxcoord = 4
            subycoord = -2
        f6.annotate('{}'.format(obj), xy=(sn4686[i], diffo[i]), xytext=(subxcoord, subycoord), ha=side, 
                    textcoords='offset points', fontsize=13)
    #ax = plt.gca() 
    #majorLocator, majorFormatter, minorLocator = find_ticks(5, 0.2)
    #ax.xaxis.set_major_locator(majorLocator)
    #ax.xaxis.set_major_formatter(majorFormatter)
    #ax.xaxis.set_minor_locator(minorLocator)
    if save_figs:
        plt.savefig(os.path.abspath(destination))
        print 'Plot ', os.path.abspath(destination), 'saved!'
    plt.show() 

else:
    use_given_lineinfo = use_given_lineinfo_list[object_number]
    object_name = objects_list[object_number]
    str_data, float_data, s2n = get_S2N(use_given_lineinfo, object_name)
    element, ion, forb, howforb = str_data
    wave, flambda, flux, ferr, fperr, intensity, ierr, iperr, ew, ewerr, ewperr = float_data
    diffs_abunds, diffs_temps = get_diff2benchmark(object_number)
    rounded_wavs = metallicity.round_wavs(wave)
    float_lines2look4 = []
    for l in lines2look4:
        if l in rounded_wavs:
            float_lines2look4.append(l)
    print '{:<11} {:>7} {:<4} {:>12} {:>10} {:>6} {:>6}'.format('Wavelength', 'Element', 'Ion', 'Forbidden', 
                                                                  'Intensity', '%error', 'S/N')
    for rw, w, el, io, f, hf, i, e, sn in zip(rounded_wavs, wave, element, ion, forb, howforb, intensity, ierr, s2n):
        #print '{:<10.2f} {:>7} {:<6} {:>4} {:>6} {:>8.2f} {:>6.2f} {:>8.2f}'.format(w, el, io, f, hf, i, e, sn)
        if rw in lines2look4:
            print '{:<10.2f} {:>7} {:<6} {:>4} {:>6} {:>8.2f} {:>6.2f} {:>8.2f}'.format(w, el, io, f, hf, i, e, sn)
    print '\nDifferences with respect to benchmark abundances for: '
    print ' {:>5} {:>5} {:>5} {:>5} {:>5} {:>5} {:>7} {:>6}'.format('He', 'O', 'C', 'N', 'Ne', 'S', 'TeO3', 'TeO2')
    for he, o, c, n, ne, s in diffs_abunds:
        print '{:>6.2f} {:>5.2f} {:>5.2f} {:>5.2f} {:>5.2f} {:>5.2f} {:>6} {:>6}'.format(he, o, c, n, ne, s, 
                                                                int(diffs_temps[0][0]), int(diffs_temps[0][1]) )

print '\n Code finished!'
