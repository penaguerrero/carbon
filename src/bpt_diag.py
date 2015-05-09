import numpy as np
import string
import os
import matplotlib.pyplot as plt

'''
This script creates a BPT plot using our line intensities AND those from Lopez-Sanchez et al.
'''
####################################################################################################################################

### SET PARAMETERS

save_plot = False
img_type = 'jpg'
show_LopSan_data = False
sampleBPTtxt = False

####################################################################################################################################

### DATA

# name of the object
#                  0            1           2          3         4        5         6         7         8
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                 9           10         11         12         13       14        15         16         17

# Were the line measurements done manually? 
#                            0     1     2     3      4       5      6      7     8
manual_measurement_list = [False, True, True, False, False, False, False, False, False, 
                           True, False, True, False, False, True, False, True, False]
#                            9     10     11     12    13     14     15    16    17

full_names_list = ['IIIZw 107', 'IRAS 08339+6517', 'MRK 1087', 'MRK 1199', 'MRK 5', 'MRK 960', 'NGC 1741', 'POX 4', 'SBS 0218+003',
               'SBS 0948+532', 'SBS 0926+606', 'SBS 1054+365', 'SBS 1319+579', 'TOL 1457-262', 'TOL 9', 'ARP 252', 'IRAS 08208+2816',
               'SBS 1415+437']

# divide the sample in objects for in which we found Te and those in which we did not
objects_with_Te = ['mrk5', 'pox4', 'sbs0948', 'sbs1054', 'tol1457', 'iras08208']
objects_Te_literature = ['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk960', 'ngc1741', 'sbs0218',
                         'sbs0926', 'sbs1319', 'tol9', 'arp252', 'sbs1415']

# File containing the SDSS outlines
bpt_contours = '../cloudy_tests/pycloudy/BPT_TEST/bpt_contours.txt'

# File with BPT line intensities from Lopez-Sanchez
lop_bpt_linefile = '../cloudy_tests/pycloudy/BPT_TEST/BPTintensitiesLop09.txt'


####################################################################################################################################

#### FUNCTIONS

def get_measured_lines(object_name, manual_measurement):
    '''This function reads the line intensity files and returns ONLY the lines of interest.'''
    linesofinterest = [4861, 5007, 6563]
    # get the benchmark measurements
    if manual_measurement:
        file2use = '_measuredLI_RedCor_Ebv.txt'   # my manual line measurements
    else:
        file2use = '_RedCor_Ebv.txt'   # my code's line measurements
    measured_lines_file = 'carbon/results/'+object_name+'/'+object_name+file2use
    measured_lines_file_path_list = string.split(os.getcwd(), sep='carbon')
    measured_lines_file_path = os.path.join(measured_lines_file_path_list[0], measured_lines_file)
    #print 'Reading: ', object_name+file2use
    # now retrieve the columns we need from the file
    meas_lineIDs = []   # list of IDs in angstroms
    meas_Isrel2Hbeta = []   # list of lines relative to Hbeta reddening corrected 
    meas_Ierr = []   # absolute error of the line intensity
    meas_Iper = []   # percentage error of the line intensity  
    meas_EW = []   # NOT BEING USED FOR THE MOMENT
    meas = open(measured_lines_file_path)
    _ = meas.readline()  # columns header
    for line in meas:
        line = line.strip()   # gets rid of \n at the end of the line
        cols = line.split()   # splits the line into a list of columns
        # the columns of interest are: ID=0, Intensity=9, Ierr=10, Iper=11, EW=12
        lineid = int(np.round(float(cols[0]), decimals=0))
        if lineid in linesofinterest:
            meas_lineIDs.append(lineid)
            meas_Isrel2Hbeta.append(float(cols[9]))
            meas_Ierr.append(float(cols[10]))
            meas_Iper.append(float(cols[11]))
            meas_EW.append(float(cols[12]))
    meas.close()
    #print 'meas_lineIDs, meas_Isrel2Hbeta', meas_lineIDs, meas_Isrel2Hbeta
    measured_lines = [meas_lineIDs, meas_Isrel2Hbeta, meas_Ierr, meas_Iper, meas_EW]
    return measured_lines

####################################################################################################################################

#### CODE

I4861, I5007, I6563 = np.array([]), np.array([]), np.array([])
I4861err, I5007err, I6563err = np.array([]), np.array([]), np.array([])

# get our STIS line intensities
for object_number, object_name in enumerate(objects_list):
    manual_measurement = manual_measurement_list[object_number]
    measured_lines = get_measured_lines(object_name, manual_measurement)
    meas_lineIDs, meas_Isrel2Hbeta, meas_Ierr, meas_Iper, meas_EW = measured_lines
    I4861 = np.append(I4861, meas_Isrel2Hbeta[0])
    I5007 = np.append(I5007, meas_Isrel2Hbeta[1])
    I6563 = np.append(I6563, meas_Isrel2Hbeta[2])
    I4861err = np.append(I4861err, meas_Ierr[0])
    I5007err = np.append(I5007err, meas_Ierr[1])
    I6563err = np.append(I6563err, meas_Ierr[2])
#print 'I5007err', I5007err
#print 'I6563err', I6563err

# gel the BPT Lopez-Sanchez line intensitites
lopdata = np.loadtxt(lop_bpt_linefile, skiprows=5, usecols=(1,2,3,4,5,6,7,8,9,10), unpack=True)
lop4861, lop4861err, lop5007, lop5007err, lop6548, lop6548err, lop6563, lop6563err, lop6584, lop6584err = lopdata

# get the SDSS outlines, columns are:
# log([NII]/Halpha), median of log([OIII]/Halpha) in each bin of log([NII]/Halpha), and 3*sigma deviation around the median, 
# where sigma is the standard deviation of log([OIII]/Halpha) values in each bin of log([NII]/Halpha) in following file:
logNHa, medOHaminus3sig, medOHa, medOHaplus3sig = np.loadtxt(bpt_contours, skiprows=1, unpack=True)

# Make text file of BPT lines for our sample
if sampleBPTtxt:
    samp_bpt_txt = '../cloudy_tests/pycloudy/BPT_TEST/sampleBPTlines.txt'
    f = open(samp_bpt_txt, 'w+')
    print >> f, '# {:<18} {:<17} {:<17} {:<18} {:<16} {:<18}'.format('Galaxy', '4861+-err', '5007+-err', 
                                                                   '6548+-err', '6563+-err', '6584+-err')
    for idx, obj in enumerate(objects_list):
        print >> f, '{:<15} {:>6.2f} {:>6.2f} {:>10.2f} {:>6.2f} {:>10.2f} {:>6.2f} {:>10.2f} {:>6.2f} {:>10.2f} {:>6.2f}'.format(
                                            obj, I4861[idx], I4861err[idx], I5007[idx], I5007err[idx], lop6548[idx], lop6548err[idx], 
                                                                            I6563[idx], I6563err[idx], lop6584[idx], lop6584err[idx] ) 
    f.close()
    print 'File  ', samp_bpt_txt, '  written!'


# PLOT

fig = plt.figure(1, figsize=(12, 10))
plt.xlim(-2.3, 0.5)
plt.ylim(-1.5, 1.5)
xlab = r'log ([NII]6584/H$\alpha$)'
ylab = r'log ([OIII]5007/H$\beta$)'
plt.xlabel(xlab)
plt.ylabel(ylab)
#  Line intensities from Lopez-Sanchez et al.
#plt.title('Data from Lopez-Sanchez et al.')
x = np.log10(lop6584/lop6563)
y = np.log10(lop5007/100.0)
if show_LopSan_data:
    plt.plot(x, y, 'ro')
# STIS data
x0 = np.log10(lop6584/I6563)
y0 = np.log10(I5007/100.0)
#plt.plot(x0, y0, 'k*')
e1 = lop6584/I6563 * np.sqrt( (lop6584err/lop6584)**2 + (I6563err/I6563)**2 )
x0err = e1 / (2.303 * (lop6584/I6563) )
e1 = I5007err/100.0
y0err = e1 / (2.303 * (I5007/100.0) )
# SDSS data
plt.plot(logNHa, medOHa, 'b')
plt.plot(logNHa, medOHaminus3sig, 'b--')
plt.plot(logNHa, medOHaplus3sig, 'b--')
# Insert name of object in the plot
i = 0
for x0i, y0i, z in zip(x0, y0, objects_list):
    if z in objects_with_Te:
        fmt='ko'
    elif z in objects_Te_literature:
        fmt='wo'
    #fmt = 'go'
    #plt.plot(x0i, y0i, fmt)   # without error bars
    plt.errorbar(x0i, y0i, xerr=x0err[i], yerr=y0err[i], fmt=fmt, ecolor='k')
    # Annotate the points 5 _points_ above and to the left of the vertex
    subxcoord = -3
    subycoord = 5
    side = 'right'
    if (z == 'sbs1054') or (z == 'sbs0948') or (z == 'tol9') or (z == 'sbs0926'):
        subxcoord = 5
        side = 'left'
    if (z == 'mrk960') or (z == 'sbs0218'):
        subxcoord = 5
        subycoord = -12
        side = 'left'
    if (z == 'tol1457') or (z == 'mrk5') or (z == 'ngc1741') or (z == 'arp252') or (z == 'iiizw107'):
        subycoord = -12
    plt.annotate('{}'.format(z), xy=(x0i,y0i), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')
    # ONLY Lopez-Sanchez data
    if show_LopSan_data:
        subxcoord = 5
        subycoord = -10
        side = 'left'
        plt.annotate('{}'.format(z), xy=(x[i],y[i]), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points', color='r')
    i = i+1

if save_plot:
    bptfig = '../cloudy_tests/pycloudy/BPT_TEST/sample_BPT.'+img_type
    fig.savefig(bptfig)
    print 'Figure ', bptfig, ' saved!'
    
plt.show()

print '\n Done!'
