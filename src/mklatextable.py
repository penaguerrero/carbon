import numpy as np
import string
import os


'''
This script creates the LaTeX tables of the intensities with errors for the paper. 
    The output text must be copied and pasted into the .tex file
'''

# choose the set of object to make each table: 1= from obj 0-5, 2= from obj 6 to 11, 3= from obj 12 to 17
obj_set = 3
# do you want to use line measurements from analysis with E(B-V) or with CHbeta?
CHbeta = True

#########################################

txt_out = '../results/IntensitiesTable'

#                  0          1           2         3         4        5           6             7         8
objects_list = ['mrk960', 'sbs0218', 'mrk1087', 'ngc1741', 'mrk5', 'mrk1199', 'iras08208', 'iras08339', 'sbs0926',
                'arp252', 'sbs0948', 'tol9',  'sbs1054', 'pox4', 'sbs1319', 'sbs1415', 'tol1457', 'iiizw107']
#                  9         10        11        12        13       14         15         16         17

#                            0     1     2     3      4       5      6      7     8
use_given_lineinfo_list = [False, False, True, False, False, False, True, True, False,
                           False, True, True, True, False, False, False, False, False]
#                            9     10    11    12     13     14     15    16    17
if CHbeta:
    #                            0     1     2     3      4       5      6      7     8
    use_given_lineinfo_list = [False, False, False, False, False, False, False, True, False,
                               False, True, True, False, False, False, False, False, False]
    #                            9     10     11     12    13     14     15    16    17
    
if obj_set == 1:
    objects_list = objects_list[0:6]
    use_given_lineinfo_list = use_given_lineinfo_list[0:6]
elif obj_set == 2:
    objects_list = objects_list[6:12]
    use_given_lineinfo_list = use_given_lineinfo_list[6:12]
elif obj_set == 3:
    objects_list = objects_list[12:18]
    use_given_lineinfo_list = use_given_lineinfo_list[12:18]
print objects_list


lines4paper = '../results/lineIntensities4paper.txt'

#### FUNCTIONS
def read_lines4paper(text_file):
    refwav = []
    IDs = []
    tf = open(text_file, 'r')
    list_rows_of_file = tf.readlines()
    tf.close()
    for row in list_rows_of_file:
        # Split each row into columns
        line_data = row.split()
        refwav.append(float(line_data[0]))
        kk = string.split(line_data[1], sep='&')
        kk2 = kk[1]  
        IDs.append(kk2)
    return refwav, IDs

def read_reasuredIs(text_file):
    wav = []
    flambda = []
    ele = []
    ion =[]
    forb = []
    howforb = []
    F = []
    Ferr = []
    Ferrper = []
    I = [] 
    Ierr = [] 
    Ierrper = []
    EW = []
    EWerr = []
    EWerrper = []
    cols_in_file = [wav, flambda, ele, ion, forb, howforb, F, Ferr, Ferrper, I, Ierr, Ierrper, EW, EWerr, EWerrper]
    tf = open(text_file, 'r')
    list_rows_of_file = tf.readlines()
    tf.close()
    for row in list_rows_of_file:
        # Disregard comment symbol
        if '#'  not in row:
            # Split each row into columns
            line_data = row.split()
            # append the element into each column in the cols_in_file
            for item, col in zip(line_data, cols_in_file):
                if '.' in item:
                    item = float(item)
                col.append(item)
    return cols_in_file

def find_mcnd(intensity, error, obj_set, last=None):
    line = ''
    if obj_set == 1:
        andpercent = '&&'
    else:
        andpercent = '&'
    if intensity <= 0.0:
        intensity = '\mcnd'
        error = ''
        line = '{:>15} {:<2}'.format(intensity, andpercent)
        if last is not None:
            line = '{:>15} \\\\'.format(intensity)
    else:
        err_perc = (error *100.0) / intensity 
        if err_perc >= 41.0:
            line = ' $<$&{:<10} {:<2}'.format(intensity, andpercent)
        else:  
            line = '{:<6}$\pm$&{:<4} {:<2}'.format(intensity, error, andpercent)
        if last is not None:
            if err_perc >= 41.0:
                line = ' $<$&{:<10} \\\\'.format(intensity)
            else:
                line = '{:<6}$\pm$&{:<4} \\\\'.format(intensity, error)          
    return line


#### CODE

refwav, IDs = read_lines4paper(lines4paper)
all_wavs = []
all_Is = []
all_Iers = []
 
for object_number in range(len(objects_list)):
    # get the correct intensities file
    use_given_lineinfo = use_given_lineinfo_list[object_number]
    if use_given_lineinfo:
        file2use = '_measuredLI_RedCor_Ebv.txt'   # my manual line measurements
    else:
        file2use = '_RedCor_Ebv.txt'   # my code's line measurements
    object_name = objects_list[object_number]
    measured_lines_file = 'carbon/results/'+object_name+'/'+object_name+file2use
    measured_lines_file_path_list = string.split(os.getcwd(), sep='carbon')
    measured_lines_file_path = os.path.join(measured_lines_file_path_list[0], measured_lines_file)
    print 'using: ', measured_lines_file_path
    cols_in_file = read_reasuredIs(measured_lines_file_path)
    wav, flambda, ele, ion, forb, howforb, F, Ferr, Ferrper, I, Ierr, Ierrper, EW, EWerr, EWerrper = cols_in_file
    # find lines for paper
    wav = np.round(wav, decimals=0)
    found_wav = []
    found_I = []
    found_Ier = []
    for rw, rid in zip(refwav, IDs):
        i = -1.0
        ier = 0.0
        if rw in wav:
            idx = np.where( wav==rw )[0][0]
            if ele[idx] in rid:
                i = np.round(I[idx], decimals=1)
                ier = np.round(Ierr[idx], decimals=1)
        # append the intensities and errors to a list for only this object
        found_wav.append(rw)
        found_I.append(i)
        found_Ier.append(ier)
        #print '{:<8} & {:<25} & {:>6}$\pm${:<4} &&'.format(int(rw), rid, i, ier)
    
    #print len(refwav), len(wav)
    
    # append the list of intensities of this object to the list of all objects
    all_Is.append(found_I)
    all_Iers.append(found_Ier)      

#print len(refwav), len(all_Is[0]), len(all_Is[1]), len(all_Is[2]), len(all_Is[3]), len(all_Is[4]), len(all_Is[5])

# Now save it to a text file
tout = txt_out+'_set'+repr(obj_set)+'.txt'
tf = open(tout, 'w+')
for idx, rw in enumerate(refwav):
    i1, i2, i3, i4, i5, i6 = all_Is[0][idx], all_Is[1][idx], all_Is[2][idx], all_Is[3][idx], all_Is[4][idx], all_Is[5][idx]
    ier1, ier2, ier3, ier4, ier5, ier6 = all_Iers[0][idx], all_Iers[1][idx], all_Iers[2][idx], all_Iers[3][idx], all_Iers[4][idx], all_Iers[5][idx]
    line1 = find_mcnd(i1, ier1, obj_set)
    line2 = find_mcnd(i2, ier2, obj_set)
    line3 = find_mcnd(i3, ier3, obj_set)
    line4 = find_mcnd(i4, ier4, obj_set)
    line5 = find_mcnd(i5, ier5, obj_set)
    line6 = find_mcnd(i6, ier6, obj_set, last=True)
    rwid = '{:<8} & {:<25} &'.format(int(rw), IDs[idx])
    print rwid, line1, line2, line3, line4, line5, line6
    print >> tf, rwid, line1, line2, line3, line4, line5, line6
tf.close()
 
print 'This file contains the LaTeX table: ', tout
print 'Done!'
