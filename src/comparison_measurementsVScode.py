import numpy
import os
import metallicity
#from matplotlib import pyplot as plt

############################################################################################################################################
'''
This scripts compares my line manual measuremets with the measurements obtained with the code. 
The script returns:
 * table of comparison values of all available abundances
'''

# This comparison is for the results of the sample objects in the following order:
#                objects in the sample 
#                  0            1           2          3         4        5         6         7         8
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                 9           10         11         12         13       14        15         16         17
# if wanting to compare all objects type all, else type in a list the objects to make the comparison
object_numbers = 'all'

# Do you want to use C_Hbeta to correct for extinction?   (if set to false the values of A_V and A_B will be used)
use_Chbeta = False

# Go into the directory of each object. Files assumed to be in: /Users/home_direcotry/Documents/AptanaStudio3/src/
results_path = "../results/"

# but just to make sure we are in the right place, get and work with full paths
full_results_path = os.path.abspath(results_path)

objects2compare = []
if object_numbers != 'all':
    for obj_number in object_numbers:
        objects2compare.append(objects_list[obj_number])
else:
    objects2compare = objects_list
    
# the results are to be saved in a list of dictionaries with the same order as the objects2compare list:
manual_measurements_list = []
code_measurements_list = []

for obj_name in objects2compare:
    # Go into the object folder
    results4object_path = os.path.join(full_results_path, obj_name)
        
    # get the two files of interest in each object folder
    manual = obj_name+'_measuredLI_IonicTotAbundances_Ebv.txt'
    code = obj_name+'_IonicTotAbundances_Ebv.txt'
    if use_Chbeta:
        manual = obj_name+'_measuredLI_IonicTotAbundances_CHbeta.txt'
        code = obj_name+'_IonicTotAbundances_CHbeta.txt'
    manual_measurements_file = os.path.join(results4object_path, manual)
    code_measurements_file = os.path.join(results4object_path, code)
    files2compare = [manual_measurements_file, code_measurements_file]
        
    # read the results file
    for resfile in files2compare:
        # creating the dictionary of abundances for each object
        elem_abun = metallicity.read_IonTotAbs(resfile)
        if 'measured' in resfile:     
            manual_measurements_list.append(elem_abun)
        else:
            code_measurements_list.append(elem_abun)

# Make the comparisons
Odiffs = []
Ndiffs = []
Nediffs = []
Sdiffs = []
Cldiffs = []
Ardiffs = []
Fediffs = []
Cdiffs = []

print '{:<10} {:>25}'.format('Object', '12 + log O/H')
print '{:<10} {:>13} {:>12} {:>7}'.format(' ', 'manual', 'code', 'diff')
for i in range(len(objects2compare)):
    O_diff = manual_measurements_list[i]['O'][3] - code_measurements_list[i]['O'][3]
    Odiffs.append(O_diff)
    N_diff = manual_measurements_list[i]['N'][3] - code_measurements_list[i]['N'][3]
    Ndiffs.append(N_diff)
    Ne_diff = manual_measurements_list[i]['Ne'][3] - code_measurements_list[i]['Ne'][3]
    Nediffs.append(Ne_diff)
    S_diff = manual_measurements_list[i]['S'][3] - code_measurements_list[i]['S'][3]
    Sdiffs.append(S_diff)
    Cl_diff = manual_measurements_list[i]['Cl'][3] - code_measurements_list[i]['Cl'][3]
    Cldiffs.append(Cl_diff)
    Ar_diff = manual_measurements_list[i]['Ar'][3] - code_measurements_list[i]['Ar'][3]
    Ardiffs.append(Ar_diff)
    Fe_diff = manual_measurements_list[i]['Fe'][3] - code_measurements_list[i]['Fe'][3]
    Fediffs.append(Fe_diff)
    C_diff = manual_measurements_list[i]['C'][3] - code_measurements_list[i]['C'][3]
    Cdiffs.append(C_diff)
    
    print '{:<10} {:>8.2f} {:>5.2f} {:>6.2f} {:>5.2f} {:>6.2f} {:>8.2f}'.format(objects2compare[i], manual_measurements_list[i]['O'][3], manual_measurements_list[i]['O'][4], 
                                                                                    code_measurements_list[i]['O'][3], code_measurements_list[i]['O'][4], 
                                                                                    manual_measurements_list[i]['N'][3])
        
#print '{:<12} {:>60}'.format('Object', 'Diff = manual_measurements - code_measuremets')
#print '{:<12} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}'.format(' ', 'O_diff', 'N_diff', 'Ne_diff', 'S_diff', 'Cl_diff', 'Ar_diff', 'Fe_diff', 'C_diff')
#for obj, ox, ni, ne, su, cl, ar, fe, ca in zip(objects2compare, Odiffs, Ndiffs, Nediffs, Sdiffs, Cldiffs, Ardiffs, Fediffs, Cdiffs):
#    print '{:<12} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f}'.format(obj, ox, ni, ne, su, cl, ar, fe, ca)


out_file = results4object_path = os.path.join(full_results_path, 'diffs_manualVScode.txt')
of = open(out_file, 'w+')
print >> of, '{:<12} {:>60}'.format('Object', 'Diff = manual_measurements - code_measuremets')
print >> of, '{:<12} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}'.format(' ', 'O_diff', 'N_diff', 'Ne_diff', 'S_diff', 'Cl_diff', 'Ar_diff', 'Fe_diff', 'C_diff')
for obj, ox, ni, ne, su, cl, ar, fe, ca in zip(objects2compare, Odiffs, Ndiffs, Nediffs, Sdiffs, Cldiffs, Ardiffs, Fediffs, Cdiffs):
    print >> of, '{:<12} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f}'.format(obj, ox, ni, ne, su, cl, ar, fe, ca)
of.close()

print '\n Text file written: ', out_file
   
print '\n ... Code finished!'
