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
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
# if wanting to compare all objects type all, else type in a list the objects to make the comparison
object_numbers = [0,1]

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

print code_measurements_list[0]['O'][0]


print '\n ... Code finished!'
