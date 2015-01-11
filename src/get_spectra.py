import os
import glob
import metallicity
from PIL import Image

'''
This code uses the 1d fits extraction files from all three filters (g230l, g430l, and g750l), 
in order to let me choose which one has the best S/N. 
'''

#####################################################################################################################################

# short name of the object
#                 0           1           2            3         4        5          6        7         8
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                 9           10         11         12         13       14        15         16         17

''' CHOOSE THE OBJECT NUMBER '''
object_number = 4

''' Choose the expected extension of the spectra. '''
img_format = '.eps'

#####################################################################################################################################

# reference IDs for sample objects:
obqn_num_list = ['180', '080', '030', '060', '050', '010', '040', '140', '020', 
                 '110', '190', '130', '150', '170', '120', '100', '070', '160']

full_names_list = ['IIIZw 107', 'IRAS 08339+6517', 'MRK 1087', 'MRK 1199', 'MRK 5', 'MRK 960', 
                   'NGC 1741', 'POX 4', 'SBS 0218+003', 'SBS 0948+532', 'SBS 0926+606', 'SBS 1054+365', 
                   'SBS 1319+579', 'TOL 1457-262', 'TOL 9', 'ARP 252', 'IRAS 08208+2816', 'SBS 1415+437']

obj_name = objects_list[object_number]
obqn_num = obqn_num_list[object_number]
full_name = full_names_list[object_number]

path_oneDspecs = '../HSTdata/'
path_oneDspecs_notMain = '../HSTdata/not_main_specs/'
path_results = '../results/1Dspecs/'
path_plots = '../results/plots/'

# In case of wanting to study the non main spectra, change to True
main_specs = True
if main_specs == False:
    oneDspecs = glob.glob(path_oneDspecs_notMain+'obqn'+obqn_num+'*x*.fits')    

# Read the spectra...
# This is assuming that the fits files have been cleaned from cosmic rays with the cosmic_rays script
oneDspecs = glob.glob(path_oneDspecs+'obqn'+obqn_num+'*x*.fits')    

# Create the name for the plot
selectedspecs_file = obj_name+'_selectedspecs'
txtf = path_results+selectedspecs_file
plot_name = os.path.join(path_plots, obj_name)
# plot the spectra for the input object
used_specs = metallicity.OneDspecs(oneDspecs, txtf, plot_name, img_format, full_name)

print("I'm done and I live.")

