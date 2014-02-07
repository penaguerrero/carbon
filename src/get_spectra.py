import os
import glob
import metallicity
from PIL import Image

'''
This code uses the 1d fits extraction files from all three filters (g230l, g430l, and g750l), 
in order to let me choose which one has the best S/N. 

# reference IDs for sample objects:
0. iiizw107 = 180              9. sbs0948 = 110
1. iras08339 = 080            10. sbs0926 = 190
2. mrk1087 = 030              11. sbs1054 = 130
3. mrk1199 = 060              12. sbs1319 = 150
4. mrk5 = 050                 13. tol1457 = 170
5. mrk960 = 010               14. tol9 = 120
6. ngc1741 = 040              15. arp252 = 100
7. pox4 = 140                 16. iras08208 = 070
8. sbs0218 = 020              17. sbs1415 = 160
'''

path_oneDspecs = '../HSTdata/'
path_oneDspecs_notMain = '../HSTdata/not_main_specs/'
path_results = '../results/1Dspecs/'
path_plots = '../results/plots/'

print('Please type first 3 digits of obqn number for the object of interest (i.e. obqn18010_... type 180)')
print('    OR ')
print('type the full path to a .txt file containing the HST obs number (obqn number) and the name of the object separated by a coma.')
print('      i.e.  180, iiizw107')
print('            030, mrk1087')

obqn_num = raw_input()

# In case of wanting to study the non main spectra, change to True
main_specs = True
if main_specs == False:
    oneDspecs = glob.glob(path_oneDspecs_notMain+'obqn'+obqn_num+'*x*.fits')    

# Read the spectra...
# This is assuming that the fits files have been cleaned from cosmic rays with the cosmic_rays script
oneDspecs = glob.glob(path_oneDspecs+'obqn'+obqn_num+'*x*.fits')    

print('Please type name of object or some identifier (i.e. ngc6822)')
obj_name = raw_input()
# Create the name for the plot
selectedspecs_file = obj_name+'_selectedspecs'
txtf = path_results+selectedspecs_file
img_format = '.jpg'
plot_name = os.path.join(path_plots, obj_name+img_format)
# plot the spectra for the input object
used_specs = metallicity.OneDspecs(oneDspecs, txtf, plot_name)

print("I'm done and I live.")

