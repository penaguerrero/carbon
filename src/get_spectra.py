import os
import glob
import metallicity


'''
This code uses the 1d fits extraction files from all three filters (g230l, g430l, and g750l), 
in order to let me choose which one has the best S/N. 
'''

path_oneDspecs = '../HSTdata/'
path_results = '../results/1Dspecs/'
path_plots = '../results/plots/'

print('Please type first 3 digits of obqn number for the object of interest (i.e. obqn18010_... type 180)')
print('    OR ')
print('type the full path to a .txt file containing the HST obs number (obqn number) and the name of the object separated by a coma.')
print('      i.e.  180, iiizw108')
print('            030, mrk1087')

obqn_num = raw_input()
oneDspecs = glob.glob(path_oneDspecs+'obqn'+obqn_num+'*x*.fits')

print('Please type name of object or some identifier (i.e. ngc6822)')
obj_name = raw_input()
# Create the name for the plot
selectedspecs_file = obj_name+'_selectedspecs'
txtf = path_results+selectedspecs_file
plot_name = os.path.join(path_plots, obj_name+'.eps')
# plot the spectra for the input object
used_specs = metallicity.OneDspecs(oneDspecs, txtf, plot_name)

print("I'm done and I live.")

