import os

'''
This script converts from .jpg format to .eps for all objects at once.
'''

objects_list = ['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 
                'ngc1741', 'pox4', 'sbs0218', 'sbs0948', 'sbs0926', 'sbs1054', 
                'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']

# Go into the directory of each object. Files assumed to be in: /Users/home_direcotry/Documents/AptanaStudio3/src/
plots_path = "../results/plots/"
# but just to make sure we are in the right place, get and work with full paths
full_results_path = os.path.abspath(plots_path)

for obj in objects_list:
    join_path = os.path.join(full_results_path, obj)
    os.system('convert '+join_path+'_suplots.jpg '+join_path+'_suplots.eps')

print ' Script finished!'
