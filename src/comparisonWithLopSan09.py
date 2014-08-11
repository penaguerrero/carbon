import numpy
import os
#from matplotlib import pyplot as plt

############################################################################################################################################
'''
This scripts takes the found metallicities with my measurements and those by the code and produces:
 * table of comparison values of all available abundances
'''

# Go into the directory of each object. Files assumed to be in: /Users/home_direcotry/Documents/AptanaStudio3/src/
results_path = "../results/"
# but just to make sure we are in the right place, get and work with full paths
full_results_path = os.path.abspath(results_path)
# read the results file
resfile = os.path.join(full_results_path, 'comparison_abunds.txt')
rf = numpy.loadtxt(resfile, skiprows=7, usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28), unpack=True)
O, Oerr, Oother, Oerrother, N, Nerr, Nother, Nerrother, Ne, Neerr, Neother, Neerrother, S, Serr, Sother, Serrother, Ar, Arerr, Arother, Arerrother, Fe, Feerr, Feother, Feerrother, T, Terr, Tother, Terrother = rf

# These are the results are for the sample objects in the following order:
#                objects in the sample 
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']

