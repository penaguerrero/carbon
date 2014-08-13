import numpy
import os
#from matplotlib import pyplot as plt

############################################################################################################################################
'''
This scripts finds the difference between our adopted metallicity values and those presented in Lopez-Sanchez & Esteban (2009).
The results are:
 * table of comparison values of all available abundances
'''

# Go into the directory of each object. Files assumed to be in: /Users/home_direcotry/Documents/AptanaStudio3/src/
results_path = "../results/"
# but just to make sure we are in the right place, get and work with full paths
full_results_path = os.path.abspath(results_path)
# read the results file
resfile = os.path.join(full_results_path, 'comparison_abunds.txt')
rf = numpy.loadtxt(resfile, skiprows=6, usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28), unpack=True)
O, Oerr, Oother, Oerrother, N, Nerr, Nother, Nerrother, Ne, Neerr, Neother, Neerrother, S, Serr, Sother, Serrother, Ar, Arerr, Arother, Arerrother, Fe, Feerr, Feother, Feerrother, T, Terr, Tother, Terrother = rf

# These are the results are for the sample objects in the following order:
#                objects in the sample 
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']

default_cols = ['Object', 'log_O/H', 'err', 'log_N/H', 'err', 'log_Ne/H', 'err', 'log_S/H', 'err', 'log_Ar/H', 'err', 'log_Fe/H', 'err', 'Te_[OIII]', 'err']

diff_O = []
diff_N = []
diff_Ne = []
diff_S = []
diff_Ar = []
diff_Fe = []
diff_T = []

colhdr1 = '{:<10} {:>17} {:>4} {:>29} {:>4} {:>29} {:>4} {:>29} {:>4} {:>29} {:>4} {:>29} {:>4} {:>32} {:>4}'
print colhdr1.format(default_cols[0], default_cols[1], default_cols[2], default_cols[3], default_cols[4], default_cols[5], default_cols[6], default_cols[7], default_cols[8], default_cols[9], default_cols[10], default_cols[11], default_cols[12], default_cols[13], default_cols[14])
colhdr2 = '{:<10} {:>13} {:>12} {:>7} {:>13} {:>12} {:>7} {:>13} {:>12} {:>7} {:>13} {:>12} {:>7} {:>13} {:>12} {:>7} {:>13} {:>12} {:>7} {:>15} {:>15} {:>9}'
print colhdr2.format(' ', 'this work', 'Lop-San09', 'diff', 'this work', 'Lop-San09', 'diff', 'this work', 'Lop-San09', 'diff', 'this work', 'Lop-San09', 'diff', 'this work', 'Lop-San09', 'diff', 'this work', 'Lop-San09', 'diff', 'this work', 'Lop-San09', 'diff')
for i in range(len(O)): 
    odiff = O[i] - Oother[i]
    diff_O.append(odiff)
    ndiff = N[i] - Nother[i]
    diff_N.append(ndiff)
    nediff = Ne[i] - Neother[i]
    diff_Ne.append(nediff)
    sdiff = S[i] - Sother[i]
    diff_S.append(sdiff)
    ardiff = Ar[i] - Arother[i]
    diff_Ar.append(ardiff)
    fediff = Fe[i] - Feother[i]
    diff_Fe.append(fediff)
    tdiff = T[i] - Tother[i]
    diff_T.append(tdiff)

    table_format = '{:<10} {:>8.2f} {:>5.2f} {:>6.2f} {:>5.2f} {:>6.2f} {:>8.2f} {:>5.2f} {:>6.2f} {:>5.2f} {:>6.2f} {:>8.2f} {:>5.2f} {:>6.2f} {:>5.2f} {:>6.2f} {:>8.2f} {:>5.2f} {:>6.2f} {:>5.2f} {:>6.2f} {:>8.2f} {:>5.2f} {:>6.2f} {:>5.2f} {:>6.2f} {:>8.2f} {:>5.2f} {:>6.2f} {:>5.2f} {:>6.2f} {:>10} {:<6} {:>8} {:<6} {:>7}'
    print table_format.format(objects_list[i], O[i], Oerr[i], Oother[i], Oerrother[i], odiff, 
                              N[i], Nerr[i], Nother[i], Nerrother[i], ndiff,
                              Ne[i], Neerr[i], Neother[i], Neerrother[i], nediff,
                              S[i], Serr[i], Sother[i], Serrother[i], sdiff,
                              Ar[i], Arerr[i], Arother[i], Arerrother[i], ardiff,
                              Fe[i], Feerr[i], Feother[i], Feerrother[i], fediff,
                              T[i], Terr[i], Tother[i], Terrother[i], tdiff)

    