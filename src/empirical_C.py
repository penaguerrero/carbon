import numpy as np
import string
import os

'''
This script calculates the differences between the carbon abundances determined from observations and Cloudy to those with the 
relation provided in the C paper.
'''

out_file = '../results/Cdifferences.txt'

#                  0          1           2         3         4        5           6             7         8
objects_list = ['mrk960', 'sbs0218', 'mrk1087', 'ngc1741', 'mrk5', 'mrk1199', 'iras08208', 'iras08339', 'sbs0926',
                'arp252', 'sbs0948', 'tol9',  'sbs1054', 'pox4', 'sbs1319', 'sbs1415', 'tol1457', 'iiizw107']
#                  9         10        11        12        13       14         15         16         17

#                            0     1     2     3      4       5      6      7     8
use_given_lineinfo_list = [False, False, False, False, False, False, True, True, False,
                           False, True, True, True, False, False, False, False, False]
#                            9     10     11     12    13     14     15    16    17


#### Functions

def get_abunds(abunds_file, objects_list):
    #                      0            1           2          3         4        5         6         7         8
    initial_obj_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
                       'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
    #                      9           10         11         12         13       14        15         16         17
    # get the data, which is in the initial order 
    he, o, c, n, ne, s, to3, to2 = np.loadtxt(abunds_file, skiprows=1, usecols=(2,3,4,5,6,7,8,9), unpack=True)
    # order the data according to objects_list
    He, O, C, N, Ne, S, TO3, TO2 = [], [], [], [], [], [], [], [] 
    for obj in objects_list:
        if obj in initial_obj_list:
            idx = initial_obj_list.index(obj)
            He.append(he[idx])
            O.append(o[idx])
            C.append(c[idx])
            N.append(n[idx])
            Ne.append(ne[idx])
            S.append(s[idx])
            TO3.append(to3[idx])
            TO2.append(to2[idx])
    He = np.array(He)
    O = np.array(O)
    C = np.array(C)
    N = np.array(N)
    Ne = np.array(Ne)
    S = np.array(S)
    TO3 = np.array(TO3)
    TO2 = np.array(TO2)
    return He, O, C, N, Ne, S, TO3, TO2

#### Code

# Benchmark abundances
Final_meas_abunds = 'carbon/results/Final_meas_abunds.txt'
file_path_list = string.split(os.getcwd(), sep='carbon')
full_file_path = os.path.join(file_path_list[0], Final_meas_abunds)
#print 'Reading: ', full_file_path
He, O, C, N, Ne, S, TO3, TO2 = get_abunds(full_file_path, objects_list)

# Cloudy abundances
Cldy_best_abunds = 'carbon/results/Cloudy_best_abunds.txt'
Cldy_file_path = os.path.join(file_path_list[0], Cldy_best_abunds)
CldyHe, CldyO, CldyC, CldyN, CldyNe, CldyS, CldyTO3, CldyTO2 = get_abunds(Cldy_file_path, objects_list)

#approxC = 0.5 - 0.18 * (N - O) + N   # previous
#approxC = 0.35 - 0.3 * (N - O) + N   # MCMC - paper
approxC = 0.30 - 0.18 * (N - O) + N   # MCMC without Sun, IZw18, Orion, and 30 Dor

tf = open(out_file, 'w+')
print '{:<15} {:>3} {:>5} {:>10} {:>7} {:>10}'.format('Object', 'DiffObs', 'C/H', 'Approx C', 'CldyC', 'DiffCldy')
print >> tf, '{:<15} {:>3} {:>5} {:>10} {:>7} {:>10}'.format('Object', 'DiffObs', 'C/H', 'Approx C', 'CldyC', 'DiffCldy')
diffs = []
difCldy = []
for obj, c, apc, cldyc in zip(objects_list, C, approxC, CldyC):
    d = apc - c
    diffs.append(np.abs(d))
    dc = apc - cldyc
    difCldy.append(np.abs(dc))
    print '{:<15} {:>5.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f}'.format(obj, d, c, apc, cldyc, dc)
    print >> tf, '{:<15} {:>5.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f}'.format(obj, d, c, apc, cldyc, dc)
avg_diff = sum(diffs) / len(diffs)
avg_diffCldy = sum(difCldy) / len(diffs)

print '\n                        Obs   Cloudy'        
print 'Max difference:       ', max(diffs), max(difCldy)
print 'Min difference:       ', min(diffs), min(difCldy)
print 'Average difference is:', np.round(avg_diff, decimals=4), np.round(avg_diffCldy, decimals=4)
print >> tf,  '\n                        Obs   Cloudy'        
print >> tf,  'Max difference:       ', max(diffs), max(difCldy)
print >> tf,  'Min difference:       ', min(diffs), min(difCldy)
print >> tf,  'Average difference is:', np.round(avg_diff, decimals=4), np.round(avg_diffCldy, decimals=4)
tf.close()
