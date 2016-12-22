import matplotlib
from matplotlib import pyplot as plt


"""
This script compares the temperatures we measured with those presented in LS08.
"""

number_objects = range(1,19)
#                  0          1           2         3         4        5           6             7         8
objects_list = ['mrk960', 'sbs0218', 'mrk1087', 'ngc1741', 'mrk5', 'mrk1199', 'iras08208', 'iras08339', 'sbs0926',
                'arp252', 'sbs0948', 'tol9',  'sbs1054', 'pox4', 'sbs1319', 'sbs1415', 'tol1457', 'iiizw107']
#                  9         10        11        12        13       14         15         16         17

#              0       1       2      3      4      5       6      7     8  
TeO3_paper = [9600., 12600., 7100., 9600., 12500., 6850., 9700., 8400., 13700. ,
              8450., 13200, 7600., 12450., 12500., 13200., 15500., 14300., 10900.]
#               9      10     11     12      13     14      15       16      17

#                  0       1      2      3      4      5      6      7     8  
TeO3_paper_err = [1800., 1000,  1000.,  600., 1200., 1800., 1100., 3000., 900.,
                  1000., 2000., 1000., 1000., 1200., 2000., 1000., 1000., 2300.]
#                  9      10     11     12      13     14      15       16      17
                  
#             0       1       2      3      4      5       6      7     8  
temp_LSE = [9500., 13200., 7100., 9600., 12400., 6950., 10100., 8700., 13600.,
            8700., 13100., 7600., 13700., 14000., 13400., 15500., 14000., 10900.]
#             9      10     11     12      13     14      15       16      17

#                0     1     2     3     4     5     6      7     8  
temp_LSE_err = [800., 600., 900., 600., 700., 800., 700., 1000., 700., 
                900., 600., 1000., 900., 500., 500., 700., 700., 900.]
#                9    10     11     12   13     14    15    16    17

# which ion did the temp come from
s0, s1, s2, s3 = 'LSE', '[OIII]', '[SIII]', '[OIII]+[SIII]'
#               0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17
temp_source = [s2, s3, s0, s0, s3, s2, s1, s1, s1, s2, s1, s0, s1, s3, s3, s1, s1, s1]

# calculate the percentage of error in the temperature
for object_number, obj_name in enumerate(objects_list):
    perc_err_paper = TeO3_paper_err[object_number]*100. / TeO3_paper[object_number]
    perc_err_LSE = temp_LSE_err[object_number]*100. / temp_LSE[object_number]
    print obj_name, TeO3_paper[object_number], temp_LSE[object_number]
    print 'Percentage error: LSE = %0.2f  Paper = %0.2f  -> source: %s' % (perc_err_LSE, perc_err_paper, temp_source[object_number])

# PLOT
font = {#'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}
matplotlib.rc('font', **font)
fig = plt.figure(1, figsize=(12, 10))
plt.title('High Ionization Degree Te Comparison')
plt.xlim(-2, 19)
#plt.ylim(-1.5, 1.5)
xlab = 'Galaxy Number Within Sample'
ylab = 'Temperature [K]'
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.minorticks_on()
for x, y, ye, z, y1, ye1 in zip(number_objects, TeO3_paper, TeO3_paper_err, objects_list, temp_LSE, temp_LSE_err):
    print x, y, ye, z
    if (z=='mrk1087') or (z=='ngc1741') or (z=='tol9'):
        #plt.plot(x, y, 'g*', markersize=10)
        plt.errorbar(x, y, yerr=ye, fmt='g*', ecolor='g', markersize=10)
    else:
        #plt.plot(x, y, 'rd', markersize=9)
        #plt.plot(x, y1, 'bo', markersize=6)
        plt.errorbar(x, y,yerr=ye, fmt='rd', ecolor='r', markersize=9)
        plt.errorbar(x, y1, yerr=ye1, fmt='bo', ecolor='b', markersize=6)
    # Annotate the points 5 _points_ above and to the left of the vertex
    subxcoord = -3
    subycoord = 5
    side = 'right'
    if (z == 'pox4') or (z == 'sbs1319'):
        subxcoord = 5
        side = 'left'
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')

save_plot = False
if save_plot:
    destination = '../results/plots/temps_comparison.jpg'
    fig.savefig(destination)
    print 'Figure ', destination, ' saved!'
plt.show()
    
    
    