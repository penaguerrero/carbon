import numpy as np
import matplotlib
import matplotlib.pyplot as plt

'''
This scripts creates a plot of the number of flat slopes of C/O vs O/H diagram.
'''

####################################################################################################################################

### Set parameters

save_plot = False
img_type = 'jpg'

####################################################################################################################################


#flat_COslope = ['Mrk 5', 'Mrk 960', 'NGC 1741', 'SBS 0948+532', 'SBS 0926+606A',
#                'SBS 1319+579', 'Tol 1457-262', 'IRAS 08208+2816']
#increasing_COslope = ['III Zw 107', 'IRAS 08339+6517', 'Mrk 1087', 'Mrk 1199', 'POX 4', 
#                      'SBS 0218+003', 'SBS 1054+365', 'Arp 252', 'Tol 9', 'SBS 1415+437']

# Line equation: m           b            m_err                  b_err
mbmrk5 = [-0.60451719,  4.47074333, 0.012138200677538969, 0.10136121195546967]
mbmrk960 = [ 0.71005708, -6.3905215, 0.008531309604514677, 0.066198697796159611 ]
mbngc1741 = [ 0.49493364, -4.50215331, 0.0071237308982452533, 0.056456682416435748]
mbsbs0948 = [-0.18366901,  0.71177388, 0.0087244782012557091, 0.070385970721644583]
mbsbs0926 = [0.44421684, -3.99120974, 0.0076071577964996664, 0.060978080707455716]
mbsbs1319 = [0.36054584, -3.48130119, 0.012437443664508467, 0.10215894075316137]
mbtol1457 = [0.30303092, -3.31014405, 0.010403330668201182, 0.083833133850903505]
mbiras08208 = [-0.21022537,  0.7441178, 0.0092989733482604954, 0.075986519332696428]

mbiiizw107 = [0.7259597,  -6.43744034, 0.01009782073166971, 0.080189510468095887]
mbiras08339 = [1.04368186, -8.6068914, 0.0099873450439584809, 0.078719797349716125]
mbmrk1087 = [0.73216924, -6.54880978, 0.01253610328034193, 0.10042472661401323]
mbmrk1199 = [0.81410212, -6.81234656, 0.0074432180590034511, 0.058185782542267736]
mbpox4 = [-0.06654894,  0.13408691, 0.02059001513073503, 0.17060611321416666]
mb0218 = [0.70832901, -6.24874429, 0.0081624377822432804, 0.064123722170562267]
mbsbs1054 = [1.00461765, -8.65622485, 0.010575343752986337, 0.082293726506359563]
mbarp252 = [0.58879739, -4.81748649, 0.0086043423471880622, 0.066451611967226404]
mbtol9 = [0.76768818, -6.68925723, 0.010004139590394539, 0.078553062544415686]
mbsbs1415 = [0.45205676, -4.75039861, 0.0086447512512986911, 0.068962493868314281]

full_names_list = ['MRK 960', 'SBS 0218+003', 'MRK 1087', 'NGC 1741', 'MRK 5', 'MRK 1199', 
                   'IRAS 08208+2816', 'IRAS 08339+6517', 'SBS 0926+606', 'ARP 252', 'SBS 0948+532', 'TOL 9', 
                   'SBS 1054+365', 'POX 4', 'SBS 1319+579', 'SBS 1415+437', 'TOL 1457-262', 'IIIZw 107']

objects_list = ['mrk960', 'sbs0218', 'mrk1087', 'ngc1741', 'mrk5', 'mrk1199', 'iras08208', 'iras08339', 'sbs0926',
                'arp252', 'sbs0948', 'tol9',  'sbs1054', 'pox4', 'sbs1319', 'sbs1415', 'tol1457', 'iiizw107']

all_slopes = np.array([mbmrk960[0], mb0218[0], mbmrk1087[0], mbngc1741[0], mbmrk5[0], mbmrk1199[0],
                       mbiras08208[0], mbiras08339[0], mbsbs0926[0], mbarp252[0], mbsbs0948[0], mbtol9[0], 
                       mbsbs1054[0], mbpox4[0], mbsbs1319[0], mbsbs1415[0], mbtol1457[0], mbiiizw107[0]])
all_slopes_errs = np.array([mbmrk960[2], mb0218[2], mbmrk1087[2], mbngc1741[2], mbmrk5[2], mbmrk1199[2],
                       mbiras08208[2], mbiras08339[2], mbsbs0926[2], mbarp252[2], mbsbs0948[2], mbtol9[2], 
                       mbsbs1054[2], mbpox4[2], mbsbs1319[2], mbsbs1415[2], mbtol1457[2], mbiiizw107[2]])

mean_slope = np.mean(all_slopes)
median_slope = np.median(all_slopes)
print 'mean_slope = ', mean_slope # =0.449179206667
print 'median_slope = ', median_slope # =0.541865515

# PLOT
font = {#'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}
matplotlib.rc('font', **font)

fig = plt.figure(1, figsize=(10, 10))
#yup = 1.1
#ylo = -1.9
plt.xlim(0, 19)
#plt.ylim(ylo, yup)
#plt.xticks(np.arange(0, 19, 1))
#plt.yticks(np.arange(-1.8, yup, 0.2))
plt.minorticks_on()
plt.xlabel('Galaxy Number Within Sample')
plt.ylabel("Slope of log(C/O) vs 12+log(O/H)")
number_objects = range(1,19)
plt.plot(number_objects, all_slopes, 'rD', markersize=9)
#plt.errorbar(number_objects, all_slopes, yerr=all_slopes_errs, fmt='rD', ecolor='k', markersize=9)
for x, y, ye, z in zip(number_objects, all_slopes, all_slopes_errs, objects_list):
    subxcoord = -3
    subycoord = 5
    side = 'right'
    if (z == 'sbs0948') or (z == 'mrk5') or (z == 'tol1457'):
        subxcoord = -4
        subycoord = -12
    if (z == 'mrk1087') or (z == 'ngc1741'):
        subxcoord = 5
        side = 'left'
    if (z == 'sbs1054') or (z == 'pox4') or (z == 'mrk1199'):
        subxcoord = 4
        subycoord = -12
        side = 'left'
    if (z == 'mrk960'):
        subxcoord = -17
        subycoord = 11
        side = 'left'
    if (z == 'sbs0218'):
        subxcoord = -17
        subycoord = -14
        side = 'left'
    plt.annotate('{}'.format(z), xy=(x,y), xytext=(subxcoord, subycoord), ha=side, textcoords='offset points')


'''
objects = ['Flat slope', 'Increasing slope']
y_pos = np.arange(len(objects))
performance = [len(flat_COslope), len(increasing_COslope)]
plt.bar(y_pos, performance, align='center', alpha=0.5)
plt.xticks(y_pos, objects)
plt.ylabel('Number of objects')
plt.title('Slope type of C/O vs O/H diagram')
'''
if save_plot:
    figname = '../results/plots/COslopes.'+img_type
    fig.savefig(figname)
    print 'Figure ', figname, ' saved!'

plt.show()
