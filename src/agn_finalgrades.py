import numpy as np
from matplotlib import pyplot as plt
import pylab 
import matplotlib.mlab as mlab
import math

fname = 'agn_day1.txt'

propnum, means, devs = np.loadtxt(fname, skiprows=2, usecols=(0,11,12), unpack=True)

props = np.array([])
for i, p in enumerate(propnum):
    props = np.append(props, i)

MEAN = sum(means)/len(means)
var = (means - MEAN)**2
stdev = np.sqrt( 1.0/(len(means)-1.0) * sum(var))
print 'mean =', MEAN, '    standard dev =', stdev

stmeans = (means - MEAN)/stdev
#print stmeans

#histo = np.histogram(stmeans, bins=5, range=(-3.,3.))
#print histo
#plt.plot(histo[0])

# the histogram od the data
n, bins, patches = pylab.hist(stmeans, bins=10, histtype='stepfilled')
pylab.setp(patches, 'facecolor', 'b', 'alpha', 0.75)

# add a line showing the expected distribution
yh = pylab.normpdf(bins, 0.0, stdev)
lh = pylab.plot(bins, yh, 'k--', linewidth=1.5)

m = 0
variance = 1
sigma = math.sqrt(variance)
x = np.linspace(-3,3,100)
height = len(means)*100.0 / 40.0
y = mlab.normpdf(x,m,sigma) * height
#plt.plot(x,y)

plt.title('Partial Final Grading ($\mu$=2.97)') 
plt.xlabel('$\sigma$ [1 $\sigma$=0.68]')
plt.ylabel('Proposals graded')
#plt.plot(stmeans, props, 'ko')
fig = 'partial_final_grading.jpg'
#plt.savefig(fig)
plt.show()