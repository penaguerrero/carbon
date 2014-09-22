import numpy as np
from matplotlib import pyplot as plt
import emcee
import triangle
import scipy.optimize as op

'''
This script exemplifies how to fit a line to noisy data (where we do not believe the errors).  
'''

# true parameters for the line
m = 1.2
b = 5.5
pp = [m, b]

# the true line
def line_eq(theta, x):
    m, b = theta
    return m * x + b

# generate the points
n = 100   # number of total points
x = 10*np.random.rand(n)   # generate random numbers between 0.0 and 1.0 and multiply it times 10
np.sort(x)
y = line_eq(pp, x)
# generate the scatter to the data points for the "observed" data
dy = [np.random.uniform(-2., 3.) for _ in range(n)]
yerr = np.array([np.random.uniform(-2., 2.) for _ in range(n)])   # error in the scattered values

# likelihood function
def lnlike(theta, xobs, yobs, yerrobs):
    model = line_eq(theta, xobs)
    chi2 = (yobs - model)**2 / yerrobs**2
    chi2 = chi2.sum()
    return - chi2/2.0   # we are returning the log of the likelihood function

# define the priors
# lets say that we knew  0.0 < m < 5.5  and that  0.0 < b < 10.0
def lnprior(theta):
    m, b = theta
    if 0.0 < m < 5.5 and 0.0 < b < 10.0:
        return 0.0
    return -np.inf

# then the probability function will be
def lnprob(theta, x, y, yerr):
    #print 'theta = ', theta
    #print 'x =', x
    #print 'y =', y
    #print 'yerr =', yerr
    #raw_input()
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)

# now set stuff for the mcmc, start with number of dimensions, walkers, and number of runs 
ndim, nwalkers, nruns = 2, 100, 50
# now we need a starting point for each walker, a ndim-vector to get a nwalkers-by-ndim array.
# there are 2 ways to initialize the walkers:
'''
# a) a small Gaussian ball around the maximum likelihood result, for which we use optimize
nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [m, b], args=(x, y, yerr))
m_ml, b_ml = result["x"]
pos = [result["x"] + .1*np.random.rand(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
'''
# b) randomly and letting the first few walks to be "burn-in" in order to explore the parameter space
p0 = np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[x, y, yerr])
pos, prob, state = sampler.run_mcmc(p0, 50)   # this allows for the first 50 steps to be "burn-in" type

# To store the chain....
f = open("line_chain.dat", "w")
f.close()
count = 1
for posn, prob, state in sampler.sample( pos, iterations=20, storechain=True ):
    print "COUNT", count
    if count % 1 == 0:
        f = open("line_chain.dat", "a")
        for k in range( posn.shape[0] ):
            strout = ""
            for p in pos[k]: strout += "{:8.3f} ".format( p )
            strout += "{:20.3f}".format( prob[k] )
            print strout
            f.write(strout+"\n")
        f.close()
    count += 1

sampler.reset()   # then restart the mcmc at the final position of the "burn-in", pos

pos, prob, state = sampler.run_mcmc(pos, nruns)   # do the mcmc which nruns steps
#print 'pos=', pos
#print 'prob=', prob
#print 'state=', state

# plot the true line and the "observed" data
plt.figure()
plt.plot( x, y, "b", lw=5, alpha=0.4 )   # the true line
plt.plot( x, y+dy, "ko", alpha=0.5 )   # the "observed" data
# best model
wh = np.where( prob == prob.max() )[0][0]
p = pos[ wh, : ]
#print 'parameters that best fit the data are: [m  b] =', p
plt.plot( x, line_eq( p, x ), "r", lw=5, alpha=0.4 )
line1 ='values of the %i dimensions that best fit the data in %i runs are the following:' % (ndim, nruns)
line2 = 'm = %0.3f   b = %0.3f' % (p[0], p[1])
print line1
print line2
f = open("line_chain.dat", "a")
f.write(line1+"\n")
f.write(line2)
f.close()

for p in pos:
    plt.plot( x, line_eq(p, x), "r", alpha=0.1 )
plt.show()

samples = sampler.chain[:, nruns*0.2, :].reshape((-1, ndim))
fig = triangle.corner(samples, labels=["$m$", "$b$"], truths=[m, b])
fig.show()
fig.savefig("triangle_test.jpg")

# Calculate the uncertainties based on the 16th, 50th and 84th percentiles
samples[:, 1] = np.exp(samples[:, 1])
m_mcmc, b_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
print 'm_mcmc, b_mcmc:', m_mcmc, b_mcmc
m_mcmc2, b_mcmc2 = map(lambda v: (v), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
print 'm_mcmc2, b_mcmc2:', m_mcmc2, b_mcmc2
