import numpy as np
import emcee
import matplotlib.pyplot as pl

'''
This is a test code to demonstrate how to draw samples from the multivariate Gaussian density given by:
            p(x) propto exp[ -1/2 * (x-mu)^t * Sum^-1 * (x-mu) ],
where mu is an N-dimensional vector position of the mean of the density and Sum is the square N-by-N covariance matrix. 
'''

# function that returns the density p(x) for specific values of x, mu, and Sum^-1.
# NOTE THAT: emcee actually requires the natural logarithm of p.
def lnprob(x, mu, icov):
    diff = x-mu   # position of the SINGLE walker 
    return -np.dot(diff, np.dot(icov,diff))/2.0

# set specific values of the "hyperoarameters" in 50 dimensions
ndim = 50 

# Create an array of the given shape and propagate it with random samples from a uniform distribution over [0, 1).
means = np.random.rand(ndim)

# cov = Sum
cov = 0.5 -np.random.rand(ndim ** 2).reshape((ndim, ndim))  # gives a new shape to an array without changing its data
cov = np.triu(cov)   # Return a copy of a matrix with the elements below the k-th diagonal zeroed
cov += cov.T - np.diag(cov.diagonal())   # Return specified diagonals
cov = np.dot(cov, cov)   # For N dimensions it is a sum product over the last axis of a and the second-to-last of b: 
#                                                dot(a, b)[i,j,k,m] = sum(a[i,j,:] * b[k,:,m])

# but we need the inverse of cov
icov = np.linalg.inv(cov)

# for this example we will use 250 walkers, but we need to guess a starting point for each walker SO
# the initial guess should be a 250-by-50 array, or a list of 250 arrays that each have 50 elements.
# Though it is not a very good guess, we will guess a random number between 0 and 1 for each component.
nwalkers = 250
p0 = np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))

# set up the ensemble sampler object
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[means, icov]) 
# this is calling lnprob(p, means, icov), where p is the position of a single walker. Without the args parameter,
# the function would be called as lnprob(p).

# Allow the first 100 steps to be "burn-in" steps for the MCMC chain (100 is made up, does not have to be that but 
# is a good idea to have the first few steps to let the walkers explore the parameter space a bit and get settled into
# the maximum of the density. This starts with p0:
pos, prob, state = sampler.run_mcmc(p0, 100)   # save the final position in pos and restart from there 
sampler.reset()
sampler.run_mcmc(pos, 1000)   # take 1000 steps staring from pos
# The result is a numpy array with the shape (250, 1000, 50). This object can be flattened to a shape of (250000, 50)
# using the EnsembleSampler.flatchain object. This produced  250,000 unbiased samples of the density p(x), which
# can now be easily plotted in a histogram:
for i in range(ndim):
    if i==0 or i==24 or i==49:
        pl.figure()
        pl.hist(sampler.flatchain[:,i], 100, color='k', histtype='step')
        pl.title('Dimension {0:d}'.format(i))
pl.show()

# To check if everything went as planed, print the mean acceptance fraction of the ensemble using
# EnsebleSampler.acceptance_fraction(), which should be between 0.25 and 0.5 for this example.
print ('Mean acceptance fraction {0:.2f}'.format(np.mean(sampler.acceptance_fraction)))

