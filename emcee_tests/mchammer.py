#!/usr/bin/env python
# -----------------------------------------------------------------------------
# MCHAMMER
#
#   Learning the emcee package.  This simple example generates points from a
#   cubic equation (with some scatter) and then attempts to recover the
#   original parameters of the cubic using MCMC.
#
#   v1.0 : Laura L Watkins [lauralwatkins@gmail.com] - MPIA, 2012/03/15
# -----------------------------------------------------------------------------

import numpy as np
import emcee
import matplotlib.pyplot as plt

def myline( x, p ):
    return p[0] + p[1] * x + p[2] * x**2 + p[3] * x**3

def lnlike( p, x, y, e ):
    m = myline( x, p )
    csq = ( y - m )**2 / e**2
    csq = csq.sum()
    return -csq

def lnprior(theta):
    a, b, c, d = theta
    if 10. < a < 15. and -10. < b < 10.0 and -10. < c < 10.0 and -10. < d < 10.0:
        return 0.0
    return -np.inf    
    
def lnprob( theta, x, y, e ):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, e)

# original parameters
aa = 12.
bb = -5.
cc = -3.
dd = 2.
pp = [ aa, bb, cc, dd ]

# generate data points
n = 100
x = np.random.rand( n ) * 10. - 5.
x.sort()
y = myline( x, pp )                     # perfect model
dy = np.random.rand( n ) * 50. - 25.   # add scatter to data points
e = np.random.rand( n ) * 50. - 25.    # errors given to mcmc

# plot model and data
plt.figure()
plt.plot( x, y, "b", lw=5, alpha=0.4 )
plt.plot( x, y + dy, "ko", alpha=0.5 )

# mcmc
ndim, nwalkers, nruns = 4, 50, 100
pos = [np.random.rand(ndim) * 30. - 15. for i in xrange(nwalkers)]
sampler = emcee.EnsembleSampler( nwalkers, ndim, lnprob, args=[x, y+dy, e] )

pos, lnp, rstate = sampler.run_mcmc( pos, nruns )
'''
f = open( "chain.dat", "wb" )
f.write( "# Omega Centauri MCMC\n"
    + "# walkers: {:}\n".format( nwalkers )
    + "# dims: {:}\n".format( ndim )
    + "# runs: {:}\n".format( nruns )
    + "# " + "-" * 77 + "\n" );
f.close()
# for i in range( nwalkers ):
#     p = pos[i]
#     fmt_model = "{:6.3f}  {:6.3f}  ".format( p[0], p[1] ) \
#         + "{:6.3f}  {:6.3f}  {:9.3f}\n".format( p[2], p[3], lnp[i] )
#     
#     # write_model to output file
#     f = open( "chain.dat", "ab" )
#     f.write( fmt_model )
#     f.close()
'''
# To store the chain....
f = open("mchammer_chain.dat", "w")
f.close()
count = 1
for posn, prob, state in sampler.sample( pos, iterations=20, storechain=True ):
    print "COUNT", count
    if count % 1 == 0:
        f = open("mchammer_chain.dat", "a")
        for k in range( posn.shape[0] ):
            strout = ""
            for p in pos[k]: strout += "{:8.3f} ".format( p )
            strout += "{:20.3f}".format( prob[k] )
            print strout
            f.write(strout+"\n")
        f.close()
    count += 1


# reset the chain to remove the burn-in samples
sampler.reset()

# starting from the final position in the burn-in chain
pos, lnp, rstate = sampler.run_mcmc( pos, 20, rstate0=rstate )

# best model
wh = np.where( lnp == lnp.max() )[0][0]
p = pos[ wh, : ]
plt.plot( x, myline( x, p ), "r", lw=5, alpha=0.4 )

for p in pos:
    plt.plot( x, myline( x, p ), "r", alpha=0.1 )

plt.show()
