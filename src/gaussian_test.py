import numpy as np
#import os
#import metallicity
from science import spectrum
from matplotlib import pyplot as plt
from scipy import optimize
from scipy import stats

region = 2

if region == 1:
    region = 'nuv'
elif region == 2:
    region = 'opt'
elif region == 3:
    region = 'nir'

spec = '../results/iiizw107/iiizw107_'+region+'spec_NORMandREBIN.txt'


### Functions
    
def gauss(x, a, meanX, sig):
    return a*stats.norm.pdf(x, loc=meanX, scale=sig)

def fit_gauss(x, y):
    sig, mean = spectrum.find_std(x)
    # p0 is the initial guess for the fitting coefficients (mean and sigma)
    a = 0.2 # initial guess
    p0 = [a, mean, sig]
    coeff, _ = optimize.curve_fit(gauss, x, y, p0=p0)
    gf = gauss(x, *coeff)
    print 'these are the coefficients:  norm=', coeff[0], '  mean=', coeff[1], '  sigma=', coeff[2]
    return gf, sig

def multi_gauss(x, *args):
    m1, m2, s1, s2, a1, a2 = args
    ret = gauss(x, a1, m1, s1)#a1 * stats.norm.pdf(x, loc=m1, scale=s1)
    ret += gauss(x, a2, m2, s2)
    return ret

def fit_multi_gauss(x, y, m1, m2):
    s1, s2, a1, a2 = 1.2, 1.2, 3.0e-15, 3.0e-15
    p0 = [m1, m2, s1, s2, a1, a2]
    coeffs, _ = optimize.curve_fit(multi_gauss, x, y, p0=p0)
    mg = multi_gauss(x, *coeffs)
    print 'these are the coefficients:  ', coeffs
    return mg, coeffs



### Code

wav, flx = np.loadtxt(spec, unpack=True)
#data_arr = [wav, flx]
cont_arr0 = spectrum.theo_cont(wav, scale_factor=0.0)   # this is an array of wavelength and continuum

lolim = 3710.0#4990.0
uplim = 3740.0#5024.5
linew = wav[(wav>=lolim) & (wav<=uplim)]
linef = flx[(wav>=lolim) & (wav<=uplim)]
data_arr = [linew, linef]
cont_arr = spectrum.theo_cont(linew, scale_factor=0.0)
elements = 100
wavelength, flux, _, flux_cont = spectrum.fill_EWarr(data_arr, cont_arr, lolim, uplim, elements)

# fit the line with a gaussian
gf, sig = fit_gauss(wavelength, flux)
fwhm = spectrum.FWHM(sig)

#meanX = 3726.30454854
mean1 = 3722.0
mean2 = 3729.0
norm1 = 3.0e-15
norm2 = 3.0e-15
sig1 = 2.2
sig2 = 2.2
gsum, coeffs = fit_multi_gauss(wavelength, flux, mean1, mean2)
g1 = gauss(wavelength, coeffs[4], coeffs[0], coeffs[2])
g2 = gauss(wavelength, coeffs[5], coeffs[1], coeffs[3])

plt.xlabel('Wavelength  [$\AA$]')
plt.ylabel('Flux  [ergs/s/cm$^2$/$\AA$]')
plt.xlim(lolim-100.0, uplim+100.0)
plt.plot(wav, flx, 'k')
plt.plot(wavelength, flux, 'b')
plt.plot(wavelength, gf, 'r')
plt.plot(wavelength, gsum, 'y--')
plt.plot(wavelength, g1, 'c--')
plt.plot(wavelength, g2, 'c--')
plt.plot(cont_arr0[0], cont_arr0[1], 'g')
plt.show()

