import os
import pyfits
import string
import pylab as plt
from matplotlib import pyplot
#from PIL import Image

import cosmics

'''
This script corrects 2d spectra from cosmic rays and saves it into *_clean.fits, the new file have the same header and the 
same fits extensions so that extractions can be done with the x1d task.
'''

#############################################################################################################################

# name of the object
#                 0           1           2            3         4        5          6        7        8
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                 9           10          11       12          13       14       15         16         17

''' Choose parameters to run script'''
object_number = 0#8#11
# choose the region to correct: opt=1  or  nir=2
region = 1

# paremeters for contrast in plotted image
vmin = 0
vmax = 65
# useful cmap options: Blues, Greens, Greys, hot, Purples, binary, rainbow,  for default type None
cmap = 'Greys'

# Choose the parameters for running the code
# suggested code parameters: data, gain = 2.2, readnoise = 5.0, sigclip = 10.0, sigfrac = 0.3, objlim = 5.0)
# suggested VanDok parameters: data, gain = 7., readnoise = 5.0, sigclip = 5.0, sigfrac = 0.3, objlim = 5.0)
gain = 0.09
readnoise = 0.09
sigclip = 10.0
sigfrac = 0.5
objlim = 5.0
# number of maximum iterations, suggested=4
maxiter = 10

# Does the file already cosmic ray corrected exists?
cosmic_corr_file_exist = True

#############################################################################################################################

# this is assuming we are in the src directory
HSTdata_path = os.path.abspath('../HSTdata')
# Reading the fits file, 2d spectra
'''
reference IDs for sample objects:
0. iiizw107 = 180              9. sbs0948 = 110
1. iras08339 = 080            10. sbs0926 = 190
2. mrk1087 = 030              11. sbs1054 = 130
3. mrk1199 = 060              12. sbs1319 = 150
4. mrk5 = 050                 13. tol1457 = 170
5. mrk960 = 010               14. tol9 = 120
6. ngc1741 = 040              15. arp252 = 100
7. pox4 = 140                 16. iras08208 = 070
8. sbs0218 = 020              17. sbs1415 = 160
'''
# optical region: 3000 to 5500 Angstroms
opt_list = ['obqn18020_crj.fits', 'obqn08020_crj.fits', 'obqn03020_crj.fits', 'obqn06020_crj.fits', 'obqn05020_crj.fits', 
            'obqn01020_crj.fits', 'obqn04020_crj.fits', 'obqn14020_crj.fits', 'obqn02030_crj.fits', 'obqn11040_crj.fits',
            'obqn19030_crj.fits', 'obqn13030_crj.fits', 'obqn15040_crj.fits', 'obqn17020_crj.fits', 'obqn12030_crj.fits',
            'obqn10020_crj.fits', 'obqn07020_crj.fits', 'obqn16030_crj.fits']
# near IR: 5500 adn 10,000 Anstroms
nir_list = ['obqn18030_crj.fits', 'obqn08030_crj.fits', 'obqn03030_crj.fits', 'obqn06030_crj.fits', 'obqn05030_crj.fits', 
            'obqn01030_crj.fits', 'obqn04030_crj.fits', 'obqn14030_crj.fits', 'obqn02040_crj.fits', 'obqn11050_crj.fits',
            'obqn19040_crj.fits', 'obqn13040_crj.fits', 'obqn15050_crj.fits', 'obqn17030_crj.fits', 'obqn12040_crj.fits',
            'obqn10030_crj.fits', 'obqn07030_crj.fits', 'obqn16040_crj.fits']
if region == 1:
    obj2dim = opt_list[object_number]
elif region == 2:
    obj2dim = nir_list[object_number]
img1 = os.path.join(HSTdata_path, obj2dim)

ori_img = pyfits.open(img1)
ori_img.info()
hdr = ori_img[0].header
ori_data = ori_img[1].data
err = ori_img[2].data
dq = ori_img[3].data
ori_img.close()
h1 = ori_img[1].header
h2 = ori_img[2].header
h3 = ori_img[3].header
real_gain = hdr['ccdgain']
print 'real gain =', real_gain, '    used gain = ', gain

# get the object ID, obqn number and unite it with the path to create the temporaty name
obj_id = string.split(obj2dim, sep='.')

if cosmic_corr_file_exist == False:    
    data, h = cosmics.fromfits(img1, hdu=1)
    # building the object
    c = cosmics.cosmicsimage(data, gain=gain, readnoise=readnoise, sigclip=sigclip, sigfrac=sigfrac, objlim=objlim)
    
    # running the script, max iterations suggested = 4
    c.run(maxiter = maxiter)
    
    # unite obqn number with the path to create the temporaty name
    clean_fits_temp = os.path.join(HSTdata_path, obj_id[0]+'_temporary.fits')
    # writting the cleaned image into a new FITS file, conserving the original header:
    cosmics.tofits(clean_fits_temp, c.cleanarray, h)
    
    # mask used
    #cosmics.tofits(obj_id+"mask.fits", c.mask, h)
    
    # open the temporary file to add the original extensions from the multifits file
    img = pyfits.open(clean_fits_temp, mode='update')
    img.info()
    clean_data = img[0].data
    new_fits = os.path.join(HSTdata_path, obj_id[0]+'_clean.fits')
    hdu0 = pyfits.PrimaryHDU(header=ori_img[0].header)
    hdu1 = pyfits.ImageHDU(clean_data, header=h1, name='SCI')
    hdu2 = pyfits.ImageHDU(err, header=h2, name='ERR')
    hdu3 = pyfits.ImageHDU(dq, header=h3, name='DQ')
    hdulist = pyfits.HDUList([hdu0, hdu1, hdu2, hdu3])
    hdulist.writeto(new_fits)
    # Remove the temporary file
    os.remove(clean_fits_temp)
    print 'File  %s  has been erased.' % (clean_fits_temp)
    
elif cosmic_corr_file_exist == True:
    clean_fits = os.path.join(HSTdata_path, obj_id[0]+'_clean.fits')
    img = pyfits.open(clean_fits)
    img.info()
    clean_data = img[1].data

# plot the spectrum before cosmic ray correction as a 2-d array
#ori_data = ori_data[300:750, :]
# Before the cosmic ray correcting
f, (ax1, ax2) = plt.subplots(1, 2, figsize=(20,10))
im1 = ax1.imshow(ori_data, vmin=vmin, vmax=vmax, origin='lower', cmap=cmap)
ax1.set_title(obj_id[0]+'_before')
#im1.colorbar(im1, orientation='vertical')
# After correcting for cosmic rays
#clean_data = clean_data[300:750, :]
im2 = ax2.imshow(clean_data, vmin=vmin, vmax=vmax, origin='lower', cmap=cmap)
ax2.set_title(obj_id[0]+'_after')
cbar_ax = f.add_axes([0.92, 0.15, 0.05, 0.7])
f.colorbar(im1, cax=cbar_ax)
plt.show()

'''
x1d_fits = os.path.abspath('/Users/pena/Documents/AptanaStudio3/carbon/HSTdata/obqn18020clean05_16_x1d.fits')
img = pyfits.open(x1d_fits)
img.info()
data = img[1].data
w1 = data.WAVELENGTH.flatten()
f1 = data.FLUX.flatten()
pyplot.plot(w1, f1, 'r')
x1d_fits2 = os.path.abspath('/Users/pena/Documents/AptanaStudio3/carbon/HSTdata/obqn18020cleaniter3_16_x1d.fits')
img2 = pyfits.open(x1d_fits2)
img2.info()
data2 = img2[1].data
w2 = data2.WAVELENGTH.flatten()
f2 = data2.FLUX.flatten()
pyplot.plot(w2, f2, 'b')
#x1d_fits2 = os.path.abspath('/Users/pena/Documents/STIS/prop12472/obqn15040_sx1.fits')
x1d_fits2 = os.path.abspath('/Users/pena/Documents/AptanaStudio3/carbon/HSTdata/obqn18020_opt16_x1d.fits')
img2 = pyfits.open(x1d_fits2)
img2.info()
data2 = img2[1].data
w2 = data2.WAVELENGTH.flatten()
f2 = data2.FLUX.flatten()
pyplot.plot(w2, f2, 'g')
x1d_fits3 = os.path.abspath('/Users/pena/Documents/AptanaStudio3/carbon/HSTdata/obqn18020cleaniter10_16_x1d.fits')
img3 = pyfits.open(x1d_fits3)
img3.info()
data3 = img3[1].data
w3 = data3.WAVELENGTH.flatten()
f3 = data3.FLUX.flatten()
pyplot.plot(w3, f3, 'k')
pyplot.show()
'''


print '\n Code has finished!'

