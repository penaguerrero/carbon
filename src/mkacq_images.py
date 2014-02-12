import os
import pyfits
import pylab as plt
import numpy
import PIL.Image as Image
#from PIL import Image
from scipy import ndimage

'''
This code gets the acquisition raw fits file for the given object and extracts the science part to save it into 
2 jpg images: 1 before and 1 after the recentering of the telescope.

This script can also convert the jpg images into eps files if necessary.
'''

############################################################################################################################################

# name of the object
#                 0           1           2            3         4        5          6        7        8
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                 9           10          11       12          13       14       15         16         17

''' Choose parameters to run script'''
object_number = 11
vmin = 10
vmax = 60
#vmin = 10
#vmax = 90

# want to save images?  (type 'y' for yes or 'n' for no)
save_plt = 'n'
# useful cmap options: Blues, Greens, Greys, hot, Purples, binary, rainbow, for default type None
cmap = None#'Greys'

# do you want to change the format of the image?
change_format = True
new_filetype = 'eps'

# Do you want to see the plots separate or together? (choose separate=True for creating the jpg images)
separate = True

# Do you want to see the rotation angle of the slit in the 2d spectra?
rotate = False

############################################################################################################################################

object_name = objects_list[object_number]
#                     0                       1                    2                    3                    4 
img_list = ['obqn18wbq_raw.fits', 'obqn08a6q_raw.fits', 'obqn03ogq_raw.fits', 'obqn06wqq_raw.fits', 'obqn05pcq_raw.fits',
            #         5                      6                     7                    8                    9
            'obqn01aaq_raw.fits', 'obqn04a5q_raw.fits', 'obqn14x5q_raw.fits', 'obqn02knq_raw.fits', 'obqn11g9q_raw.fits',
            #        10                     11                    12                   13                   14
            'obqn19dzq_raw.fits', 'obqn13ufq_raw.fits', 'obqn15r2q_raw.fits', 'obqn17lyq_raw.fits', 'obqn12q8q_raw.fits',
            #        15                     16                    17 
            'obqn10c4q_raw.fits', 'obqn07pcq_raw.fits', 'obqn16plq_raw.fits']
img_name = img_list[object_number]

full_names_list = ['IIIZw 107', 'IRAS 08339+6517', 'MRK 1087', 'MRK 1199', 'MRK 5', 'MRK 960', 'NGC 1741', 'POX 4', 'SBS 0218+003',
               'SBS 0948+532', 'SBS 0926+606', 'SBS 1054+365', 'SBS 1319+579', 'TOL 1457-262', 'TOL 9', 'ARP 252', 'IRAS 08208+2816',
               'SBS 1415+437']
full_name = full_names_list[object_number]

#fits_path = '/Users/pena/Documents/STIS/prop12472/raw_data'
fits_path = '../HSTdata/raw_acq_imgs'
jpgs_path = '../results/plots/object_images'

fits_img = os.path.join(fits_path, img_name)
fits_raw = pyfits.open(fits_img)
fits_raw.info()
img = fits_raw[1].data      # Intensity data
img_recenter = fits_raw[4].data      # Intensity data
h1 = fits_raw[1].header
PA_APER = h1['PA_APER']
ORIENT = h1['ORIENTAT']
PA = ORIENT-45.  # the 45.35 came from the STIS data Handbook for a 52x0.2 slit
print 'PA = ', PA


def add_slit_line(lolim, uplim, teta):
    if (teta <= 90.0) and (teta >= -90.0):
        alfa = -1*(teta) * numpy.pi / 180.0
        print 'quadrant 1 or 2'
    elif (teta > 90.0) and (teta < 180.0):
        alfa = (180.0 - teta) * numpy.pi / 180.0
        print 'quadrant 3'
    elif (teta > 180.0) and (teta < 270.0):
        alfa = (-1)*(180.0 - teta) * numpy.pi / 180.0
        print 'quadrant 4'
    elif (teta < -90.0) and (teta > -180.):
        alfa = (180.0 - teta) * numpy.pi / 180.0
        print 'quadrant 4'
    elif (teta < -180.0) and (teta > -270.):
        alfa = ((-1)*teta - 180.0) * numpy.pi / 180.0
        print 'quadrant 3'
    x = [uplim/2, uplim/2]
    y = [lolim, uplim]
    ca = uplim/2.0
    co = (ca / numpy.cos(alfa)) * numpy.sin(alfa)
    z0 = x[0] - co
    z1 = x[0] + co
    z = [z0, z1]
    # draw plane
    #plt.plot(x, y, 'k-', lw=2)
    #plt.plot(y, x, 'k-', lw=2)
    # now draw slit
    plt.plot(z, y, 'w-', lw=2)
    plt.xlim(lolim, uplim)
    plt.ylim(lolim, uplim)


jpg_img_name = object_name+'_before'+'.jpg'
destination = os.path.join(jpgs_path, jpg_img_name)
jpg_img_name_after = object_name+'_after'+'.jpg'
destination_after = os.path.join(jpgs_path, jpg_img_name_after)
lolim = 0
uplim = 95    
pangle = round(PA, 2)
if separate:
    # before
    f1 = plt.figure(1, figsize=(10, 10))
    if rotate == True:
        rotate_im = ndimage.rotate(img, PA)#, reshape=False)
        im = plt.imshow(rotate_im, vmin=vmin, vmax=vmax, origin='lower', cmap=cmap)
    else:
        im = plt.imshow(img, vmin=vmin, vmax=vmax, origin='lower', cmap=cmap)
    plt.title(object_name+'_before')
    plt.colorbar(im, orientation='vertical')
    if save_plt == 'y':
        plt.savefig(destination)
        print('Plot %s was saved!' % destination)
    elif save_plt == 'n':
        print('Plot not saved.')
    plt.show()
    # after
    f2 = plt.figure(1, figsize=(10, 10))
    if rotate == True:
        rotate_im2 = ndimage.rotate(img_recenter, PA)#, reshape=False)
        im2 = plt.imshow(rotate_im2, vmin=vmin, vmax=vmax, origin='lower', cmap=cmap)
    else:
        im2 = plt.imshow(img_recenter, vmin=vmin, vmax=vmax, origin='lower', cmap=cmap)
    #plt.title(object_name+'_after')
    #plt.title(object_name+'_PA'+repr(pangle))
    color = 'w'
    size = 'x-large'
    plt.text(60, 10, full_name, size=size, color=color)
    plt.text(60, 6, 'P.A. = '+repr(pangle), size=size, color=color)
    add_slit_line(lolim, uplim, PA)
    cbar_ax = f2.add_axes([0.92, 0.11, 0.03, 0.776])
    plt.colorbar(im2, cax=cbar_ax)#, orientation='vertical')
    if save_plt == 'y':
        plt.savefig(destination_after)
        print('Plot %s was saved!' % destination_after)
    elif save_plt == 'n':
        print('Plot not saved.')
    plt.show()
else:
    # plot the images before and after recentering
    # Before
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(20,10))
    if rotate == True:
        rotate_im = ndimage.rotate(img, PA)#, reshape=False)
        im1 = ax1.imshow(rotate_im, vmin=vmin, vmax=vmax, origin='lower', cmap=cmap)
        rotate_im2 = ndimage.rotate(img_recenter, PA)#, reshape=False)
        im2 = ax2.imshow(rotate_im2, vmin=vmin, vmax=vmax, origin='lower', cmap=cmap)
    else:
        im1 = ax1.imshow(img, vmin=vmin, vmax=vmax, origin='lower', cmap=cmap)
        im2 = ax2.imshow(img_recenter, vmin=vmin, vmax=vmax, origin='lower', cmap=cmap)
    ax1.set_title(object_name+'_before')
    # After
    add_slit_line(lolim, uplim, PA)
    ax2.set_title(object_name+'_after_recenter'+'_PA'+repr(pangle))
    cbar_ax = f.add_axes([0.92, 0.15, 0.02, 0.7])
    f.colorbar(im1, cax=cbar_ax)
    plt.show()


# Change the format
if change_format:
    # get the after jpg image and convert it to eps
    original = Image.open(destination_after)
    #original.show()
    destination_path = os.path.join(jpgs_path, object_name)
    new_img = destination_path+'.'+new_filetype
    original.save(new_img)
    print('Plot %s was saved!' % new_img)

print '\n Code has finished!'

