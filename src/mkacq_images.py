import os
import pyfits
import pylab as plt
from PIL import Image
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
object_number = 1
vmin = 18
vmax = 65

# want to save images?  (type 'y' for yes or 'n' for no)
save_plt = 'n'
# useful cmap options: Blues, Greens, Greys, hot, Purples, binary, rainbow, for default type None
cmap = None#'Greys'

# do you want to change the format of the image?
change_format = True
new_filetype = 'eps'

# Do you want to see the plots separate or together? (choose separate=True for creating the jpg images)
separate = False

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

#fits_path = '/Users/pena/Documents/STIS/prop12472/raw_data'
fits_path = '/Users/pena/Documents/AptanaStudio3/carbon/HSTdata/raw_acq_imgs'
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

jpg_img_name = object_name+'_before'+'.jpg'
destination = os.path.join(jpgs_path, jpg_img_name)
jpg_img_name_after = object_name+'_after'+'.jpg'
destination_after = os.path.join(jpgs_path, jpg_img_name_after)
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
    plt.title(object_name+'_after')
    plt.colorbar(im2, orientation='vertical')
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
    ax2.set_title(object_name+'_after')
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

