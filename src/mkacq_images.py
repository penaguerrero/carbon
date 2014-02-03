import os
import pyfits
import pylab as plt
from PIL import Image

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
object_number = 0#8#11
vmin = 18
vmax = 65

# want to save images?  (type 'y' for yes or 'n' for no)
save_plt = 'n'
# useful cmap options: Blues, Greens, Greys, hot, Purples, binary, rainbow, for default type None
cmap = 'Greys'

# do you want to change the format of the image?
change_format = False
new_filetype = 'eps'

# Do you want to see the plots separate or together? (choose separate=True for creating the jpg images)
separate = False

############################################################################################################################################

object_name = objects_list[object_number]
img_list = ['obqn18wbq_raw.fits', 'obqn08a6q_raw.fits', 'obqn03ogq_raw.fits', 'obqn06wqq_raw.fits', 'obqn05pcq_raw.fits',
            'obqn01aaq_raw.fits', 'obqn04a5q_raw.fits', 'obqn14x5q_raw.fits', 'obqn02knq_raw.fits', 'obqn11g9q_raw.fits',
            'obqn19dzq_raw.fits', 'obqn13ufq_raw.fits', 'obqn15r2q_raw.fits', 'obqn17lyq_raw.fits', 'obqn12q8q_raw.fits',
            'obqn10c4q_raw.fits', 'obqn07pcq_raw.fits', 'obqn16plq_raw.fits']
img_name = img_list[object_number]

fits_path = '/Users/pena/Documents/STIS/prop12472/raw_data'
jpgs_path = '../results/plots/object_images'

fits_img = os.path.join(fits_path, img_name)
fits_raw = pyfits.open(fits_img)
fits_raw.info()
img = fits_raw[1].data      # Intensity data
img_recenter = fits_raw[4].data      # Intensity data

jpg_img_name = object_name+'_before'+'.jpg'
destination = os.path.join(jpgs_path, jpg_img_name)
jpg_img_name_after = object_name+'_after'+'.jpg'
destination_after = os.path.join(jpgs_path, jpg_img_name_after)
if separate:
    # before
    im = plt.imshow(img, vmin=vmin, vmax=vmax, origin='lower')#, cmap=cmap)
    plt.title(object_name+'_before')
    plt.colorbar(im, orientation='vertical')
    if save_plt == 'y':
        plt.savefig(destination)
        print('Plot %s was saved!' % destination)
    elif save_plt == 'n':
        print('Plot not saved.')
    plt.show()
    # after
    im = plt.imshow(img_recenter, vmin=vmin, vmax=vmax, origin='lower')#, cmap='Greys')
    plt.title(object_name+'_after')
    plt.colorbar(im, orientation='vertical')
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
    im1 = ax1.imshow(img, vmin=vmin, vmax=vmax, origin='lower', cmap=cmap)
    ax1.set_title(object_name+'_before')
    #im1.colorbar(im1, orientation='vertical')
    # After
    im2 = ax2.imshow(img_recenter, vmin=vmin, vmax=vmax, origin='lower', cmap=cmap)
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

