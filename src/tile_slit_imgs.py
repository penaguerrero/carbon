import os
import string
import matplotlib
from pylab import *
import PIL.Image as Image
from glob import glob


#####################################################################################################################

'''

This script crates a tile of the slit images of the objects in the sample.

Choose parameters to run script  
'''


# Specify the type of image to be saved
save_fig = True
img_format = '.jpg'


#####################################################################################################################

##### LISTS OF INPUT INFORMATION

# Number of plots per tile 
npt = 9

# name of the object in order of RA (same as paper)
#                  0          1           2         3         4        5           6             7         8
objects_list = ['mrk960', 'sbs0218', 'mrk1087', 'ngc1741', 'mrk5', 'mrk1199', 'iras08208', 'iras08339', 'sbs0926',
                'arp252', 'sbs0948', 'tol9',  'sbs1054', 'pox4', 'sbs1319', 'sbs1415', 'tol1457', 'iiizw107']
#                  9         10        11        12        13       14         15         16         17
number_of_tiles = len(objects_list)/npt

full_names_list = ['MRK 960', 'SBS 0218+003', 'MRK 1087', 'NGC 1741', 'MRK 5', 'MRK 1199', 
                   'IRAS 08208+2816', 'IRAS 08339+6517', 'SBS 0926+606', 'ARP 252', 'SBS 0948+532', 'TOL 9', 
                   'SBS 1054+365', 'POX 4', 'SBS 1319+579', 'SBS 1415+437', 'TOL 1457-262', 'IIIZw 107']

##### FUNCTIONS
def show_img_at_loc(data, obj, row, column, number_of_plot):
    img = Image.open(data)
    # These commands make the image bigger
    #dpi = rcParams['figure.dpi']
    #figsize = img.size[0]/dpi, img.size[1]/dpi
    #figs = figure(figsize=figsize)
    #ax = axes([0,0,1,1], frameon=False)
    #ax.set_axis_off()
    loc = str(row)+str(column)+str(number_of_plot)
    subplot(int(loc))
    axis('off')
    title(obj)
    imshow(img)

#### CODE

path_list = string.split(os.getcwd(), sep='carbon')
carbon_dir = 'carbon/results/plots/object_images/'
# Name of tile image
font = {#'family' : 'Vera Sans',
        'weight' : 'regular',
        'size'   : 14}
matplotlib.rc('font', **font)

# get list of images
list_of_imgs = glob(os.path.join(path_list[0], carbon_dir+'*_after.jpg'))
 
# Create the actual tile or tiles
img_number = 1
fig = figure(1, figsize=(15, 15))
fig.subplots_adjust(wspace=0.1, hspace=0.15)
# First plot images from object number 0 through 8
for obj in objects_list:
    for data in list_of_imgs:
        if obj in data:
            show_img_at_loc(data, full_names_list[img_number-1], 3, 3, img_number)
    img_number = img_number + 1
    if img_number > npt:
        break
plt.annotate('E', xy=(1.22, -0.2), xycoords='axes fraction', xytext=(0, -0.22), 
            arrowprops=dict(arrowstyle="<-", color='k'))
plt.annotate('N', xy=(1.21, -0.21), xycoords='axes fraction', xytext=(1.19, 1), 
            arrowprops=dict(arrowstyle="<-", color='k'))
if save_fig:
    path_from_carbon_dir = carbon_dir+'imgs_tile_1'+img_format
    destination = os.path.join(path_list[0], path_from_carbon_dir)
    savefig(os.path.abspath(destination))
    print 'Plot ', os.path.abspath(destination), 'saved!'
show()

# Now plot images from object 9 through 17
img_number = 1
fig = figure(1, figsize=(15, 15))
fig.subplots_adjust(wspace=0.1, hspace=0.15)
for i in range(npt, len(objects_list)):
    for data in list_of_imgs:
        if objects_list[i] in data:
            show_img_at_loc(data, full_names_list[i], 3, 3, img_number)
    img_number = img_number + 1
plt.annotate('E', xy=(1.22, -0.2), xycoords='axes fraction', xytext=(0, -0.22), 
            arrowprops=dict(arrowstyle="<-", color='k'))
plt.annotate('N', xy=(1.21, -0.21), xycoords='axes fraction', xytext=(1.19, 1), 
            arrowprops=dict(arrowstyle="<-", color='k'))
if save_fig:
    path_from_carbon_dir = carbon_dir+'imgs_tile_2'+img_format
    destination = os.path.join(path_list[0], path_from_carbon_dir)
    savefig(os.path.abspath(destination))
    print 'Plot ', os.path.abspath(destination), 'saved!'
show()

print 'Code finished!'
