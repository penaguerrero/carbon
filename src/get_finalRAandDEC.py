import glob
from astropy.io import fits

"""
This script gets the final RA and DEC for each object.
"""

# short name of the object
#                 0           1           2            3         4        5          6        7         8
objects_list =['iiizw107', 'iras08339', 'mrk1087', 'mrk1199', 'mrk5', 'mrk960', 'ngc1741', 'pox4', 'sbs0218',
               'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol1457', 'tol9', 'arp252', 'iras08208', 'sbs1415']
#                 9           10         11         12         13       14        15         16         17

#####################################################################################################################

#### SET PARAMETERS FOR SCRIPT
 
object_number = 'all'    # choice from 0 through 17 or 'all'
save_txt_file = False    # save the text file with the resulting RA and Dec: True or False 

#####################################################################################################################

#### FUNCTIONS

def deg2HMS(ra='', dec='', rounding=False):
    # function copied from: http://www.bdnyc.org/?p=584
    RA, DEC, rs, ds = '', '', '', ''
    if dec:
        if str(dec)[0] == '-':
            ds, dec = '-', abs(dec)
        deg = int(dec)
        decM = abs(int((dec-deg)*60))
        if rounding:
            decS = int((abs((dec-deg)*60)-decM)*60)
        else:
            decS = (abs((dec-deg)*60)-decM)*60
        DEC = '{0}{1} {2} {3}'.format(ds, deg, decM, decS)
    
    if ra:
        if str(ra)[0] == '-':
            rs, ra = '-', abs(ra)
        raH = int(ra/15)
        raM = int(((ra/15)-raH)*60)
        if rounding:
            raS = int(((((ra/15)-raH)*60)-raM)*60)
        else:
            raS = ((((ra/15)-raH)*60)-raM)*60
        RA = '{0}{1} {2} {3}'.format(rs, raH, raM, raS)
    
    if ra and dec:
        return (RA, DEC)
    else:
        return RA or DEC


#####################################################################################################################

#### CODE

# reference IDs for sample objects:
obqn_num_list = ['180', '080', '030', '060', '050', '010', '040', '140', '020', 
                 '110', '190', '130', '150', '170', '120', '100', '070', '160']

full_names_list = ['IIIZw 107', 'IRAS 08339+6517', 'MRK 1087', 'MRK 1199', 'MRK 5', 'MRK 960', 
                   'NGC 1741', 'POX 4', 'SBS 0218+003', 'SBS 0948+532', 'SBS 0926+606', 'SBS 1054+365', 
                   'SBS 1319+579', 'TOL 1457-262', 'TOL 9', 'ARP 252', 'IRAS 08208+2816', 'SBS 1415+437']

path_oneDspecs = '../HSTdata/'

RAs = []
DECs = []

for i, obj in enumerate(objects_list):
    if object_number != 'all':
        obj = objects_list[object_number]
        obqn_num = obqn_num_list[object_number]
        full_name = full_names_list[object_number]
    else:
        obqn_num = obqn_num_list[i]
        full_name = full_names_list[i]
    
    # get all the fits files of the same object into a list
    oneDspecs = glob.glob(path_oneDspecs+'obqn'+obqn_num+'*x*.fits')
    #print 'Object: ', obj
    for spec in oneDspecs:
        hdr = fits.getheader(spec, 0)
        RA_deg = fits.getval(spec, "RA_TARG", 0)
        DEC_deg = fits.getval(spec, "DEC_TARG", 0)
        RA = deg2HMS(ra=RA_deg)#, rounding=True)
        DEC = deg2HMS(dec=DEC_deg)#, rounding=True)
        #print 'RA  = ', RA_deg, '  --> ', RA
        #print 'DEC = ', DEC_deg, ' --> ', DEC
        #raw_input()
    RAs.append(RA)
    DECs.append(DEC)
    
    if object_number != 'all':
        print 'Object: ', full_name
        print '    RA  = ', RA_deg, ' --> ', RA
        print '    DEC = ', DEC_deg, ' --> ', DEC
        exit()


# Save the RAs and DECs into a text file or show it on screen
out_file = "../results/Final_RA_Dec_sample.txt"
line0 = "{:<18} {:<20} {:<16}".format("Galaxy name", "RA", "Dec")
print line0
if save_txt_file:
    f = open(out_file, "w+")
    f.write(line0+"\n")
for i, fn in enumerate(full_names_list):
    line1 = "{:<18} {:<20} {:<16}".format(fn, RAs[i], DECs[i])
    print line1
    if save_txt_file:
        f.write(line1+"\n")
if save_txt_file:
    f.close()
    print "\n RA and Dec saved in file: ", out_file

print "\n   Script finished!"

    