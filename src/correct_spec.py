import pyfits
import os

'''
This script creates the text files of wavelength and flux in case there are some corrections to make by hand to the spectra.
        * these text files can later be runned by the task rspectext in order to create a fits file that splot can read
        AFTER RUNNING THIS SCRIPT:
        1. in splot correct the fits file
        2. reconvert the corrected fits file to text with wspectext
'''

obj = 'sbs1415'
img = 'obqn16012_nuv16_512_x1d'
part = 'nuv'

spec_path = os.path.abspath('../HSTdata/')
spec = os.path.join(spec_path,img+".fits")

t1 = pyfits.open(spec)
t1.info()
data1 = t1[1].data
wav = data1.WAVELENGTH.flatten()
flx = data1.FLUX.flatten()

results_path = os.path.abspath('../results')
out_path = os.path.join(results_path, obj)
outfile = os.path.join(out_path, obj+"_"+part+".txt")
fout = open(outfile, 'w+')
for w, f in zip(wav, flx):
    fout.write('{:<10.16f} {:26.10e}\n'.format(w, float(f)))
fout.close()
print 'Printed wavelengths and fluxes in:', outfile
print ''
print 'Code finished!'