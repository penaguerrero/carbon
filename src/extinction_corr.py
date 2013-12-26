import numpy as np
import matplotlib.pyplot as plt
import pyneb as pn 
import os
import copy
import science

############################################################################################################################################

'''  Choose parameters to run script  '''

# 1) Select a number from objects_list, i = :
objects_list =['arp252', 'iiizw107', 'iras08208', 'iras08339', 'mrk5', 'mrk960', 'mrk1087', 'mrk1199', 'ngc1741', 
               'pox4', 'sbs0218', 'sbs0948', 'sbs0926', 'sbs1054', 'sbs1319', 'tol9', 'tol1457']
#       arp252 = 0,  iiizw107 = 1,  iras08208 = 2,  iras08339 = 3,  mrk5 = 4,  mrk960 = 5, mrk1087 = 6,  mrk1199 = 7,  ngc1741 = 8,  
#       pox4 =9,  sbs0218 = 10,  sbs0948 = 11, sbs0926 = 12,  sbs1054 = 13,  sbs1319 = 14,  tol9 =15,  tol1457 = 16
object_number = 14
object_name = objects_list[object_number]

# 2) Do you want to create a unified lines text file?
create_txt = True

# Set theoretical Halpha/Hbeta ratio
I_theo_HaHb = 2.86 

# Set initial value of EWabsHbeta (this is a guessed value taken from HII regions)
# for HII region type objects typical values are 2.0-4.0 
EWabsHbeta = 0.2

# Set value for extinction
# for HII region type objects there is no restriction to max but values MUST be positive
C_Hbeta = 0.03

############################################################################################################################################

#### Read the observed lines from the table of lines_info.txt and normalize to Hbeta
### Path of the text files of wavelengths and fluxes for the objects. 
### This works as long as the folder structure is always the same. This file is assumed to be in
###           /Users/home_direcotry/Documents/AptanaStudio3/src/
### so I want to go one back to find the results folder
results_path = "../results/"
### but just to make sure we are in the right place, get and work with full paths
full_results_path = os.path.abspath(results_path)
text_files_path = os.path.join(full_results_path, "1Dspecs/")
### Go into the object folder
results4object_path = os.path.join(full_results_path, object_name)
### Use nuv, opt, and nir files
specs = [0, 1, 2]
add_str = "_lineinfo"
object_file = os.path.join(results4object_path, object_name)
text_file_list, _ = science.spectrum.get_obj_files2use(object_file, add_str, specs)

### Define the file name for for all the lines
name_out_file = os.path.join(results4object_path, object_name+"_linesNUV2NIR.txt")

cols_in_file, all_err_cont_fit = science.spectrum.gather_specs(text_file_list, name_out_file, reject=50.0, start_w=None, create_txt=create_txt)
catalog_wavelength, observed_wavelength, element, ion, forbidden, how_forbidden, width, flux, continuum, EW = cols_in_file

### Create an array of the numbers
data = np.array([catalog_wavelength, observed_wavelength, flux, continuum, EW])
# data: 0=catalog_wavelength, 1=observed_wavelength, 2=flux, 3=continuum, 4=equivalent_widths

### Step 1 of first iteration of reddening correction: Assume that there is no collisional excitation
### get the EW_abs of the H and He lines with respect to EW_abs(Hbeta)
Hline_and_EWs, Heline_and_EWs = science.spectrum.readlines_EWabsRelHbeta()
### Both H and He lines are going to be used to correct for underlying absorption
undabs_wav = []
undabs_EW = []
for i in range(len(Hline_and_EWs[0])):
    undabs_wav.append(Hline_and_EWs[0][i])
    undabs_EW.append(Hline_and_EWs[1][i])
for i in range(len(Heline_and_EWs[0])):
    undabs_wav.append(Heline_and_EWs[0][i])
    undabs_EW.append(Heline_and_EWs[1][i])
lines_undabs_and_EW = [undabs_wav, undabs_EW]
### Now add a 0.000 to all other observed lines
corr_undelyingAbs_EWs = []
for w in catalog_wavelength:
    if int(w) in undabs_wav:
        w_idx_undabs = undabs_wav.index(int(w))
        e = undabs_EW[w_idx_undabs]
    else:
        e = 0.000
    corr_undelyingAbs_EWs.append(e)

### Round all catalog lines to make it easier to find lines
rounded_catalog_wavelength = []
for item in catalog_wavelength:
    rw = np.round(item)
    rounded_catalog_wavelength.append(rw)
'''
### Calculate fluxes of 3726 and 3739   -- this works as long as the sum of the individual lines add up to better than 85% of the total 3727 measurement
idx3726 = rounded_catalog_wavelength.index(3726)
idx3729 = rounded_catalog_wavelength.index(3729)
idx3727 = rounded_catalog_wavelength.index(3727)
# assume that real3726 + real3729 = real3727 and that measured3726 + measured3729 = measured3727
# then we need a constant K that allows the fluxes of real3726 and real3729 to be comparable with measured3726 and measured3729
measured3727 = flux[idx3727]
measured3726 = flux[idx3726]
measured3729 = flux[idx3729]
K = measured3727 / (measured3726 + measured3729)
real3726 = K * measured3726
real3729 = K * measured3729
# insert these fluxes and new equivalent wids in the corresponding places
#print 'PREVIOUS flux of 3726 =', flux[idx3726], ' and 3729 =', flux[idx3729], '    sum =', flux[idx3726]+flux[idx3729]
#print '          EWs of 3726 =', EW[idx3726], '   and 3729 =', EW[idx3729], '      sum =', EW[idx3726]+EW[idx3729]
flux[idx3726] = real3726
flux[idx3729] = real3729
EW[idx3726] = real3726 / continuum[idx3726]
EW[idx3729] = real3729 / continuum[idx3729]
#print ' NEW     flux of 3726 =', flux[idx3726], ' and 3729 =', flux[idx3729], '    sum =', flux[idx3726]+flux[idx3729]
#print '          EWs of 3726 =', EW[idx3726], '   and 3729 =', EW[idx3729], '      sum =', EW[idx3726]+EW[idx3729]
'''
### Remove UNDERLYING ABSORPTION for optical lines to get Intensities
intensities = []
for EWabsLine,cont,flx in zip(corr_undelyingAbs_EWs, continuum, flux):
    I = EWabsLine * EWabsHbeta * cont + flx
    intensities.append(I)

### Step 2 of first iteration of reddening correction: Using Seaton
Is_corr = []
cHbeta = 0.434*C_Hbeta
### f_lambda is the reddening law, in this case I got it from Seaton 1979 based on Tol 2146-319
f_lambda = science.spectrum.readreddCorr(rounded_catalog_wavelength)
for I, fl in zip(intensities,f_lambda):
    I_c = I * 10**(cHbeta*(1+fl)) 
    Is_corr.append(I_c)
    #print I, fl, I_c

### Step 3: Normalize observed line ratios to Hbeta
Hb_idx = rounded_catalog_wavelength.index(4861.0)
Hbeta = Is_corr[Hb_idx]
#print 'Corrected intensity of Hbeta = ', Hbeta
norm_Icor = []
for I, w in zip(Is_corr, observed_wavelength):
    norm_f = (I / Hbeta) * 100.00
    norm_Icor.append(norm_f)
    #print w, I, Hbeta, norm_f

### ERROR Calculation
### STIS Data Handbook states this depends on many factors so I set it up to 2% because
### they suggest 1% to increase considerably for faint objects.
abs_flux_calibration_err = 2.0
### Find the faintest detected emission line: get rid of negative fluxes
catalog_emission_lines = []
element_emission_lines = []
ion_emission_lines = []
forbidden_emission_lines = []
how_forbidden_emission_lines = []
wavs_emission_lines = []
final_f_lambda = []
positive_normfluxes = []
positive_norm_intensities = []
positive_norm_Icorr = []
EW_emission_lines = []
pos_calc_cont = []
for i in range(len(norm_Icor)):
    if flux[i] > 0.00000:
        catalog_emission_lines.append(rounded_catalog_wavelength[i])
        element_emission_lines.append(element[i])
        ion_emission_lines.append(ion[i])
        forbidden_emission_lines.append(forbidden[i])
        how_forbidden_emission_lines.append(how_forbidden[i])
        wavs_emission_lines.append(observed_wavelength[i])
        final_f_lambda.append(f_lambda[i])
        norm_flux = flux[i] / flux[Hb_idx] * 100.
        positive_normfluxes.append(norm_flux)
        normI = intensities[i] / intensities[Hb_idx] * 100.
        positive_norm_intensities.append(normI)
        positive_norm_Icorr.append(norm_Icor[i])
        EW_emission_lines.append(EW[i])
        pos_calc_cont.append(continuum[i])
print '  ***  There are ', len(positive_norm_intensities), ' emission lines in this object!'

### I am defining the chosen fainest line with a S/N=3 or error=33%
faintest_line = min(positive_normfluxes)
print 'first faintest_line = ', faintest_line
kk = copy.deepcopy(positive_normfluxes)
min_counts = 1
end_loop = False
while end_loop != True:
    if min_counts == 5:
        end_loop = True
    else:
        for f in kk:
            if f == faintest_line:
                min_counts = min_counts + 1 
                kk.pop(kk.index(faintest_line))
                faintest_line = min(kk)
                print 'new faintest_line = ', faintest_line
chosen_faintest_line = faintest_line
print 'my definition of S/N~3 after min_counts =', min_counts
print 'chosen_faintest_line =', chosen_faintest_line, '+- 33%'

### In order to account for the error in the continuum fitting, I realized that the 
### polynomial order does not really change the intensity, however, the window width does!
### To estimate the error on the polynomial fitting I added 1/sum(err), 
### were err=1/sqroot(points in window width)
#print 'all_err_cont_fit', all_err_cont_fit
### Add errors and present them in percentage of line intensities
percent_err_I = []
absolute_err_I = []
emission_lines_file = os.path.join(results4object_path, object_name+'_emlines.txt')
emfile = open(emission_lines_file, 'w+')
print >> emfile, '# Positive EW = emission        Negative EW = absorption' 
print >> emfile, '# C_Hbeta = %0.3f' % C_Hbeta
print >> emfile, '# EWabs_Hbeta = %0.3f' % EWabsHbeta
print >> emfile, '# I_theo_HaHb = %0.3f' % I_theo_HaHb
print >> emfile, ('# {:<13} {:<8} {:>6} {:<12} {:<12} {:<16} {:<18} {:<18} {:<16} {:<16} {:<16}'.format('Wavelength', 'Element', 'Ion', 'Forbidden', 'How much', 'Flux [cgs]', 'Intensity [cgs]', '% Err', 'Abs err', 'EW [A]', 'S/N'))
for w, el, i, f, h, F_norm, I_norm, eqw in zip(catalog_emission_lines, element_emission_lines, ion_emission_lines, forbidden_emission_lines, 
                             how_forbidden_emission_lines,positive_normfluxes, positive_norm_intensities, EW_emission_lines):
    if w <= 2000.0:
        e = all_err_cont_fit[0]
    elif w > 2000.0 and w < 5000.0:
        e = all_err_cont_fit[1]
    elif w >= 5000.0:
        e = all_err_cont_fit[2]
    per_err_I = (e*e + 10000*(chosen_faintest_line/9)/(F_norm) + abs_flux_calibration_err*abs_flux_calibration_err)**0.5
    #per_err_I = (e*e + 10000*(chosen_faintest_line/1)/F_norm + abs_flux_calibration_err*abs_flux_calibration_err)**0.5
    #per_err_I = (100*(1/(e*e*F_norm)) + abs_flux_calibration_err*abs_flux_calibration_err)**0.5
    abs_err = (per_err_I * F_norm) / 100.
    sn = F_norm/per_err_I*100.
    print w, '   F_norm = %0.2f   I_norm = %0.2f   err_porcent = %0.2f   abs_err = %0.2f' % (F_norm, I_norm, per_err_I, abs_err), '   S/N = ', sn
    print >> emfile, ('{:<15.0f} {:<8} {:>6} {:<12} {:<12} {:<16.3f} {:<18.3f} {:<18.3f} {:<16.3f} {:<16.3f} {:<16.3f}'.format(w, el, i, f, h, F_norm, I_norm, per_err_I, abs_err, eqw, sn))
    percent_err_I.append(per_err_I)
    absolute_err_I.append(abs_err)
emfile.close()


##### FROM THIS POINT, THIS IS FROM THE ORIGINAL PYNEB SCRIPT
# Convert wavelength to x
def x(wave):
    return 10000. / wave

# Define an extinction law (to be used below)
def my_X(wave):
    x = 10000. / wave
    Rv = 3.1
    X_lin = x/2. # linear part of the extinction law
    X_bump = 0.5*x**2. -6*x + 20. # bump part of the extinction law
    return Rv*np.where(x<5., X_lin, X_bump)

# Define a reddening correction object
RC = pn.RedCorr()

# List the available laws
#RC.printLaws()

# Plot the available laws
#RC.plot(laws='all')
#plt.show()

# Choose the one we intend to use 
#RC.law = 'S 79 H 83'
RC.law = 'CCM 89'
# or define a new one
#RC.UserFunction = my_X
#RC.law = 'user'

'''
# Plot the selected law as a function of x
# Define an array in lambda to do the plot
wave= np.logspace(2.5, 5, 100)
# Plot commands
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylim([0, 15])
ax.set_xlim([0, 10])
ax.plot(x(wave), my_X(wave), label='%s' % (RC.law))
plt.xlabel('1/$\lambda$ ($\mu^{-1}$)')
plt.ylabel('A$_V$/E(B-V)')
plt.legend(loc='upper left')
plt.show()
'''   
##### UP TO HERE: THIS IS FROM THE ORIGINAL PYNEB SCRIPT

### For my intensities and using Seaton's law
pyneb_Idered = []
I_dered_norCorUndAbs = []
for w, nF, nI, Ic, e  in zip(catalog_emission_lines, positive_normfluxes, positive_norm_intensities, positive_norm_Icorr, EW_emission_lines):    
    # Correct based on the given law and c(Hb)
    RC = pn.RedCorr(law= 'CCM 89', cHbeta=cHbeta)
    I_dered = nI * RC.getCorrHb(w)
    pyneb_Idered.append(I_dered)
    #print '\nCorrect based on the given law and the observed Ha/Hb ratio:'
    #print str(w) + ': F_obs =', nF, 'I_obs =', nI, ' I_dered =', I_dered, '  my_I_dered =', Ic, '  EW =', e
    IdnUA = nF * RC.getCorrHb(w)
    I_dered_norCorUndAbs.append(IdnUA)
    
### Find observed Halpha/Hbeta ratio
## 
Halpha_idx = catalog_emission_lines.index(6563)
Halpha = positive_norm_Icorr[Halpha_idx]
pynebHalpha = pyneb_Idered[Halpha_idx]
Hbeta_idx = catalog_emission_lines.index(4861)
Hbeta = positive_norm_Icorr[Hbeta_idx]
pynebHbeta = pyneb_Idered[Hbeta_idx]
I_obs_HaHb = Halpha/Hbeta
pynebI_obs_HaHb = pynebHalpha/pynebHbeta
print ''
print 'cHbeta = %0.5f' % cHbeta
print '            Using Seaton                    Using ', RC.law
print 'Corrected for reddening and underlying absorption'
print catalog_emission_lines[Halpha_idx], '    ', Halpha, '                  ', pynebHalpha
print catalog_emission_lines[Hbeta_idx], '    ', Hbeta, '                          ', pynebHbeta
print 'theoretical ratio Ha/Hb = %0.2f' % (I_theo_HaHb)
print '      observed Ha/Hb = %0.2f           observed Ha/Hb = %0.2f' % (I_obs_HaHb, pynebI_obs_HaHb)
w1 = 6563
w2 = 4861
raw_w1 = flux[rounded_catalog_wavelength.index(w1)]
raw_w2 = flux[rounded_catalog_wavelength.index(w2)]
raw_contw1 = continuum[rounded_catalog_wavelength.index(w1)]
raw_contw2 = continuum[rounded_catalog_wavelength.index(w2)]
print 'Before extinction correction'
print 'FLUXES:    %i = %e    %i = %e    ratio = %0.2f' % (w1, raw_w1, w2, raw_w2, raw_w1/raw_w2)
print 'continuum: %i = %e    %i = %e' % (w1, raw_contw1, w2, raw_contw2)

# Finding the f_lambda values:
pyneb_flambda = []
for w,Icor,Iobs in zip(catalog_emission_lines, I_dered_norCorUndAbs, positive_normfluxes):#pyneb_Idered, positive_norm_intensities):
    f12 = (np.log10(Icor) - np.log10(Iobs)) / cHbeta
    #print w, Iobs, Icor, f12
    pyneb_flambda.append(f12)
    
# Write the first round of reddeding correction in pyneb readable format
tfile1stRedCor = os.path.join(results4object_path, object_name+"_1stRedCor.txt")
tf = open(tfile1stRedCor, 'w+')
#print >> tf,  ('{:<15} {:>15} {:>15}'.format('Ion_line', 'Intensity', 'Abs Error'))
#print >> tf, 'cHbeta  %0.3f' % cHbeta
lines_pyneb_matches = []
for cw, el, io, Ic, er in zip(catalog_emission_lines, element_emission_lines, ion_emission_lines, positive_norm_Icorr, absolute_err_I):
    cw = int(cw)
    cw = str(cw)
    el = str(el)
    io = str(io)
    wavid = cw+'A'
    pynebid = el+io
    #print 'pynebid =', pynebid, '    wavid =', wavid
    lineID = pynebid + '_' + wavid
    matching_line = [cw, Ic, er]
    if pynebid in pn.LINE_LABEL_LIST:
        if wavid in pn.LINE_LABEL_LIST[pynebid]:
            print >> tf,  ('{:<15} {:>15.3f} {:>15.3f}'.format(lineID, Ic, er))
            lines_pyneb_matches.append(matching_line)
        else:
            pynebid = el+io+'_'+wavid+'+'
            if pynebid in pn.BLEND_LIST:
                print >> tf,  ('{:<15} {:>15.3f} {:>15.3f}'.format(pynebid, Ic, er))
                lines_pyneb_matches.append(matching_line)
            else:
                continue
tf.close()
print 'File   %s   writen!' % tfile1stRedCor

# Determine first approximation of temperatures and densities
# OBSERVATIONS
# Define an Observation object and assign it to name 'obs'
obs = pn.Observation()
# read data from file created specifically for pyneb reading
obs.readData(tfile1stRedCor, fileFormat='lines_in_rows', corrected=True, errIsRelative=False)
# Intensities
#print 'len(lines_pyneb_matches)', len(lines_pyneb_matches)
for line in obs.lines:
    iline = line.corrIntens
    for matching_line in lines_pyneb_matches:
        if matching_line[0] in line.label:
            #print 'found it!', matching_line[0]
            #print line.wave
            #print line.corrIntens
            #iline = line.corrIntens
            if line.wave == 4363:
                I1 = line.corrIntens
                print '4363 has an intensity of', I1[0]
            elif line.wave == 5007:
                I2 = line.corrIntens
                print '5007 has an intensity of', I2[0]
            elif line.wave == 3726:
                IO21 = line.corrIntens
                print '3726 has an intensity of', IO21[0]
            elif line.wave == 3729:
                IO22 = line.corrIntens
                print '3729 has an intensity of', IO22[0]
            elif line.wave == 6312:
                IS31 = line.corrIntens
                print '6312 has an intensity of', IS31[0]
            elif line.wave == 9069:
                IS32 = line.corrIntens
                print '9069 has an intensity of', IS32[0]
            elif line.wave == 9531:
                IS33 = line.corrIntens
                print '9531 has an intensity of', IS33[0]
            elif line.wave == 5518:
                ICl31 = line.corrIntens
                print '5518 has an intensity of', ICl31[0]
            elif line.wave == 5538:
                ICl32 = line.corrIntens
                print '5538 has an intensity of', ICl32[0]


# Define all atoms to make calculations
all_atoms = pn.getAtomDict()
 
# simultaneously compute temperature and density from pairs of line ratios
# First of all, a Diagnostics object must be created and initialized with the relevant diagnostics.
diags = pn.Diagnostics()   # Instantiate the Diagnostics class
diags.getAllDiags()  # see what Diagnostics exist
# temperature determination from an intensity ratio
# explore some specific atom in the atoms collection
O3 = pn.Atom("O", "3")
O3ratio = I1[0] / I2[0]
print 'ratio of O3 = ', O3ratio
TO3 = O3.getTemDen(O3ratio, den=100., wave1=4363, wave2=5007)
print 'guess 1 of temp of O3 = ', TO3 
### Get temperature of OII from OIII. 
# Using equation of Peimbert, Peimbert, & Luridiana. 2002, ApJ, 565, 668 (Based from data of Stasinska's models.)
TO2pei = 2430. + TO3 * (1.031 - TO3/54350.)
print 'This is the theoretically obtained temperature of O2 from Peimbert etal 2002 = ', TO2pei
# Using equation of Garnett, D. R. 1992, AJ, 103, 1330
TO2gar = 0.7 * TO3 + 3000.
print 'Theoretically obtained temperature of O2 from Garnet 1992 = ', TO2gar
O2 = pn.Atom("O", "2")
O2ratio = IO22[0] / IO21[0]
print 'ratio of O2 = ', O2ratio
denO2 = O2.getTemDen(O2ratio, den=10000., wave1=3729, wave2=3726) 
print 'guess 1 of density of O2 = ', denO2
S3 = pn.Atom("S", "3")
S3ratio = IS31[0] / (IS32[0]) 
print 'ratio of S3 = ', S3ratio
TS3 = S3.getTemDen(S3ratio, den=100., wave1=6312, wave2=9532)
print 'guess 1 of temp of S3 = ', TS3 
#Cl3 = pn.Atom("Cl", "3")
#Cl3ratio = ICl32[0] / (ICl31[0]) 
#dCl3 = Cl3.getTemDen(S3ratio, temp=11000., wave1=5538, wave2=5518)
#print 'guess 1 of density of Cl3 = ', dCl3 

# Simultaneously determine temps and densities
try:
    tem_Ar3, den_S2 = diags.getCrossTemDen('[ArIII] 5192/7300+', '[SII] 6731/6716', obs=obs)
    tem_O3, den_Cl3 = diags.getCrossTemDen('[OIII] 4363/5007', '[ClIII] 5538/5518', obs=obs)
    tem_O3, den_Ar4 = diags.getCrossTemDen('[OIII] 4363/5007', '[ArIV] 4740/4711', obs=obs)
    tem_O3, den_O2 = diags.getCrossTemDen('[OIII] 4363/5007', '[OII] 3926/3929', obs=obs)
except:
    tem_Ar3 = 'NA'
    den_S2 = 'NA'
    tem_O3 = 'NA'
    den_O2 = 'NA'
    den_Ar4 = 'NA'
    den_Cl3 = 'NA'
    
# Printout of physical conditions
print 'den_O2: ', den_O2
print 'den_S2: ', den_S2
print 'tem_O3: ', tem_O3
print 'den_Ar4: ', den_Ar4

'''
# Alternate way of computing T(OIII)
if tem_O3 == 'NA':
    tem_O3 = all_atoms['O3'].getTemDen(i5007/i4363, den=100., wave1=5007, wave2=4363)
# Include in diags the relevant line ratios
diags.addDiag([
                ## temperatures
                #'[NII] 5755/6548',
                '[OII] 7320/3737+',
                '[OIII] 4363/5007',
                '[ArIII] 5192/7136',
                '[ArIII] 5192/7300+',
                '[ArIV] 7230+/4720+',
                #'[SIII] 6312/9532',
                ## densities
                #'[NI] 5198/5200',
                #'[OII] 3729/3736',
                '[ArIV] 4740/4711',
                #'[SII] 4072+/6725+',
                '[SII] 6731/6716',
                ])
'''



