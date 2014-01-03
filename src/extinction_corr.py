import numpy as np
#import matplotlib.pyplot as plt
import pyneb as pn 
import os
import copy
import science
#import uncertainties

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

# Choose case
case = 'B'

# Set theoretical Halpha/Hbeta ratio
I_theo_HaHb = 2.86 

# Set initial value of EWabsHbeta (this is a guessed value taken from HII regions)
# for HII region type objects typical values are 2.0-4.0 
EWabsHbeta = 0.2

# Set value for extinction
# for HII region type objects there is no restriction to max but values MUST be positive
C_Hbeta = 0.024

# Do you want to do the second iteration for reddening correction (that is collisional excitation)?
reddeningCorrection2 = True

############################################################################################################################################

### FUNCTIONS
def underlyingAbsCorr(EWabsHbeta, corr_undelyingAbs_EWs, continuum, flux):
    intensities = []
    for EWabsLine,cont,flx in zip(corr_undelyingAbs_EWs, continuum, flux):
        I = EWabsLine * EWabsHbeta * cont + flx
        intensities.append(I)
    return intensities

def find_emission_lines(rounded_catalog_wavelength, element, ion, forbidden,
                        how_forbidden, observed_wavelength, flux, intensities, EW, continuum):
    Hb_idx = rounded_catalog_wavelength.index(4861.0)
    catalog_emission_lines = []
    element_emission_lines = []
    ion_emission_lines = []
    forbidden_emission_lines = []
    how_forbidden_emission_lines = []
    wavs_emission_lines = []
    positive_normfluxes = []
    positive_norm_intensities = []
    EW_emission_lines = []
    pos_calc_cont = []
    for i in range(len(flux)):
        if flux[i] > 0.00000:
            catalog_emission_lines.append(rounded_catalog_wavelength[i])
            element_emission_lines.append(element[i])
            ion_emission_lines.append(ion[i])
            forbidden_emission_lines.append(forbidden[i])
            how_forbidden_emission_lines.append(how_forbidden[i])
            wavs_emission_lines.append(observed_wavelength[i])
            norm_flux = flux[i] / flux[Hb_idx] * 100.
            positive_normfluxes.append(norm_flux)
            normI = intensities[i] / intensities[Hb_idx] * 100.
            positive_norm_intensities.append(normI)
            EW_emission_lines.append(EW[i])
            pos_calc_cont.append(continuum[i])
    #print '  ***  There are ', len(positive_norm_intensities), ' emission lines in this object!'
    return (catalog_emission_lines, wavs_emission_lines, element_emission_lines, ion_emission_lines, forbidden_emission_lines, 
            how_forbidden_emission_lines, positive_normfluxes, pos_calc_cont, positive_norm_intensities, EW_emission_lines)
    
def Halpha2Hbeta_dered(I_theo_HaHb, cHbeta, catalog_emission_lines, positive_normfluxes, positive_norm_intensities):
    ''' Function to dered and obtain the Halpha/Hbeta ratio nicely printed along with the unreddend values to compare. '''
    # Normalized fluxes are the raw flux normalized to Hbeta. 
    # Intensities are already corrected for underlying absorption and are also normalized to Hbeta.
    Idered = []
    I_dered_norCorUndAbs = []
    for w, nF, nI in zip(catalog_emission_lines, positive_normfluxes, positive_norm_intensities):    
        # Obtain the reddening corrected intensities based on the given law and c(Hb)
        RC = pn.RedCorr(law= 'CCM 89', cHbeta=cHbeta)
        I_dered = nI * RC.getCorrHb(w)
        Idered.append(I_dered)
        # Obtain the reddening corrected intensities WITHOUT the correction due to 
        # underlying absorption based on the given law and c(Hb)
        IdnUA = nF * RC.getCorrHb(w)
        I_dered_norCorUndAbs.append(IdnUA)
    # Find observed Halpha/Hbeta ratio
    Halpha_idx = catalog_emission_lines.index(6563)
    Hbeta_idx = catalog_emission_lines.index(4861)
    Halpha = positive_normfluxes[Halpha_idx]
    Hbeta = positive_normfluxes[Hbeta_idx]
    raw_ratio = Halpha/Hbeta
    I_Halpha = Idered[Halpha_idx]
    I_Hbeta = Idered[Hbeta_idx]
    I_obs_HaHb = I_Halpha/I_Hbeta
    print ''
    print 'cHbeta = %0.5f' % cHbeta
    print '            Using', RC.law, '                   Normalized fluxes before extinction correction'
    print 'Corrected for reddening and underlying absorption'
    print catalog_emission_lines[Halpha_idx], '    ', I_Halpha, '                  ', Halpha
    print catalog_emission_lines[Hbeta_idx], '    ', I_Hbeta, '                          ', Hbeta
    print 'theoretical ratio Ha/Hb = %0.3f' % (I_theo_HaHb)
    print '      observed Ha/Hb = %0.3f           raw Ha/Hb = %0.3f' % (I_obs_HaHb, raw_ratio)
    return Idered, I_dered_norCorUndAbs

def det_I_uncertainties(catalog_emission_lines, positive_normfluxes, all_err_cont_fit):
    ''' Uncertainty Calculation '''
    ### STIS Data Handbook states this depends on many factors so I set it up to 2% because
    ### they suggest 1% to increase considerably for faint objects.
    abs_flux_calibration_err = 2.0
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
                    #print 'new faintest_line = ', faintest_line
    chosen_faintest_line = faintest_line
    print 'my definition of S/N~3 after min_counts =', min_counts
    print 'chosen_faintest_line =', chosen_faintest_line, '+- 33%'
    ### In order to account for the error in the continuum fitting, I realized that the 
    ### polynomial order does not really change the intensity, however, the window width does!
    ### To estimate the error on the polynomial fitting I added 1/sum(err), 
    ### were err=1/sqroot(points in window width)
    #print 'all_err_cont_fit', all_err_cont_fit
    ### Add errors and present them in percentage of line intensities
    percent_Iuncert = []
    absolute_Iuncert = []
    S2N = []
    for w, F_norm in zip(catalog_emission_lines, positive_normfluxes):
        if w <= 2000.0:
            e = all_err_cont_fit[0]
        elif w > 2000.0 and w < 5000.0:
            e = all_err_cont_fit[1]
        elif w >= 5000.0:
            e = all_err_cont_fit[2]
        per_Iu = (e*e + 10000*(chosen_faintest_line/9)/(F_norm) + abs_flux_calibration_err*abs_flux_calibration_err)**0.5
        abs_Iu = (per_Iu * F_norm) / 100.
        sn = F_norm/per_Iu*100.
        percent_Iuncert.append(per_Iu)
        absolute_Iuncert.append(abs_Iu)
        S2N.append(sn)
        #print w, '   F_norm = %0.2f   I_norm = %0.2f   uncert_percent = %0.2f   abs_uncert = %0.2f' % (F_norm, I_norm, per_Iu, abs_Iu), '   S/N = ', sn
    return percent_Iuncert, absolute_Iuncert, S2N

def find_flambdas(cHbeta, I_dered_norCorUndAbs, positive_normfluxes):
    # Finding the f_lambda values
    flambdas = []
    for Icor, Iobs in zip(I_dered_norCorUndAbs, positive_normfluxes):
        f12 = (np.log10(Icor) - np.log10(Iobs)) / cHbeta
        flambdas.append(f12)
    return flambdas

def corr_ColExcit(TO2gar, catalog_emission_lines, element_emission_lines, Idered, Hlines):
    '''
    This function is the second iteration of reddening correction. It fits the observed H lines given to the theoretical
    ones found with INTRAT by Storey & Hummer (1995).
    # Idered = first iteration of reddening correction
    # Hlines = list of the hydrogen wavelengths to look for (typically from Halpha to H12).
    # theoCE = theoretical hydrogen intensities corrected for collisional excitation.
    '''
    ### For collisional excitation, interpolate from Table 1 of Peimbert, Luridiana, Peimbert (2007, ApJ, 666, 636)
    # Table 1
    # Objects = NGC346, NGC2363, Haro29, SBS0335-052, IZw18
    TeOII = [12600.0, 13800.0, 14000.0, 15600.0, 15400]
    xalpha = [0.011, 0.037, 0.033, 0.086, 0.070]
    #xbeta = [0.007, 0.027, 0.021, 0.066, 0.053]
    # x_lambda = I_col/I_tot
    # now do the interpolation according to the [OII] temperature
    xL = []
    xalpha_interp = np.interp(TO2gar, TeOII, xalpha)
    xL.append(xalpha_interp)
    xbeta = xalpha_interp * 0.67
    xL.append(xbeta)
    L = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]
    for l in L:
        xl = xalpha_interp / ( 2**((l-2.0)/3.0) )
        xL.append(xl)
    # For the following hydrogen lines recalculate the intensity correcting for collisional exitation
    norm_IcorrCE = [] #intensities corrected for collisional excitation
    obs_ratios = []
    found_Hlines = []
    for w, el, I in zip(catalog_emission_lines, element_emission_lines, Idered):
        for h, l in zip(Hlines, xL):
            if (w == h) and (el == 'H'):
                found_Hlines.append(h)
                newI = I * (1-l)
                normI = newI/100.
                print w, 'before', I, '  after collisionally excited corrected', newI, '  ratio2Hbeta', normI
                norm_IcorrCE.append(newI)
                obs_ratios.append(normI)
                if w == 4102:
                    norm_H6theo = normI
    return norm_IcorrCE, obs_ratios, found_Hlines, norm_H6theo

def find_Chi_of_CE(TO2gar, catalog_emission_lines, element_emission_lines, Idered, Hlines, theoCE, percent_Iuncert):
    # Correct for collisional excitation
    IcorrCE, obs_ratios, found_Hlines, norm_H6theo = corr_ColExcit(TO2gar, catalog_emission_lines, element_emission_lines, Idered, Hlines)
    # Recalculate the intensities of the most prominent hydrogen lines (Halpha through H12) to match them 
    # with the theoretical ratios given by INTRAT (Storey & Hummer, 1995, MNRAS, 272, 41).
    uncert = []
    for H, Ic, obsr in zip(found_Hlines, IcorrCE, obs_ratios):
        idx = Hlines.index(H)
        if H in catalog_emission_lines:
            H_index = catalog_emission_lines.index(H)
            u = (1 - (obsr / theoCE[idx])) * 100.0 / percent_Iuncert[H_index]
            I = Idered[H_index]
            print H, 'theo_ratio =', theoCE[idx], 'obs_ratio', obsr, '   Icorr =', Ic, '   Idered = ', I#, '   percent_Iuncert[H_index]=', percent_Iuncert[H_index]
            print '   error de excitacion colisional = ', u 
            uncert.append(u)
    # In order to find the best combination of C_Hbeta and EWabsHbeta determine chi squared
    sqs = []
    for u in uncert:
        nu = u * u
        sqs.append(nu)
    Chi_sq = sum(sqs)
    return Chi_sq, norm_H6theo
    
def redcor2(I_theo_HaHb, theoCE, Hlines, TO2gar, Idered, C_Hbeta, EWabsHbeta, catalog_emission_lines, corr_undelyingAbs_EWs, 
            rounded_catalog_wavelength, element, ion, forbidden, how_forbidden, observed_wavelength, flux, 
            intensities, EW, continuum, all_err_cont_fit):
    number_iterations = 14 #this number must be even
    EWabsHbeta_increase = 0.1
    C_Hbeta_increase = 0.01
    Chi_sq_models = []
    EWabsHbeta_values = []
    EWabsHbeta_values.append(EWabsHbeta)
    C_Hbeta_values = []
    C_Hbeta_values.append(C_Hbeta)
    Halpha_idx = catalog_emission_lines.index(6563.)
    Hbeta_idx = catalog_emission_lines.index(4861.)
    H6_idx = catalog_emission_lines.index(4102)
    dif_TheoObs_H6Hb_values = []
    # Find the faintest detected emission line: get rid of negative fluxes
    catalog_emission_lines, _, element_emission_lines, _, _, _, positive_normfluxes, _, positive_norm_intensities, _ = find_emission_lines(rounded_catalog_wavelength, element, ion, forbidden, how_forbidden, observed_wavelength, flux, intensities, EW, continuum)
    I_obs_H6Hb = catalog_emission_lines[H6_idx] / catalog_emission_lines[Hbeta_idx]
    # Determine uncertainties
    percent_Iuncert, _, _ = det_I_uncertainties(catalog_emission_lines, positive_normfluxes, all_err_cont_fit)
    Chi_sq, I_theo_H6Hb = find_Chi_of_CE(TO2gar, catalog_emission_lines, element_emission_lines, Idered, Hlines, theoCE, percent_Iuncert)
    Chi_sq_models.append(Chi_sq)
    dif_TheoObs_H6Hb = np.fabs(I_theo_H6Hb - I_obs_H6Hb) 
    dif_TheoObs_H6Hb_values.append(dif_TheoObs_H6Hb)
    I_obs_HaHb = Idered[Halpha_idx] / Idered[Hbeta_idx]
    print ' ***    I_theo_HaHb =', I_theo_HaHb, '   I_obs_HaHb =', I_obs_HaHb
    diff_HaHb_values = []
    diff_HaHb = np.fabs(I_theo_HaHb - I_obs_HaHb)
    diff_HaHb_values.append(diff_HaHb)
    # First, variate EWabsHbeta with C_Hbeta fixed
    for EWabsHbeta_iterations in range(0, number_iterations):
        print 'EWabsHbeta_iterations', EWabsHbeta_iterations
        if I_theo_HaHb < I_obs_HaHb:
            EWabsHbeta = EWabsHbeta + EWabsHbeta_increase
        elif I_theo_HaHb > I_obs_HaHb:
            EWabsHbeta = EWabsHbeta - EWabsHbeta_increase
            if EWabsHbeta < 0.0:
                EWabsHbeta = 0.00001
        EWabsHbeta_values.append(EWabsHbeta)
        intensities = underlyingAbsCorr(EWabsHbeta, corr_undelyingAbs_EWs, continuum, flux)
        # Find the faintest detected emission line: get rid of negative fluxes
        catalog_emission_lines, _, element_emission_lines, _, _, _, positive_normfluxes, _, positive_norm_intensities, _ = find_emission_lines(rounded_catalog_wavelength, element, ion, forbidden, how_forbidden, observed_wavelength, flux, intensities, EW, continuum)
        # Determine uncertainties
        percent_Iuncert, _, _ = det_I_uncertainties(catalog_emission_lines, positive_normfluxes, all_err_cont_fit)
        # Dered again and find the Chi_squared of that model
        cHbeta = 0.434*C_Hbeta
        Idered, _ = Halpha2Hbeta_dered(I_theo_HaHb, cHbeta, catalog_emission_lines, positive_normfluxes, positive_norm_intensities)
        I_obs_H6Hb = catalog_emission_lines[H6_idx] / catalog_emission_lines[Hbeta_idx]
        Chi_sq, norm_H6theo = find_Chi_of_CE(TO2gar, catalog_emission_lines, element_emission_lines, Idered, Hlines, theoCE, percent_Iuncert)
        Chi_sq_models.append(Chi_sq)
        dif_TheoObs_HaHb = np.fabs(norm_H6theo - I_obs_H6Hb) 
        dif_TheoObs_H6Hb_values.append(dif_TheoObs_HaHb)
        I_obs_HaHb = Idered[Halpha_idx] / Idered[Hbeta_idx]
        print ' ***    I_theo_HaHb =', I_theo_HaHb, '   I_obs_HaHb =', I_obs_HaHb
        diff_HaHb = np.fabs(I_theo_HaHb - I_obs_HaHb)
        diff_HaHb_values.append(diff_HaHb)
        EWabsHbeta_iterations = EWabsHbeta_iterations + 1
    # Second, variate C_Hbeta with EWabsHbeta fixed
    EWabsHbeta = EWabsHbeta_values[0]
    for C_Hbeta_iterations in range(0, number_iterations):
        print 'C_Hbeta_iterations =', C_Hbeta_iterations
        if I_theo_HaHb < I_obs_HaHb:
            C_Hbeta = C_Hbeta + C_Hbeta_increase
        elif I_theo_HaHb > I_obs_HaHb:
            C_Hbeta = C_Hbeta - C_Hbeta_increase
            if C_Hbeta < 0.0:
                C_Hbeta = 0.00001            
        C_Hbeta_values.append(C_Hbeta)
        cHbeta = 0.434*C_Hbeta
        intensities = underlyingAbsCorr(EWabsHbeta_values[0], corr_undelyingAbs_EWs, continuum, flux)
        # Find the faintest detected emission line: get rid of negative fluxes
        catalog_emission_lines, _, element_emission_lines, _, _, _, positive_normfluxes, _, positive_norm_intensities, _ = find_emission_lines(rounded_catalog_wavelength, element, ion, forbidden,
                        how_forbidden, observed_wavelength, flux, intensities, EW, continuum)
        # Determine uncertainties
        percent_Iuncert, _, _ = det_I_uncertainties(catalog_emission_lines, positive_normfluxes, all_err_cont_fit)
        # Dered again and find the Chi_squared of that model
        Idered, _ = Halpha2Hbeta_dered(I_theo_HaHb, cHbeta, catalog_emission_lines, positive_normfluxes, positive_norm_intensities)
        I_obs_H6Hb = catalog_emission_lines[H6_idx] / catalog_emission_lines[Hbeta_idx]
        Chi_sq, norm_H6theo = find_Chi_of_CE(TO2gar, catalog_emission_lines, element_emission_lines, Idered, Hlines, theoCE, percent_Iuncert)
        Chi_sq_models.append(Chi_sq)
        dif_TheoObs_HaHb = np.fabs(norm_H6theo - I_obs_H6Hb) 
        dif_TheoObs_H6Hb_values.append(dif_TheoObs_HaHb)
        I_obs_HaHb = Idered[Halpha_idx] / Idered[Hbeta_idx]
        print ' ***    I_theo_HaHb =',I_theo_HaHb, '   I_obs_HaHb =', I_obs_HaHb
        diff_HaHb = np.fabs(I_theo_HaHb - I_obs_HaHb)
        diff_HaHb_values.append(diff_HaHb)
        C_Hbeta_iterations = C_Hbeta_iterations + 1
    # With all 41 models find the one that has the smallest Chi_sq
    print 'Chi_sq_models:', Chi_sq_models
    minChi = min(Chi_sq_models)
    minChi_idx = Chi_sq_models.index(minChi)
    print 'minChi', minChi, 'minChi_idx', minChi_idx
    # Now find the model that has the closest observed Hdelta/Hbeta ratio to the theoretical one
    min_dif_TheoObs_H6Hb = min(dif_TheoObs_H6Hb_values)
    min_dif_TheoObs_H6Hb_idx = dif_TheoObs_H6Hb_values.index(min_dif_TheoObs_H6Hb)
    print 'min_dif_TheoObs_HaHb = ', min_dif_TheoObs_H6Hb, 'min_dif_TheoObs_HaHb_idx', min_dif_TheoObs_H6Hb_idx
    # Calculate the final dereddend values but keep in mind that model 0 is the first reddening iteration, 
    # if there were number_iterations = 10
    # model 5 through 9 is where EWabsHbeta varied and C_Hbeta was fixed at model 0, and
    # model 10 though 14 is where C_Hbeta varied and EWabsHbeta was fixed at model 0.
    print 'LENGTHS OF CHbeta, EWabsHbeta, and diff_HaHb_values lists: ', len(C_Hbeta_values), len(EWabsHbeta_values), len(diff_HaHb_values)
    if minChi_idx == min_dif_TheoObs_H6Hb_idx:
        print 'min indeces are the same!'
    tolerance = 0.005
    if (I_obs_HaHb < I_theo_HaHb+tolerance) and (I_obs_HaHb > I_theo_HaHb-tolerance):
        print ' VALUE WITHIN TOLERANCE!'
        EWabsHbeta = EWabsHbeta_values[0]
        C_Hbeta = C_Hbeta_values[minChi_idx - number_iterations]
    elif (I_obs_HaHb > I_theo_HaHb+tolerance) or (I_obs_HaHb < I_theo_HaHb-tolerance):
        print ' VALUE STILL FAR... LOOKING FOR ALTERNATIVE...'
        min_diff_HaHb = min(diff_HaHb_values)
        min_diff_HaHb_idx = diff_HaHb_values.index(min_diff_HaHb)
        if min_diff_HaHb_idx >= number_iterations:
            idx = min_diff_HaHb_idx - number_iterations
            EWabsHbeta = EWabsHbeta_values[0]
            C_Hbeta = C_Hbeta_values[idx]
        else:
            EWabsHbeta = EWabsHbeta_values[min_diff_HaHb_idx]
            C_Hbeta = C_Hbeta_values[0]
    print 'Chi_sq_models', Chi_sq_models
    print 'dif_TheoObs_H6Hb_values', dif_TheoObs_H6Hb_values
    cHbeta = 0.434*C_Hbeta
    intensities = underlyingAbsCorr(EWabsHbeta, corr_undelyingAbs_EWs, continuum, flux)
    # Find the faintest detected emission line: get rid of negative fluxes
    catalog_emission_lines, wavs_emission_lines, element_emission_lines, ion_emission_lines, forbidden_emission_lines, how_forbidden_emission_lines, positive_normfluxes, pos_calc_cont, positive_norm_intensities, EW_emission_lines = find_emission_lines(rounded_catalog_wavelength, element, ion, forbidden,
                        how_forbidden, observed_wavelength, flux, intensities, EW, continuum)
    emission_lines_info = [catalog_emission_lines, wavs_emission_lines, element_emission_lines, ion_emission_lines, forbidden_emission_lines, how_forbidden_emission_lines, positive_normfluxes, pos_calc_cont, positive_norm_intensities, EW_emission_lines]
    # Dered again and find the Chi_squared of that model
    norm_Idered, I_dered_norCorUndAbs = Halpha2Hbeta_dered(I_theo_HaHb, cHbeta, catalog_emission_lines, positive_normfluxes, positive_norm_intensities)
    flambdas = find_flambdas(cHbeta, I_dered_norCorUndAbs, positive_normfluxes)
    dereddening_info = [EWabsHbeta, C_Hbeta, norm_Idered, I_dered_norCorUndAbs, flambdas]
    # Determine uncertainties    
    percent_Iuncert, absolute_Iuncert, S2N = det_I_uncertainties(catalog_emission_lines, positive_normfluxes, all_err_cont_fit)
    uncertainties_info = [percent_Iuncert, absolute_Iuncert, S2N]
    I_obs_HaHb = norm_Idered[Halpha_idx] / norm_Idered[Hbeta_idx]
    print ' ***    I_theo_HaHb =', I_theo_HaHb, '   I_obs_HaHb =', I_obs_HaHb
    print''
    print 'First iteration of reddening correction:   EWabsHbeta = %0.3f  C_Hbeta = %0.3f' % (EWabsHbeta_values[0], C_Hbeta_values[0])
    print '    The best combination was:              EWabsHbeta = %0.3f  C_Hbeta = %0.3f' % (EWabsHbeta, C_Hbeta)
    print '                                                     this means cHbeta = %0.3f' % (cHbeta)
    return (emission_lines_info, dereddening_info, uncertainties_info)


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
intensities = underlyingAbsCorr(EWabsHbeta, corr_undelyingAbs_EWs, continuum, flux)

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

### Uncertainty Calculation
# Find the faintest detected emission line: get rid of negative fluxes
catalog_emission_lines, wavs_emission_lines, element_emission_lines, ion_emission_lines, forbidden_emission_lines, how_forbidden_emission_lines, positive_normfluxes, pos_calc_cont, positive_norm_intensities, EW_emission_lines = find_emission_lines(rounded_catalog_wavelength, element, ion, forbidden,
                        how_forbidden, observed_wavelength, flux, intensities, EW, continuum)
### Determine uncertainties
percent_Iuncert, absolute_Iuncert, S2N = det_I_uncertainties(catalog_emission_lines, positive_normfluxes, all_err_cont_fit)

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

# Find dered intensities with and without underlying absorption
norm_Idered, I_dered_norCorUndAbs = Halpha2Hbeta_dered(I_theo_HaHb, cHbeta, catalog_emission_lines, positive_normfluxes, positive_norm_intensities)

flambdas = find_flambdas(cHbeta, I_dered_norCorUndAbs, positive_normfluxes)

percent_Iuncert, absolute_Iuncert, S2N = det_I_uncertainties(catalog_emission_lines, positive_normfluxes, all_err_cont_fit)

# Write the first round of reddeding correction in pyneb readable format
tfile1stRedCor = os.path.join(results4object_path, object_name+"_Case"+case+"_1stRedCor.txt")
tf = open(tfile1stRedCor, 'w+')
#print >> tf,  ('{:<15} {:>15} {:>15}'.format('Ion_line', 'Intensity', 'Abs Error'))
#print >> tf, 'cHbeta  %0.3f' % cHbeta
lines_pyneb_matches = []
#for cw, el, io, Ic, er in zip(catalog_emission_lines, element_emission_lines, ion_emission_lines, positive_norm_Icorr, absolute_Iuncert):
for cw, el, io, Ic, er in zip(catalog_emission_lines, element_emission_lines, ion_emission_lines, norm_Idered, absolute_Iuncert):
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

### Density measurement from [Fe III] lines -- taken from Peimbert, Pena-Guerrero, Peimbert (2012, ApJ, 753, 39)
I4986 = 0.0
I4987 = 0.0
I4658 = 0.0
for w, i in zip(catalog_emission_lines, norm_Idered):
    if int(w) == 4986:
        I4986 = i
    elif int(w) == 4987:
        I4987 = i
    elif int(w) == 4658:
        I4658 = i
if (I4986 != 0.0) and (I4987 != 0) and (I4658 !=0):
    log_Fe3den = 2 - ( (np.log10((I4986+I4987)/I4658) - 0.05 - 0.25*(np.log10(TO3-4))) / (0.66 - 0.18*(np.log10(TO3-4))) )
    Fe3den = 10**(log_Fe3den)
    print 'Density measured from [Fe3] lines:', Fe3den
else:
    print 'No [Fe3] density available.'

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


### Second iteration of extinction correction: Collisional Excitation
### Following analysis presented in Pena-Guerrero, Peimbert, Peimbert, Ruiz (2012, ApJ, 746, 115) & Peimbert, Pena-Guerrero, Peimbert (2012, ApJ, 753, 39).
# If the [O II] temperature was not obtained directly from observations, get an estimate
### Get temperature of OII from OIII. 
# Using equation of Peimbert, Peimbert, & Luridiana (2002, ApJ, 565, 668) - Based on data of Stasinska's models.
TO2pei = 2430. + TO3 * (1.031 - TO3/54350.)
print 'This is the theoretically obtained temperature of O2 from Peimbert etal 2002 = ', TO2pei
# Using equation of Garnett, D. R. 1992, AJ, 103, 1330
TO2gar = 0.7 * TO3 + 3000.
print 'Theoretically obtained temperature of O2 from Garnet 1992 = ', TO2gar

if reddeningCorrection2 == True:
    ### Correct for collisional excitation
    # Hydrogen lines to be considered for correction
    #        Halpha,  H5,   H6,   H7,   H8,   H9,  H10,   H11,  H12
    #Hlines = [6563, 4340, 4101, 3967, 3889, 3835, 3798, 3771, 3750]
    # HOWERVER H7 and H8 are contaminated by HeI (H8 is also contaminated with [NeIII] 
    Hlines = [6563, 4340, 4102, 3967, 3798, 3771, 3750]
    # Recalculate the intensities of the most prominent hydrogen lines (Halpha through H12) to match them with the
    # theoretical ratios given by INTRAT (Storey & Hummer, 1995, MNRAS, 272, 41).
    # This has to be done separately object by object for both Case A and B with its respective [OII] temperature.
    #               Halpha, H5,   H6,   H7,    H8,     H9,     H10,    H11,    H12
    theoCE_caseA = [2.80, 0.47, 0.265, 0.164, 0.109, 0.0760, 0.0553, 0.0415, 0.0320]
    theoCE_caseB = [2.85, 0.469, 0.260, 0.160, 0.105, 0.733, 0.0532, 0.0398, 0.0306]
    # write the file with these theoretical ratios
    file_theoreticalHratios = os.path.join(results4object_path, object_name+"_Case"+case+'_HtheoRatios.txt')
    Htheofile = open(file_theoreticalHratios, 'w+')
    print >> Htheofile, '# Theoretical ratios obtained with INTRAT by Storey & Hummer (1995, MNRAS, 272, 41).' 
    print >> Htheofile, '# Temperature of [O II] used : ', TO2gar
    print >> Htheofile, ('# {:<13} {:<10} '.format('Wavelength', 'H_X / H_beta'))
    # Determine which theoretical ratios to use
    if case == 'A':
        theoCE = theoCE_caseA
    elif case == 'B':
        theoCE = theoCE_caseB
    for w, r in zip(Hlines, theoCE):
        print >> Htheofile, ('{:<15} {:<10.4f} '.format(w, r))
    Htheofile.close()
    # Do the second iteration of correction: collisional excitation
    I_theo_HaHb = theoCE[0]
    emission_lines_info, dereddening_info, uncertainties_info = redcor2(I_theo_HaHb, theoCE, Hlines, TO2gar, norm_Idered, C_Hbeta, EWabsHbeta, catalog_emission_lines, corr_undelyingAbs_EWs, 
                rounded_catalog_wavelength, element, ion, forbidden, how_forbidden, observed_wavelength, flux, 
                intensities, EW, continuum, all_err_cont_fit)
    # separate the results
    catalog_emission_lines, wavs_emission_lines, element_emission_lines, ion_emission_lines, forbidden_emission_lines, how_forbidden_emission_lines, positive_normfluxes, pos_calc_cont, positive_norm_intensities, EW_emission_lines = emission_lines_info
    EWabsHbeta, C_Hbeta, norm_Idered, I_dered_norCorUndAbs, flambdas = dereddening_info
    cHbeta = 0.434*C_Hbeta
    percent_Iuncert, absolute_Iuncert, S2N = uncertainties_info
    
    ### In order to account for the error in the continuum fitting, I realized that the 
    ### polynomial order does not really change the intensity, however, the window width does!
    ### To estimate the error on the polynomial fitting I added 1/sum(err), 
    ### were err=1/sqroot(points in window width)
    #print 'all_err_cont_fit', all_err_cont_fit
    ### Add errors and present them in percentage of line intensities
    emission_lines_file = os.path.join(results4object_path, object_name+"_Case"+case+'_emlines.txt')
    emfile = open(emission_lines_file, 'w+')
    print >> emfile, '# Positive EW = emission        Negative EW = absorption' 
    print >> emfile, '# C_Hbeta = %0.3f   or   cHbeta = %0.3f' % (C_Hbeta, cHbeta)
    print >> emfile, '# EWabs_Hbeta = %0.3f' % EWabsHbeta
    print >> emfile, '# I_theo_HaHb = %0.3f' % I_theo_HaHb
    print >> emfile, ('# {:<13} {:<8} {:>6} {:<12} {:<12} {:<16} {:<16} {:<18} {:<18} {:<16} {:<16} {:<16}'.format('Wavelength', 'Element', 'Ion', 'Forbidden', 'How much', 'f_lambda', 'Flux [cgs]', 'Intensity [cgs]', '% Err', 'Abs err', 'EW [A]', 'S/N'))
    for w, el, i, f, h, ls, F_norm, I_norm, pu, au, eqw in zip(catalog_emission_lines, element_emission_lines, ion_emission_lines, forbidden_emission_lines, 
                                 how_forbidden_emission_lines, flambdas, positive_normfluxes, positive_norm_intensities, percent_Iuncert, absolute_Iuncert, EW_emission_lines):
        sn = F_norm/pu*100.
        #print w, '   F_norm = %0.2f   I_norm = %0.2f   uncert_percent = %0.2f   abs_uncert = %0.2f' % (F_norm, I_norm, pu, au), '   S/N = ', sn
        print >> emfile, ('{:<15.0f} {:<8} {:>6} {:<12} {:<12} {:<16.3f} {:<16.3f} {:<18.3f} {:<18.3f} {:<16.3f} {:<16.3f} {:<16.3f}'.format(w, el, i, f, h, ls, F_norm, I_norm, pu, au, eqw, sn))
    emfile.close()

print '\n Code finished for Case', case
