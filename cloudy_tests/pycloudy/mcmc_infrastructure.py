from __future__ import division
import os
import numpy as np
import matplotlib.pyplot as plt
import pyCloudy as pc
import time
import copy
import string
import emcee
import triangle
from glob import glob
#from scipy import stats
from science import spectrum

''' 
This script contains all the functions to serve as infrastructure for running a MCMC chain. 
'''

def find_cloudyexe(cloudyexe_path):
    ''' Changing the location and version of the cloudy executable. '''
    full_cloudyexe_path = os.path.abspath(cloudyexe_path)
    pc.config.cloudy_exe = full_cloudyexe_path
    
def lngauss(x, x0, s):
    return -0.5*((x-x0)/s)**2

def lntophat(x, a, b):
    # b has to be grater than a
    if a > b:
        bb = a
        aa = b
    elif a == b:
        print 'Boundaries are the same in lntophat function... change them please.'
        exit()
    else:
        aa = a
        bb = b
    if (aa < x) and (x < bb):
        return np.log(1/(bb-aa))
    else:
        return -1.e10

class PyCloudy_model:
    def __init__(self, model_name, dens, emis_tab, abunds, stb99_table, age, dir, verbosity, options, 
                 iterations, keep_files):
        '''    This set of initial conditions is specific for this particular project.
        model_name = just the name of the model (without path)
        dens = Hydrogen density in units of log cm^-3
        emis_tab = lines to find (the integral over the volume of the emissivity gives the intensity of the line)
        abunds = this is the set of target abundances for Cloudy to match, in units of log(X/H)
        stb99_table = ionizing source, should be in format: 'table star "starburst99.mod"'
        age = look for this age in the stb99 table, should be in format: 'age=4.0' 
        dir = where to put the output files
        verbosity = how many out files do you want:
             -1 = Quiet and no log into log_ variable
             0 = Quiet
             1 = Only Errors
             2 = Errors and Warnings
             3 = very verbose, print Errors, Warnings and Messages 
             4 = debug
        options = This only to make the code run faster
        iterations = integer -- total number of iterations, remember that line intensities are not reliable at 1st iteration
        keep_files = the files that we want to keep (the others will be erased)
        '''
        self.model_name = model_name
        self.dens = dens
        self.emis_tab = emis_tab
        self.abunds = abunds
        self.stb99_table = stb99_table
        self.age = age
        #if dir == None:
        #    self.dir = './'
        #else:
        self.dir = dir
        #if verbosity == None:
        #    pc.log_.level = 3
        #else:
        self.verbosity = verbosity
        #if options != None:
        self.options = options
        #else:
        #    self.options = None
        #if iterations != None:
        self.iterations = iterations
        #else:
        #    self.iterations = None
        #if keep_files == None:
        #    self.keep_files = None
        #else:
        self.keep_files = keep_files

    def runCldy_with_initconds(self):
        '''    
        Read this specific set of initial conditions for this specific for this particular project.
        '''
        self.full_model_name = '{0}{1}'.format(self.dir, self.model_name)
        # Define verbosity to high level (will print errors, warnings and messages)
        # Defining the object that will manage the input file for Cloudy
        c_input = pc.CloudyInput(self.full_model_name)
        # Filling the object with the parameters
        # Defining the ionizing SED: Effective temperature and luminosity.
        # it is set for the total flux of Halpha
        c_input.set_star(SED = self.stb99_table, SED_params = self.age, lumi_unit='f(nu)', lumi_value=-12.3316)
        # Defining the density. You may also use set_dlaw(parameters) if you have a density law defined in dense_fabden.cpp.
        c_input.set_cste_density(self.dens)
        c_input.set_abund(ab_dict = self.abunds, nograins = True)
        if self.options != None:
            c_input.set_other(self.options)
        c_input.set_emis_tab(self.emis_tab) # better use read_emis_file(file) for long list of lines, where file is an external file.
        if self.iterations != None:
            c_input.set_iterate(self.iterations) # (0) for no iteration, () for one iteration, (N) for N iterations.
        # Writing the Cloudy inputs. to_file for writing to a file (named by full_model_name). verbose to print on the screen.
        c_input.print_input(to_file = True, verbose = False)
            
        # Printing some message to the screen
        pc.log_.message('Running {0}'.format(self.model_name), calling = 'stb99_test1')
        
        # Running Cloudy with a timer. Here we reset it to 0.
        pc.log_.timer('Starting Cloudy', quiet = True, calling = 'stb99_test1')
        c_input.run_cloudy()
        pc.log_.timer('Cloudy ended after seconds:', calling = 'stb99_test1')
        
    def read_Cldyouts(self):
        '''
        This function will read the output files of Cloudy.
        full_model_name = is the full path of the model 
        '''
        # Reading the Cloudy outputs in the Mod CloudyModel object
        Mod = pc.CloudyModel(self.full_model_name)
        dir(Mod) # This is the online answering way
        Mod.print_stats()
        #Mod.print_lines()
        Mod.get_ab_ion_vol_ne('O',2)
        Mod.get_T0_ion_vol_ne('O', 2)
        Mod.log_U_mean
        Mod.log_U_mean_ne
        print('T0 = {0:7.1f}K, t2 = {1:6.4f}'.format(Mod.T0, Mod.t2))
        print('Hbeta Equivalent width = {0:6.1f}, Hbeta Surface Brightness = {1:4.2e}'.format(Mod.get_Hb_EW(), Mod.get_Hb_SB()))
        # printing line intensities into a file and on the screen
        line_file = os.path.join(self.dir, self.model_name+'_lines.txt')
        of = open(line_file, 'w+')
        print >> of, '{:<5} = {:>8}'.format('T(O3)', Mod.get_T0_ion_vol_ne('O', 2))
        print >> of, '{:<5} = {:>8}'.format('T(O2)', Mod.get_T0_ion_vol_ne('O', 1))
        print >> of, 'T0 = {0:7.1f}K, t2 = {1:6.4f}'.format(Mod.T0, Mod.t2)
        print >> of, 'Hbeta Equivalent width = {0:6.1f}'.format(Mod.get_Hb_EW())
        print >> of, '{:<10} {:>11} {:>9}'.format('Line', 'Intensity', 'I/Hbeta')
        for line in Mod.emis_labels:
            #print('{0} {1:10.3e} {2:7.2f}'.format(line, Mod.get_emis_vol(line), Mod.get_emis_vol(line) / Mod.get_emis_vol('H__1__4861A') * 100.))
            print >> of, '{0} {1:10.3e} {2:7.2f}'.format(line, Mod.get_emis_vol(line), 
                                                         Mod.get_emis_vol(line) / Mod.get_emis_vol('H__1__4861A') * 100.)
        of.close()
        return line_file
        
    def clean_Cldyoutfiles(self):
        '''  Clean the output files you want... Default are:
             SAVE_LIST = [['radius', '.rad'],
                          ['continuum', '.cont'],
                          ['physical conditions', '.phy'],
                          ['overview', '.ovr'],
                          ['heating', '.heat'],
                          ['cooling', '.cool'],
                          ['optical depth', '.opd']
                          ]
        and 
             SAVE_LIST_ELEMS = [['hydrogen','.ele_H'],
                                ['helium','.ele_He'],
                                ['carbon','.ele_C'],
                                ['nitrogen','.ele_N'],
                                ['oxygen','.ele_O'],
                                ['argon','.ele_Ar'],
                                ['neon','.ele_Ne'],
                                ['sulphur','.ele_S'],
                                ['chlorin','.ele_Cl'],
                                ['iron','.ele_Fe'],
                                ['silicon','.ele_Si']]
    
        '''
        # I tried to change these but it needs them in order to read output.... 
        #pc.config.SAVE_LIST = SAVE_LIST
        #pc.config.SAVE_LIST_ELEMS = SAVE_LIST_ELEMS
        ''' except_files = is the list of files we actually want to keep...''' 
        default_extensions_list = ['.emis', # just to avoid overload with the huge run
                                   '.rad', '.cont', '.phy', '.ovr', '.heat', '.cool', '.opd',  # files pyCloudy needs to read outputs
                                   '.ele_H', '.ele_He', '.ele_C', '.ele_N', '.ele_O', '.ele_Ar', # out elements
                                   '.ele_Ne', '.ele_S', '.ele_Cl', '.ele_Fe', '.ele_Si']
        final_default_extensions_list = copy.deepcopy(default_extensions_list)
        if self.keep_files != None: 
            for df in default_extensions_list:
                if df in self.keep_files:
                    final_default_extensions_list.pop(df)
        # Do the cleaning
        for fdf in final_default_extensions_list:
            file2beerased = glob(self.dir+self.model_name+'*'+fdf)
            #print file2beerased[0]
            os.remove(file2beerased[0])
        
    def mk_model(self):
        self.runCldy_with_initconds()
        print'***  got initial conditions and ran Cloudy!'
        lines_file = self.read_Cldyouts()
        print '***  reading Cloudy outputs!'
        self.clean_Cldyoutfiles()
        print'***  cleaning Cloudy output files...'
        return lines_file
    

class Get_Chi2_info:
    ''' This part of the code will determine the Chi squared from by comparing the measured line intensities with those 
    modeled by Cloudy.'''
    def __init__(self, measured_lines, modeled_lines):
        '''
        measured_lines = full path of the file with the benchmark/measured lines
        modeled_lines = full path of the file with the Cloudy modeled lines
        '''
        self.measured_lines = measured_lines
        self.modeled_lines = modeled_lines 
        self.read_files()

    def read_files(self):
        # MEASURED lines
        self.meas_lineIDs = []   # list of IDs in angstroms
        self.meas_Isrel2Hbeta = []   # list of lines relative to Hbeta reddening corrected 
        self.meas_Ierr = []   # absolute error of the line intensity
        self.meas_Iper = []   # percentage error of the line intensity  
        self.meas_EW = []   # NOT BEING USED FOR THE MOMENT
        meas = open(self.measured_lines)
        _ = meas.readline()  # columns header
        for line in meas:
            line = line.strip()   # gets rid of \n at the end of the line
            cols = line.split()   # splits the line into a list of columns
            # the columns of interest are: ID=0, Intensity=9, Ierr=10, Iper=11, EW=12
            self.meas_lineIDs.append(int(np.round(float(cols[0]), decimals=0)))
            self.meas_Isrel2Hbeta.append(float(cols[9]))
            self.meas_Ierr.append(float(cols[10]))
            self.meas_Iper.append(float(cols[11]))
            self.meas_EW.append(float(cols[12]))
        meas.close()
        
        # MODELED lines
        self.mod_lineIDs = []
        self.ions = []
        #self.mod_intensities = []   # NOT BEING USED FOR THE MOMENT
        self.mod_Isrel2Hbeta = []
        mod = open(self.modeled_lines, 'r')
        # Save the temperatures, t2, HbetaEW, and the column headers into special headers
        mod_TO3 = mod.readline()
        mod_TO2 = mod.readline()
        self.mod_T0t2 = mod.readline()   # NOT BEING USED FOR THE MOMENT
        self.HbetaEW = mod.readline()   # NOT BEING USED FOR THE MOMENT
        _ = mod.readline()  # columns header
        for line in mod:   # read the rest of the file line by line as a string
            line = line.strip()   # gets rid of \n at the end of the line
            cols = line.split()   # splits the line into a list of columns
            kk1 = string.split(cols[0], sep='A')
            kk2 = string.split(kk1[0], sep='__')
            ion = kk2[0]   # means that we have a case type 'HE_2__4686A'
            if len(kk2) > 2:   # means that we have a case type 'H__1__6563A'
                ion = kk2[0]+'_'+kk2[1]
            lid = int(kk2[-1])
            #print ion, lid
            self.ions.append(ion)
            self.mod_lineIDs.append(lid)
            #self.mod_intensities.append(float(cols[1]))   # NOT BEING USED FOR THE MOMENT
            self.mod_Isrel2Hbeta.append(float(cols[2]))
        mod.close()
        kk = string.split(mod_TO3, sep='=')
        self.mod_TO3 = kk[1]
        kk = string.split(mod_TO2, sep='=')
        self.mod_TO2 = kk[1]
        print 'model_Te_O3 =', self.mod_TO3, 'model_Te_O2 =', self.mod_TO2
    
    def calculateChi2(self):
        # Find Chi2 using:  Chi2 = Sum( {[I/Hbeta]_obs - [I/Hbeta]_theo}^2 / errI_obs^2 )
        comp = []
        for idtheo in self.mod_lineIDs:
            Itheo_idx = self.mod_lineIDs.index(idtheo)
            Itheo = self.mod_Isrel2Hbeta[Itheo_idx]
            if idtheo in self.meas_lineIDs:
                idx = self.meas_lineIDs.index(idtheo)
                Iobs = self.meas_Isrel2Hbeta[idx]
                Ierr = self.meas_Ierr[idx]
                print self.meas_lineIDs[idx], 'Itheo =', Itheo, '  Iobs =', Iobs, '+-', Ierr
                c = (Iobs - Itheo)**2 / Ierr**2
                print 'individual', c
                comp.append(c)
        Chi2 = sum(comp)
        return Chi2
    
    def getChi2(self):
        self.read_files()
        Chi2 = self.calculateChi2()
        print '\nChi^2 for this model = ', Chi2
        return Chi2
    
class MCMC:
    ''' This class does the actual Markov Chain Monte Carlo (MCMC). This is done to marginalize over some parameters and 
    estimate of the posterior probability function, that is the distribution of parameters that is consistent with the
    data set. The posterior probability function can be described as follows:
                            p(m,b,f|x,y,sigma)  proportional2  p(m,b,f) * p(y|x,sigma,m,b,f),
    where p(y|x,sigma,m,b,f) is the likelihood function and p(m,b,f) is the prior function.   
    '''
    def __init__(self, object_name, manual_measurement, initial_Cloudy_conditions): 
        self.object_name = object_name   # name of the object
        self.manual_measurement = manual_measurement   # manual (manual_measurement=True) or code's measurements (manual_measurement=False) 
        self.get_measured_lines()   # this is part of the initialization because I only want it to locate and read the file once
        self.initial_Cloudy_conditions = initial_Cloudy_conditions
    
    def get_measured_lines(self):
        # get the benchmark measurements
        if self.manual_measurement:
            file2use = '_measuredLI_RedCor_Ebv.txt'   # my manual line measurements
        else:
            file2use = '_RedCor_Ebv.txt'   # my code's line measurements
        measured_lines_file_path = os.path.abspath('../../results/'+self.object_name+'/'+self.object_name+file2use)
        # now retrieve the columns we need from the file
        meas_lineIDs = []   # list of IDs in angstroms
        meas_Isrel2Hbeta = []   # list of lines relative to Hbeta reddening corrected 
        meas_Ierr = []   # absolute error of the line intensity
        meas_Iper = []   # percentage error of the line intensity  
        meas_EW = []   # NOT BEING USED FOR THE MOMENT
        meas = open(measured_lines_file_path)
        _ = meas.readline()  # columns header
        for line in meas:
            line = line.strip()   # gets rid of \n at the end of the line
            cols = line.split()   # splits the line into a list of columns
            # the columns of interest are: ID=0, Intensity=9, Ierr=10, Iper=11, EW=12
            meas_lineIDs.append(int(np.round(float(cols[0]), decimals=0)))
            meas_Isrel2Hbeta.append(float(cols[9]))
            meas_Ierr.append(float(cols[10]))
            meas_Iper.append(float(cols[11]))
            meas_EW.append(float(cols[12]))
        meas.close()
        self.measured_lines = [meas_lineIDs, meas_Isrel2Hbeta, meas_Ierr, meas_Iper, meas_EW]
        return self.measured_lines

    def get_modeled_lines(self, theta):
        # get the Cloudy model
        He, C, N, O, Ne, S = theta   # convert the list into a value and a dictionary
        abunds = {}
        abunds['He'] = He-12.
        abunds['C'] = C-12.
        abunds['N'] = N-12.
        abunds['O'] = O-12.
        abunds['Ne'] = Ne-12.
        abunds['S'] = S-12.
        self.model_name, dens, emis_tab, _, stb99_table, age, self.dir, verbosity, options, iterations, keep_files = self.initial_Cloudy_conditions
        initiate_PyCmodel = PyCloudy_model(self.model_name, dens, emis_tab, abunds, stb99_table, age, self.dir, verbosity, options, iterations, keep_files)
        cldymod = initiate_PyCmodel.mk_model()
        modeled_lines_file_path = os.path.abspath(cldymod)
        mod_lineIDs = []
        ions = []
        #mod_intensities = []   # NOT BEING USED FOR THE MOMENT
        mod_Isrel2Hbeta = []
        mod = open(modeled_lines_file_path, 'r')
        # Save the temperatures, t2, HbetaEW, and the column headers into special headers
        mod_TO3 = mod.readline()
        mod_TO2 = mod.readline()
        mod_T0t2 = mod.readline()   # NOT BEING USED FOR THE MOMENT
        HbetaEW = mod.readline()   # NOT BEING USED FOR THE MOMENT
        _ = mod.readline()  # columns header
        for line in mod:   # read the rest of the file line by line as a string
            line = line.strip()   # gets rid of \n at the end of the line
            cols = line.split()   # splits the line into a list of columns
            kk1 = string.split(cols[0], sep='A')
            kk2 = string.split(kk1[0], sep='__')
            ion = kk2[0]   # means that we have a case type 'HE_2__4686A'
            if len(kk2) > 2:   # means that we have a case type 'H__1__6563A'
                ion = kk2[0]+'_'+kk2[1]
            lid = int(kk2[-1])
            #print ion, lid
            ions.append(ion)
            mod_lineIDs.append(lid)
            #mod_intensities.append(float(cols[1]))   # NOT BEING USED FOR THE MOMENT
            mod_Isrel2Hbeta.append(float(cols[2]))
        mod.close()
        kk = string.split(mod_TO3, sep='=')
        mod_TO3 = kk[1]   # NOT BEING USED FOR THE MOMENT
        kk = string.split(mod_TO2, sep='=')
        mod_TO2 = kk[1]   # NOT BEING USED FOR THE MOMENT
        #print 'model_Te_O3 =', self.mod_TO3, 'model_Te_O2 =', self.mod_TO2
        modeled_lines = [mod_lineIDs, ions, mod_Isrel2Hbeta]
        # Do the cleaning
        final_default_extensions_list = ['.in', '.out', '.txt']
        for fdf in final_default_extensions_list:
            file2beerased = glob(self.dir+self.model_name+'*'+fdf)
            #print file2beerased[0]
            os.remove(file2beerased[0])
        return modeled_lines
        
    def find_Itheo_in_Iobs(self, IDobs, Iobs, Iobserr, IDmod, Imod):
        new_Iobs = []
        new_Iobserr = []
        new_Imod = []
        for idtheo in IDmod:
            Itheo_idx = IDmod.index(idtheo)
            Itheo = Imod[Itheo_idx]
            if idtheo in IDobs:
                idx = IDobs.index(idtheo)
                iobs= Iobs[idx]
                ierr = Iobserr[idx]
                new_Iobs.append(iobs)
                new_Iobserr.append(ierr)
                new_Imod.append(Itheo)
        return new_Iobs, new_Iobserr, new_Imod

    def lnlikehd(self, theta, IDobs, Iobs, Iobserr):
        modeled_lines = self.get_modeled_lines(theta)
        mod_lineIDs, ions, mod_Isrel2Hbeta = modeled_lines
        new_Iobs, new_Iobserr, new_Imod = self.find_Itheo_in_Iobs(IDobs, Iobs, Iobserr, mod_lineIDs, mod_Isrel2Hbeta)
        model = np.array(new_Imod)
        y = np.array(new_Iobs)
        e = np.array(new_Iobserr)
        csq = ( y - model )**2 / e**2
        Chi2 = csq.sum()
        return -Chi2 / 2.0
    
    def lnpriors(self, theta):
        ''' This function encodes any previous knowledge that we have about the parameters: results from other experiments, 
        physically acceptable ranges, etc. It is necessary that you write down priors if you are going to use MCMC because 
        all that MCMC does is draw samples from a probability distribution and you want that to be a probability distribution
        for your parameters. We will use 2 functions: Gaussian and top-hat. 
        '''
        He, C, N, O, Ne, S = theta
        # Set the abundances set to that of the observed values. Since the boundary conditions are already built-in Cloudy,
        # we will allow them to vary in a wide range. The abundances dictionary is expected to have 6 elements:
        he = lntophat(He, 9.5, 11.0) 
        #c = lntophat(C, 6.0, 8.5)
        #n = lntophat(N, 6.0, 8.8) 
        o = lntophat(O, 7.2, 8.9) 
        #ne = lntophat(Ne, 7.0, 8.7)  
        #s = lntophat(S, 4.0, 7.0) 
        #print 'top hat results:', he , c , n , o , ne , s
        # additional constraints
        ClowerthanO = lntophat(C, 6.0, O)
        NlowerthanO = lntophat(N, 6.0, O)
        NelowerthanO = lntophat(Ne, 7.0, O)
        SlowerthanO = lntophat(S, 4.0, O)
        print 'top hat results:', he , ClowerthanO , NlowerthanO , o , NelowerthanO , SlowerthanO
        # check that all conditions are met
        #if he != -1.e10 and c != -1.e10 and n != -1.e10 and o != -1.e10 and ne != -1.e10 and s != -1.e10:
            #return he + c + n + o + ne + s + ClowerthanO + NlowerthanO + NelowerthanO + SlowerthanO
        if he != -1.e10 and ClowerthanO != -1.e10 and NlowerthanO != -1.e10 and o != -1.e10 and NelowerthanO != -1.e10 and SlowerthanO != -1.e10:
            return he + o + ClowerthanO + NlowerthanO + NelowerthanO + SlowerthanO
        else:
            return -1.e10

    def lnprob(self, theta, IDobs, Iobs, Iobserr):
        print 'theta =', theta
        lp = self.lnpriors(theta)
        print 'probabilities: ', lp
        if lp != -1.e10:
            return lp + self.lnlikehd(theta, IDobs, Iobs, Iobserr)
        else:
            return lp
        
    def run_chain(self):
        meas_lineIDs, meas_Isrel2Hbeta, meas_Ierr, meas_Iper, meas_EW = self.measured_lines
        ndim, nwalkers, nruns = 6, 20, 50
        p0 = [[np.random.uniform(9.5, 11.),np.random.uniform(7., 8.1),np.random.uniform(7., 8.7),np.random.uniform(7., 8.9),np.random.uniform(7., 8.2),np.random.uniform(5., 6.7)] for _ in range(nwalkers)]
        #p0 = [[np.random.uniform(8., 11.),np.random.uniform(6., 9.),np.random.uniform(6., 9.),np.random.uniform(6., 9.),np.random.uniform(6., 9.),np.random.uniform(4., 8.)] for _ in range(nwalkers)]
        sampler = emcee.EnsembleSampler(nwalkers, ndim, self.lnprob, args=(meas_lineIDs, meas_Isrel2Hbeta, meas_Ierr))
        #p0, prob, state = sampler.run_mcmc(p0, 10)   # this allows for the first few steps to be "burn-in" type
        #sampler.reset()   # then restart the mcmc at the final position of the "burn-in", p0
        pos, prob, rstate = sampler.run_mcmc(p0, nruns)   # do the mcmc starting at p0 for nruns steps
        # To store the chain....
        chain_file = os.path.abspath(self.dir+self.model_name+"_chain.dat")
        f = open(chain_file, "w")
        f.close()
        count = 1
        for posn, prob, state in sampler.sample( pos, iterations=20, storechain=True ):
            print "COUNT", count
            if count % 1 == 0:
                f = open(chain_file, "a")
                for k in range( posn.shape[0] ):
                    strout = ""
                    for p in pos[k]: strout += "{:8.3f} ".format( p )
                    strout += "{:20.3f}".format( prob[k] )
                    print strout
                    f.write(strout+"\n")
                f.close()
            count += 1
        # best model
        wh = np.where( prob == prob.max() )[0][0]
        p = pos[ wh, : ]
        print 'Values of the', ndim,'dimensions that best fit the data in', nruns, 'runs, are the following:'
        #print '   density of H =', p[0]
        print '   abundances of: He =', p[0], ' C =', p[1], ' N =', p[2], 'O =', p[3], 'Ne =', p[4], 'S =', p[5]
        # samples has an attribute called chain, which is an array with shape (nwalkers, nruns, ndim)
        samples = sampler.chain[:, nruns*0.2:, :].reshape((-1, ndim)) # this discards the first 20% of the runs
        # but for now we'll use all the runs:
        #samples = sampler.chain[:, :, :].reshape((-1, ndim))
        fig = triangle.corner(samples, labels=["$He$", "$C$", "$N$", "$O$", "$Ne$", "$S$", "$\ln\,f$"], truths=[m_true, b_true, np.log(f_true)])
        fig.savefig(os.path.abspath(self.dir+"triangle_test1.jpg"))
        fig = triangle.corner(samples)
        fig.savefig(os.path.abspath(self.dir+"triangle_test2.jpg"))
        # Calculate the uncertainties based on the 16th, 50th and 84th percentiles
        samples[:, ndim-1] = np.exp(samples[:, ndim-1])
        p_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
        print 'mcmc values and uncertainties according to 16th, 50th, and 84th percentiles:'
        print p_mcmc

        
        