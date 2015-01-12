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
import scipy.optimize as op
import copy_reg, pickle
import types
#from emcee.utils import MPIPool
import uuid
import shlex
import subprocess

''' 
This script contains all the functions to serve as infrastructure for running a MCMC chain. 
'''

def _lnprob(*args):
    obj = args[1]
    args = list(args)
    del args[1]
    x = obj.lnprob(*args)
    return x


def find_cloudyexe(cloudyexe_path):
    ''' Changing the location and version of the cloudy executable. '''
    cloudyexe = 'Addons/cloudy/c13.03all/c13.03/source/cloudy.exe'
    cloudyexe_path_list = string.split(os.getcwd(), sep='AptanaStudio3')
    full_cloudyexe_path = os.path.join(cloudyexe_path_list[0], cloudyexe)
    #full_cloudyexe_path = os.path.abspath(cloudyexe_path)
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
        return -np.inf#1.e10

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
        self.lines_file = self.mk_model()


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
        #dir(Mod) # This is the online answering way
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
        self.lines_file = self.read_Cldyouts()
        print '***  reading Cloudy outputs!'
        self.clean_Cldyoutfiles()
        print'***  cleaning Cloudy output files...'
        return self.lines_file
    

    
### MCMC methodology

def get_measured_lines(object_name, manual_measurement):
    # get the benchmark measurements
    if manual_measurement:
        file2use = '_measuredLI_RedCor_Ebv.txt'   # my manual line measurements
    else:
        file2use = '_RedCor_Ebv.txt'   # my code's line measurements
    measured_lines_file = 'carbon/results/'+object_name+'/'+object_name+file2use
    measured_lines_file_path_list = string.split(os.getcwd(), sep='carbon')
    measured_lines_file_path = os.path.join(measured_lines_file_path_list[0], measured_lines_file)
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
    measured_lines = [meas_lineIDs, meas_Isrel2Hbeta, meas_Ierr, meas_Iper, meas_EW]
    return measured_lines

def get_modeled_lines(theta, initial_Cloudy_conditions):   # get the Cloudy model
    # Turn the ratios into abundances in the form log(X/H) and save into a dictionary for Cloudy use
    He, O, CoverO, NoverO, NeoverO, SoverO = theta
    abunds = {}
    abunds['He'] = He-12.
    abunds['O'] = O-12.
    abunds['C'] = CoverO+O-12
    abunds['N'] = NoverO+O-12.
    abunds['Ne'] = NeoverO+O-12.
    abunds['S'] = SoverO+O-12.
    model_name, dens, emis_tab, _, stb99_table, age, dir, verbosity, options, iterations, keep_files = initial_Cloudy_conditions
    print "I received the following initial conditions:"
    print 'model_name, dens, emis_tab, abunds, stb99_table, age, dir, verbosity, options, iterations, keep_files ='
    print model_name, dens, emis_tab, abunds, stb99_table, age, dir, verbosity, options, iterations, keep_files
    cldymod = PyCloudy_model(model_name, dens, emis_tab, abunds, stb99_table, age, dir, verbosity, options, iterations, keep_files)
    #cldymod = initiate_PyCmodel.mk_model()
    modeled_lines_file_path = os.path.abspath(cldymod.lines_file)
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
    mod_TO3 = kk[1] 
    kk = string.split(mod_TO2, sep='=')
    mod_TO2 = kk[1]  
    #print 'model_Te_O3 =', mod_TO3, 'model_Te_O2 =', mod_TO2
    modeled_Otemperatures = [mod_TO3, mod_TO2]
    modeled_lines = [mod_lineIDs, ions, mod_Isrel2Hbeta]
    # Do the cleaning
    final_default_extensions_list = ['.in', '.out', '.txt']
    for fdf in final_default_extensions_list:
        file2beerased = glob(dir+model_name+'*'+fdf)
        #print file2beerased[0]
        os.remove(file2beerased[0])
    return modeled_Otemperatures, modeled_lines

def find_Itheo_in_Iobs(IDobs, Iobs, Iobserr, IDmod, Imod):
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

def lnlikehd(new_Iobs, new_Iobserr, new_Imod):
    model = np.array(new_Imod)
    y = np.array(new_Iobs)
    e = np.array(new_Iobserr)
    csq = ( y - model )**2 / e**2
    Chi2 = csq.sum()
    return -Chi2 / 2.0

def lnpriors(theta):
    ''' 
    This function encodes any previous knowledge that we have about the parameters: results from other experiments, 
    physically acceptable ranges, etc. It is necessary that you write down priors if you are going to use MCMC because 
    all that MCMC does is draw samples from a probability distribution and you want that to be a probability distribution
    for your parameters. We will use 2 functions: Gaussian and top-hat.
    '''
    He, O, CoverO, NoverO, NeoverO, SoverO = theta
    print 'theta = He, O, C/O, N/O, Ne/O, S/O = ', theta
    # Set the abundances set to that of the observed values. Since the boundary conditions are already built-in Cloudy,
    # we will allow them to vary in a wide range. The abundances dictionary is expected to have 6 elements:
    he = lntophat(He, 9.5, 12.0) 
    o = lntophat(O, 7.5, 8.7) 
    c = lntophat(CoverO, -1.6, 1.7)
    n = lntophat(NoverO, -1.7, -0.6) 
    ne = lntophat(NeoverO, -1.0, 0.01) 
    s = lntophat(SoverO, -2.3, -1.4) 
    print 'top hat results:', he , c , n , o , ne , s
    # check that all conditions are met
    if he != -np.inf and c != -np.inf and n != -np.inf and o != -np.inf and ne != -np.inf and s != -np.inf:
        return he + c + n + o + ne + s 
    else:
        return -np.inf

def lnprob(theta, measured_lines, initial_Cloudy_conditions, mod_temps):
    lp = lnpriors(theta)
    if lp != -np.inf:
        # Get the benchmark intensity lines
        IDobs, Iobs, Iobserr, _, _ = measured_lines
        print ' ---> Running model number: ', len(mod_temps)+1
        # Get the Cloudy modeled lines
        modeled_Otemperatures, modeled_lines = get_modeled_lines(theta, initial_Cloudy_conditions)
        mod_temps.append(modeled_Otemperatures)
        mod_TO3, mod_TO2 = modeled_Otemperatures
        print 'model_Te_O3 =', mod_TO3, 'model_Te_O2 =', mod_TO2
        IDmod, _, Imod = modeled_lines
        # Get the relevant lines from the measurements and the model
        new_Iobs, new_Iobserr, new_Imod = find_Itheo_in_Iobs(IDobs, Iobs, Iobserr, IDmod, Imod)
        # Get the Chi2 and add it to the priors' probability 
        return lp + lnlikehd(new_Iobs, new_Iobserr, new_Imod)
    else:
        return lp
    

def new_lnprob(theta, object_name, manual_measurement, init_Cldy_conds_file, mod_temps, threads, chain_file):
    #print '*** theta = ', theta
    #print '*** init_Cldy_conds_file = ', init_Cldy_conds_file
    #print '*** mod_temps = ', mod_temps
    lp = lnpriors(theta)
    if lp != -np.inf:
        print ' * Running model for these abundances: ', theta
        # Read the files with the calculated Chi2 
        lnChi2, TO3, TO2 = get_lnlike_from_file(theta, object_name)
        # if that theta is not in the files already created make one
        if lnChi2 == 'create_job':
            print 'I AM CREATING A JOB...'
            lnChi2, TO3, TO2 = get_new_lnChi2(theta, object_name, manual_measurement, init_Cldy_conds_file)
            modeled_Otemperatures = [TO3, TO2]
            mod_temps.append(modeled_Otemperatures)
            # Store the chain.... 
            f = open(chain_file, "a")
            he, o, co, no, neo, so = theta
            print >> f, "{:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<15} {:<15} {:<20.3f}".format(he, o, co, no, neo, so, TO3, TO2, lnChi2)
            f.close()
            if threads > 1:
                file_with_main_counter = '/home/pena/Documents/AptanaStudio3/carbon/cloudy_tests/pycloudy/'+object_name+'_chaincounter.txt'
                fcount = open(file_with_main_counter, 'a')
                print >> fcount, '{:<15} {:<15} {:<20}'.format(TO3, TO2, 1)
                fcount.close()
                print ' ---> Ran test walker.'
            else:
                print ' ---> Ran model number: ', len(mod_temps)+1
            return lp + lnChi2
        else:
            modeled_Otemperatures = [TO3, TO2]
            mod_temps.append(modeled_Otemperatures)
            # Store the chain.... 
            f = open(chain_file, "a")
            he, o, co, no, neo, so = theta
            print >> f, "{:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<15} {:<15} {:<20.3f}".format(he, o, co, no, neo, so, TO3, TO2, lnChi2)
            f.close()
            print ' ---> Ran model number: ', len(mod_temps)+1
            return lp + lnChi2
    return lp

#### functions needed for the new_lnprob function
def create_HTConjob(jobname, positions, unique_names_list, object_name, manual_measurement, init_Cldy_conds_file, single_job=None, rint=None):
    j = open(jobname, 'w+')
    #print positions, unique_names_list
    print >> j, '# CCC test'
    if single_job != None:
        print >> j, 'Name = Chi2_'+object_name+'_'+rint
    print >> j, 'Universe =  vanilla'
    #print >> j, 'Notification = NEVER'
    print >> j, 'Requirements = (Machine != "science1.stsci.edu")'
    print >> j, 'Getenv = True'
    print >> j, 'Executable = ./likelihood.py'
    print >> j, 'Log = logs/$(Name).log'#_$(Cluster).log'
    print >> j, 'Error = logs/$(Name)_$(Cluster)_error.log'
    print >> j, 'Output = logs/$(Name)_$(Process).out'
    print >> j, 'Should_transfer_files = IF_NEEDED'
    print >> j, 'When_to_transfer_output = ON_EXIT'
    for theta, unique_filename in zip(positions, unique_names_list):
        he, o, co, no, neo, so = theta
        #print 'THESE ARE THE ARGUMENTS FOR CONDOR:'
        #print 'Arguments = '+str(he)+' '+str(o)+' '+str(co)+' '+str(no)+' '+str(neo)+' '+str(so)+' '+str(unique_filename)+' '+object_name+' '+str(manual_measurement)+' '+init_Cldy_conds_file
        if single_job == None:
            print >> j, 'Name = '+object_name+'_'+str(unique_filename)
        print >> j, 'Arguments = '+str(he)+' '+str(o)+' '+str(co)+' '+str(no)+' '+str(neo)+' '+str(so)+' '+str(unique_filename)+' '+object_name+' '+str(manual_measurement)+' '+init_Cldy_conds_file
        print >> j, 'Queue\n'
    j.close()

#def create_END_HTCondor_DAG():
#    jobfile = os.path.abspath('line_run.dag.halt')
#    j = open(jobfile, 'w+')
#    j.close()

def mk_uniquenames_list(nwalkers):
    newunique_names_list = []
    for i in range(nwalkers):
        unique_filename = uuid.uuid4()
        newunique_names_list.append(unique_filename)
    return newunique_names_list

def send_condor_job(job_name, object_name, rint):
    '''This function is used ONLY for individual jobs, i.e. test walkers.'''
    args = shlex.split('condor_submit '+job_name)
    subproc = subprocess.Popen(args)
    subproc.communicate()
    subproc.wait()
    #args2 = str('condor_wait -echo -debug logs/Chi2_'+object_name+'_'+rint+'.log').split()
    args2 = str('condor_wait logs/Chi2_'+object_name+'_'+rint+'.log').split()
    #print(args2)
    subproc2 = subprocess.Popen(args2, stdout=subprocess.PIPE)
    sb2_out, sb2_err = subproc2.communicate()
    #print(sb2_out)
    #print(sb2_err)
    sb2_retval = subproc2.wait()
    #print("condor_wait exited with value {}".format(sb2_retval))

def get_lnlike_from_file(theta, object_name):
    he, o, co, no, neo, so = theta
    unique_names_list = glob(os.path.abspath('Chi2_'+object_name+'*.txt'))
    He_files = []
    O_files = []
    CO_files = []
    NO_files = []
    NeO_files = []
    SO_files = []
    TO3_files = []
    TO2_files = []
    Chi2list = []
    for unique_filename in unique_names_list:
        #print 'unique_filename = ', unique_filename
        try:
            hef, of, cof, nof, neof, sof, lnChi2f, to3f, to2f = np.loadtxt(str(unique_filename), unpack=True)
        except ValueError:
            uf = open(str(unique_filename), 'r')
            line_idx = 1
            for line in uf.readlines():
                kk = string.split(line)
                if line_idx == 1:
                    hef = float(kk[0])
                    of = float(kk[1])
                    cof = float(kk[2])
                    nof = float(kk[3])
                    neof = float(kk[4])
                    sof = float(kk[5])
                    lnChi2f = float(kk[6]) 
                    to3f = float(kk[7]) 
                elif line_idx == 2:
                    to2f = float(kk[0])
                line_idx = line_idx + 1 
            uf.close()
        He_files.append(hef)
        O_files.append(of)
        CO_files.append(cof)
        NO_files.append(nof)
        NeO_files.append(neof)
        SO_files.append(sof)
        Chi2list.append(lnChi2f)
        TO3_files.append(to3f)
        TO2_files.append(to2f)
    if (he in He_files) and (o in O_files) and (co in CO_files) and (no in NO_files):
        idx = He_files.index(he)
        lnChi2 = Chi2list[idx]
        TO3 = TO3_files[idx]
        TO2 = TO3_files[idx]
        #print 'theta = ', theta, '  corresponding Chi2 = ', lnChi2
        #print '  with temperatures of:   TO3 =', TO3, '  TO2 =', TO2
    else:
        lnChi2 = 'create_job'
        TO3 = 0
        TO2 = 0
    return lnChi2, TO3, TO2

def send_and_wait4models(pos_matrix, unique_names_list, object_name, manual_measurement, init_Cldy_conds_file):
    jobname = 'stepCldymods_'+object_name+'.job'
    create_HTConjob(jobname, pos_matrix, unique_names_list, object_name, manual_measurement, init_Cldy_conds_file, 
                    single_job=None, rint=None)
    # Send job and wait for it to come back...
    args = shlex.split('condor_submit '+jobname)
    subproc = subprocess.Popen(args)
    subproc.communicate()
    subproc.wait()
    for uname in unique_names_list:
        args2 = str('condor_wait logs/'+object_name+'_'+str(uname)+'.log').split()
        #print(args2)
        subproc2 = subprocess.Popen(args2, stdout=subprocess.PIPE)
        sb2_out, sb2_err = subproc2.communicate()
        #print(sb2_out)
        #print(sb2_err)
        sb2_retval = subproc2.wait()

def run_position_matrix_models(nwalkers, pos_matrix, object_name, manual_measurement, init_Cldy_conds_file):
    # make a new unique file name list so that files can be erased without loosing any other info
    unique_names_list = []
    for i in range(nwalkers):
        unique_filename = uuid.uuid4()
        unique_names_list.append(unique_filename)       
    # Write the job file for the next set of parallel model runs in groups of 50 jobs max
    if nwalkers > 50:
        # divide the unique names list into groups of 50
        groups_of_uniquienames = []
        njobs = int(nwalkers / 50)
        iniat = 0
        endat = 50
        for nj in range(njobs):
            g = unique_names_list[iniat:endat]
            groups_of_uniquienames.append(g)
            iniat = iniat + 50
            endat = endat + 50 
        for groupof50 in groups_of_uniquienames:
            send_and_wait4models(pos_matrix, groupof50, object_name, manual_measurement, init_Cldy_conds_file)
    else:
        send_and_wait4models(pos_matrix, unique_names_list, object_name, manual_measurement, init_Cldy_conds_file)

def get_new_lnChi2(theta, object_name, manual_measurement, init_Cldy_conds_file):
    #print 'ok, got to obtain new Chi2... Now I will create the job, send it and wait for it...'
    ndim = len(theta)
    unique_name = mk_uniquenames_list(1)
    next_pos = np.array([]).reshape(0, ndim)
    next_pos = np.vstack((next_pos, theta))
    str_randint = str(np.random.randint(0, 1000))
    #print 'this is the random integer string:', str_randint
    jobname = os.path.abspath('test_walker_'+object_name+'_'+str_randint+'.job')
    create_HTConjob(jobname, next_pos, unique_name, object_name, manual_measurement, init_Cldy_conds_file, single_job=True, rint=str_randint)
    send_condor_job(jobname, object_name, str_randint)
    got_file_back_from_condor = False
    Chifile = 'Chi2_'+object_name+'_'+str(unique_name[0])+'.txt'
    while got_file_back_from_condor != True:
       time.sleep(0.025)
       if os.path.isfile(Chifile):
           got_file_back_from_condor = True
           #print 'ok, condor job done, reading the Chi2 file...', Chifile
    he, o, co, no, neo, so, lnChi2, to3, to2 = np.loadtxt(Chifile, unpack=True)
    #print theta, '... from Chi2 file temperatures are:  TO3 =', to3, '  TO2 =', to2, '  lnChi2 = ', lnChi2
    return lnChi2, to3, to2

def clean_directory(object_name):
    unique_names_list = glob(os.path.abspath('Chi2_'+object_name+'*.txt'))
    for unique_filename in unique_names_list:
        os.remove(unique_filename)
    os.system('rm -f test_walker*'+object_name+'*.job')
    os.system('rm -f logs/*'+object_name+'*_error.log')
    os.system('rm -f logs/*'+object_name+'*.out')

def run_chain41step(chain_file, pos_matrix, nwalkers, ndim, new_lnprob, object_name, manual_measurement, init_Cldy_conds_file, mod_temps, threads):
    # With the new probability function that allows to use HTCondor
    sampler = emcee.EnsembleSampler(nwalkers, ndim, new_lnprob, args=[object_name, manual_measurement, init_Cldy_conds_file, 
                                                                      mod_temps, threads, chain_file],
                                    threads=threads)#pool=pool)
    next_pos_matrix, prob, rstate = sampler.run_mcmc(pos_matrix, 1)
    '''
    # Store the chain....
    #chain_file = os.path.abspath(dir+model_name+"_chain0.dat")
    f = open(chain_file, "w+")
    # make sure there is space for the header information to be added at the end
    lines2skip = 16
    for l in range(lines2skip):
        print >> f, '#'
    f.close()
    count = 1
    for posn, prob, state in sampler.sample( pos, iterations=1, storechain=False ):
        #print "COUNT", count
        if count % 1 == 0:
            f = open(chain_file, "w+")
            for k in range( posn.shape[0] ):
                strout = ""
                for p in pos[k]: strout += "{:<8.3f} ".format( p )
                TO3, TO2 = mod_temps[k]
                strout += "{:<15} {:<15}".format( float(TO3), float(TO2) )
                strout += "{:<20.3f}".format( prob[k] )
                #print strout
                f.write(strout+"\n")
            f.close()
        count += 1
    #print 'The chain was saved in:', chain_file
    '''
    return next_pos_matrix, mod_temps
    
####

def matrix_lastpos(He, O, CO, NO, NeO, SO, nwalkers):
    last_pos_idx = len(He)-nwalkers
    He = He[last_pos_idx:]
    O = O[last_pos_idx:]
    CO = CO[last_pos_idx:]
    NO = NO[last_pos_idx:]
    NeO = NeO[last_pos_idx:]
    SO = SO[last_pos_idx:]
    pos_matrix = np.array([]).reshape(0, 6)
    for he, o, co, no, neo, so in zip(He, O, CO, NO, NeO, SO):
        theta = np.array([he, o, co, no, neo, so])
        pos_matrix = np.vstack((pos_matrix, theta))
    return pos_matrix

def read_chain4starting_posmatrix(chainfile, object_name, nwalkers):
    #chainfile = os.path.abspath('mcmc_'+object_name+'_chain0.dat')
    f = open(chainfile, 'r+')
    TO3_list = []
    TO2_list = []
    Chi2_list = []
    i = 0
    for line in f.readlines():
        ''' skip the first 16 lines of the file, which contain the previous running time, benchmark abundances, best
        chi2 fitted model, and the uncertainties associated with that particular model.'''
        if i == 0:
            prev_run_time_string = line.strip()
            prev_run_time_list = string.split(prev_run_time_string)#, sep=" ")
            prev_run_time = float(prev_run_time_list[3])
        # skip the rest of the info until the actual chain
        elif i > 16:
        # Check that if the file has strings in the main body
            if '[' or '\n' in line:
                abunds_temps_chi2_string = line.strip()
                abunds_temps_chi2_list = string.split(abunds_temps_chi2_string, sep="[")
                temps_chi2_list = string.split(abunds_temps_chi2_list[1], sep="\\n")
                to3 = float(string.replace(temps_chi2_list[0], "' ", ""))
                to2 = float(string.replace(temps_chi2_list[1], "', ' ", ""))
                chi = float(string.replace(temps_chi2_list[2], "']", ""))
                TO3_list.append(to3)
                TO2_list.append(to2)
                Chi2_list.append(chi)
            else:
                He, O, CO, NO, NeO, SO, TO3, TO2, Chi2 = np.loadtxt(chainfile, skiprows=17, unpack=True)
                # last positions
                pos_matrix = matrix_lastpos(He, O, CO, NO, NeO, SO, nwalkers)
                return prev_run_time, pos_matrix, He, O, CO, NO, NeO, SO, TO3, TO2, Chi2
        i = i + 1
    f.close()
    # Turn lists into numpy arrays
    TO3 = np.array(TO3_list)
    TO2 = np.array(TO2_list)
    Chi2 = np.array(Chi2_list)
    # Now read the abundances part of the file
    He, O, CO, NO, NeO, SO = np.loadtxt(chainfile, skiprows=17, usecols=(0,1,2,3,4,5), unpack=True)
    # last positions
    pos_matrix = matrix_lastpos(He, O, CO, NO, NeO, SO, nwalkers)
    return prev_run_time, pos_matrix, He, O, CO, NO, NeO, SO, TO3, TO2, Chi2


# with the regular probability function
#def run_chain_and_plot(model_name, dir, true_abunds, theta, nwalkers, nruns, measured_lines, initial_Cloudy_conditions, mod_temps):
# With the new probability function that allows to use HTCondor
def run_chain_and_plot(model_name, dir, true_abunds, theta, nwalkers, nruns, object_name, manual_measurement, init_Cldy_conds_file, mod_temps, threads=1, recover=False):
    # start the timer to compute the whole running time
    start_time = time.time()
    
    #ndim, nwalkers, nruns = 6, 100, 100
    ndim = len(theta)
    randadd2point = lambda x: x+np.random.rand(1)
    pos_matrix = [[float(randadd2point(x)) for x in true_abunds] for i in range(nwalkers)] 
    '''
    pool = MPIPool()
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
    '''
    # Determine object name
    object_name_list = string.split(model_name, sep="_")
    object_name = object_name_list[1]
    
    # Name of the file to store the chain
    chain_file = os.path.abspath(dir+model_name+"_chain0.dat")
    
    # If there was a previous chain ran, start from there
    if recover:
        print '\nUsing previous position matrix for next Cloudy models...'
        prev_run_time, pos_matrix, He, O, CO, NO, NeO, SO, TO3, TO2, Chi2 = read_chain4starting_posmatrix(chain_file, object_name, nwalkers)
        # create a new chain file and write the previous info
        f = open(chain_file, "w+")
        # make sure there is space for the header information to be added at the end
        lines2skip = 16
        for l in range(lines2skip):
            print >> f, '#'
        # Store the chain.... print all the previous steps 
        for he, o, co, no, neo, so, to3, to2, chi2 in zip(He, O, CO, NO, NeO, SO, TO3, TO2, Chi2):
            print >> f, "{:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<15} {:<15} {:<20.3f}".format(he, o, co, no, neo, so, to3, to2, chi2)
        f.close()
    else:
        if len(mod_temps) == 0:   # so that it creates the file only once, at the beginning. 
            f = open(chain_file, "w+")
            # make sure there is space for the header information to be added at the end
            lines2skip = 16
            for l in range(lines2skip):
                print >> f, '#'
            f.close()
       
    # Run MCMC with HTCondor...
    # Create file to keep track of runs if threads is more than 1
    if threads > 1:
        file_with_main_counter = '/home/pena/Documents/AptanaStudio3/carbon/cloudy_tests/pycloudy/'+object_name+'_chaincounter.txt'
        fcount = open(file_with_main_counter, 'w')
        fcount.close()
    
    # Create the jobs, send them, and wait for them to come back
    run_position_matrix_models(nwalkers, pos_matrix, object_name, manual_measurement, init_Cldy_conds_file)
    for main_counter in range(nruns):
        print '\n --->  Starting step number: ', main_counter
        # read Chi2 from text file and generate new position matrix 
        pos_matrix, mod_temps = run_chain41step(chain_file, pos_matrix, nwalkers, ndim, new_lnprob, object_name, manual_measurement, 
                                     init_Cldy_conds_file, mod_temps, threads)
        # create the jobs for next step, send them, and wait for them to come back
        run_position_matrix_models(nwalkers, pos_matrix, object_name, manual_measurement, init_Cldy_conds_file)
        # Clean...
        clean_directory(object_name)
        
    '''
    # with the regular probability function
    #sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[measured_lines, initial_Cloudy_conditions, mod_temps])
                                    #threads=2)#pool=pool)
    pos, prob, rstate = sampler.run_mcmc(p0, 1)
    '''
    # Store the chain....
    time2run = 'Chain finished! Took  %s  seconds to finish.' % (time.time() - start_time)
    meas_lineIDs, _, _, _, _ = measured_lines
    lines4chi = 'Used %i lines to determine Chi2.' % len(meas_lineIDs)
    # best model
    wh = np.where( prob == prob.max() )[0][0]
    p = pos[ wh, : ]
    TOs = mod_temps[wh]
    temps = 'model_Te_O3 = %s     model_Te_O2 = %s' % (TOs[0], TOs[1])
    
    trueabs0 = 'Values of the BENCHMARK abundances:\n'
    trueabs1 = 'He = %0.2f   O = %0.2f   C/O = %0.2f   N/O = %0.2f   Ne/O = %0.2f   S/O = %0.2f' % (true_abunds[0], true_abunds[1], true_abunds[2], 
                                                                                                    true_abunds[3], true_abunds[4], true_abunds[5])
    trueabs2= 'He = %0.2f   O = %0.2f   C = %0.2f   N = %0.2f   Ne = %0.2f   S = %0.2f' % (true_abunds[0], true_abunds[1], true_abunds[2]+true_abunds[1],
                                                                                           true_abunds[3]+true_abunds[1], true_abunds[4]+true_abunds[1], 
                                                                                           true_abunds[5]+true_abunds[1])
    line1 = 'Values of the %i dimensions that best fit the data in %i runs with %i walkers, are the following:' % (ndim, nruns, nwalkers)
    line2 = 'He = %0.2f   O = %0.2f   C/O = %0.2f   N/O = %0.2f   Ne/O = %0.2f   S/O = %0.2f' % (p[0], p[1], p[2], p[3], p[4], p[5])
    line3 = 'He = %0.2f   O = %0.2f   C = %0.2f   N = %0.2f   Ne = %0.2f   S = %0.2f' % (p[0], p[1], p[2]+true_abunds[1], 
                                                                                         p[3]+true_abunds[1], p[4]+true_abunds[1],
                                                                                         p[5]+true_abunds[1])
    # reorder the sampler chain in order to plot
    #samples = sampler.chain[:, nruns*0.2:, :].reshape((-1, ndim)) # this discards the first 20% of the runs
    # samples has an attribute called chain, which is an array with shape (nwalkers, nruns, ndim)
    # but for now we'll use all the runs:
    samples = sampler.chain[:, :, :].reshape((-1, ndim))
    # Calculate the uncertainties based on the 16th, 50th and 84th percentiles
    samples[:, ndim-1] = np.exp(samples[:, ndim-1])
    percentiles = 'mcmc values and uncertainties according to 16th, 50th, and 84th percentiles:'
    p_mcmc1 = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    p_mcmc2 = map(lambda v: (v), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    
    # write statistics in final chain file...
    final_chain_file = os.path.abspath(dir+model_name+"_chain.dat")    
    info = [time2run, trueabs0, trueabs1, trueabs2, line1, line2, line3, temps, lines4chi, percentiles, p_mcmc2, p_mcmc1]
    f0 = open(chain_file, "r")
    f1 = open(final_chain_file, "w+")
    linfo = len(info)
    idx = 0
    for line in f1.readlines():
        if idx < lena:
            print >> f0, str(info[idx])
        else:
            print >> f1, line.strip()     
        idx = idx+1      
    f0.close()
    f1.close()
    os.remove(chain_file)

    print '\n'
    print time2run
    print temps
    print lines4chi
    print trueabs0
    print trueabs1
    print trueabs2
    print line1
    print line2
    print line3
    print percentiles
    print p_mcmc1
    print p_mcmc2
    
    # PLOTS
    # plot of abundance ratios including the benchmark abundances
    fig = triangle.corner(samples, labels=["$He$", "$O$", "$C/O$", "$N/O$", "$Ne/O$", "$S/O$"], 
                          truths=[true_abunds[0], true_abunds[1], true_abunds[2], true_abunds[3],
                                  true_abunds[4], true_abunds[5]])
    fig.savefig(os.path.abspath(dir+model_name+"_ratios_"+repr(nwalkers)+"w"+repr(nruns)+"r"+"_initb.jpg"))
    # plot of the ratios without the benchmark abundances
    fig = triangle.corner(samples, labels=["$He$", "$O$", "$C/O$", "$N/O$", "$Ne/O$", "$S/O$"])
    fig.savefig(os.path.abspath(dir+model_name+"_ratios2_"+repr(nwalkers)+"w"+repr(nruns)+"r"+"_initb.jpg"))

        
