#!/usr/bin/env python
from __future__ import division
import os
import numpy as np
import matplotlib.pyplot as plt
import time
import copy
import string
import emcee
import triangle
from glob import glob
from science import spectrum
import scipy.optimize as op
import types
import uuid
import shlex
import subprocess

''' 
This script contains all the functions to serve as infrastructure for running a MCMC chain. 
'''

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
        return -np.inf


def get_measured_lines(object_name, manual_measurement, debug=False):
    if debug:
        print 'Entering get_measured_lines function...'
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
    if debug:
        print 'get_measured_lines function excited ok!'
    return measured_lines


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


def lnpriors(theta, debug=False):
    ''' 
    This function encodes any previous knowledge that we have about the parameters: results from other experiments, 
    physically acceptable ranges, etc. It is necessary that you write down priors if you are going to use MCMC because 
    all that MCMC does is draw samples from a probability distribution and you want that to be a probability distribution
    for your parameters. We will use 2 functions: Gaussian and top-hat.
    '''
    if debug:
        print 'Entering lnpriors function...'
    He, O, CoverO, NoverO, NeoverO, SoverO = theta
    # Set the abundances set to that of the observed values. Since the boundary conditions are already built-in Cloudy,
    # we will allow them to vary in a wide range. The abundances dictionary is expected to have 6 elements:
    he = lntophat(He, 9.5, 12.0) 
    o = lntophat(O, 7.5, 8.7) 
    c = lntophat(CoverO, -1.6, 1.7)
    n = lntophat(NoverO, -1.7, -0.4) 
    ne = lntophat(NeoverO, -1.0, 0.01) 
    s = lntophat(SoverO, -2.3, -1.4) 
    if debug:
        print 'theta = ', theta
        print 'top hat results:', he , c , n , o , ne , s
    # check that all conditions are met
    if he != -np.inf and c != -np.inf and n != -np.inf and o != -np.inf and ne != -np.inf and s != -np.inf:
        return he + c + n + o + ne + s 
    else:
        return -np.inf
    
def new_lnprob(theta, object_name, manual_measurement, init_Cldy_conds_file, debug):
    if debug:
        print 'Entering new_lnprob function...'
        print '*** theta = ', theta
        print '*** init_Cldy_conds_file = ', init_Cldy_conds_file
    lp = lnpriors(theta)
    if lp != -np.inf:
        print ' * Running model for these abundances: ', theta
        # Read the files with the calculated Chi2 
        lnChi2, TO3, TO2 = get_lnlike_from_file(theta, object_name, debug)
        # if that theta is not in the files already created make one
        if lnChi2 == 'create_job':
            if debug:
                print 'I AM CREATING A JOB...'
            lnChi2, TO3, TO2 = get_new_lnChi2(theta, object_name, manual_measurement, init_Cldy_conds_file, debug)
            print ' ---> Ran test walker, probability recorded. '
        return lp + lnChi2
    return lp

'''
def store_chain(object_name, chain_file, debug=False):
    # Store the chain....
    if debug:
        print 'Entering store_chain function...'
    He_files, O_files, CO_files, NO_files, NeO_files, SO_files, TO3_files, TO2_files, Chi2list = get_abstempschis(object_name, debug=debug)
    if debug:
        print 'lengths of lists to append to chain_file:'
        print len(He_files), len(O_files), len(CO_files), len(NO_files), len(NeO_files), len(SO_files), len(TO3_files), len(TO2_files), len(Chi2list)
    f = open(chain_file, "a")
    for he, o, co, no, neo, so, to3, to2, chi2 in zip(He_files, O_files, CO_files, NO_files, NeO_files, SO_files, TO3_files, TO2_files, Chi2list):
        if to3 and to2 != 0:
            print >> f, "{:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<15} {:<15} {:<20.3f}".format(he, o, co, no, neo, so, to3, to2, chi2)
        #strout = "{:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<15} {:<15} {:<20.3f}".format(he, o, co, no, neo, so, to3, to2, chi2)
        #f.write(strout+'\n')
    f.close()
'''
def store_chain(object_name, pos_matrix, chain_file, debug=False):
    # Store the chain.... 
    if debug:
        print 'Entering store_chain function...'
    He_files, O_files, CO_files, NO_files, NeO_files, SO_files = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
    for set in pos_matrix:
        He_files = np.append(He_files, set[0])
        O_files = np.append(O_files, set[1])
        CO_files = np.append(CO_files, set[2])
        NO_files = np.append(NO_files, set[3])
        NeO_files = np.append(NeO_files, set[4])
        SO_files = np.append(SO_files, set[5])
    TO3_files, TO2_files, Chi2list = get_from_files(object_name, pos_matrix, debug=debug)
    if debug:
        print 'lengths of lists to append to chain_file:'
        print len(He_files), len(O_files), len(CO_files), len(NO_files), len(NeO_files), len(SO_files), len(TO3_files), len(TO2_files), len(Chi2list)
    f = open(chain_file, "a")
    for he, o, co, no, neo, so, to3, to2, chi2 in zip(He_files, O_files, CO_files, NO_files, NeO_files, SO_files, TO3_files, TO2_files, Chi2list):
        if to3 and to2 != 0:
            print >> f, "{:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<15} {:<15} {:<20.3f}".format(he, o, co, no, neo, so, to3, to2, chi2)
    f.close()
    
'''def get_from_files(object_name, debug=False):
    TO2, TO3 = [], []
    He_files, O_files, CO_files, NO_files, NeO_files, SO_files, TO3_files, TO2_files, Chi2list = get_abstempschis(object_name, debug=debug)
    for he, o, co, no, neo, so, to3, to2, chi in zip(He_files, O_files, CO_files, NO_files, NeO_files, SO_files, TO3_files, TO2_files, Chi2list):
        set = [he, o, co, no, neo, so]
        TO2 = np.append(TO2, to2)
        TO3 = np.append(TO3, to3)
    return TO3, TO2'''
def get_from_files(object_name, pos_matrix, debug=False):
    TO2, TO3, probs = np.array([]), np.array([]), np.array([])
    He_files, O_files, CO_files, NO_files, NeO_files, SO_files, TO3_files, TO2_files, Chi2list = get_abstempschis(object_name, debug=debug)
    i = 0
    for he, o, co, no, neo, so, to3, to2, chi in zip(He_files, O_files, CO_files, NO_files, NeO_files, SO_files, TO3_files, TO2_files, Chi2list):
        set = [he, o, co, no, neo, so]
        if set in pos_matrix:
            TO2 = np.append(TO2, TO2_files[i])
            TO3 = np.append(TO3, TO3_files[i])
            probs = np.append(probs, Chi2list[i])
        i = i + 1
    print 'Lengths of temperatures and likelihood arrays: ', len(TO3), len(TO2), len(probs)
    return TO3, TO2, probs

def create_HTConjob(jobname, positions, unique_names_list, object_name, manual_measurement, init_Cldy_conds_file, single_job=None, rint=None, debug=False):
    if debug:
        print 'Entering createHTConjob function...'
    j = open(jobname, 'w+')
    print >> j, '# CCC test'
    if single_job != None:
        print >> j, 'Name = Chi2_'+object_name+'_'+rint
    print >> j, 'Universe =  vanilla'
    print >> j, 'Notification = NEVER'
    print >> j, 'Requirements =  ( (Machine != "science3.stsci.edu") || (Machine != "science1.stsci.edu") )'
    #print >> j, 'Requirements =  (Machine == "science5.stsci.edu" )'
    print >> j, 'Getenv = True'
    print >> j, 'Executable = ./likelihood.py'
    print >> j, 'Log = logs/$(Name).log'#_$(Cluster).log'
    print >> j, 'Error = logs/$(Name)_$(Cluster)_error.log'
    print >> j, 'Output = logs/$(Name)_$(Process).out'
    print >> j, 'Should_transfer_files = IF_NEEDED'
    print >> j, 'When_to_transfer_output = ON_EXIT'
    for theta, unique_filename in zip(positions, unique_names_list):
        he, o, co, no, neo, so = theta
        if debug:
            print 'THESE ARE THE ARGUMENTS FOR CONDOR for job: ', object_name+'_'+str(unique_filename)
            print 'Arguments = '+str(he)+' '+str(o)+' '+str(co)+' '+str(no)+' '+str(neo)+' '+str(so)+' '+str(unique_filename)+' '+object_name+' '+str(manual_measurement)+' '+init_Cldy_conds_file
        if single_job == None:
            print >> j, 'Name = '+object_name+'_'+str(unique_filename)
        print >> j, 'Arguments = '+str(he)+' '+str(o)+' '+str(co)+' '+str(no)+' '+str(neo)+' '+str(so)+' '+str(unique_filename)+' '+object_name+' '+str(manual_measurement)+' '+init_Cldy_conds_file
        print >> j, 'Queue\n'
    j.close()


def mk_uniquenames_list(nwalkers):
    newunique_names_list = []
    for i in range(nwalkers):
        unique_filename = uuid.uuid4()
        newunique_names_list.append(unique_filename)
    return newunique_names_list


def send_condor_job(job_name, object_name, rint, debug=False):
    '''This function is used ONLY for individual jobs, i.e. test walkers.'''
    if debug:
        print 'Entering send_condor_job function ...'
        print 'will send HTCondor job and wait for text file'
    #subproc = os.system('condor_submit '+job_name)
    args = shlex.split('condor_submit '+job_name)
    subproc = subprocess.call(args)
    args2 = str('condor_wait logs/Chi2_'+object_name+'_'+rint+'.log').split()
    subproc2 = subprocess.call(args2, stdout=subprocess.PIPE)
    

def get_abstempschis(object_name, debug=False):
    ''' This function gets all the abundances, probabilities, and temperatures from the text files. '''
    if debug:
        print 'Entering get_abstempschis function...'
        print 'obtaining information from all of the Chi2 text files for object: ', object_name
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
    if debug:
        print 'lengths of He_files, O_files, CO_files, NO_files, NeO_files, SO_files, TO3_files, TO2_files, Chi2list: '
        print len(He_files), len(O_files), len(CO_files), len(NO_files), len(NeO_files), len(SO_files), len(TO3_files), len(TO2_files), len(Chi2list)
    return He_files, O_files, CO_files, NO_files, NeO_files, SO_files, TO3_files, TO2_files, Chi2list


def get_lnlike_from_file(theta, object_name, debug=False):
    ''' This function gets the probability from text file if it exists, if not then orders a job. '''
    if debug:
        print 'Entering get_lnlike_from_file function...'
    he, o, co, no, neo, so = theta
    He_files, O_files, CO_files, NO_files, NeO_files, SO_files, TO3_files, TO2_files, Chi2list = get_abstempschis(object_name)
    if (he in He_files) and (o in O_files) and (co in CO_files) and (no in NO_files):
        idx = He_files.index(he)
        lnChi2 = Chi2list[idx]
        TO3 = TO3_files[idx]
        TO2 = TO3_files[idx]
    else:
        lnChi2 = 'create_job'
        TO3 = 0
        TO2 = 0
    if debug:
        print 'obtained lnChi2, TO3, and TO2 from text file: ', lnChi2, TO3, TO2 
    return lnChi2, TO3, TO2


def send_and_wait4models(pos_matrix, unique_names_list, object_name, manual_measurement, init_Cldy_conds_file, debug=False):
    ''' This function sends the job and waits for files for the multiple models at the beginning of the run. '''
    if debug:
        print 'Entering send_and_wait4models function...'
    jobname = 'stepCldymods_'+object_name+'.job'
    create_HTConjob(jobname, pos_matrix, unique_names_list, object_name, manual_measurement, init_Cldy_conds_file, 
                    single_job=None, rint=None, debug=debug)
    # Send job and wait for it to come back...
    args = shlex.split('condor_submit '+jobname)
    subproc = subprocess.call(args)
    for uname in unique_names_list:
        args2 = str('condor_wait logs/'+object_name+'_'+str(uname)+'.log').split()
        subproc2 = subprocess.call(args2, stdout=subprocess.PIPE)
    if debug:
        print 'Created and sent job for multiple models, and waited for corresponding text files.'


def run_position_matrix_models(nwalkers, pos_matrix, object_name, manual_measurement, init_Cldy_conds_file, max_parallel_jobs=25, debug=False):
    ''' This function creates the unique names for the multiple jobs at the beginning of the run, sends the jobs,
    and waits for all the text files expected. All of this in groups of a set maximum number of jobs to prevent
    from using up all resources in the Cluster. '''
    if debug:
        print 'Entering run_position_matrix_models function...'
    # make a new unique file name list so that files can be erased without loosing any other info
    unique_names_list = []
    for i in range(nwalkers):
        unique_filename = uuid.uuid4()
        unique_names_list.append(unique_filename)       
    # Write the job file for the next set of parallel model runs in groups of max_parallel_runs jobs max
    if nwalkers >= max_parallel_jobs:
        # divide the unique names list into groups of multiples of 50
        groups_of_uniquienames = []
        njobs = int(nwalkers / max_parallel_jobs)
        iniat = 0
        endat = max_parallel_jobs
        for nj in range(njobs):
            g = unique_names_list[iniat:endat]
            groups_of_uniquienames.append(g)
            iniat = iniat + max_parallel_jobs
            endat = endat + max_parallel_jobs
        for groupofmaxjobs in groups_of_uniquienames:
            send_and_wait4models(pos_matrix, groupofmaxjobs, object_name, manual_measurement, init_Cldy_conds_file, debug=debug)
    else:
        send_and_wait4models(pos_matrix, unique_names_list, object_name, manual_measurement, init_Cldy_conds_file, debug=debug)


def get_new_lnChi2(theta, object_name, manual_measurement, init_Cldy_conds_file, debug):
    ''' This function is used only for individual jobs for the test walkers. It creates a unique name for the text file, 
    creates a job, sends the job, and waits for corresponding text out file.'''
    if debug:
        print 'Entering get_new_lnChi2 function...'
        print 'ok, got to obtain new Chi2... Now I will create the job, send it and wait for the out text file...'
    ndim = len(theta)
    unique_name = mk_uniquenames_list(1)
    next_pos = np.array([]).reshape(0, ndim)
    next_pos = np.vstack((next_pos, theta))
    str_randint = str(unique_name[0])
    jobname = os.path.abspath('test_walker_'+object_name+'_'+str_randint+'.job')
    create_HTConjob(jobname, next_pos, unique_name, object_name, manual_measurement, init_Cldy_conds_file, 
                    single_job=True, rint=str_randint, debug=debug)
    send_condor_job(jobname, object_name, str_randint, debug=debug)
    path4file = '/home/pena/Documents/AptanaStudio3/carbon/cloudy_tests/pycloudy/'+object_name
    Chifiles_list = glob(path4file+'/Chi2_'+object_name+'*.txt')
    Chifile = 'Chi2_'+object_name+'_'+str(unique_name[0])+'.txt'
    full_path_Chifile = os.path.join(path4file, Chifile)
    time.sleep(5)
    if debug:
       print 'ok, got the file! reading from...', full_path_Chifile
    try:
        he, o, co, no, neo, so, lnChi2, to3, to2 = np.loadtxt(full_path_Chifile, unpack=True)
    except IOError:
        print 'This model failed...'
        lnChi2 = -np.inf  # changed value to infinity from 'create_job'
        to3, to2 = 0, 0
    if debug:
        print theta, '... from Chi2 file temperatures are:  TO3 =', to3, '  TO2 =', to2, '  lnChi2 = ', lnChi2
    return lnChi2, to3, to2


def clean_directory(object_name, debug=False):
    ''' This function makes sure that there is no overflow of files. It is used at the end of each run to clean both the
    pycloudy and logs directories.'''
    if debug:
        print 'Entering clean_directory function...'
    unique_names_list = glob(os.path.abspath('Chi2_'+object_name+'*.txt'))
    for unique_filename in unique_names_list:
        os.remove(unique_filename)
    os.system('rm -f test_walker*'+object_name+'*.job')
    os.system('rm -f logs/*'+object_name+'*.out')
    os.system('rm -f logs/*'+object_name+'*.log')
    # If the Cloudy out file is empty erase it, else move it to the out files directory.
    files_path = "/home/pena/Documents/AptanaStudio3/carbon/cloudy_tests/pycloudy/"+object_name
    files_list = glob(files_path+"/mcmc_"+object_name+"*.out")
    for f in files_list:
        empty_file = os.stat(f).st_size == 0
        if empty_file:
            os.system("rm "+f)
        else:
            os.system("mv "+f+" /user/pena/pycloudy_stuff/"+object_name+"/outfiles/")
    
    
def matrix_lastpos(He, O, CO, NO, NeO, SO, nwalkers):
    ''' This function creates the position matrix of the correct shape when using a previous chain. '''
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
    ''' This function obtains the data to create the previous position matrix to initialize the chain from
    last step of previous run. '''
    try:
        He, O, CO, NO, NeO, SO, TO3, TO2, Chi2 = np.loadtxt(chainfile, skiprows=16, unpack=True)
        # last positions
        pos_matrix = matrix_lastpos(He, O, CO, NO, NeO, SO, nwalkers)
        prev_run_time = 0.0
        
    except ValueError:
        f = open(chainfile, 'r')
        TO3_list = []
        TO2_list = []
        Chi2_list = []
        i = 0
        for line in f.readlines():
            ''' skip the first 16 lines of the file, which contain the previous running time, benchmark abundances, best
            chi2 fitted model, and the uncertainties associated with that particular model.'''
            if i == 0:
                if '#' not in line:
                    prev_run_time_string = line.strip()
                    prev_run_time_list = string.split(prev_run_time_string)
                    prev_run_time = float(prev_run_time_list[3])
            # skip the rest of the info until the actual chain
            elif i > 15:
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
            i = i + 1
        f.close()
        
        # Turn lists into numpy arrays
        TO3 = np.array(TO3_list)
        TO2 = np.array(TO2_list)
        Chi2 = np.array(Chi2_list)
        # Now read the abundances part of the file
        He, O, CO, NO, NeO, SO = np.loadtxt(chainfile, skiprows=16, usecols=(0,1,2,3,4,5), unpack=True)
        # last positions
        pos_matrix = matrix_lastpos(He, O, CO, NO, NeO, SO, nwalkers)
        
    return prev_run_time, pos_matrix, He, O, CO, NO, NeO, SO, TO3, TO2, Chi2


# With the new probability function that allows to use HTCondor
def run_chain_and_plot(model_name, dir, true_abunds, theta, nwalkers, nruns, object_name, manual_measurement, init_Cldy_conds_file, mod_temps, threads=1, recover=False, debug=False):
    ''' This is the function that controls the entire chain. 
    '''
    if debug:
        print 'Entering run_chain_and_plot function...'
        
    # start the timer to compute the whole running time
    start_time = time.time()
    
    # Define number of dimensions
    ndim = len(theta)
    
    # INITIALIZATION OF POSITION MATRIX
    # OPTION 1: initialize with a random addition to the previously know values 
    #randadd2point = lambda x: x+np.random.rand(1)
    #pos_matrix = [[float(randadd2point(x)) for x in true_abunds] for i in range(nwalkers)] 
    # OPTION 2: initialize randomly within allowed values
    pos_matrix = [[np.random.uniform(9.6, 11.9), np.random.uniform(7.6, 8.6), 
                   np.random.uniform(-1.5, 1.6), np.random.uniform(-1.6, -0.5),
                   np.random.uniform(-0.9, 0.009), np.random.uniform(-2.2, -1.3)] for i in range(nwalkers)]

    # Determine object name
    object_name_list = string.split(model_name, sep="_")
    object_name = object_name_list[1]
    
    # Name of the file to store the chain
    chain_file = '/home/pena/Documents/AptanaStudio3/carbon/cloudy_tests/pycloudy/'+object_name+'/'+model_name+'_chain0.dat'
    
    # If there was a previous chain ran, start from there
    restart_from = 0
    if recover:
        if debug:
            print '\nUsing previous position matrix for next Cloudy models...'
        string_step = raw_input('If wanting to restart the chain from a certain step type the number, else press Enter   ')
        if string_step == " ":
            restart_from = 0
        else:
            restart_from = int(string_step) - 1
        prev_run_time, pos_matrix, He, O, CO, NO, NeO, SO, TO3, TO2, Chi2 = read_chain4starting_posmatrix(chain_file, object_name, nwalkers)
        # create a new chain file and write the previous info
        f = open(chain_file, "w")
        # make sure there is space for the header information to be added at the end
        lines2skip = 16
        for l in range(lines2skip):
            print >> f, '#'
        # Store the chain.... print all the previous steps 
        for he, o, co, no, neo, so, to3, to2, chi2 in zip(He, O, CO, NO, NeO, SO, TO3, TO2, Chi2):
            print >> f, "{:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<8.3f} {:<15} {:<15} {:<20.3f}".format(he, o, co, no, neo, so, to3, to2, chi2)
        f.close()
    else:
        f = open(chain_file, "w")
        f.close()
        if debug:
            print 'Created initial position matrix... staring chain...'
        if len(mod_temps) == 0:   # so that it creates the file only once, at the beginning. 
            f = open(chain_file, "a")
            # make sure there is space for the header information to be added at the end
            lines2skip = 16
            for l in range(lines2skip):
                print >> f, '#'
            f.close()
       
    # Run MCMC with HTCondor...
    
    # Create the jobs, send them, and wait for them to come back
    max_parallel_jobs = 25#threads
    run_position_matrix_models(nwalkers, pos_matrix, object_name, manual_measurement, init_Cldy_conds_file, 
                               max_parallel_jobs=max_parallel_jobs, debug=debug)
    
    # Initialization of emcee Ensemble Sampler object
    sampler = emcee.EnsembleSampler(nwalkers, ndim, new_lnprob, args=[object_name, manual_measurement, init_Cldy_conds_file, debug],
                                    threads=threads)

    for main_counter in range(restart_from, nruns):
        print '\n*** --->  Starting step number: ', main_counter+1, '\n'
        # If restarting the chain set counter to desired step, else start from 0
        if debug:
            print 'obtaining next position matrix...'
        pos_matrix, prob, rstate = sampler.run_mcmc(pos_matrix, 1)
        print 'got it! shape of NEXT position matrix: ', np.shape(pos_matrix)
        # Store the chain
        if debug:
            print 'got it! shape of NEXT position matrix: ', np.shape(pos_matrix)
            #print pos_matrix
            print 'ok! Now storing the chain for step number: ', main_counter+1, '\n'
        #store_chain(object_name, chain_file, debug=debug)
        store_chain(object_name, pos_matrix, chain_file, debug=debug)
        
        # Clean...
        if debug:
            print 'Done! Proceed to clean current and logs directories...'
        clean_directory(object_name, debug=debug)
        print 'Directories clean! '

        if debug:
            print 'Now run models for next position matrix of step number: ', main_counter+1, '\n'
        # create the jobs for next step, send them, and wait for them to come back but ONLY if it is not the last run
        if main_counter < nruns:
            run_position_matrix_models(nwalkers, pos_matrix, object_name, manual_measurement, init_Cldy_conds_file, 
                                       max_parallel_jobs=max_parallel_jobs, debug=debug)
            print '\nModels ready for next step... '
        else:
            #clean_directory(object_name, debug=debug)
            #print 'Directories clean! '
            print 'Finished step number ', main_counter+1, '\n'
            
        
    # Store the statistical parameters of the chain....
    time2run = 'Chain finished! Took  %s  days to finish the chain.' % ( (time.time() - start_time) / 86400.0)
    print time2run
        
