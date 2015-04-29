#!/usr/bin/env python
import os
import numpy as np
import string
import copy
import pyCloudy as pc
from sys import argv
from glob import glob


'''
This script will run the cloudy model and calculate the likelihood of such model
 INPUT: 
    - theta = the set of abundances to be set for running Cloudy (currently it is 6 numbers read independently).
    - unique_filename = a string to make the text file unique.
    - object_name = name of the object to read the correct benchmark intensities.
    - manual_measurement = were the measurements made manually? (True or False). This is to read the 
                         correct benchmark intensities.
    - initial_Cloudy_conditions = name of the file that sets the rest of the parameters for running Cloudy.
    * Total inputs expected: 10
 OUTPUT:
    - A text file with a unique name containing the modeled abundances, the calculated lnChi2, and 
     the modeled temperatures of O3 and O2.
'''


#################################################################################################################

# Changing the location and version of the cloudy executable. Will turn the relative path into an absolute path.
#cloudyexe = 'Addons/cloudy/c13.03all/c13.03/source/cloudy.exe'
#cloudyexe_path_list = string.split(os.getcwd(), sep='AptanaStudio3')
#cloudyexe_path = os.path.join(cloudyexe_path_list[0], cloudyexe)
# Path HARDCODED FOR HTCONDOR...
cloudyexe_path = '/home/pena/Documents/Addons/cloudy/c13.03all/c13.03/source/cloudy.exe'
print 'Path for the Cloudy executable: ', cloudyexe_path
pc.config.cloudy_exe = cloudyexe_path


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
        self.dir = dir
        self.verbosity = verbosity
        self.options = options
        self.iterations = iterations
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
        # Filling the object with the parameters. Defining the ionizing SED: Effective temperature and luminosity,
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
        pc.log_.message('Running {0}'.format(self.model_name), calling = self.model_name)
        # Running Cloudy with a timer. Here we reset it to 0.
        pc.log_.timer('Starting Cloudy', quiet = True, calling = self.model_name)
        c_input.run_cloudy()
        pc.log_.timer('Cloudy ended after seconds:', calling = self.model_name)
        
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
        # Path HARDCODED FOR HTCONDOR...
        cloudyexe_path = '/home/pena/Documents/Addons/cloudy/c13.03all/c13.03/source/cloudy.exe'
        print 'Path for the Cloudy executable: ', cloudyexe_path
        print 'executable exists? ', os.path.exists(cloudyexe_path)
        pc.config.cloudy_exe = cloudyexe_path
        # Now run the model
        self.runCldy_with_initconds()
        print'***  got initial conditions and ran Cloudy!'
        self.lines_file = self.read_Cldyouts()
        print '***  reading Cloudy outputs!'
        self.clean_Cldyoutfiles()
        print'***  cleaning Cloudy output files...'
        return self.lines_file


# Functions needed for independent runs
def lnlikehd(object_name, unique_filename, new_Iobs, new_Iobserr, new_Imod, modeled_Otemperatures, bptid, bptImod):
    model = np.array(new_Imod)
    y = np.array(new_Iobs)
    e = np.array(new_Iobserr)
    csq = ( y - model )**2 / e**2
    Chi2 = csq.sum()
    lnprobChi2 = -Chi2 / 2.0
    # print probability and model temperatures to a file
    path2file = '/home/pena/Documents/AptanaStudio3/carbon/cloudy_tests/pycloudy/'+object_name+'/'
    textfile4Chi2 = path2file+'Chi2_'+object_name+'_'+unique_filename+'.txt'
    tf = open(textfile4Chi2, 'w+')
    mod_TO3, mod_TO2 = modeled_Otemperatures
    mod_TO3 = float(mod_TO3)
    mod_TO2 = float(mod_TO2.strip())
    he, o, co, no, neo, so = theta
    print >> tf, '{:<6} {:<6} {:<6} {:<6} {:<6} {:<6} {:<16} {:<20} {:<20}'.format(he, o, co, no, neo, so, lnprobChi2, mod_TO3, mod_TO2)
    tf.close()
    # create file with BPT lines: 4861, 5007, 6548, 6563, 6584
    bptfile = path2file+'bpt_'+object_name+'_'+unique_filename+'.txt'
    tf = open(bptfile, 'w+')
    print >> tf, '{:<6} {:<6} {:<6} {:<6} {:<6}'.format(bptid[0], bptid[1], bptid[2], bptid[3], bptid[4])
    print >> tf, '{:<6} {:<6} {:<6} {:<6} {:<6}'.format(bptImod[0], bptImod[1], bptImod[2], bptImod[3], bptImod[4])
    tf.close()

def read_cldyfile(clyfile):
    f = open(clyfile, 'r')
    i = 0
    for line in f.readlines():
        if i == 0:
            model_name = line.strip()
        elif i == 1:
            dens = float(line.strip())
        elif i == 2:
            emis_tab_string = line.strip()
        elif i == 3:
            theta_string = line.strip()
        elif i == 4:
            stb99_table = line.strip()
        elif i == 5:
            age = line.strip()
        elif i == 6:
            dir  = line.strip()
        elif i == 7:
            verb = int(line.strip())
        elif i == 8:
            options_string = line.strip()
        elif i == 9:
            iterations = int(line.strip())
        elif i == 10:
            keep_files = line.strip()
        i = i + 1
    f.close()
    emis_tab_list = string.replace(emis_tab_string, "[", "")
    emis_tab_list = string.replace(emis_tab_list, "]", "")
    emis_tab_list = string.replace(emis_tab_list, "'", "")
    emis_tab = string.split(emis_tab_list, sep=",")
    theta_list = string.replace(theta_string, "[", "")
    theta_list = string.replace(theta_list, "]", "")
    theta_list = string.split(theta_list, sep=",")
    theta = []
    for t in theta_list:
        theta.append(float(t))
    options_list = string.replace(options_string, "(", "")
    options_list = string.replace(options_list, ")", "")
    options_list = string.replace(options_list, "'", "")
    options_list = string.split(options_list, sep=",")
    options = tuple(options_list)
    if keep_files == 'None':
        keep_files = None
    initial_Cloudy_conditions = [model_name, dens, emis_tab, theta, stb99_table, age, dir, verb, options, iterations, keep_files]
    return initial_Cloudy_conditions

def get_measured_lines(object_name, manual_measurement):
    # get the benchmark measurements
    if manual_measurement:
        file2use = '_measuredLI_RedCor_Ebv.txt'   # my manual line measurements
    else:
        file2use = '_RedCor_Ebv.txt'   # my code's line measurements
    # Path HARDCODED for HTCondor
    measured_lines_file_path = '/home/pena/Documents/AptanaStudio3/carbon/results/'+object_name+'/'+object_name+file2use
    #measured_lines_file = 'carbon/results/'+object_name+'/'+object_name+file2use
    #measured_lines_file_path_list = string.split(os.getcwd(), sep='carbon')
    #measured_lines_file_path = os.path.join(measured_lines_file_path_list[0], measured_lines_file)
    # now retrieve the columns we need from the file
    meas_lineIDs = []   # list of IDs in angstroms
    meas_Isrel2Hbeta = []   # list of lines relative to Hbeta reddening corrected 
    meas_Ierr = []   # absolute error of the line intensity
    meas_Iper = []   # percentage error of the line intensity  
    meas_EW = []   # NOT BEING USED FOR THE MOMENT
    meas = open(measured_lines_file_path)
    _ = meas.readline()  # columns header
    linesofinterest = [4861, 6563, 4340, 5876, 4686, 3727, 5007, 9532, 1907, 1910, 1909]
    for line in meas:
        line = line.strip()   # gets rid of \n at the end of the line
        cols = line.split()   # splits the line into a list of columns
        # the columns of interest are: ID=0, Intensity=9, Ierr=10, Iper=11, EW=12
        lineid = int(np.round(float(cols[0]), decimals=0))
        if lineid in linesofinterest:
            meas_lineIDs.append(lineid)
            meas_Isrel2Hbeta.append(float(cols[9]))
            meas_Ierr.append(float(cols[10]))
            meas_Iper.append(float(cols[11]))
            meas_EW.append(float(cols[12]))
    meas.close()
    measured_lines = [meas_lineIDs, meas_Isrel2Hbeta, meas_Ierr, meas_Iper, meas_EW]
    return measured_lines

def get_modeled_lines(theta, initial_Cloudy_conditions, unique_filename):   # get the Cloudy model
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
    model_name = model_name+"_"+unique_filename
    print "The name of this model is: ", model_name
    #print "I received the following initial conditions:"
    #print 'model_name, dens, emis_tab, abunds, stb99_table, age, dir, verbosity, options, iterations, keep_files ='
    #print model_name, dens, emis_tab, abunds, stb99_table, age, dir, verbosity, options, iterations, keep_files
    print "Ran a Cloudy model for the following abundances:"
    print abunds
    cldymod = PyCloudy_model(model_name, dens, emis_tab, abunds, stb99_table, age, dir, verbosity, options, iterations, keep_files)
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
    #final_default_extensions_list = ['.in', '.out', '.txt']
    final_default_extensions_list = ['.in', '.txt']
    for fdf in final_default_extensions_list:
        file2beerased = glob(dir+model_name+'*'+fdf)
        #print file2beerased[0]
        os.remove(file2beerased[0])
    print 'mod_lineIDs: ', mod_lineIDs
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

def get_BPTlines(IDmod, Imod):
    bptlines = [4861, 5007, 6548, 6563, 6584]
    bptid =[]
    bptImod = []
    for bl in bptlines:
        if bl in IDmod:
            idx = IDmod.index(bl)
            bptid.append(bl)
            bptImod.append(Imod[idx])
        else:
            bptid.append(bl)
            bptImod.append(0.0)
    print 'bptid: ', bptid
    return bptid, bptImod
    

#################################################################################################################


# Set the given arguments to corresponding variables
he = float(argv[1])
o = float(argv[2])
co = float(argv[3])
no = float(argv[4])
neo = float(argv[5])
so = float(argv[6])
unique_filename = argv[7]
object_name = argv[8]
manual_measurement = argv[9]
# make sure that the variable is the correct type
if manual_measurement == 'True':
    manual_measurement = True
else:
    manual_measurement = False
clyfile = argv[10]

# Abundance set
theta = [he, o, co, no, neo, so]

# Get the benchmark intensity lines
measured_lines = get_measured_lines(object_name, manual_measurement)
IDobs, Iobs, Iobserr, _, _ = measured_lines

# Get the Cloudy modeled lines
initial_Cloudy_conditions = read_cldyfile(clyfile)
#print 'THESE ARE THE CONDITIONS GIVEN TO RUN CLOUDY:'
#print 'theta = ', theta
#print 'iniCldy_conds = ', initial_Cloudy_conditions
modeled_Otemperatures, modeled_lines = get_modeled_lines(theta, initial_Cloudy_conditions, unique_filename)
IDmod, _, Imod = modeled_lines

# Get the relevant lines from the measurements and from the model
new_Iobs, new_Iobserr, new_Imod = find_Itheo_in_Iobs(IDobs, Iobs, Iobserr, IDmod, Imod)

# Get the BPT line intensities to make diagrams
bptid, bptImod = get_BPTlines(IDmod, Imod)

# Get the Chi2 and add it to the priors' probability AND save file with BPT lines
lnlikehd(object_name, unique_filename, new_Iobs, new_Iobserr, new_Imod, modeled_Otemperatures, bptid, bptImod)

