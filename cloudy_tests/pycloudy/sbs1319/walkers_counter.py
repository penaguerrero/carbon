import os
import numpy as np
import argparse
import string
from types import *

'''
This script only counts the number of walkers so far tested, considering that HTCondor has been used.
'''

parser = argparse.ArgumentParser(description='Count number of models ran so far.')
parser.add_argument("object_name",
                    action='store',
                    default=None,
                    help='The abbreviated name of the object, i.e. sbs1319.')
parser.add_argument("-n",
                    action='store_true',
                    dest="not_final",
                    default=False,
                    help='If the chain is not finished, use -n = not_final to read the correct chain file.')
parser.add_argument("-t",
                    action='store_true',
                    dest="temps_file",
                    default=False,
                    help='If wanting to use the temperatures file, use -t. For the case of multiple threads.')
args = parser.parse_args()


def read_chain4starting_posmatrix(chainfile, object_name):
    ''' This function obtains the data to create the previous position matrix to initialize the chain from
    last step of previous run. '''
    #chainfile = os.path.abspath('mcmc_'+object_name+'_chain0.dat')
    f = open(chainfile, 'r')
    TO3_list = []
    TO2_list = []
    Chi2_list = []
    for line in f.readlines():
        ''' skip the first 16 lines of the file, which contain the previous running time, benchmark abundances, best
        chi2 fitted model, and the uncertainties associated with that particular model.'''
        # Check that if the file has strings in the main body
        if '#' not in line:
            try:
                if '[' or '\n' in line:
                    abunds_temps_chi2_string = line.strip()
                    abunds_temps_chi2_list = string.split(abunds_temps_chi2_string, sep="[")
                    temps_chi2_list = string.split(abunds_temps_chi2_list[1], sep="\n")
                    to3 = float(string.replace(temps_chi2_list[0], "' ", ""))
                    to2 = float(string.replace(temps_chi2_list[1], "', ' ", ""))
                    chi = float(string.replace(temps_chi2_list[2], "']", ""))
                    TO3_list.append(to3)
                    TO2_list.append(to2)
                    Chi2_list.append(chi)
            except IndexError:
                He, O, CO, NO, NeO, SO, TO3, TO2, Chi2 = np.loadtxt(chainfile, skiprows=17, unpack=True)
                return He, O, CO, NO, NeO, SO, TO3, TO2, Chi2
                
    f.close()
    # Turn lists into numpy arrays
    TO3 = np.array(TO3_list)
    TO2 = np.array(TO2_list)
    Chi2 = np.array(Chi2_list)
    # Now read the abundances part of the file
    try:
        He, O, CO, NO, NeO, SO = np.loadtxt(chainfile, skiprows=16, usecols=(0,1,2,3,4,5), unpack=True)
    except ValueError:
        He, O, CO, NO, NeO, SO, TO3, TO2, Chi2 = 0,0,0,0,0,0,0,0,0
    return He, O, CO, NO, NeO, SO, TO3, TO2, Chi2

def get_number_of_walkers(object_name, not_final):
    path2file = '/home/pena/Documents/AptanaStudio3/carbon/cloudy_tests/pycloudy/'+object_name+'/mcmc_'
    if not_final:
        chain_counter = path2file+object_name+'_chain0.dat'
    else:
        chain_counter = path2file+object_name+'_chain.dat'
    if os.path.isfile(chain_counter):
        try:
            he, o, co, no, neo, so, to3, to2, lnChi2 = read_chain4starting_posmatrix(chain_counter, object_name)
            print '\nFrom reading: ', chain_counter
            if not isinstance(lnChi2, int):
                print 'So far, you have  %i  walkers. \n' % len(lnChi2)
            else:
                print '\nFile %s \n is still empty, wait for the first step of the chain to complete and check again. \n' % chain_counter
        except IOError:
            print '\nFile %s \n is still empty, wait for the first step of the chain to complete and check again. \n' % chain_counter
    else:
        print '\nFile ', chain_counter, ' does not exist.\n'

def get_length_of_temps(object_name):
    chain_counter = '/home/pena/Documents/AptanaStudio3/carbon/cloudy_tests/pycloudy/'+object_name+'_chaincounter.txt' 
    try:
        tempsO3, tempsO2, models_ran = np.loadtxt(chain_counter, unpack=True)
        print '\nFrom reading: ', chain_counter
        print 'So far, you have determined  %i  temperatures. \n' % len(tempsO3)
    except IOError:
        if os.path.isfile(chain_counter):
            print '\nFile %s \n is still empty, wait for the first step of the chain to complete and check again. \n' % chain_counter
        else:
            print '\nFile ', chain_counter, ' does not exist.\n'

# Set the the variables
object_name = args.object_name
not_final = args.not_final
temps_file = args.temps_file

# choose correct chain file
if not_final:
    get_number_of_walkers(object_name, not_final)
elif not_final==False and temps_file==False:
    get_number_of_walkers(object_name, not_final)
    
# Multiple threads. If the temperatures file is to be used:
if temps_file:
    get_length_of_temps(object_name)
