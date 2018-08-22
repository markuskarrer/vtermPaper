'''
this script creates ncdf data from McSnows mass2fr..dat file with the original and some postprocessed variables (like temporal averages)
'''
#from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug

#import modules
import numpy as np
import sys #for some debugging
import re #for selecting numbers from file-string
import subprocess #for using shell commands with subprocess.call()
import os #import environment variables
from netCDF4 import Dataset 
import time

#self written #these files are linked here from the pythoncode/functions directory
import __postprocess_McSnow

#read variables passed by shell script
tstep = os.environ["tstep"] #string with 4 numbers and 'min'
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"] #"more readable" string of the experiment specifications
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output

#directory of experiments
directory = "/home/mkarrer/Dokumente/McSnow/MCSNOW/experiments/"

#read timestep from string
dtc = int(re.search(r'dtc(.*?)_nrp', experiment).group(1)) #this line gets the values after dtc and before _nrp -> timestep of collision in s

#perform temporal averaging
t_after_full = 0
while t_after_full<av_tstep: #1. merge all timesteps <av_tstep in one tuple
    filestring_mass2fr = directory + experiment + "/mass2fr_" + tstep + '_' + str(t_after_full).zfill(2) + 's' + ".dat"
    print 'reading: ' + filestring_mass2fr
    if t_after_full==0: #Initialize tuple with SP dictionary for each individual timestep here
        SP_nonaveraged = (__postprocess_McSnow.read_mass2frdat(experiment,filestring_mass2fr),)
    else:
        SP_nonaveraged = SP_nonaveraged + (__postprocess_McSnow.read_mass2frdat(experiment,filestring_mass2fr),)
    t_after_full=t_after_full+dtc

#2. compute the 'average' of the SP list, which results in a longer list (approx SP_pertimestep*averaged timesteps) but approx the same number of RP
SP = __postprocess_McSnow.average_SPlists(SP_nonaveraged)
#get number of timesteps which can be averaged averaged
num_timesteps = len(SP_nonaveraged)

#create netCDF4 file
dataset = Dataset(directory + experiment + '/mass2fr_' + tstep + '_avtstep_' + str(av_tstep) + '.ncdf', 'w',  format='NETCDF4_CLASSIC')
#Global attributes
dataset.description = 'list with SP and their properties'
dataset.history = 'Created ' + time.ctime(time.time())

#create Dimensions
dim_SP_all = dataset.createDimension('dim_SP_all_av' +str((num_timesteps)*dtc),len(SP['xi']))

#create Variables and fill them
#1. grid variables
#i_SP_all = dataset.createVariable('i_SP_index_av' +str((num_timesteps)*dtc),np.int32,('dim_SP_all_av' +str((num_timesteps)*dtc),)) #"grid-variable"
#dataset.variables['i_SP_index_av' +str((num_timesteps)*dtc)][:] = range(0,len(SP['xi'])) #thats a trivial variable
#2. real variables
unit = ['kg','1','kg','1','kg m-3','kg','m','m','1','m2','m s-1','m2','kg','m']
for i,key in enumerate(SP.keys()):
    dataset.createVariable(key,np.float32,('dim_SP_all_av' +str((num_timesteps)*dtc),)); dataset.variables[key].units = unit[i]
    dataset.variables[key][:] = SP[key]

#fill variables
#i_SP_all[:] = range(0,len(SP['xi']))    
#dataset.variables["m_tot"][:] = SP["m_tot"]
#print dataset.variables['m_tot']
#print 'dim_SP_all',len(dim_SP_all),dim_SP_all.isunlimited()
print dataset.history
dataset.close()
    #print SP_averaged
#from IPython.core.debugger import Tracer ; Tracer()()