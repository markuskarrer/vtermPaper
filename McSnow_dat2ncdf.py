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
tstep = int(os.environ["tstep"]) #string with 4 numbers and 'min'
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"] #"more readable" string of the experiment specifications
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output
MC_dir = os.environ["MC"]

#directory of experiments
directory = MC_dir + "/experiments/"

#read timestep from string
dtc = int(re.search(r'dtc(.*?)_nrp', experiment).group(1)) #this line gets the values after dtc and before _nrp -> timestep of collision in s

#perform temporal averaging
t_after_full = 0 #time after the timestep with the full output in seconds
while t_after_full<av_tstep: #1. merge all timesteps <av_tstep in one tuple
    filestring_mass2fr = directory + experiment + "/mass2fr_" + str(tstep).zfill(4) + 'min_' + str(t_after_full).zfill(2) + 's' + ".dat"
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
dataset = Dataset(directory + experiment + '/mass2fr_' + str(tstep).zfill(4) + 'min_avtstep_' + str(av_tstep) + '.ncdf', 'w',  format='NETCDF4_CLASSIC')
#Global attributes
dataset.description = 'list with SP and their properties'
dataset.history = 'Created ' + time.ctime(time.time())

#create Dimensions
dim_SP_all = dataset.createDimension('dim_SP_all_av' +str((num_timesteps)*dtc),len(SP['xi']))

#create Variables and fill them
unit = ['kg','1','kg','1','kg m-3','kg','m','m','1','m2','m s-1','m2','kg','m']
for i,key in enumerate(SP.keys()):
    dataset.createVariable(key,np.float32,('dim_SP_all_av' +str((num_timesteps)*dtc),)); dataset.variables[key].units = unit[i]
    dataset.variables[key][:] = SP[key]

print "mass2fr... file created at: ",dataset.history
dataset.close()

#################################################
#also save SB-twomoment scheme output to netCDF4
#################################################

twomom_bool = True #TODO: do not hardcode this here?
if not twomom_bool:
    sys.exit(0)

#read number of vertical levels from string
nz = int(re.search(r'nz(.*?)_', experiment).group(1)) #this line gets the values after dtc and before _nrp -> timestep of collision in s

#read the twomom_d.dat file to the dictionary twomom
filestring_twomom_d = MC_dir + '/experiments/' + experiment + '/twomom_d.dat' 
twomom = __postprocess_McSnow.read_twomom_d(experiment,filestring_twomom_d,nz) #returns the twomom variables in a dictionary

#create height axis from dz (layer thickness) array
heights = np.cumsum(twomom['dz'][0,::-1])[::-1] #fix index 0 is allowed unless level thickness change with height

#create netCDF4 file
dataset = Dataset(directory + experiment + '/twomom_d.ncdf', 'w',  format='NETCDF4_CLASSIC')
#Global attributes
dataset.description = 'Output of SB twomoment-scheme (embedded in McSnow) converted to netCDF4'
dataset.history = 'Created ' + time.ctime(time.time())

#create Dimensions
dim_height = dataset.createDimension('height',nz)
dim_time = dataset.createDimension('time',twomom[twomom.keys()[0]].shape[0])

#create Variables and fill them
#1. coordinate variables
dim_times = dataset.createVariable('times',np.float64,('time',))
dim_heights = dataset.createVariable('heights',np.int32,('height'))
timestep_out = 30 #TODO: do not hardcode this
dataset.variables['times'][:] = np.arange(0,twomom[twomom.keys()[0]].shape[0])*timestep_out; dataset.variables['times'].units = 'min'
dataset.variables['heights'][:] = heights; dataset.variables['heights'].units = 'm' 
#2. real variables
#unit = ['kg','1','kg','1','kg m-3','kg','m','m','1','m2','m s-1','m2','kg','m'] #TODO: define units
for i,key in enumerate(twomom.keys()):
    dataset.createVariable(key,np.float32,('time','height')); #dataset.variables[key].units = unit[i]
    dataset.variables[key][:] = twomom[key]

print "twomom_d file created at: ",dataset.history
dataset.close()
