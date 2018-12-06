'''
plots the fluxes and the absolute values of the number density and the mass density
'''

#from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import math
import os
import subprocess
from netCDF4 import Dataset
import traceback
import sys
#self-written functions
from functions import __plotting_functions
from functions import __postprocess_McSnow
from functions import __postprocess_PAMTRA
from functions import __postprocess_SB

#read variables passed by shell script
tstep = int(os.environ["tstep"])
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"]
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output
MC_dir = os.environ["MC"]
skipMC = (os.environ["skipMC"]=="True") #allows to run the scripts also if no McSnow data is there (only 1D-SB runs)

#directory of experiments
directory = MC_dir + "/experiments/"

#read hei2massdens to get average /summed up values per height
filestring_hei2massdens = directory + experiment + "/hei2massdens.dat"
timestep = tstep/60 #TODO: do not hardcode the 30 minute output interval here
if skipMC:
    hei2massdens = __postprocess_McSnow.read_hei2massdens(filestring_hei2massdens,timestep=timestep,empty_flag=True) #get empty arrays with the same dimension as the true arrays in order to not plot the McSnow data
else:
    hei2massdens = __postprocess_McSnow.read_hei2massdens(filestring_hei2massdens,timestep=timestep)

number_of_plots = 5

figsize_height = 6.0/2.0*(number_of_plots)
fig	=	plt.figure(figsize=(8.0,figsize_height))#figsize=(4, 4))

#get timestep
i_timestep=(tstep/60)-1 #there is no output for t=0min after that there are output steps in 30 minute steps (this could vary)

#######
#read twomoment-data
#######
#define filestring
filestring = directory + experiment + "/twomom_d.dat"

#load netCDF4 file
twomom_file = Dataset(directory + experiment + '/twomom_d.ncdf',mode='r')
#create dictionary for all variables from twomom_d.ncdf
twomom = dict()

#if necessary change name of variables
varlist = twomom_file.variables
#read 2mom variables to twomom dictionary
for var in varlist:#read files and write it with different names in Data
    twomom[var] = np.squeeze(twomom_file.variables[var])

####################################
#plot mixing ratio + number density
####################################
#number density
mass_num_flag = 0 #0-> plot only number flux; 1-> plot only mass flux; 2-> plot both 

ax = plt.subplot2grid((number_of_plots, 1), (0, 0))
if mass_num_flag==2:
    ax2 = ax.twiny()
else: #in case there is no need for a second axis, just pass the first ax twice
    ax2 = ax
    
ax = __plotting_functions.plot_moments(ax,ax2,twomom,hei2massdens,i_timestep,mass_num_flag=mass_num_flag)

#mass density
mass_num_flag = 1 #0-> plot only number flux; 1-> plot only mass flux; 2-> plot both 

ax = plt.subplot2grid((number_of_plots, 1), (1, 0))
if mass_num_flag==2:
    ax2 = ax.twiny()
else: #in case there is no need for a second axis, just pass the first ax twice
    ax2 = ax
    
ax = __plotting_functions.plot_moments(ax,ax2,twomom,hei2massdens,i_timestep,mass_num_flag=mass_num_flag)

##
#plot normalized mixing ration (q/qn) labelling normq...
##
#calculate q/qn respectively Md/Nd
ax = plt.subplot2grid((number_of_plots, 1), (2, 0))

ax2 = ax
    
ax = __plotting_functions.plot_normmix(ax,ax2,twomom,hei2massdens,i_timestep)

############
#plot fluxes
############
#number flux
mass_num_flag = 0 #0-> plot only number flux; 1-> plot only mass flux; 2-> plot both 

ax = plt.subplot2grid((number_of_plots, 1), (3, 0))
if mass_num_flag==2:
    ax2 = ax.twiny()
else: #in case there is no need for a second axis, just pass the first ax twice
    ax2 = ax
    
ax = __plotting_functions.plot_fluxes(ax,ax2,twomom,hei2massdens,i_timestep,mass_num_flag=mass_num_flag)

#mass flux
mass_num_flag = 1 #0-> plot only number flux; 1-> plot only mass flux; 2-> plot both 

ax = plt.subplot2grid((number_of_plots, 1), (4, 0))
if mass_num_flag==2:
    ax2 = ax.twiny()
else: #in case there is no need for a second axis, just pass the first ax twice
    ax2 = ax
    
ax = __plotting_functions.plot_fluxes(ax,ax2,twomom,hei2massdens,i_timestep,mass_num_flag=mass_num_flag)

#save figure
plt.tight_layout()
if not os.path.exists('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment): #create directory if it does not exists
    os.makedirs('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment)
out_filestring = "/fluxes_" + testcase + "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'])