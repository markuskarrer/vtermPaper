'''
plot the bulk properties of the ICON-SB and the 1D-SB from the semi-idealized runs in one panel
'''

#from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug

import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
from netCDF4 import Dataset
from netCDF4 import chartostring
import traceback
import sys
import re #to search in string
import pandas as pd
#self-written functions
from functions import __plotting_functions
from functions import __postprocess_McSnow
from functions import __postprocess_PAMTRA
from functions import __postprocess_SB
from functions import __general_utilities
#read variables passed by shell script
tstep = int(os.environ["tstep"])
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"]
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output
MC_dir = os.environ["MC"]
adapt_version = int(os.environ["adapt_version"]) #reading the files of the appropriate adaption version
skipMC = (os.environ["skipMC"]=="True") #allows to run the scripts also if no McSnow data is there (only 1D-SB runs) #ATTENTION: not completely implemented yet

#directory of experiments
directory = MC_dir + "/experiments/"

####
#organize figure
###
number_of_plots = 3

figsize_height = 6.0/2.0*(number_of_plots)
fig	=	plt.figure(figsize=(8.0,figsize_height))#figsize=(4, 4))

#get timestep
i_timestep=(tstep/60)-1 #there is no output for t=0min after that there are output steps in 30 minute steps (this could vary)
#######
#read McSnow-data
#######
#read hei2massdens to get average /summed up values per height
filestring_hei2massdens = directory + experiment + "/hei2massdens.dat"
timestep = tstep/60 #TODO: do not hardcode the 30 minute output interval here

if skipMC:
    hei2massdens = __postprocess_McSnow.read_hei2massdens(filestring_hei2massdens,timestep=timestep,empty_flag=True) #get empty arrays with the same dimension as the true arrays in order to not plot the McSnow data
else:
    hei2massdens = __postprocess_McSnow.read_hei2massdens(filestring_hei2massdens,timestep=timestep)
    
#######
#read twomoment-data (1D)
#######
#define filestring
filestring = directory + experiment + "/twomom_d.dat"

#load netCDF4 file
twomom_file = Dataset(directory + experiment + '/twomom_d.ncdf',mode='r')
#create dictionary for all variables from twomom_d.ncdf
twomom = dict()

#if necessary change name of variables
varlist = twomom_file.variables
#read PAMTRA variables to pamData dictionary
for var in varlist:#read files and write it with different names in Data
    twomom[var] = np.squeeze(twomom_file.variables[var])

#########################
#read ICON-2mom (3D) data
##########################

#read time which should be analyzed from testcase string
date = (re.search(r'day(.*?)hour', testcase).group(1)) #this line gets the date from the testcase string
hour = int(re.search(r'hour(.*?)min', testcase).group(1)) #this line gets the hour from the testcase string
minute = int(re.search(r'min(.*?)s', testcase).group(1)) #this line gets the min from the testcase string
ssat = int(os.environ["ssat"]) #get supersaturation from governing script (this is the semi-idealized part of this setup)
#TODO also read in seconds
#convert time to s to find nearest output timestep
analyzed_time_s = hour*3600+minute*60

#read meteogram file
veras_tripexpol_simul = "/data/inscape/icon/experiments/juelich/testbed/testbed_"
if date=='20151124':
    filename = "/data/inscape/icon/experiments/tripex_220km/METEOGRAM_patch002_joyce.nc"
elif date.startswith('2018') or date.startswith('2019'): #TODO: this is a bit dirty in case you have other simulations in 2018/2019
    filename = veras_tripexpol_simul + date + "/METEOGRAM_patch001_" + date + "_joyce.nc"

#open netCDF4 file
nc = Dataset(filename)
#get formated string of date
time = pd.to_datetime(chartostring(nc.variables['date'][:]))
#get timestep in array for analysis
timestep_start = np.argmin(np.absolute(nc.variables["time"][:]-analyzed_time_s))


#initialize dictionary
twomom3D = dict()
#load some variables to dictionary
twomom3D["heights"] = nc.variables["height_2"] #mid-level height (relevant for most variables (except: W,.. )

#TODO: do not hardcode average time
average_time_s = 60*30 #120 minute average
timestep_end = np.argmin(np.absolute(nc.variables["time"][:]-analyzed_time_s-average_time_s))
#print some info in terminal
print 'analyzing time: ',time[timestep_start],' to: ',time[timestep_end]

for new_key,icon_key in zip(["temp","pres","qv","qc","qnc","qr","qnr","qi","qni","qs","qns","qg","qng","qh","qnh","rho","rhw"],
                            ["T",   "P",   "QV","QC","QNC","QR","QNR","QI","QNI","QS","QNS","QG","QNG","QG","QNH","RHO","REL_HUM"]):
    twomom3D[new_key] = np.mean(nc.variables[icon_key][timestep_start:timestep_end,:],axis=0) #ATTENTION: from here on they are temporally averaged
    if new_key in ("qc","qnc","qr","qnr","qi","qni","qs","qns","qg","qng","qh","qnh"):
        twomom3D[new_key] = twomom3D[new_key][None,:] #add singleton dimension for hydrometeors because of plotting routine
        
        twomom3D[new_key] = __general_utilities.q2abs(twomom3D[new_key],twomom3D["qv"],twomom3D["temp"],twomom3D["pres"],q_all_hydro=0) #(twomom3D["qc"]+twomom3D["qr"]+twomom3D["qi"]+twomom3D["qs"]+twomom3D["qg"]+twomom["qh"]))
        
####################################
#plot mixing ratio + number density
####################################
#number density
mass_num_flag = 0 #0-> plot only number flux; 1-> plot only mass flux; 2-> plot both 

ax11 = plt.subplot2grid((number_of_plots, 1), (0, 0))
if mass_num_flag==2:
    ax12 = ax.twiny()
else: #in case there is no need for a second axis, just pass the first ax twice
    ax12 = ax11
    
ax = __plotting_functions.plot_moments(ax11,ax12,twomom,hei2massdens,i_timestep,mass_num_flag=mass_num_flag)

#mass density
mass_num_flag = 1 #0-> plot only number flux; 1-> plot only mass flux; 2-> plot both 

ax21 = plt.subplot2grid((number_of_plots, 1), (1, 0))
if mass_num_flag==2:
    ax22 = ax.twiny()
else: #in case there is no need for a second axis, just pass the first ax twice
    ax22 = ax21
    
ax = __plotting_functions.plot_moments(ax21,ax22,twomom,hei2massdens,i_timestep,mass_num_flag=mass_num_flag)

##
#plot normalized mixing ratio (better known as mean particle mass) (q/qn) labelling normq...
##
#calculate q/qn respectively Md/Nd
ax31 = plt.subplot2grid((number_of_plots, 1), (2, 0))

ax32 = ax31
    
ax = __plotting_functions.plot_normmix(ax31,ax32,twomom,hei2massdens,i_timestep)

####################################
#add mixing ratio + number density from ICON-SB (3D-data)
####################################
#number density
mass_num_flag = 0 #0-> plot only number flux; 1-> plot only mass flux; 2-> plot both 


#plot_moments in __plotting_functions.py needs special names in the dictionary

i_timestep = 0 #dirty workaround: the variables do have just one timesteps here, because this is choose before
ax = __plotting_functions.plot_moments(ax11,ax12,twomom3D,hei2massdens,i_timestep,mass_num_flag=mass_num_flag,forced_linestyle='--')


#mass density
mass_num_flag = 1 #0-> plot only number flux; 1-> plot only mass flux; 2-> plot both 

    
ax = __plotting_functions.plot_moments(ax21,ax22,twomom3D,hei2massdens,i_timestep,mass_num_flag=mass_num_flag,forced_linestyle='--')

##
#plot normalized mixing ratio (better known as mean particle mass) (q/qn) labelling normq...
##
#calculate q/qn respectively Md/Nd
    
ax = __plotting_functions.plot_normmix(ax31,ax32,twomom3D,hei2massdens,i_timestep,forced_linestyle='--')

###set y-limits to a reasonable range
#get maximum of mass in 3D
i_maxq = int(np.argmax(twomom3D["qi"]+twomom3D["qs"]))
height_of_max_mass3D = np.ma.getdata(twomom3D["heights"])[i_maxq]
#get initialization height
init_height = twomom["heights"][0]
for ax_now in [ax11,ax21,ax31]:
    ax_now.set_ylim([height_of_max_mass3D-1000,init_height+1000])
    ax_now.axhline(y=height_of_max_mass3D,color='grey',linestyle='--')
    ax_now.axhline(y=init_height,color='grey',linestyle='--')
    ax_now.plot(0,0,color='k',linestyle='--',label='3D-SB')
    ax_now.plot(0,0,color='k',linestyle='-',label='1D-SB')
    ax_now.legend()

#save figure
plt.tight_layout()
if not os.path.exists('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment): #create direktory if it does not exists
    os.makedirs('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment)
out_filestring = "/3D1Dcomparison_" + testcase + "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'])