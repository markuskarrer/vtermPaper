'''
create an overview panel with properties of McSnow and PAMTRA output
'''

from IPython.core.debugger import Tracer ; debug=Tracer()

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
from functions import __general_utilities
#read variables passed by shell script
tstep = int(os.environ["tstep"])
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"]
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output
MC_dir = os.environ["MC"]
adapt_version = int(os.environ["adapt_version"]) #reading the files of the appropriate adaption version
skipMC = (os.environ["skipMC"]=="True") #allows to run the scripts also if no McSnow data is there (only 1D-SB runs) #ATTENTION: not completely implemented yet

if "separated_by_sensruns" in os.environ.keys():
    separated_by_sensruns = (os.environ["separated_by_sensruns"]=="True") #plot a line for different sensitivity runs
    separated_by_sensruns_onestring= os.environ["model_setup_specifier_onestring"]
else:
    print "separated_by_sensruns not found in environ: set to False"
    separated_by_sensruns = False
if "switch_off_processes" in os.environ.keys():
    switch_off_processes_str = os.environ["switch_off_processes"]
else:
    switch_off_processes_str = ''
#directory of experiments
directory = MC_dir + "/experiments/"

num_plots = 4 #number of subplots from pamtra variables and twomoment output
figsize_height = 6.0/2.0*(num_plots)
fig	=	plt.figure(figsize=(8.0,figsize_height))#figsize=(4, 4))

height_bounds = [5000,0] #set some default heigh-bounds
    
###### 
#plot SB pamtra output data
#####
if separated_by_sensruns:
    sensrun_list = separated_by_sensruns_onestring

    sensrun_list = separated_by_sensruns_onestring.split('_')
else:
    sensrun_list = [""]

#plot different lines for different sensitivity runs
shape_spec_list=['nu0mu0.5e'] #this is the default
pprop_list=[]
stick_list=['stick-2'] #-2: Lin83
sensrun_set_dic = dict() #dictionary which contains all the settings (varying particles, shapedist and stickeff)
linestyles=["-","--","-.",":"]
colors=["b","r","g"]
markers = ["","x","o"]
axKaW= plt.subplot2grid((num_plots, 1), (0, 0))
axXKa= plt.subplot2grid((num_plots, 1), (1, 0))
axvDoppler= plt.subplot2grid((num_plots, 1), (2, 0))
axZe= plt.subplot2grid((num_plots, 1), (3, 0))
#axmeanmass= plt.subplot2grid((num_plots, 1), (4, 0))
for i_sensrun, sensrun_now in enumerate(sensrun_list): #loop over different sensitivity runs to get all information about the settings
    if 'nu' in sensrun_now:
        shape_spec = 'nu' +sensrun_now.split('nu')[1] #this gets the specifier of the size dist parameters
        if shape_spec not in shape_spec_list:
            shape_spec_list.append(shape_spec)
        try:
            sensrun_set_dic[sensrun_now]["sizedist"] = shape_spec
        except: #in the first occurence of sensrun_now we need to create that entry first
            sensrun_set_dic[sensrun_now] = {} 
            sensrun_set_dic[sensrun_now]["sizedist"] = shape_spec
            
        sensrun_set_dic[sensrun_now]["pprop"] = sensrun_now.split('nu')[0] #this is the default
        if 'stick' in sensrun_now: #both: shapedist and stick are different than default
            stick_spec = 'stick' +sensrun_now.split('stick')[1] #this gets the specifier of the sticking efficiency
            if stick_spec not in stick_spec:
                stick_list.append(stick_spec)
            sensrun_set_dic[sensrun_now]["stick"] = stick_spec

        else: 
            sensrun_set_dic[sensrun_now]["stick"] = stick_list[0] #this is the default
        if not sensrun_set_dic[sensrun_now]["pprop"] in pprop_list:
            pprop_list.append(sensrun_set_dic[sensrun]["pprop"])
    elif 'stick' in sensrun_now: #default shapedist but different stick
        stick_spec = 'stick' +sensrun_now.split('stick')[1] #this gets the specifier of the sticking efficiency
        pprop_spec = sensrun_now.split('stick')[0] #this gets the specifier of the particle properties 
        try:
            sensrun_set_dic[sensrun_now]["stick"] = stick_spec
        except:
            sensrun_set_dic[sensrun_now] = {} 
            sensrun_set_dic[sensrun_now]["stick"] = stick_spec
        if stick_spec not in stick_list:
            stick_list.append(stick_spec)
        sensrun_set_dic[sensrun_now]["pprop"] = pprop_spec #this is the default
        sensrun_set_dic[sensrun_now]["sizedist"] = shape_spec_list[0] #this is the default
        if not sensrun_set_dic[sensrun_now]["pprop"] in pprop_list:
            pprop_list.append(sensrun_set_dic[sensrun_now]["pprop"])
    else: #here we have the run with default shape dist and stickeff
        pprop_list.append(sensrun_now)
        sensrun_set_dic[sensrun_now] = {} 
        sensrun_set_dic[sensrun_now]["pprop"] = sensrun_now #this is the default
        sensrun_set_dic[sensrun_now]["sizedist"] = shape_spec_list[0] #this is the default
        sensrun_set_dic[sensrun_now]["stick"] = stick_list[0] #this is the default

    if separated_by_sensruns:#modify the experiment string
        experiment_splitted = experiment.split('_',2) #this results e.g. in ['1d', 'powerlawJplate', 'xi10000000_
        experiment_splitted[1] = sensrun_now #replace experiment part
        experiment = "_".join(experiment_splitted)


sorted_sensrunnames=sorted(sensrun_set_dic.keys(), key=lambda x:x.lower())
for i_sensrun, sensrun_now in enumerate(sorted_sensrunnames): #loop over different SB sensitivity runs
    print sensrun_now
    ###modify the experiment string
    #replace the specifier (e.g. powerlawOLD AtlasJmixcolumndend 
    experiment_splitted = experiment.split('_',2) #this results e.g. in ['1d', 'powerlawJplate', 'xi10000000_
    if not ("nu" in sensrun_now): #if size dist is not default it is combined in the beginning of the run name
        print sensrun_now,sensrun_set_dic[sensrun_now].keys()
        experiment_splitted[1] = sensrun_set_dic[sensrun_now]["pprop"] #replace experiment part
    else:
        experiment_splitted[1] = sensrun_now 
    experiment = "_".join(experiment_splitted)
    #replace the stickeff number
    experiment_split0 = experiment.split('stick',1) #this results e.g. in ['1d_powerlawJplate...','_dt...'
    experiment_split1 = experiment.split('_dt1',1) #this results e.g. in ['1d_powerlawJplate...stick-1','_melt.'
    experiment = experiment_split0[0] + sensrun_set_dic[sensrun_now]["stick"] + '_dt1' + experiment_split1[1] #put together and insert the sticking efficiency
    #debug()

    #define file to read
    pam_filestring = directory + experiment + "/PAMTRA_2mom_" + testcase +  '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min.nc'
    ###
    #read pamtra output to pamData dictionary
    ###
    pamData = __postprocess_PAMTRA.read_pamtra(pam_filestring)
 
    #read atmospheric variables once
    filestring_atmo = directory + experiment + "/atmo.dat"
    atmo = __postprocess_McSnow.read_atmo(experiment,filestring_atmo)
    model_top=5000. #[m]
    n_heights=50 
    heightvec = np.linspace(model_top/n_heights,model_top,n_heights)
     
    #interpolate atmospheric variables to heightvec
    atmo_interpolated = __postprocess_McSnow.interpolate_atmo(atmo,heightvec)

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
    
    #get the indices for the settings 
    i_prop = pprop_list.index(sensrun_set_dic[sensrun_now]['pprop'])
    i_stick = stick_list.index(sensrun_set_dic[sensrun_now]['stick']) 
    i_shape = shape_spec_list.index(sensrun_set_dic[sensrun_now]['sizedist']) 
    twomom['mmeans'] = twomom['qs']/twomom['qns'] #calculate mean mass
    #calculate vars from pamtra-output
    pamData["DWR_X_Ka"] = pamData["Ze"][:,1]-pamData["Ze"][:,2]   
    pamData["DWR_Ka_W"] = pamData["Ze"][:,0]-pamData["Ze"][:,1]   
    pamData["Ze"][:,1][pamData["Ze"][:,1]<-100]=np.nan
    pamData["ZeKa"] = pamData["Ze"][:,1]
    pamData["Radar_MeanDopplerVel"][:,1][pamData["Radar_MeanDopplerVel"][:,1]<=0.0]=np.nan #remove negative VDoppler
    pamData["vDopplerKa"] = pamData["Radar_MeanDopplerVel"][:,1] 
    #axmeanmass.semilogx(twomom['mmeans'][-1],twomom["heights"],label=sensrun_now)
    #get temperature as first y-axis
    y_coordinate=atmo_interpolated["T"]-273.15
    #plot 
    for ax,varname in zip([axKaW,axXKa,axvDoppler,axZe],["DWR_X_Ka","DWR_Ka_W","vDopplerKa","ZeKa"]):
        ax.plot(pamData[varname],y_coordinate,linestyle=linestyles[i_stick],color=colors[i_prop],marker=markers[i_shape])   
        #ax2 = ax.twinx() #does not work yet
        #ax2.plot(np.nan*np.ones_like(pamData["height"]),pamData["height"])
        

    #labels for each sensrun
    if sensrun_now=="AtlasJmixcolumndend":
        sensrun_now_label="Mix2 (Cotton86)"
    elif sensrun_now=="AtlasJmixcolumndendstick3":
        sensrun_now_label="Mix2 (Connolly/2)"
    elif sensrun_now=="powerlawOLDstick3":
        sensrun_now_label="powerlaw def (Connolly/2)"
    elif sensrun_now=="powerlawOLD":
        sensrun_now_label="powerlaw def (Cotton86)"
    else:
        sensrun_now_label=sensrun_now
    for ax in [axKaW,axXKa,axvDoppler,axZe]:
        ax.plot(np.nan,np.nan,label=sensrun_now_label,linestyle=linestyles[i_stick],color=colors[i_prop],marker=markers[i_shape])

#adjust limits, add legend
for ax in [axKaW,axXKa]:
    ax.set_xlim([0,15])
for ax in [axKaW,axXKa,axvDoppler,axZe]:
    ax.set_ylim([0,-30])
    ax.set_ylabel("T [$^\circ$C]")
    ax.legend()
#axmeanmass.legend()
axKaW.set_xlabel("DWR Ka-W [dB]")
axXKa.set_xlabel("DWR X-Ka [dB]")
axvDoppler.set_xlabel("MDV [m/s]")
axZe.set_xlabel("Ze [dBz]")
#save figure
plt.tight_layout
if not os.path.exists('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment): #create direktory if it does not exists
    os.makedirs('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment)
out_filestring = "/radar_profiles_" + switch_off_processes_str + "_" + testcase + "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'])
