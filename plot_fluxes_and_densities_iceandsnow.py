'''
plots the fluxes and the absolute values of the number density and the mass density
'''

from IPython.core.debugger import Tracer ; debug=Tracer() #insert this line somewhere to debug

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
sys.path.append("/home/mkarrer/Dokumente/Doripamtra_ICON/comparison")
import calc_theoretical_DWRs_MK 

#read variables passed by shell script
tstep = int(os.environ["tstep"])
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"]
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output
MC_dir = os.environ["MC"]

if "separated_by_sensruns" in os.environ.keys():
    separated_by_sensruns = (os.environ["separated_by_sensruns"]=="True") #plot a line for different sensitivity runs
    separated_by_sensruns_onestring= os.environ["model_setup_specifier_onestring"]
else:
    print "separated_by_sensruns not found in environ: set to False"
    separated_by_sensruns = False
    separated_by_sensruns_onestring=""
if "switch_off_processes" in os.environ.keys():
    switch_off_processes_str = os.environ["switch_off_processes"]
else:
    switch_off_processes_str = ''

#directory of experiments
directory = MC_dir + "/experiments/"

number_of_plots = 5

figsize_height = 6.0/2.0*(number_of_plots)
fig	=	plt.figure(figsize=(8.0,figsize_height))#figsize=(4, 4))

#get timestep
i_timestep=-1 #(tstep/10)-1 #there is no output for t=0min after that there are output steps in 30 minute steps (this could vary)


#plot different lines for different sensitivity runs
sensrun_list = separated_by_sensruns_onestring #this list has all the runs (varying particles, shapedist and stickeff)
sensrun_list = separated_by_sensruns_onestring.split('_')
shape_spec_list=['nu0mu0.5e'] #this is the default
pprop_list=[]
stick_list=['stick-2'] #-2: Lin83
sensrun_set_dic = dict() #dictionary which contains all the settings (varying particles, shapedist and stickeff)
for sensrun in sensrun_list: #extract an array with different size dist params
    if 'nu' in sensrun:
        shape_spec = 'nu' +sensrun.split('nu')[1] #this gets the specifier of the size dist parameters
        if shape_spec not in shape_spec_list:
            shape_spec_list.append(shape_spec)
        try:
            sensrun_set_dic[sensrun]["sizedist"] = shape_spec
        except: #in the first occurence of sensrun we need to create that entry first
            sensrun_set_dic[sensrun] = {} 
            sensrun_set_dic[sensrun]["sizedist"] = shape_spec
            
        sensrun_set_dic[sensrun]["pprop"] = sensrun.split('nu')[0] #this is the default
        if 'stick' in sensrun: #both: shapedist and stick are different than default
            stick_spec = 'stick' +sensrun.split('stick')[1] #this gets the specifier of the sticking efficiency
            if stick_spec not in stick_spec:
                stick_list.append(stick_spec)
            sensrun_set_dic[sensrun]["stick"] = stick_spec

        else: 
            sensrun_set_dic[sensrun]["stick"] = stick_list[0] #this is the default
        if not sensrun_set_dic[sensrun]["pprop"] in pprop_list:
            pprop_list.append(sensrun_set_dic[sensrun]["pprop"])
    elif 'stick' in sensrun: #default shapedist but different stick
        stick_spec = 'stick' +sensrun.split('stick')[1] #this gets the specifier of the sticking efficiency
        pprop_spec = sensrun.split('stick')[0] #this gets the specifier of the particle properties 
        try:
            sensrun_set_dic[sensrun]["stick"] = stick_spec
        except:
            sensrun_set_dic[sensrun] = {} 
            sensrun_set_dic[sensrun]["stick"] = stick_spec
        if stick_spec not in stick_list:
            stick_list.append(stick_spec)
        sensrun_set_dic[sensrun]["pprop"] = pprop_spec #this is the default
        sensrun_set_dic[sensrun]["sizedist"] = shape_spec_list[0] #this is the default
        if not sensrun_set_dic[sensrun]["pprop"] in pprop_list:
            pprop_list.append(sensrun_set_dic[sensrun]["pprop"])
    else: #here we have the run with default shape dist and stickeff
        pprop_list.append(sensrun)
        sensrun_set_dic[sensrun] = {} 
        sensrun_set_dic[sensrun]["pprop"] = sensrun #this is the default
        sensrun_set_dic[sensrun]["sizedist"] = shape_spec_list[0] #this is the default
        sensrun_set_dic[sensrun]["stick"] = stick_list[0] #this is the default
axnum = plt.subplot2grid((number_of_plots, 1), (0, 0))
axmass = plt.subplot2grid((number_of_plots, 1), (1, 0))
axmean = plt.subplot2grid((number_of_plots, 1), (2, 0))
axnumflux = plt.subplot2grid((number_of_plots, 1), (3, 0))
axmassflux = plt.subplot2grid((number_of_plots, 1), (4, 0))
linestyles=["-","--"]
colors=["b","r"]
colors_ice=["cyan","orange"]
colors_snow=["navy","maroon"]
markers = ["","x","o"]
sorted_sensrunnames=sorted(sensrun_set_dic.keys(), key=lambda x:x.lower())
for i_sensrun, sensrun_now in enumerate(sorted_sensrunnames): #loop over different SB sensitivity runs
    print sensrun_now
    #read properties (for calculating Dmean)
    p = calc_theoretical_DWRs_MK.init_class()
    if "powerlawOLD" in sensrun_now :
        particle="Mix2SB" #this has m-D from snowSBB and v-D from Mix2
    elif "AtlasJmixcolumndend" in sensrun_now:
        particle="Mix2" #this has m-D from snowSBB and v-D from Mix2
    else:
        print "define particle type for convert to Dmean"
        sys.exit()

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
    #read hei2massdens to get average /summed up values per height
    filestring_hei2massdens = directory + experiment + "/hei2massdens.dat"
    timestep = -1 #tstep/10 #TODO: do not hardcode the .. minute output interval here
    hei2massdens = __postprocess_McSnow.read_hei2massdens(filestring_hei2massdens,timestep=timestep,empty_flag=True) #get empty arrays with the same dimension as the true arrays in order to not plot the McSnow data
       
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
    print sensrun_now,i_prop,pprop_list
    i_stick = stick_list.index(sensrun_set_dic[sensrun_now]['stick']) 
    i_shape = shape_spec_list.index(sensrun_set_dic[sensrun_now]['sizedist']) 
    twomom['mmean'] = (twomom['qi']+twomom['qs'])/(twomom['qni']+twomom['qns']) #calculate mean mass
    twomom['mmeani'] = (twomom['qi'])/(twomom['qni']) #calculate mean mass
    twomom['mmeans'] = (twomom['qs'])/(twomom['qns']) #calculate mean mass
    twomom['q'] = twomom['qi']+twomom['qs']
    twomom['qn'] = twomom['qni']+twomom['qns']
    twomom['f'] = twomom['fi']+twomom['fs']
    twomom['fn'] = twomom['fni']+twomom['fns']
    twomom['Dmean'] = p[particle].a_ms**(-1./p[particle].b_ms)*twomom["mmean"]**(1./p[particle].b_ms)*1e3
    twomom['Dmeani'] = p[particle].a_ms**(-1./p[particle].b_ms)*twomom["mmeani"]**(1./p[particle].b_ms)*1e3
    twomom['Dmeans'] = p[particle].a_ms**(-1./p[particle].b_ms)*twomom["mmeans"]**(1./p[particle].b_ms)*1e3
    #debug() 
    ####
    #plotting
    ####
    for var,ax,varlabel in zip(['qn','q','Dmean','fn','f'],[axnum,axmass,axmean,axnumflux,axmassflux],["number density [1/m3]","mass density [kg/m3]","mean diameter [mm]","number flux [1/m3/s]","mass flux [kg/m3/s]"]): #these two lines do all the plotting
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
        
        ax.semilogx(twomom[var][i_timestep,:],twomom['heights'],color=colors[i_prop],linestyle=linestyles[i_stick],marker=markers[i_shape],markevery=20,label=sensrun_now_label) 
        ax.semilogx(twomom[var + 'i'][i_timestep,:],twomom['heights'],color=colors_ice[i_prop],linestyle=linestyles[i_stick],marker=markers[i_shape],markevery=20,alpha=0.3,label=sensrun_now_label+ "(ice)") 
        ax.semilogx(twomom[var + 's'][i_timestep,:],twomom['heights'],color=colors_snow[i_prop],linestyle=linestyles[i_stick],marker=markers[i_shape],markevery=20,alpha=0.3,label=sensrun_now_label + "(snow)") 
        ax.set_xlabel(varlabel)
        ax.set_ylabel("height [m]")


#activate the legend
for ax in [axnum,axmass,axmean,axnumflux,axmassflux]:
    ax.legend()

#save figure
plt.tight_layout()
if not os.path.exists('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment): #create directory if it does not exists
    os.makedirs('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment)
out_filestring = "/fluxes_" + switch_off_processes_str + '_' + separated_by_sensruns_onestring + '_' + testcase + "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'])
