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
if "tstep_end" in  os.environ.keys(): 
    tstep_end = int(os.environ["tstep_end"]) #must be there if a temporal average over output steps is wanted
else:
    tstep_end = 0 

experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"]
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output
MC_dir = os.environ["MC"]
skipMC = (os.environ["skipMC"]=="True") #allows to run the scripts also if no McSnow data is there (only 1D-SB runs)
skipSB = (os.environ["skipSB"]=="True") #allows to hide plots for SB

if "McSnow_geom_specifier_onestring" in  os.environ.keys():
    McSnow_geom_specifier_onestring = os.environ["McSnow_geom_specifier_onestring"]
else:
    McSnow_geom_specifier_onestring = ''
if "separated_by_sensrunsMC" in os.environ.keys():
    separated_by_sensrunsMC = (os.environ["separated_by_sensrunsMC"]=="True") #plot a line for different sensitivity runs
    separated_by_sensruns_onestring= os.environ["model_setup_specifier_onestring"]
else:
    print "separated_by_sensruns not found in environ: set to False"
    separated_by_sensruns = False
if "separated_by_fallspeedsens" in os.environ.keys():
    separated_by_fallspeedsens = (os.environ["separated_by_fallspeedsens"]=="True") #plot a line for different sensitivity runs
else:
    print "separated_by_fallspeedsens not found in environ: set to False"
    separated_by_fallspeedsens = False
if "switch_off_processes" in os.environ.keys():
    switch_off_processes_str = os.environ["switch_off_processes"]
else:
    switch_off_processes_str = ''
if "MCtermvel_specifier_onestring" in os.environ.keys():
    MCtermvel_specifier_onestring = os.environ["MCtermvel_specifier_onestring"]
else:
    MCtermvel_specifier_onestring = ''
if "depos_from0mto" in os.environ.keys():
    depos_from0mto = os.environ["depos_from0mto"]
else:
    depos_from0mto = ''
if "riming_fromdeposto" in os.environ.keys():
    riming_fromdeposto = os.environ["riming_fromdeposto"]
else:
    riming_fromdeposto = ''
#directory of experiments
directory = MC_dir + "/experiments/"

number_of_plots = 5

figsize_height = 6.0/2.0*(number_of_plots)
fig	=	plt.figure(figsize=(8.0,figsize_height))#figsize=(4, 4))

#get timestep
i_timestep=(tstep/60)-1 #there is no output for t=0min after that there are output steps in 30 minute steps (this could vary)
i_timestep_end=(tstep_end/60)-1
#plot different lines for different sensitivity runs
if separated_by_sensrunsMC:
    sensrun_list = separated_by_sensruns_onestring
    if separated_by_sensruns_onestring=="_": #this is just a placeholder -> no sensruns -> dont split
        sensrun_list = [""]
    else:
        sensrun_list = separated_by_sensruns_onestring.split('_')
else:
    sensrun_list = [""]
#plot different lines for different McSnow geometries
McSnow_geom_list = McSnow_geom_specifier_onestring.split('_')
MCtermvel_list = MCtermvel_specifier_onestring.split('_')

axnum = plt.subplot2grid((number_of_plots, 1), (0, 0))
axmass = plt.subplot2grid((number_of_plots, 1), (1, 0))
axnumflux = plt.subplot2grid((number_of_plots, 1), (2, 0))
axmassflux = plt.subplot2grid((number_of_plots, 1), (3, 0))
axmean = plt.subplot2grid((number_of_plots, 1), (4, 0))

precip_string = "surface precip. rate\n[kg m-2 h-1]" #header of precipitation annotation
for i_sensrun, sensrun_now in enumerate(sensrun_list): #loop over different SB sensitivity runs
    for i_sensMC, sensrun_now_MC in  enumerate(McSnow_geom_list): #loop over different McSnow geometry sensitivity runs
        for i_sensMCfallspeed, sensrun_now_MC_fallspeed in enumerate(MCtermvel_list): #loop over different fall speed models (boehm, KC05,HW10,powerlaw, powerlawSB)
            if separated_by_sensrunsMC:#modify the experiment string
                #from IPython.core.debugger import Tracer ; Tracer()()
                #raw_input(sensrun_now_MC_fallspeed)
                if separated_by_fallspeedsens:
                    linestyleorder=['--','-','-.',':']
                    if sensrun_now_MC_fallspeed=="powerlawSBdefault": #"fallspeed model" is the old powerlaw
                        experiment_splitted = experiment.split('_vt',1) #this results e.g. in [_rm10_rt0_', '3_at2_
                        experiment = "_".join([experiment_splitted[0],'vt5',experiment_splitted[1][2:]]) #[1:] cuts the old velocity index plus the following '_'                
                    elif sensrun_now_MC_fallspeed=="powerlaw": ##"fallspeed model" is the new powerlaw
                        experiment_splitted = experiment.split('_vt',1) #this results e.g. in [_rm10_rt0_', '3_at2_
                        experiment = "_".join([experiment_splitted[0],'vt4',experiment_splitted[1][2:]]) #[1:] cuts the old velocity index
                        print "powerlaw",experiment
                    elif sensrun_now_MC_fallspeed=="Atlas": ##"fallspeed model" is the new powerlaw
                        experiment_splitted = experiment.split('_vt',1) #this results e.g. in [_rm10_rt0_', '3_at2_
                        experiment = "_".join([experiment_splitted[0],'vt6',experiment_splitted[1][2:]]) #[1:] cuts the old velocity index
                        print "powerlaw",experiment
                    elif sensrun_now_MC_fallspeed=="boehm": #other sensruns (not modifying the fallspeed model
                        #fix the fall speed model again
                        experiment_splitted = experiment.split('_vt',1) #this results e.g. in [_rm10_rt0_', '3_at2_
                        experiment = "_".join([experiment_splitted[0],'vt3',experiment_splitted[1][2:]]) #[1:] cuts the old velocity index #ATTENTION: vt2 and vt1 not implemented
                    McSnow_geom_list_str = MCtermvel_specifier_onestring
                else:
                    linestyleorder=['-','--','-.',':']

                    #other sensruns (not modifying the fallspeed model
                    #experiment_splitted = experiment.split('_',2) #this results e.g. in ['1d', 'powerlawJplate', 'xi10000000_
                    #experiment_splitted[1] = sensrun_now #replace experiment part
                    #experiment = "_".join(experiment_splitted)
                
            
                    if len(McSnow_geom_list)>1:#modify the experiment string to read the file with the write geometry
                        experiment_splitted = experiment.split('_',4) #this results e.g. in ['1d', 'powerlawJplate', 'xi10000000_
                        experiment_splitted[3] = sensrun_now_MC #replace experiment part
                        experiment = "_".join(experiment_splitted)
                        #from IPython.core.debugger import Tracer ; Tracer()()
                        #raw_input("sensrun_now_MC: " + sensrun_now_MC + " " + str(experiment_splitted))
                        McSnow_geom_list_str = McSnow_geom_specifier_onestring
                    else:
                    
                        McSnow_geom_list_str = str(McSnow_geom_list[0])
            else:
                linestyleorder=['-','--','-.',':']
                McSnow_geom_list_str=McSnow_geom_specifier_onestring
                
        

            #read hei2massdens to get average / summed up values per height
            filestring_hei2massdens = directory + experiment + "/hei2massdens.dat"
            print filestring_hei2massdens
            timestep = tstep/10 #TODO: do not hardcode the 10 minute output interval here
            timestep_end = tstep_end/10
            if skipMC:
                hei2massdens = __postprocess_McSnow.read_hei2massdens(filestring_hei2massdens,timestep=timestep,empty_flag=True) #get empty arrays with the same dimension as the true arrays in order to not plot the McSnow data
            else:
                hei2massdens = __postprocess_McSnow.read_hei2massdens(filestring_hei2massdens,timestep=timestep,timestep_end=timestep_end)
            model_top_height=hei2massdens["z"][-1]

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
                #do a skipSB similar to skipMC above
                if skipSB and var!="heights":
                    #create empty variables so we dont see them in the plot
                    twomom[var] = np.nan*np.ones_like(twomom_file.variables[var])
                else:
                    twomom[var] = np.squeeze(twomom_file.variables[var])

            #from IPython.core.debugger import Tracer ; Tracer()()
            ####################################
            #plot mixing ratio + number density + normalized mixing ratio
            ####################################
            #number density
            var_flag = 0
            
            axnum = __plotting_functions.plot_MC_profiles(axnum,hei2massdens,i_timestep,var_flag=var_flag,forced_linestyle=linestyleorder[i_sensrun+i_sensMC+i_sensMCfallspeed],forced_markerMC=['','','',''][i_sensrun+i_sensMC+i_sensMCfallspeed],top_height=model_top_height)
            

            #mass density
            var_flag = 1
                
            axmass = __plotting_functions.plot_MC_profiles(axmass,hei2massdens,i_timestep,var_flag=var_flag,forced_linestyle=linestyleorder[i_sensrun+i_sensMC+i_sensMCfallspeed],forced_markerMC=['','','',''][i_sensrun+i_sensMC+i_sensMCfallspeed],top_height=model_top_height)



            ############
            #plot fluxes
            ############
            #number flux
            var_flag = 2
                
            axnumflux = __plotting_functions.plot_MC_profiles(axnumflux,hei2massdens,i_timestep,var_flag=var_flag,forced_linestyle=linestyleorder[i_sensrun+i_sensMC+i_sensMCfallspeed],forced_markerMC=['','','',''][i_sensrun+i_sensMC+i_sensMCfallspeed],top_height=model_top_height)

            #mass flux
            var_flag = 3
                
            axmassflux = __plotting_functions.plot_MC_profiles(axmassflux,hei2massdens,i_timestep,var_flag=var_flag,forced_linestyle=linestyleorder[i_sensrun+i_sensMC+i_sensMCfallspeed],forced_markerMC=['','','',''][i_sensrun+i_sensMC+i_sensMCfallspeed],top_height=model_top_height)

            #calculate  Md/Nd in __plotting_functions.plot_MC_profiles
            var_flag = 4
        
            axmean = __plotting_functions.plot_MC_profiles(axmean,hei2massdens,i_timestep,var_flag=var_flag,forced_linestyle=linestyleorder[i_sensrun+i_sensMC+i_sensMCfallspeed],forced_markerMC=['','','',''][i_sensrun+i_sensMC+i_sensMCfallspeed],top_height=model_top_height)
            
            #from IPython.core.debugger import Tracer ; Tracer()()
            if separated_by_fallspeedsens:
                if i_sensMCfallspeed==0:
                    reference_precip = hei2massdens["Fm"][0]*3600. #save precip rate from reference run (always the first) to calculate relative deviations
                    rel_diff=""
                else:
                    rel_diff=" ({:+.1f}%)".format((hei2massdens["Fm"][0]*3600.-reference_precip)/reference_precip*100)
                if sensrun_now_MC_fallspeed=="boehm": 
                    sensrun_now_MC_fallspeed="binary" #change the label to improve understandability of the two Figures
                precip_string+="\n" + str(sensrun_now_MC_fallspeed) + ": " + '{:.3f}'.format(hei2massdens["Fm"][0]*3600.) + rel_diff
            else:
                if i_sensMC==0:
                    reference_precip = hei2massdens["Fm"][0]*3600. #save precip rate from reference run (always the first) to calculate relative deviations
                    rel_diff=""
                else:
                    rel_diff=" ({:+.1f}%)".format((hei2massdens["Fm"][0]*3600.-reference_precip)/reference_precip*100)
                precip_string+="\n" + str(McSnow_geom_list[i_sensMC]) + ": " + '{:.3f}'.format(hei2massdens["Fm"][0]*3600.) + rel_diff
if len(sensrun_list)>1:#add labels for the different sensruns
    for i_sensrun,sens_run in enumerate(sensrun_list):
        for ax in [axnum,axmass,axmean,axnumflux,axmassflux]:
            ax.plot(np.nan,np.nan,color='k',linestyle=linestyleorder[i_sensrun+i_sensMC+i_sensMCfallspeed],label=sens_run)
            ax.legend()
            
if len(McSnow_geom_list)>1:#add labels for the different sensruns
    for i_sensMC,sensrun_now_MC in enumerate(McSnow_geom_list):
        for ax in [axnum,axmass,axmean,axnumflux,axmassflux]:
            ax.plot(np.nan,np.nan,color='k',linestyle=linestyleorder[i_sensrun+i_sensMC+i_sensMCfallspeed],label=McSnow_geom_list[i_sensrun+i_sensMC+i_sensMCfallspeed])
            ax.legend()
if len(MCtermvel_list)>1:#add labels for the different sensruns
    for i_sensMCfallspeed, sensrun_now_MC_fallspeed in enumerate(MCtermvel_list): #loop over different fall speed models (boehm, KC05,HW10,powerlaw, powerlawSB)
        for ax in [axnum,axmass,axmean,axnumflux,axmassflux]:
            if sensrun_now_MC_fallspeed=="boehm": 
                sensrun_now_MC_fallspeed="binary" #change the label to improve understandability of the two Figures
            ax.plot(np.nan,np.nan,color='k',linestyle=linestyleorder[i_sensrun+i_sensMC+i_sensMCfallspeed],label=sensrun_now_MC_fallspeed,marker=['','','',''][i_sensrun+i_sensMC+i_sensMCfallspeed])
            ax.legend()
#add precipitation string
axmassflux.text(0.0,1.0,precip_string,fontsize=10,horizontalalignment='left',verticalalignment='top',transform = axmassflux.transAxes)
#save figure
plt.tight_layout()

if not depos_from0mto=='0': #we dont have simple SDA-runs, so we need different labels
    processes_add = "_deposz" + depos_from0mto + "_rimez" + riming_fromdeposto
else:
    processes_add = "_"
if not os.path.exists('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment): #create directory if it does not exists
    os.makedirs('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment)
out_filestring = "/fluxes_" + switch_off_processes_str + '_' + McSnow_geom_list_str + separated_by_sensruns_onestring + '_' + testcase + processes_add + "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'])

if separated_by_fallspeedsens:
    #copy to plots4paper folder
    subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf /home/mkarrer/Dokumente/plots/4paper/McSnow_run_profiles_bulkparam.pdf',shell=True)
else:
    #copy to plots4paper folder
    subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf /home/mkarrer/Dokumente/plots/4paper/McSnow_run_profiles_monodepbinaryconst.pdf',shell=True)