'''
plots height profiles of some variables (mean monomer number)
'''

from IPython.core.debugger import Tracer ; debug = Tracer() #insert this line somewhere to debug

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
from matplotlib import rc
#read variables passed by shell script
tstep = int(os.environ["tstep"])
if "tstep_end" in  os.environ.keys(): 
    tstep_end = int(os.environ["tstep_end"]) #must be there if a temporal average over output steps is wanted
else:
    tstep_end = 0 

experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"]
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output
MC_dir = os.environ["MCexp"]
skipMC = (os.environ["skipMC"]=="True") #allows to run the scripts also if no McSnow data is there (only 1D-SB runs)
skipSB = (os.environ["skipSB"]=="True") #allows to hide plots for SB
plot_totalice = (os.environ["plot_totalice"]=="True")

#0. N_mono dependency 
#read bool
if "separated_by_geom" in os.environ.keys():
    separated_by_geom= (os.environ["separated_by_geom"]=="True") #plot a line for different sensitivity runs
    #read string for sensitivity runs
    if "McSnow_geom_specifier_onestring" in os.environ.keys():
        McSnow_geom_specifier_onestring= os.environ["McSnow_geom_specifier_onestring"]
    else:
        print "error in 4paper_plot_fluxes_and_densities_MCsensruns.py: McSnow_geom_specifier_onestring not found but separated_by_geom=True" 
else:
    separated_by_geom =False
    McSnow_geom_specifier_onestring=''

#1. vterm fit
#read bool
if "separated_by_fallspeedsens" in os.environ.keys():
    separated_by_fallspeedsens = (os.environ["separated_by_fallspeedsens"]=="True") #plot a line for different sensitivity runs
    #read string for sensitivity runs
    if "MCtermvel_specifier_onestring" in os.environ.keys():
        MCtermvel_specifier_onestring = os.environ["MCtermvel_specifier_onestring"]
    else:
        print "error in 4paper_plot_fluxes_and_densities_MCsensruns.py: MCtermvel_specifier_onestring not found but separated_by_fallspeedsens=True" 
else:
    separated_by_fallspeedsens = False
    MCtermvel_specifier_onestring=''

#2. habits 
#read bool
if "separated_by_habit" in os.environ.keys():
    separated_by_habit= (os.environ["separated_by_habit"]=="True") #plot a line for different sensitivity runs
else:
    separated_by_habit=False
#read string for sensitivity runs
if "habits_onestring" in os.environ.keys():
    habits_onestring= os.environ["habits_onestring"]
else:
    print "error in 4paper_plot_fluxes_and_densities_MCsensruns.py:habits_onestring not found but separated_by_habit=True" 
#else:
#    separated_by_habit=False
#    habits_onestring=''

#switch of processes (for McSnow this is so far not used)
if "switch_off_processes" in os.environ.keys():
    switch_off_processes_str = os.environ["switch_off_processes"]
else:
    switch_off_processes_str = ''

#control microphysical processes
if "depos_from0mto" in os.environ.keys():
    depos_from0mto = os.environ["depos_from0mto"]
else:
    depos_from0mto = ''
if "riming_fromdeposto" in os.environ.keys():
    riming_fromdeposto = os.environ["riming_fromdeposto"]
else:
    riming_fromdeposto = ''

def weighted_median(data, weights):
    """
    Args:
      data (list or numpy.array): data
      weights (list or numpy.array): weights
    """
    data, weights = np.array(data).squeeze(), np.array(weights).squeeze()
    s_data, s_weights = map(np.array, zip(*sorted(zip(data, weights))))
    midpoint = 0.5 * sum(s_weights)
    if any(weights > midpoint):
        w_median = (data[weights == np.max(weights)])[0]
    else:
        cs_weights = np.cumsum(s_weights)
        idx = np.where(cs_weights <= midpoint)[0][-1]
        if cs_weights[idx] == midpoint:
            w_median = np.mean(s_data[idx:idx+2])
        else:
            w_median = s_data[idx+1]
    return w_median

#directory of experiments
directory = MC_dir + "/experiments/"

number_of_plots = 1

#fig,ax	=	plt.figure(figsize=(4.0,10.0*number_of_plots/6.))#figsize=(4, 4))

#get timestep
i_timestep=(tstep/60)-1 #there is no output for t=0min after that there are output steps in 10 minute steps (this could vary)
i_timestep_end=(tstep_end/60)-1

#plot different lines for different McSnow geometries
McSnow_geom_list = McSnow_geom_specifier_onestring.split('_')
MCtermvel_list = MCtermvel_specifier_onestring.split('_')
MCsnowhabit_list = habits_onestring.split('_')

#setup axes
#ax1   = plt.subplot2grid((3, 2), (0, 0))
fig,ax1   = plt.subplots()

precip_string = "surface precip. rate\n[kg m-2 h-1]" #header of precipitation annotation
mmean_string = "$m_{mean}$\n[$\mu g$]]" #header of precipitation annotation
for i_sensMC, sensrun_now_MC in  enumerate(McSnow_geom_list):#[McSnow_geom_list[0]]):show first only #loop over different McSnow geometry sensitivity runs
    for i_sensMCfallspeed, sensrun_now_MC_fallspeed in enumerate(MCtermvel_list): #loop over different fall speed models (boehm, KC05,HW10,powerlaw, powerlawSB)
        for i_habit, habit in enumerate(MCsnowhabit_list):
            linestyleorder=['-','--',':','-.']
            if separated_by_geom: #i_mode=0
                if len(McSnow_geom_list)>1:#modify the experiment string to read the file with the write geometry
                    experiment_splitted = experiment.split('_',4) #this results e.g. in ['1d', 'powerlawJplate', 'xi10000000_
                    experiment_splitted[3] = sensrun_now_MC #replace experiment part
                    experiment = "_".join(experiment_splitted)
                    #raw_input("sensrun_now_MC: " + sensrun_now_MC + " " + str(experiment_splitted))
                    McSnow_geom_list_str = McSnow_geom_specifier_onestring
                else:
                    McSnow_geom_list_str = str(McSnow_geom_list[0])

            if separated_by_fallspeedsens:#i_mode=1
                if sensrun_now_MC_fallspeed=="powerlawlimit": #"fallspeed model" is the old powerlaw
                    experiment_splitted = experiment.split('_vt',1) #this results e.g. in [_rm10_rt0_', '3_at2_
                    experiment = "_".join([experiment_splitted[0],'vt7',experiment_splitted[1][2:]]) #[1:] cuts the old velocity index plus the following '_'                
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
            if separated_by_habit:#i_mode=2
                #split string to insert the habit information
                experiment_splitted = experiment.split('_',4)
                exp_after_xi= experiment.split('_xi',1)
                if habit=="plate": #default and does not appear in string
                    experiment = "_".join(experiment_splitted[0:4]) + "_xi" + exp_after_xi[1]
                else:
                    experiment = "_".join(experiment_splitted[0:4]) + "_" + habit + "_xi" + exp_after_xi[1]
            #read hei2massdens to get average / summed up values per height
            filestring_hei2massdens = directory + experiment + "/hei2massdens.dat"
            print habit,filestring_hei2massdens
            timestep = tstep/10 #TODO: do not hardcode the 10 minute output interval here
            timestep_end = tstep_end/10
            hei2massdens = __postprocess_McSnow.read_hei2massdens(filestring_hei2massdens,timestep=timestep,timestep_end=timestep_end)
            model_top_height=hei2massdens["z"][-1]

            #read mass2fr.dat
            filestring_SP = 'mass2fr_' + str(tstep).zfill(4) + '-' + str(tstep_end).zfill(4) + 'min_avtstep_' + str(av_tstep) + '.ncdf'
            SP = __postprocess_McSnow.read_mass2frncdf(MC_dir + 'experiments/' + experiment,filestring_SP)

            height_bounds=np.linspace(0,5000,100)
            profiles = dict()
            profiles["mm_median"] = np.nan*np.ones_like(height_bounds)
            profiles["mm_mean"] = np.nan*np.ones_like(height_bounds)
            var = "mm" #monomer number

            for i_height,dummy in enumerate(height_bounds[:-1]):
                in_bin = ((SP['height']>height_bounds[i_height]) & (SP['height']<height_bounds[i_height+1]))
                #profiles["mm_median"][i_height] = np.median(SP[var][in_bin]) 
                profiles["mm_mean"][i_height] = np.average(SP[var][in_bin],weights=SP["xi"][in_bin]) 
                profiles["mm_median"][i_height] = weighted_median(SP[var][in_bin],SP["xi"][in_bin]) 
                print profiles["mm_mean"][i_height],profiles["mm_median"][i_height] #; debug()


            ############
            #plot mean Nmono 
            ############
            var_flag = 2
                
            #ax1 = __plotting_functions.plot_MC_profiles(ax1,hei2massdens,i_timestep,var_flag=var_flag,forced_linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],forced_markerMC=['','','',''][i_sensMC+i_sensMCfallspeed+i_habit],top_height=model_top_height,plot_totalice=plot_totalice,catonly='_mm1',show_label=False)
            ax1.semilogx(profiles["mm_median"],height_bounds,linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],color="orange")
            ax1.semilogx(profiles["mm_mean"],height_bounds,linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],color="r")

            

            
if len(McSnow_geom_list)>1:#add labels for the different sensruns
    for i_sensMC,sensrun_now_MC in enumerate(McSnow_geom_list):
        for ax in [ax1]:
            labellist=["monodep/CTRL",McSnow_geom_list[1],McSnow_geom_list[2]]
            ax.plot(np.nan,np.nan,color='k',linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],label=labellist[i_sensMC+i_sensMCfallspeed+i_habit])
    #labels for mean and median
    ax.plot(np.nan,np.nan,color='orange',linestyle="-",label="median")
    ax.plot(np.nan,np.nan,color='r',linestyle="-",label="mean")
    #axes labels
    ax1.set_xlabel("Number of monomers")
    ax1.set_ylabel("Height [m]")
    ax1.legend()
    ax1.grid(b=True,which="both")
    ax1.set_xlim([0.9,1e3])
    ax1.set_ylim([0,5e3])
if len(MCtermvel_list)>1:#add labels for the different sensruns #here only on the upper left plot
    for i_sensMCfallspeed, sensrun_now_MC_fallspeed in enumerate(MCtermvel_list): #loop over different fall speed models (boehm, KC05,HW10,powerlaw, powerlawSB)
        for ax in [ax1]: #[axmean,axnumflux,axmassflux]:
            if sensrun_now_MC_fallspeed=="boehm": 
                sensrun_now_MC_fallspeed="binary" #change the label to improve understandability of the two Figures
            ax.plot(np.nan,np.nan,color='k',linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],label=sensrun_now_MC_fallspeed) #,marker=['','','',''][i_sensMC+i_sensMCfallspeed+i_habit])
            ax.legend(loc="lower right")
#make headers
#for i_ax,ax in enumerate([ax1]):
#    ax.text(0.02,0.98,["a)","b)","c)","d)","e)","f)"][i_ax],fontsize=12,fontweight='bold',horizontalalignment='left',
#             verticalalignment='top',
#             transform = ax.transAxes)
    #remove legend in all subplots but a)
    if i_ax!=0:
        ax.get_legend().remove()


#save figure
plt.tight_layout()

if not depos_from0mto=='0': #we dont have simple SDA-runs, so we need different labels
    processes_add = "_deposz" + depos_from0mto + "_rimez" + riming_fromdeposto
else:
    processes_add = "_"
if not os.path.exists('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment): #create directory if it does not exists
    os.makedirs('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment)
    
    
out_filestring = "/additional_" + switch_off_processes_str + '_' + McSnow_geom_specifier_onestring + MCtermvel_specifier_onestring + habits_onestring + '_' + testcase + processes_add + "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'])

if separated_by_fallspeedsens:
    #copy to plots4paper folder
    subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf /home/mkarrer/Dokumente/plots/4paper/McSnow_run_profiles_bulkparam_add.pdf',shell=True)
    subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png /home/mkarrer/Dokumente/plots/4paper/McSnow_run_profiles_bulkparam_add.png',shell=True)
elif separated_by_geom:
    #copy to plots4paper folder
    subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png /home/mkarrer/Dokumente/plots/4paper/McSnow_run_profiles_monodepbinaryconst_add.png',shell=True)
    subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf /home/mkarrer/Dokumente/plots/4paper/McSnow_run_profiles_monodepbinaryconst_add.pdf',shell=True)
elif separated_by_habit:
    #copy to plots4paper folder
    subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf /home/mkarrer/Dokumente/plots/4paper/McSnow_run_profiles_habits_add.pdf',shell=True)
    subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png /home/mkarrer/Dokumente/plots/4paper/McSnow_run_profiles_habits_add.png',shell=True)
