'''
plots the fluxes and the absolute values of the number density and the mass density
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

#directory of experiments
directory = MC_dir + "/experiments/"

#number_of_plots = 6

#figsize_height = 12.0/3.0*(number_of_plots)
fig	=	plt.figure(figsize=(8.0,10.0))#figsize=(4, 4))

#get timestep
i_timestep=(tstep/60)-1 #there is no output for t=0min after that there are output steps in 10 minute steps (this could vary)
i_timestep_end=(tstep_end/60)-1

#plot different lines for different McSnow geometries
McSnow_geom_list = McSnow_geom_specifier_onestring.split('_')
MCtermvel_list = MCtermvel_specifier_onestring.split('_')
MCsnowhabit_list = habits_onestring.split('_')

#setup axes
axnumflux_ice   = plt.subplot2grid((3, 2), (0, 0))
axnumflux_snow  = plt.subplot2grid((3, 2), (0, 1))
axmassflux_ice  = plt.subplot2grid((3, 2), (1, 0))
axmassflux_snow = plt.subplot2grid((3, 2), (1, 1))
axmean_ice      = plt.subplot2grid((3, 2), (2, 0))
axmean_snow     = plt.subplot2grid((3, 2), (2, 1))

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


            ############
            #plot fluxes
            ############
            #number flux
            var_flag = 2
                
            axnumflux_ice = __plotting_functions.plot_MC_profiles(axnumflux_ice,hei2massdens,i_timestep,var_flag=var_flag,forced_linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],forced_markerMC=['','','',''][i_sensMC+i_sensMCfallspeed+i_habit],top_height=model_top_height,plot_totalice=plot_totalice,catonly='_mm1',show_label=False)

            axnumflux_snow = __plotting_functions.plot_MC_profiles(axnumflux_snow,hei2massdens,i_timestep,var_flag=var_flag,forced_linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],forced_markerMC=['','','',''][i_sensMC+i_sensMCfallspeed+i_habit],top_height=model_top_height,plot_totalice=plot_totalice,catonly='_unr',show_label=False)
            #mass flux
            var_flag = 3
                
            axmassflux_ice  = __plotting_functions.plot_MC_profiles(axmassflux_ice,hei2massdens,i_timestep,var_flag=var_flag,forced_linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],forced_markerMC=['','','',''][i_sensMC+i_sensMCfallspeed+i_habit],top_height=model_top_height,plot_totalice=plot_totalice,catonly='_mm1',show_label=False)

            axmassflux_snow  = __plotting_functions.plot_MC_profiles(axmassflux_snow,hei2massdens,i_timestep,var_flag=var_flag,forced_linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],forced_markerMC=['','','',''][i_sensMC+i_sensMCfallspeed+i_habit],top_height=model_top_height,plot_totalice=plot_totalice,catonly='_unr',show_label=False)
            #calculate  Md/Nd in __plotting_functions.plot_MC_profiles
            var_flag = 4
        
            axmean_ice = __plotting_functions.plot_MC_profiles(axmean_ice,hei2massdens,i_timestep,var_flag=var_flag,forced_linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],forced_markerMC=['','','',''][i_sensMC+i_sensMCfallspeed+i_habit],top_height=model_top_height,plot_totalice=plot_totalice,catonly='_mm1',show_label=False)
            axmean_snow = __plotting_functions.plot_MC_profiles(axmean_snow,hei2massdens,i_timestep,var_flag=var_flag,forced_linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],forced_markerMC=['','','',''][i_sensMC+i_sensMCfallspeed+i_habit],top_height=model_top_height,plot_totalice=plot_totalice,catonly='_unr',show_label=False)
            
            #mean mass at ground [mu g]
            meanmassground=hei2massdens["Md"][0]/hei2massdens["Nd"][0]*1e6
            if experiment.startswith("1d___monodep"):#this should be the control run
                fulldic={"0": hei2massdens}
                np.savez("CTRL", **fulldic)
                rel_diff=""
                rel_diff_mmean=""
                factor_mmean=""
            else:
                #read full dictionary to plot CTRL profiles
                CTRLdic = np.load("CTRL.npz")
                hei2massdensCTRL = CTRLdic['0'].item()
                reference_precip = precip_rate=hei2massdensCTRL["Fm"][0]*3600. 
                rel_diff=" ({:+.1f}%)".format((hei2massdens["Fm"][0]*3600.-reference_precip)/reference_precip*100)
                
                reference_mmean = hei2massdensCTRL["Md"][0]/hei2massdensCTRL["Nd"][0]*1e6 
                rel_diff_mmean=" ({:+.1f}%)".format((meanmassground-reference_mmean)/reference_mmean*100)
                factor_mmean=" (x{:.1f})".format(meanmassground/reference_mmean)
                
                #debug() 
                if "vt6" in experiment: #vt6 is Atlas type #overlay CTRL run for fallspeedsens
                    ############
                    #plot fluxes of CTRL run
                    ############
                    #number flux
                    var_flag = 2
                    axnumflux_ice = __plotting_functions.plot_MC_profiles(axnumflux_ice,hei2massdensCTRL,i_timestep,var_flag=var_flag,forced_linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],forced_markerMC=['','','',''][i_sensMC+i_sensMCfallspeed+i_habit],top_height=model_top_height,plot_totalice=plot_totalice,catonly='_mm1',show_label=False,force_gray=True)

                    axnumflux_snow = __plotting_functions.plot_MC_profiles(axnumflux_snow,hei2massdensCTRL,i_timestep,var_flag=var_flag,forced_linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],forced_markerMC=['','','',''][i_sensMC+i_sensMCfallspeed+i_habit],top_height=model_top_height,plot_totalice=plot_totalice,catonly='_unr',show_label=False,force_gray=True)
                    #mass flux
                    var_flag = 3
                        
                    axmassflux_ice  = __plotting_functions.plot_MC_profiles(axmassflux_ice,hei2massdensCTRL,i_timestep,var_flag=var_flag,forced_linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],forced_markerMC=['','','',''][i_sensMC+i_sensMCfallspeed+i_habit],top_height=model_top_height,plot_totalice=plot_totalice,catonly='_mm1',show_label=False,force_gray=True)

                    axmassflux_snow  = __plotting_functions.plot_MC_profiles(axmassflux_snow,hei2massdensCTRL,i_timestep,var_flag=var_flag,forced_linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],forced_markerMC=['','','',''][i_sensMC+i_sensMCfallspeed+i_habit],top_height=model_top_height,plot_totalice=plot_totalice,catonly='_unr',show_label=False,force_gray=True)
                    #calculate  Md/Nd in __plotting_functions.plot_MC_profiles
                    var_flag = 4
                
                    axmean_ice = __plotting_functions.plot_MC_profiles(axmean_ice,hei2massdensCTRL,i_timestep,var_flag=var_flag,forced_linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],forced_markerMC=['','','',''][i_sensMC+i_sensMCfallspeed+i_habit],top_height=model_top_height,plot_totalice=plot_totalice,catonly='_mm1',show_label=False,force_gray=True)
                    axmean_snow = __plotting_functions.plot_MC_profiles(axmean_snow,hei2massdensCTRL,i_timestep,var_flag=var_flag,forced_linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],forced_markerMC=['','','',''][i_sensMC+i_sensMCfallspeed+i_habit],top_height=model_top_height,plot_totalice=plot_totalice,catonly='_unr',show_label=False,force_gray=True)
                    
                    if (i_sensMC+i_sensMCfallspeed+i_habit)==0:
                        axnumflux_ice.plot(np.nan,np.nan,color='gray',label="CTRL")

            if separated_by_fallspeedsens:
                '''
                if i_sensMCfallspeed==0:
                    #reference_precip = hei2massdens["Fm"][0]*3600. #save precip rate from reference run (always the first) to calculate relative deviations
                    rel_diff=""
                else:
                    rel_diff=" ({:+.1f}%)".format((hei2massdens["Fm"][0]*3600.-reference_precip)/reference_precip*100)
                if sensrun_now_MC_fallspeed=="boehm": 
                    sensrun_now_MC_fallspeed="binary" #change the label to improve understandability of the two Figures
                '''
                precip_string+="\n" + str(sensrun_now_MC_fallspeed) + ": " + '{:.3f}'.format(hei2massdens["Fm"][0]*3600.) + rel_diff
                mmean_string+="\n" + str(sensrun_now_MC_fallspeed) + ": " + '{:.3f}'.format(meanmassground) +  factor_mmean
                factor_mmean=" (x{:.1f})".format(meanmassground/reference_mmean)

            else:
                '''
                if i_sensMC==0:
                    reference_precip = hei2massdens["Fm"][0]*3600. #save precip rate from reference run (always the first) to calculate relative deviations
                    rel_diff=""
                else:
                    rel_diff=" ({:+.1f}%)".format((hei2massdens["Fm"][0]*3600.-reference_precip)/reference_precip*100)
                '''
                precip_string+="\n" + str(McSnow_geom_list[i_sensMC]) + ": " + '{:.3f}'.format(hei2massdens["Fm"][0]*3600.) + rel_diff
                mmean_string+= "\n" + str(McSnow_geom_list[i_sensMC]) + ": " + '{:.3f}'.format(meanmassground) +  factor_mmean
            
if len(McSnow_geom_list)>1:#add labels for the different sensruns
    for i_sensMC,sensrun_now_MC in enumerate(McSnow_geom_list):
        for ax in [axnumflux_ice]:
            labellist=["monodep/CTRL",McSnow_geom_list[1],McSnow_geom_list[2]]
            ax.plot(np.nan,np.nan,color='k',linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],label=labellist[i_sensMC+i_sensMCfallspeed+i_habit])
            #ax.legend(loc="lower right")
            ax.legend(loc="center left")
            ax.set_title("$N_{mono}$=1")
        for ax in [axnumflux_snow]:
            ax.set_title("$N_{mono}$>1")
if len(MCtermvel_list)>1:#add labels for the different sensruns #here only on the upper left plot
    for i_sensMCfallspeed, sensrun_now_MC_fallspeed in enumerate(MCtermvel_list): #loop over different fall speed models (boehm, KC05,HW10,powerlaw, powerlawSB)
        for ax in [axnumflux_ice]: #[axmean,axnumflux,axmassflux]:
            if sensrun_now_MC_fallspeed=="boehm": 
                sensrun_now_MC_fallspeed="binary" #change the label to improve understandability of the two Figures
            ax.plot(np.nan,np.nan,color='k',linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],label=sensrun_now_MC_fallspeed) #,marker=['','','',''][i_sensMC+i_sensMCfallspeed+i_habit])
            ax.legend(loc="lower right")
            ax.set_title("$N_{mono}$=1")
        for ax in [axnumflux_snow]:
            ax.set_title("$N_{mono}$>1")
if len(MCsnowhabit_list)>1:#add labels for the different sensruns
    for i_habit, habit in enumerate(MCsnowhabit_list):
        for ax in [axnumflux_ice]:
            labellist=MCsnowhabit_list
            if habit=="mixcolumndend":
                labelhabit="mix2"
            else:
                labelhabit=habit
            ax.plot(np.nan,np.nan,color='k',linestyle=linestyleorder[i_sensMC+i_sensMCfallspeed+i_habit],label=labelhabit)
            ax.legend(loc="lower left")
            ax.set_title("$N_{mono}$=1")
        for ax in [axnumflux_snow]:
            ax.set_title("$N_{mono}$>1")
#remove yticks from right plots
for ax in [axnumflux_snow,axmassflux_snow,axmean_snow]:
    ax.set_yticks([])
    ax.set_ylabel("")

#make headers
for i_ax,ax in enumerate([axnumflux_ice,axnumflux_snow,axmassflux_ice,axmassflux_snow,axmean_ice,axmean_snow]):
    ax.text(0.02,0.98,["a)","b)","c)","d)","e)","f)"][i_ax],fontsize=12,fontweight='bold',horizontalalignment='left',
             verticalalignment='top',
             transform = ax.transAxes)
    #remove legend in all subplots but a)
    if i_ax!=0:
        ax.get_legend().remove()
#add precipitation string
#axmassflux_ice.text(0.0,1.0,precip_string,fontsize=10,horizontalalignment='left',verticalalignment='top',transform = axmassflux_ice.transAxes)
print '\n\n'
print precip_string 
print '\n\n'
print '\n\n'
print mmean_string 
print '\n\n'
#save figure
plt.tight_layout()

if not depos_from0mto=='0': #we dont have simple SDA-runs, so we need different labels
    processes_add = "_deposz" + depos_from0mto + "_rimez" + riming_fromdeposto
else:
    processes_add = "_"
if not os.path.exists('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment): #create directory if it does not exists
    os.makedirs('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment)
    
'''
out_filestring = "/fluxes_" + switch_off_processes_str + '_' + McSnow_geom_specifier_onestring + MCtermvel_specifier_onestring + habits_onestring + '_' + testcase + processes_add + "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'])
sys.exit(1)
if separated_by_fallspeedsens:
    #copy to plots4paper folder
    subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf /home/mkarrer/Dokumente/plots/4paper/McSnow_run_profiles_bulkparam.pdf',shell=True)
    subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png /home/mkarrer/Dokumente/plots/4paper/McSnow_run_profiles_bulkparam.png',shell=True)
elif separated_by_geom:
    plt.tight_layout()

if not depos_from0mto=='0': #we dont have simple SDA-runs, so we need different labels
    processes_add = "_deposz" + depos_from0mto + "_rimez" + riming_fromdeposto
else:
    processes_add = "_"
if not os.path.exists('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment): #create directory if it does not exists
    os.makedirs('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment)
'''    
    
    
out_filestring = "/fluxes_" + switch_off_processes_str + '_' + McSnow_geom_specifier_onestring + MCtermvel_specifier_onestring + habits_onestring + '_' + testcase + processes_add + "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'])

if separated_by_fallspeedsens:
    #copy to plots4paper folder
    subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf /home/mkarrer/Dokumente/plots/4paper/McSnow_run_profiles_bulkparam.pdf',shell=True)
    subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png /home/mkarrer/Dokumente/plots/4paper/McSnow_run_profiles_bulkparam.png',shell=True)
elif separated_by_geom:
    #copy to plots4paper folder
    subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png /home/mkarrer/Dokumente/plots/4paper/McSnow_run_profiles_monodepbinaryconst.png',shell=True)
    subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf /home/mkarrer/Dokumente/plots/4paper/McSnow_run_profiles_monodepbinaryconst.pdf',shell=True)
elif separated_by_habit:
    #copy to plots4paper folder
    subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf /home/mkarrer/Dokumente/plots/4paper/McSnow_run_profiles_habits.pdf',shell=True)
    subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png /home/mkarrer/Dokumente/plots/4paper/McSnow_run_profiles_habits.png',shell=True)
