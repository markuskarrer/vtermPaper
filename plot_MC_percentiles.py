'''
plot the percentiles of diameter D and terminal velocity v from McSnow simulations resolved by height
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
from functions import __general_utilities
#read variables passed by shell script
tstep = int(os.environ["tstep"]) #start of the averaging period
tstep_end = int(os.environ["tstep_end"]) #string with 4 numbers and 'min'
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"] 
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output
MC_dir = os.environ["MC"]
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
    #read string for sensitivity runs
    if "habits_onestring" in os.environ.keys():
        habits_onestring= os.environ["habits_onestring"]
    else:
        print "error in 4paper_plot_fluxes_and_densities_MCsensruns.py:habits_onestring not found but separated_by_habit=True" 
else:
    separated_by_habit=False
    habits_onestring=''

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

#from https://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy
def weighted_quantile(values, quantiles, sample_weight=None, 
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)



#directory of experiments
directory = MC_dir + "/experiments/"

##general settings for the plotting
number_of_plots = 4

#optimize the appearance of the plot (figure size, fonts)
[fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots,legend_fontsize='medium')

#plot different lines for different McSnow geometries
McSnow_geom_list = McSnow_geom_specifier_onestring.split('_')
MCtermvel_list = MCtermvel_specifier_onestring.split('_')
MCsnowhabit_list = habits_onestring.split('_')

for i_sensMC, sensrun_now_MC in  enumerate(McSnow_geom_list): #loop over different McSnow geometry sensitivity runs
    for i_sensMCfallspeed, sensrun_now_MC_fallspeed in enumerate(MCtermvel_list): #loop over different fall speed models (boehm, KC05,HW10,powerlaw, powerlawSB)
        for i_habit, habit in enumerate(MCsnowhabit_list):
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
            #load file
            SP_file = Dataset(directory + experiment + '/mass2fr_' + str(tstep).zfill(4) + '-' + str(tstep_end).zfill(4) + 'min_avtstep_' + str(av_tstep) + '.ncdf',mode='r')
            #create dictionary for all superparticle variables
            SP = dict()

            #if necessary change name of variables
            varlist = SP_file.variables
            #read PAMTRA variables to pamData dictionary
            for var in varlist:#read files and write it with different names in Data
                SP[var] = np.squeeze(SP_file.variables[var])
            #get maximum height for plotting
            model_top = np.nanmax(SP['height'])

            #define height array
            height_bound_vec= np.linspace(0,model_top,101)
            height_cent_vec = (height_bound_vec[:-1]+height_bound_vec[1:])/2
            #initialize arrays of D and v for profiles
            SP_percentiles = dict()
            quantiles=[0.1,0.5,0.9] #[0.25,0.5,0.75]
            SP_percentiles["diam"] = np.zeros([len(quantiles),height_bound_vec.shape[0]-1]) #three different percentiles at each height #max. dimension
            SP_percentiles["vt"] = np.zeros([len(quantiles),height_bound_vec.shape[0]-1]) #three different percentiles at each height #terminal velocity
            
            debug_upper_bound=False
            if debug_upper_bound:
                prop="diam"
                in_bin = ((SP['height']>model_top-0.01*model_top))
                prop_upper_boundary = SP[prop][in_bin]
                mult_of_prop_upper_boundary = SP["xi"][in_bin]
                if prop_upper_boundary.shape[0]>1:
                    SP_percentiles[prop] =  weighted_quantile(prop_upper_boundary,quantiles,mult_of_prop_upper_boundary) #TODO: that is percentiles of SP not RP!!
                else:
                    SP_percentiles[prop] =  np.nan
                print sensrun_now_MC_fallspeed,SP_percentiles[prop]
                raw_input(); continue
                
            #loop over the above defined height vector
            for i_height,height in enumerate(height_bound_vec[:-1]):
                for key in SP_percentiles.keys():
                    #select only all SP in the corresponding height
                    in_bin = ((SP['height']>height_bound_vec[i_height]) & (SP['height']<height_bound_vec[i_height+1]))
                    prop_in_bin = SP[key][in_bin] #[((SP['height']>height_bound_vec[i_height]) & (SP['height']<height_bound_vec[i_height+1]))]
                    mult_of_prop_in_bin = SP["xi"][in_bin] #[((SP['height']>height_bound_vec[i_height]) & (SP['height']<height_bound_vec[i_height+1]))]
                    #for i_SP_in_bin,in_bin_bool in enumerate(in_bin): #loop over all SP to check whether they are in the right heidht and get their position in the SP dictionary
                    #    if in_bin_bool: #if this
                    #        for mult_now in mult_of_prop_in_bin:
                    #            prop_in_bin = np.append(prop_in_bin,SP[key][i_SP_in_bin])
                    if prop_in_bin.shape[0]>1:
                        SP_percentiles[key][:,i_height] =  weighted_quantile(prop_in_bin,quantiles,mult_of_prop_in_bin) #TODO: that is percentiles of SP not RP!!
                    else:
                        SP_percentiles[key][:,i_height] =  np.nan
                    
                print i_height,height,prop_in_bin.shape, SP_percentiles[key][:,i_height]
                
                

            ######
            ##  plot the percentiles of the simulated diameters
            #####
            i_ax = 0
            for key in ["diam","vt"]:

                #axes[i_ax].semilogx(SP_percentiles[key][2,:],height_cent_vec,color=np.array(['b','r','g'])[i_sensMC+i_sensMCfallspeed+i_habit])
                for i_quantile,quantile_now in enumerate(quantiles):
                    axes[i_ax].plot(SP_percentiles[key][i_quantile,:],height_cent_vec,color=np.array(['r','b','g','y'])[i_sensMC+i_sensMCfallspeed+i_habit],linestyle=np.array(["--","-","--"])[i_quantile],linewidth=np.array([0.5,1.5,0.5])[i_quantile])
                    if key=="diam":
                        axes[i_ax].set_xscale("log")
                    elif key=="vt":
                        axes[i_ax].set_xscale("linear")

                if key=="diam":
                    #make labels
                    axes[i_ax].set_xlabel("diameter D / m")

                    #change the axis
                    axes[i_ax].set_xlim([1e-4,4e-2]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                    axes[i_ax].set_ylim([0,model_top]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                    axes[i_ax].grid(b=True,which="both",axis="both")
                elif key=="vt":
                    #make labels
                    axes[i_ax].set_xlabel("velocity m s-1 / m")
                    axes[i_ax].set_ylabel("height / m" ) #TODO: plot also the untis of these properties

                    #change the axis
                    axes[i_ax].set_xlim([0,2.5]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                    axes[i_ax].set_ylim([0,model_top]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                    axes[i_ax].grid(b=True,which="both",axis="both")
            
                i_ax+=1
                
            ####
            #plot the difference between the percentiles
            ####
            i_ax=2
            for key in ["diam","vt"]:

                #axes[i_ax].semilogx(SP_percentiles[key][2,:],height_cent_vec,color=np.array(['b','r','g'])[i_sensMC+i_sensMCfallspeed+i_habit])
                axes[i_ax].plot(SP_percentiles[key][2,:]-SP_percentiles[key][0,:],height_cent_vec,color=np.array(['r','b','g','y'])[i_sensMC+i_sensMCfallspeed+i_habit],linestyle="-",linewidth=1.5)
                if key=="diam":
                    axes[i_ax].set_xscale("linear")
                elif key=="vt":
                    axes[i_ax].set_xscale("linear")

                if key=="diam":
                    #make labels
                    axes[i_ax].set_xlabel(r"$\Delta$ diameter D / m")

                    #change the axis
                    #axes[i_ax].set_xlim([1e-4,4e-2]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                    axes[i_ax].set_ylim([0,model_top]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                    axes[i_ax].grid(b=True,which="both",axis="both")
                elif key=="vt":
                    #make labels
                    axes[i_ax].set_xlabel(r"$\Delta$ velocity m s-1 / m")
                    axes[i_ax].set_ylabel("height / m" ) #TODO: plot also the untis of these properties

                    #change the axis
                    #axes[i_ax].set_xlim([0,1.0]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                    axes[i_ax].set_ylim([0,model_top]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                    axes[i_ax].grid(b=True,which="both",axis="both")
            
                i_ax+=1

#add labels and legend
for i_ax,ax_now in enumerate(axes):
    #generate some dummy labels indicating mean and quantiles
    if i_ax<2: #absolute values
        ax_now.plot(np.nan,np.nan,color="k",linestyle="-",linewidth=2.0,label="median")
        ax_now.plot(np.nan,np.nan,color="k",linestyle="--",linewidth=0.5,label="quantiles: " + str(quantiles[0]) + " and " + str(quantiles[-1]))
    else: #differences
        ax_now.plot(np.nan,np.nan,color="k",linestyle="-",linewidth=2.0,label="diff of quantiles: " + str(quantiles[0]) + " and " + str(quantiles[-1]))
        
    if separated_by_fallspeedsens:
        ax_now.plot(np.nan,np.nan,color="r",linestyle="-",linewidth=2.0,label="Boehm")
        ax_now.plot(np.nan,np.nan,color="b",linestyle="-",linewidth=2.0,label="Atlas")
        ax_now.plot(np.nan,np.nan,color="g",linestyle="-",linewidth=2.0,label="power law")
        ax_now.plot(np.nan,np.nan,color="y",linestyle="-",linewidth=2.0,label="power law SB")

    elif separated_by_geom:
        ax_now.plot(np.nan,np.nan,color="b",linestyle="-",linewidth=2.0,label="monodep")
        ax_now.plot(np.nan,np.nan,color="r",linestyle="-",linewidth=2.0,label="binary")
        ax_now.plot(np.nan,np.nan,color="g",linestyle="-",linewidth=2.0,label="constant")
    elif separated_by_habit:
        for i_habit, habit in enumerate(MCsnowhabit_list):
            if habit=="mixcolumndend":
                labelhabit="mix2\n(column\n+dendrite)"
            else:
                labelhabit=habit
            ax_now.plot(np.nan,np.nan,color=np.array(['r','b','g'])[i_sensMC+i_sensMCfallspeed+i_habit],label=labelhabit)
    
    #show legend
    ax_now.legend() 
    ax_now.set_ylabel("height / m" ) #TODO: plot also the untis of these properties  
############
#save figure
############
plt.tight_layout()
output_string_folder = '/home/mkarrer/Dokumente/plots/MC_percentiles/' + experiment
if not os.path.exists(output_string_folder): #create direktory if it does not exists
    os.makedirs(output_string_folder)
out_filestring = "Dv_perc_" + McSnow_geom_specifier_onestring + MCtermvel_specifier_onestring + habits_onestring + '_' + testcase
plt.savefig(output_string_folder + out_filestring + '.pdf', dpi=400)
plt.savefig(output_string_folder + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + output_string_folder + out_filestring + '.pdf'
subprocess.Popen(['evince',output_string_folder + out_filestring + '.pdf'])

if separated_by_fallspeedsens:
    #copy to plots4paper folder
    subprocess.call('cp ' + output_string_folder + out_filestring + '.pdf /home/mkarrer/Dokumente/plots/4paper/McSnow_run_percentiles_bulkparam.pdf',shell=True)
elif separated_by_geom:
    #copy to plots4paper folder
    subprocess.call('cp ' + output_string_folder + out_filestring + '.pdf /home/mkarrer/Dokumente/plots/4paper/McSnow_run_percentiles_monodepbinaryconst.pdf',shell=True)
elif separated_by_habit:
    #copy to plots4paper folder
    subprocess.call('cp ' + output_string_folder + out_filestring + '.pdf /home/mkarrer/Dokumente/plots/4paper/McSnow_run_percentiles_habits.pdf',shell=True)
