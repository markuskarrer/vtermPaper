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
#self-written functions
from functions import __plotting_functions
from functions import __postprocess_McSnow
from functions import __postprocess_PAMTRA
from matplotlib import rc


def plot_histograms(ax,SP,var="diam",color="r",linestyle="-",run_name=""):
    '''
    plot histogram at selected heights
    INPUT:
        ax: axes handle
        SP: super-particle list
    OPTIONAL:
        var: variable from SP dictionary which should be plotted
        linestyle: linestyle for the plot (to distinguish different simulations)
        run_name: name of the run/simulation
    '''
    #select heights
    #h1_bounds = [[4900,5000],[2400,2600],[0,100]]
    h1_bounds = [[2950,3050],[1950,2050],[0,100]]

    colors=["r","g","b"]
    ylabel="normalized counts" 
    if var=="diam":
        bins=np.linspace(0,0.05,100) 
        xlabel=r"$D_{max}$ [m]"
    elif var=="mm":
        bins=np.linspace(1,1000,100)
        xlabel=r"$N_{mono}$ [m]"
    elif var=="vt":
        bins=np.linspace(0.0,2.5,16)
        xlabel=r"$v_{term}$ [m s-1]"
    else:
        raise ValueError(var + " not implemented") 
    
    for i_height,height_bound in enumerate(h1_bounds):
        #select all SP in the corresponding height boundaries
        in_bin = ((SP['height']>height_bound[0]) & (SP['height']<height_bound[1]))
        prop_in_bin = SP[var][in_bin] 
        
           
        #compute and plot the histogram
        hist,bin_edges = np.histogram(prop_in_bin, bins=bins, density=True)
        bin_center=bin_edges[:-1]+np.diff(bin_edges)
        b=0
        ax.plot(bin_center,hist,c=colors[i_height],linestyle=linestyle,
            label=run_name + "height=[" + str(height_bound[0])+ "-" + str(height_bound[1]) + "] max=" + str(np.max(prop_in_bin)))
        a=0
        ax.set_xlabel(xlabel) 
        ax.set_ylabel(ylabel) 
        ax.legend()
    return ax

def get_coll_eff(SP1,SP2):
    coll_eff=1 #TODO:
    return coll_eff

def plot_agg_rates(ax,SP,which_coll="all"):
    '''
    calculate and plot aggregation rates resolved over sizes considering eq. (21) in Brdar&Seifert (2018)
    INPUT:
        ax: axes handle
        SP: list of superparticles
    OPTIONAL:
        which_coll=[all,self_ice,ice_snow,self_snow] which collections should be displayed? all->all; self_ice-> just Nmono=1 with Nmono=1;  ...
    '''

    prefactor=np.pi/4. #see eq. (21) in Brdar&Seifert (2018)

    #define diameter array (used for D_i and D_j)
    D_array=np.linspace(0,0.05,100) 


    stick_eff=1.0 #ATTENTION: sticking efficiency is not considered yet

    SP_cat=dict() #dictionary reduced to certain particles (categorized by monomers and aggregates; ice and snow)
    if which_coll=="all":
        SP_cat = SP #just make a copy of SP
    elif which_coll=="self_ice":
        
        for key in SP.keys():
            SP_cat[key] = SP[key][SP["mm"]==1] 



    for i_SP,diam in enumerate(SP_cat["diam"]):
        SP_smaller = dict() #dictionary reduced to particles smaller than SP[i_SP]
        #select all particles smaller than i_SP
        for key in SP.keys():
            SP_smaller[key] =  SP_cat[key][SP_cat["diam"]<SP_cat["diam"][i_SP]]

        for j_SP,diam in enumerate(SP_smaller["diam"]):
            pass
            #coll_eff = get_coll_eff(SP],SP[j_SP]) #this is very hard
        

        debug()

    
    return ax



if __name__ == "__main__":
    #prepare plotting          
    number_of_rows = 3
    #optimize the appearance of the plot (figure size, fonts)
    [fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_rows,aspect_ratio=1./7.,legend_fontsize='medium')
     
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


    #get timestep
    i_timestep=(tstep/60)-1 #there is no output for t=0min after that there are output steps in 10 minute steps (this could vary)
    i_timestep_end=(tstep_end/60)-1

    #plot different lines for different McSnow geometries
    McSnow_geom_list = McSnow_geom_specifier_onestring.split('_')
    MCtermvel_list = MCtermvel_specifier_onestring.split('_')
    MCsnowhabit_list = habits_onestring.split('_')


    linestyles=["-","--","-."]

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
     

       
                #read SP list (mass2fr.dat) 
                filestring_SP = 'mass2fr_' + str(tstep).zfill(4) + '-' + str(tstep_end).zfill(4) + 'min_avtstep_' + str(av_tstep) + '.ncdf'
                SP = __postprocess_McSnow.read_mass2frncdf(MC_dir + 'experiments/' + experiment,filestring_SP)
                
                ####plot some histograms
                #monomer number
                print i_habit,i_sensMC,i_sensMCfallspeed
                i_run = i_habit+i_sensMC+i_sensMCfallspeed
                if len(McSnow_geom_list)>1:#add labels for the different sensruns
                    labellist=["CTRL: ","binary: ","constant: "]
                elif len(MCtermvel_list)>1:
                    labellist=["Atlas: ","power law: "]
                else:
                    labellist=["","",""]
                plot_histograms(axes[0],SP,var="mm",linestyle=linestyles[i_run],run_name=labellist[i_run])
                #D_max
                plot_histograms(axes[1],SP,var="diam",linestyle=linestyles[i_run],run_name=labellist[i_run])

                #vterm
                plot_histograms(axes[2],SP,var="vt",linestyle=linestyles[i_run],run_name=labellist[i_run])

       

    ###plot agg_rates 
    #agg_matrix = plot_agg_rates(axes[3],SP,which_coll="self_ice")

    #formatting
    plt.tight_layout()
    #save and display plot
    out_filestring = "/add_diagnostics_" + switch_off_processes_str + '_' + McSnow_geom_specifier_onestring + MCtermvel_specifier_onestring + habits_onestring + '_' + testcase +  "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
    plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf', dpi=400)
    plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png', dpi=400)
    print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'
    subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'])
