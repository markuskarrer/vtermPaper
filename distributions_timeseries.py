'''
create panels which depict the spectra of different quantities of McSnow (including 1d-SB) and PAMTRA output
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
import time

#self-written functions
from functions import __plotting_functions
from functions import __postprocess_McSnow
from functions import __postprocess_PAMTRA
from functions import __postprocess_SB
#some hardcoded switches
plot_spectra = False #plot also the radar spectra?

#read variables passed by shell script
tstep1 = int(os.environ["tstep1"]) #start of the averaging period
tstep2 = int(os.environ["tstep2"]) #string with 4 numbers and 'min'
tstep_int = int(os.environ["tstep_int"]) #interval of timesteps
tstep_aver = int(os.environ["tstep_aver"]) #average timesteps for one distribution

experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"]
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output
MC_dir = os.environ["MC"]
adapt_version = int(os.environ["adapt_version"]) #reading the files of the appropriate adaption version
skipMC = (os.environ["skipMC"]=="True") #allows to run the scripts also if no McSnow data is there (only 1D-SB runs)
skipSB = (os.environ["skipSB"]=="True") #ATTENTION: check the matching of the heights between SB and MC again
#directory of experiments
directory = MC_dir + "/experiments/"

#number of heights in plot
number_of_heights=5

#small self-defined function
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx]


#####
##organising the whole figure (including subplots)
#####
#get number of plot automatically
raw_plots=  (tstep2-tstep1)/tstep_int+1
#choose variables to plot as histograms #the variables are not relevant if skipMC is True but the length of the arrays isused later
#plot_vars = ["d_counts_normalized","mass_in_bin_norm"]
#full_name = ["number density","mass density"]
#var_units = ["m-3","kg"]
#organize the limits of the number and mass density plots
nmin=1e2; lmin=1e-3 #define minimum of plot in y-dimension (n->number concentration; l-> mass mixing ratio)
ylims_low = [nmin,lmin]
ylims_high = [1e10,1e1]
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
number_of_plots = raw_plots #+len(plot_vars)+pam_plots+1 #+1 is for the v-D scatterplot
figsize_height = 6.0/2.0*(number_of_plots)
fig, axes	=	plt.subplots(nrows=number_of_plots, ncols=1, figsize=(8.0,figsize_height))

i_ax=-1
for tstep in range(tstep1,tstep2+1,tstep_int):
    i_ax+=1
    tstep_end=tstep+tstep_aver

    if not skipMC:
        ###
        #load SP-file (keys are: "m_tot","Frim","height","d_rime","vt","xi",    "rhor","a_rime","mr_crit","diam",    "proj_A",   "mm",         "m_rime",   "m_wat")
        ###
        SP_file = Dataset(directory + experiment + '/mass2fr_' + str(tstep).zfill(4) + '-' + str(tstep_end).zfill(4) + 'min_avtstep_' + str(av_tstep) + '.ncdf',mode='r')
        #create dictionary for all variables from PAMTRA
        SP = dict()
        #if necessary change name of variables (not done yet)
        varlist = SP_file.variables
        #read Sp variables to SP dictionary
        for var in varlist:#read files and write it with different names in Data
            SP[var] = np.squeeze(SP_file.variables[var])

        #get maximum height for plotting
        model_top = np.nanmax(SP['height'])
        if plot_spectra:
            #read pamtra output to pamData dictionary (this is done here already in order to plot the size distributions at the same heights)
            #define file to read
            if adapt_version==1:
                filename = "/adaptv1_t" + str(tstep) #+ ".nc"
            elif adapt_version==2:
                filename = "/adaptv2_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min'
            elif adapt_version==3:
                filename = "/adaptv3_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min'
            pam_file = Dataset(directory + experiment + filename + '.nc')
            nheights = pam_file["heightbins"].shape[0]+1 #this is used in the next lines for determining the binning of the heights


        else:#we need to define the number of heights if we dont have PAMTRA output
            nheights=100
        #####
        ##calculate general properties of the run
        #####
        #calculate volume of box
        zres = model_top/(nheights-1) #vertical resolution
        Vbox = __postprocess_McSnow.calculate_Vbox(experiment,zres)    
    #get relevant timestep
    i_timestep=(tstep_end/10)-1 #there is no output for t=0min after that there are output steps in 30 minute steps (this could vary)

    ###
    #plot N(m): size distribution as a function of mass
    ###
    #define array for mass
    nmassbins = 100
    m_bound_ds = np.logspace(-13,-5,nmassbins+1) #array from 1nm to 1m
    m_ds = m_bound_ds[:-1] + 0.5*np.diff(m_bound_ds) #diameter at center of bins
    #get colorbar from  matplotlib internal color order to cycle through them by your own


    if not skipMC:
        #McSnow
        #for i_height in range(0,binned_val['d_counts'].shape[0],10): #
        i_lines = 0 #counter: used to get the right color of the lines
        heightvec = np.linspace(0,model_top,nheights) #has to be calculated before separate_by_height_and_diam to fasten calculation with calconly but is overwritten by seperate_by_height_and_diam
        heightstep = nheights/number_of_heights #plot 5 different heights #del_z is model_top/n_heights
        for i_height in range(heightstep-1,heightvec.shape[0],heightstep): #
            if not plot_spectra:
                i_height-=1 #reduce i_height because of index problem ATTENTION: no final solution
            in_height_bin = np.logical_and(heightvec[i_height]<=SP["height"],SP["height"]<heightvec[i_height+1])
            [N_m,__] = np.histogram(SP["m_tot"][in_height_bin],bins=m_bound_ds,weights=SP["xi"][in_height_bin],density=False)
            binwidths = m_bound_ds[1::] - m_bound_ds[:-1]
            N_m_normed = N_m/binwidths/Vbox #get from [#]/[kg]/[m3]
            #from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug

            #plot a line to the corresponding height
            axes[i_ax] == __plotting_functions.plot1Dhistline(axes[i_ax],m_ds,N_m_normed,xlabel='mass / kg',ylabel=('number density / m-3 kg-1'),linelabel="[{:5.0f}m,{:5.0f}m]".format(heightvec[i_height],heightvec[i_height+1]),logflag=3,linestyle='--',color=colors[i_lines])
            i_lines += 1
    #SB
    #load netCDF4 file from McSnows-embedded-2mom-scheme
    twomom_file = Dataset(directory + experiment + '/twomom_d.ncdf',mode='r')
    #create dictionary for all variables from twomom_d.ncdf
    twomom = dict()
    #if necessary change name of variables
    varlist = twomom_file.variables
    #read twomom variables to twomom dictionary
    for var in varlist:#read files and write it with different names in Data
        twomom[var] = np.squeeze(twomom_file.variables[var])

    heightstep = twomom["heights"].shape[0]/number_of_heights #100m for dz=20m
    #loop over certain heights (calculation and plotting)
    #initialize arrays for distribution
    twomom["i_mdist"] = np.zeros((twomom["heights"].shape[0],m_ds.shape[0]))
    twomom["s_mdist"] = np.zeros((twomom["heights"].shape[0],m_ds.shape[0]))
    twomom["all_mdist"] = np.zeros((twomom["heights"].shape[0],m_ds.shape[0]))


    i_lines = 0 #counter: used to get the right color of the lines
    for i_height in range(twomom["heights"].shape[0]-heightstep,-1,-heightstep): #ATTENTION: if you do changes here uncomment linelabel for checking consistent height labelling with McSnow data...
        #calculate SB distribution
        try:
            twomom["i_mdist"][i_height,:],dummy = __postprocess_SB.calc_fmass_distribution_from_moments(twomom,'icecosmo5',m_ds,i_height=i_height,i_time=i_timestep)
            twomom["s_mdist"][i_height,:],dummy = __postprocess_SB.calc_fmass_distribution_from_moments(twomom,'snowjplatesnonsphere',m_ds,i_height=i_height,i_time=i_timestep)
            #sum up distributions #ATTENTION: if there are more categories used update this here (include them)
            twomom["all_mdist"][i_height,:] = np.nan_to_num(twomom["i_mdist"][i_height,:])+np.nan_to_num(twomom["s_mdist"][i_height,:]) #nan_to_num converts nan to 0 (this is helpful if one of the category doesnt exist)
        except:
            twomom["i_mdist"][i_height,:] = np.nan
            twomom["s_mdist"][i_height,:] = np.nan
            #sum up distributions #ATTENTION: if there are more categories used update this here (include them)
            twomom["all_mdist"][i_height,:] = np.nan            
        #plot distribution
        axes[i_ax] = __plotting_functions.plot1Dhistline(axes[i_ax],m_ds,twomom["all_mdist"][i_height,:],xlabel='mass / kg',ylabel='number density / m-3 kg-1',ylims=[10**8,10**19],logflag=3,color=colors[i_lines], linelabel="__none") #>">{:6.0f}m".format(twomom["heights"][i_height])) #linelabel="_none")
        i_lines += 1
    #add labels to distinguish the linestyle
    if not skipMC:
        axes[i_ax].plot(0,0,color='k',linestyle='--',label='McSnow')
    axes[i_ax].plot(0,0,color='k',linestyle='-',label='SB')
    axes[i_ax].legend()
    
    #add text for identifying time
    axes[i_ax].text(0.1, 0.9,'t='+ str(tstep) + 'min',
     horizontalalignment='center',
     verticalalignment='center',
     transform = axes[i_ax].transAxes)

    
    
#save figure
plt.tight_layout()
dir_distributions = '/home/mkarrer/Dokumente/plots/distributions/'
if not os.path.exists(dir_distributions + experiment): #create direktory if it does not exists
    os.makedirs(dir_distributions + experiment)
out_filestring = "/dist_timeseries_" + testcase + "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
plt.savefig(dir_distributions + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig(dir_distributions + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + dir_distributions + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince',dir_distributions + experiment + out_filestring + '.pdf'])