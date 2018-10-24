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
adapt_version = int(os.environ["adapt_version"]) #reading the files of the appropriate adaption version
skipMC = (os.environ["skipMC"]=="True") #allows to run the scripts also if no McSnow data is there (only 1D-SB runs)
#directory of experiments
directory = MC_dir + "/experiments/"

if not skipMC:
    ###
    #load SP-file (keys are: "m_tot","Frim","height","d_rime","vt","xi",    "rhor","a_rime","mr_crit","diam",    "proj_A",   "mm",         "m_rime",   "m_wat")
    ###
    SP_file = Dataset(directory + experiment + '/mass2fr_' + str(tstep).zfill(4) + 'min_avtstep_' + str(av_tstep) + '.ncdf',mode='r')
    #create dictionary for all variables from PAMTRA
    SP = dict()
    #if necessary change name of variables (not done yet)
    varlist = SP_file.variables
    #read Sp variables to SP dictionary
    for var in varlist:#read files and write it with different names in Data
        SP[var] = np.squeeze(SP_file.variables[var])

    #get maximum height for plotting
    model_top = np.nanmax(SP['height'])

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

    #####
    ##calculate general properties of the run
    #####
    #calculate volume of box
    zres = model_top/(nheights-1) #vertical resolution
    Vbox = __postprocess_McSnow.calculate_Vbox(experiment,zres)
    
#get relevant timestep
i_timestep=(tstep/60)-1 #there is no output for t=0min after that there are output steps in 30 minute steps (this could vary)

#####
##organising the whole figure (including subplots)
#####
#get number of plot automatically
raw_plots=  1
pam_plots = 1
#choose variables to plot as histograms #the variables are not relevant if skipMC is True but the length of the arrays isused later
plot_vars = ["d_counts_normalized","mass_in_bin_norm"]
full_name = ["number density","mass density"]
var_units = ["m-3","kg"]
#organize the limits of the number and mass density plots
nmin=1e2; lmin=1e-3 #define minimum of plot in y-dimension (n->number concentration; l-> mass mixing ratio)
ylims_low = [nmin,lmin]
ylims_high = [1e10,1e1]

###
#plot N(m): size distribution as a function of mass
###
#define array for mass
nmassbins = 100
m_bound_ds = np.logspace(-13,-5,nmassbins+1) #array from 1nm to 1m
m_ds = m_bound_ds[:-1] + 0.5*np.diff(m_bound_ds) #diameter at center of bins
#get colorbar from  matplotlib internal color order to cycle through them by your own
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
number_of_plots = raw_plots+len(plot_vars)+pam_plots
figsize_height = 6.0/2.0*(number_of_plots)
fig, axes	=	plt.subplots(nrows=number_of_plots, ncols=1, figsize=(8.0,figsize_height))
if not skipMC:
    #McSnow
    #for i_height in range(0,binned_val['d_counts'].shape[0],10): #
    i_lines = 0 #counter: used to get the right color of the lines
    heightvec = np.linspace(0,model_top,nheights) #has to be calculated before separate_by_height_and_diam to fasten calculation with calconly but is overwritten by seperate_by_height_and_diam
    heightstep=nheights/5 #plot 5 different heights #del_z is model_top/n_heights
    for i_height in range(heightstep-1,heightvec.shape[0],heightstep): #
        in_height_bin = np.logical_and(heightvec[i_height]<=SP["height"],SP["height"]<heightvec[i_height+1])
        [N_m,__] = np.histogram(SP["m_tot"][in_height_bin],bins=m_bound_ds,weights=SP["xi"][in_height_bin],density=False)
        binwidths = m_bound_ds[1::] - m_bound_ds[:-1]
        N_m_normed = N_m/binwidths/Vbox #get from [#]/[kg]/[m3]
        #from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug

        #plot a line to the corresponding height
        axes[0] == __plotting_functions.plot1Dhistline(axes[0],m_ds,N_m_normed,xlabel='mass / kg',ylabel=('number density / m-3 kg-1'),linelabel="[{:5.0f}m,{:5.0f}m]".format(heightvec[i_height],heightvec[i_height+1]),logflag=3,linestyle='--',color=colors[i_lines])
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
heightstep = twomom["heights"].shape[0]/5 #100m for dz=20m
#loop over certain heights (calculation and plotting)
#initialize arrays for distribution
twomom["i_mdist"] = np.zeros((twomom["heights"].shape[0],m_ds.shape[0]))
twomom["s_mdist"] = np.zeros((twomom["heights"].shape[0],m_ds.shape[0]))
twomom["all_mdist"] = np.zeros((twomom["heights"].shape[0],m_ds.shape[0]))


i_lines = 0 #counter: used to get the right color of the lines
for i_height in range(twomom["heights"].shape[0]-heightstep,-1,-heightstep): #ATTENTION: if you do changes here uncomment linelabel for checking consistent height labelling with McSnow data...
    #calculate SB distribution
    twomom["i_mdist"][i_height,:] = __postprocess_SB.calc_fmass_distribution_from_moments(twomom,'icecosmo5',m_ds,i_height=i_height,i_time=i_timestep)
    twomom["s_mdist"][i_height,:] = __postprocess_SB.calc_fmass_distribution_from_moments(twomom,'snowSBB',m_ds,i_height=i_height,i_time=i_timestep)
    #sum up distributions #ATTENTION: if there are more categories used update this here (include them)
    twomom["all_mdist"][i_height,:] = twomom["i_mdist"][i_height,:]+twomom["s_mdist"][i_height,:]
    #plot distribution
    axes[0] = __plotting_functions.plot1Dhistline(axes[0],m_ds,twomom["all_mdist"][i_height,:],xlabel='mass / kg',ylabel='number density / m-3 kg-1',ylims=[10**8,10**19],logflag=3,color=colors[i_lines], linelabel="__none") #>">{:6.0f}m".format(twomom["heights"][i_height])) #linelabel="_none")
    i_lines += 1
#add labels to distinguish the linestyle
if not skipMC:
    axes[0].plot(0,0,color='k',linestyle='--',label='McSnow')
axes[0].plot(0,0,color='k',linestyle='-',label='SB')
axes[0].legend()

###
#end of N(m) plot
###
#increase font size
params = {'legend.fontsize': 'large',
    'figure.figsize': (15, 5),
    'axes.labelsize': 'x-large', #size: Either a relative value of 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large' or an absolute font size, e.g., 12
    'axes.titlesize':'x-large',
    'xtick.labelsize':'x-large',
    'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
#calculate values binned to here defined h-D bins
nbins=100;diamrange=[-6,0]
if not skipMC:
    #calculate binned values (histograms)
    heightvec = np.linspace(0,model_top,nheights) #has to be calculated before seperate_by_height_and_diam to fasten calculation with calconly but is overwritten by seperate_by_height_and_diam
    heightstep = nheights/5
    calconly = heightvec[range(heightstep-1,heightvec.shape[0],heightstep)] #skip heights in seperate_by_height_and_diam calculation which are not plotted later
    binned_val,heightvec,d_bound_ds,d_ds,__ = __postprocess_McSnow.separate_by_height_and_diam(SP,nbins=nbins,diamrange=diamrange,nheights=nheights,model_top=model_top,calconly=calconly)


    #divide by box volume to achieve [#]->[#/m3]
    binned_val["d_counts"] = __postprocess_McSnow.conv_num2numpm3(binned_val["d_counts"],Vbox)
    binned_val["mass_in_bin"] = __postprocess_McSnow.conv_num2numpm3(binned_val["mass_in_bin"],Vbox)
    binwidths = d_bound_ds[1::] - d_bound_ds[:-1]
    binned_val["d_counts_normalized"] = binned_val["d_counts"] / binwidths
    binned_val["mass_in_bin_norm"] = binned_val["mass_in_bin"] / binwidths
    binned_val["d_counts_no_mult"] = __postprocess_McSnow.conv_num2numpm3(binned_val["d_counts"],Vbox)

    #do several plots as defined above in plot_vars ...
    i_plots = 0
    #initialize variable for upper x-limit of plot (diameter in m)
    dmaxn = 1e-12; dmax_newm=1e-12; dmaxm=1e-12; dmax_newn=1e-12 
    dmax_same_bool=True #if true maximum of dmaxn and dmaxn is choosen
    #plot number density and mass density for McSnow
    for ax, var, name,units,ylim_low_now,ylim_high_now,color in zip(axes[raw_plots::].flat, plot_vars, full_name,var_units,ylims_low,ylims_high,colors):
        #for i_plots,varname in enumerate(plot_vars): 
        print '####################'
        print 'plot: ' + var
        print '####################'
        print ax, var, name,units,ylim_low_now,ylim_high_now,color

        
        #for i_height in range(0,binned_val['d_counts'].shape[0],10): #
        i_lines = 0 #counter: used to get the right color of the lines
        heightstep=nheights/5 #plot 5 different heights #del_z is model_top/n_heights
        for i_height in range(heightstep-1,heightvec.shape[0],heightstep): #
            #binned_val[var] = __postprocess_McSnow.kernel_estimate(binned_val[var][i_height,:],number_of_SP) #TODO:? a kernel kernel_estimate for the distribution here

            #plot a line to the corresponding height
            #dmax_new = d_ds[max(np.argwhere(binned_val["d_counts_normalized"][i_height,:]> nmin)[-1],np.argwhere(binned_val["mass_in_bin_norm"][i_height,:]> lmin)[-1])]*1.1
            #try:
            if var=="d_counts_normalized":
                if np.argwhere(binned_val["d_counts_normalized"][i_height,:]>nmin).shape[0]>=1:
                    dmax_newn = d_ds[np.argwhere(binned_val["d_counts_normalized"][i_height,:]>nmin)[-1]]
            elif var=="mass_in_bin_norm":
                if np.argwhere(binned_val["mass_in_bin_norm"][i_height,:]>lmin).shape[0]>=1:
                    dmax_newm = d_ds[np.argwhere(binned_val["mass_in_bin_norm"][i_height,:]>lmin)[-1]]

            #find the minimum diameter in which the lines fit in
            if dmax_newn>dmaxn: #update dmax if a bigger one in this height is found
                dmaxn=dmax_newn
            if dmax_newm>dmaxm: #update dmax if a bigger one in this height is found
                dmaxm=dmax_newm
            if dmax_same_bool:
                dmaxn=max(dmaxn,dmaxm)
                dmaxm=max(dmaxn,dmaxm)  
            #plot N(D) (for i_plots=0) and M(D) (for i_plots=1)
            axes[i_plots] == __plotting_functions.plot1Dhistline(ax,d_ds,binned_val[var][i_height,:],xlabel='diameter / m',ylabel=(name + ' / ' + units),linelabel="[{:5.0f}m,{:5.0f}m]".format(heightvec[i_height],heightvec[i_height+1]),logflag=2,xlims=[0,dmaxn],linestyle='--',color=colors[i_lines])
            i_lines += 1
        i_plots+=1

else: #we have to create diameter arrays also if there is no McSnow plot
    #define arrays with sp_diam
    d_bound_ds = np.logspace(diamrange[0],diamrange[1],nbins+1) #array from 1nm to 1m
    d_ds = d_bound_ds[:-1] + 0.5*np.diff(d_bound_ds) #diameter at center of bins
    #initialize variable for upper x-limit of plot (diameter in m)
    dmaxn = 1e-12; dmax_newm=1e-12; dmaxm=1e-12; dmax_newn=1e-12 
    dmax_same_bool=True #if true maximum of dmaxn and dmaxn is choosen
    
print "#################"
print "#calculate and plot distribution of SBs ice (and snow) in same panel as McSnows distribution (or alone if skipMC==True)"
print "#################"

#add SB lines to first plot
#initialize arrays for distribution
#N(D)
twomom["i_dist"] = np.zeros((twomom["heights"].shape[0],d_ds.shape[0]))
twomom["s_dist"] = np.zeros((twomom["heights"].shape[0],d_ds.shape[0]))
twomom["all_dist"] = np.zeros((twomom["heights"].shape[0],d_ds.shape[0]))
#M(D)
twomom["i_dist_mass"] = np.zeros((twomom["heights"].shape[0],d_ds.shape[0]))
twomom["s_dist_mass"] = np.zeros((twomom["heights"].shape[0],d_ds.shape[0]))
twomom["all_dist_mass"] = np.zeros((twomom["heights"].shape[0],d_ds.shape[0]))

#calculate sum of distributions of SB cloud ice + snow (based on number & mass density) and plot them
i_lines = 0 #counter: used to get the right color of the lines
heightstep = twomom["heights"].shape[0]/5 #100m for dz=20m
for i_height in range(twomom["heights"].shape[0]-heightstep,-1,-heightstep): #ATTENTION: if you do changes here uncomment linelabel for checking consistent height labelling with McSnow data...
    twomom["i_dist"][i_height,:],twomom["i_dist_mass"][i_height,:] = __postprocess_SB.calc_distribution_from_moments(twomom,'icecosmo5',d_ds,i_height=i_height,i_time=i_timestep)
    twomom["s_dist"][i_height,:],twomom["s_dist_mass"][i_height,:] = __postprocess_SB.calc_distribution_from_moments(twomom,'snowSBB',d_ds,i_height=i_height,i_time=i_timestep)
    twomom["all_dist"][i_height,:] = twomom["i_dist"][i_height,:]+twomom["s_dist"][i_height,:]
    twomom["all_dist_mass"][i_height,:] = twomom["i_dist_mass"][i_height,:]+twomom["s_dist_mass"][i_height,:]

    if np.argwhere(twomom["all_dist"][i_height,:]>nmin).shape[0]>=1:
        dmax_newn = d_ds[np.argwhere(twomom["all_dist"][i_height,:]> nmin)[-1]]
    if np.argwhere(twomom["all_dist_mass"][i_height,:]>lmin).shape[0]>=1:
        dmax_newm = d_ds[np.argwhere(twomom["all_dist_mass"][i_height,:]> lmin)[-1]]

    if dmax_newn>dmaxn: #update dmax if a bigger one in this height is found
        dmaxn=dmax_newn
    if dmax_newm>dmaxm: #update dmax if a bigger one in this height is found
        dmaxm=dmax_newm
    if dmax_same_bool:
        dmaxn=max(dmaxn,dmaxm)
        dmaxm=max(dmaxn,dmaxm)    

    #plot a line corresponding to a specific height; N(D) on axes[1] and M(D) on axes[2]
    axes[1] = __plotting_functions.plot1Dhistline(axes[1],d_ds,twomom["all_dist"][i_height,:],xlabel='diameter / m',ylabel='number density / m-4',ylims=[ylims_low[0],ylims_high[0]],logflag=2,xlims=[0,dmaxn],color=colors[i_lines], linelabel="__none") #>">{:6.0f}m".format(twomom["heights"][i_height])) #linelabel="_none")
    axes[2] = __plotting_functions.plot1Dhistline(axes[2],d_ds,twomom["all_dist_mass"][i_height,:],xlabel='diameter / m',ylabel='mass density / kg m-4',ylims=[ylims_low[1],ylims_high[1]],logflag=2,xlims=[0,dmaxm],color=colors[i_lines], linelabel="__none") #>{:6.0f}m".format(twomom["heights"][i_height]))
    i_lines += 1
#add labels to distinguish the linestyle
if not skipMC:
    axes[1].plot(0,0,color='k',linestyle='--',label='McSnow')
    axes[2].plot(0,0,color='k',linestyle='--',label='McSnow')
axes[1].plot(0,0,color='k',linestyle='-',label='SB')
axes[2].plot(0,0,color='k',linestyle='-',label='SB')
axes[1].legend(ncol=3)  
axes[2].legend(ncol=3)  

###
#plot spectra
###
#get axis of subplot
#ax = plt.subplot2grid((number_of_plots, 1), (i_plots+1, 0))
try:
    print "##################################"
    print "plot PAMTRA output below (McSnow)"
    print "#################################"
    
    #define string of pamtraData
    pam_filestring = directory + experiment + filename + '.nc'

    #read pamtra output to pamData dictionary
    pamData = __postprocess_PAMTRA.read_pamtra(pam_filestring)
    #plot spectrogram
    freqindex = 1 #TODO: query frequency
    i_lines = 0 #TODO loop over this
    heightstep = 10 #1000m if dz=100
    for i_height in range(heightstep-1,pamData["height"].shape[0],heightstep):
        #print pamData["Radar_Velocity"][freqindex,:],pamData["Radar_Spectrum"][i_height,freqindex,:]
        axes[len(plot_vars)+raw_plots] = __plotting_functions.plot1Dhistline(axes[len(plot_vars)+raw_plots],-pamData["Radar_Velocity"][freqindex,:],pamData["Radar_Spectrum"][i_height,freqindex,:],xlabel="Doppler velocity / m s-1",ylabel="spectral power / dB",xlims=[-4,0],ylims=[-50,30],logflag=0,color=colors[i_lines],linelabel="[{:5.0f}m,{:5.0f}m]".format(pamData["height"][i_height]-0.5*(pamData["height"][1]-pamData["height"][0]),pamData["height"][i_height]+0.5*(pamData["height"][1]-pamData["height"][0])),linestyle='--') #the labelling assumes that the pamtra-heightbins are equally spaced
        i_lines+=1

    
except Exception:
	print ' \n \n in except: \n \n'
	s = traceback.format_exc()
   	serr = "there were errors:\n%s\n" % (s)
    	sys.stderr.write(serr) 

#plot SB pamtra output data
try:
    print "#############################"
    print "plot PAMTRA output below (SB)"
    print "#############################"

    #define file to read
    pam_filestring = directory + experiment + "/PAMTRA_2mom_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min.nc'
    #ax = plt.subplot2grid((len(plot_vars)+num_pam_SB_plots, 1), (i+4, 0))

    #read pamtra output to pamData dictionary
    pamData = __postprocess_PAMTRA.read_pamtra(pam_filestring)
    #plot spectrogram
    #for i_height in range(twomom["heights"].shape[0]-1,0,-50):
    freqindex = 1 #TODO: query frequency
    i_lines = 0 #TODO loop over this
    heightstep = 10 #1000m if dz=100
    for i_height in range(heightstep-1,pamData["height"].shape[0],heightstep):
        xlims=[-4,0];ylims=[-60,30]
        axes[len(plot_vars)+raw_plots] = __plotting_functions.plot1Dhistline(axes[len(plot_vars)+raw_plots],-pamData["Radar_Velocity"][freqindex,:],pamData["Radar_Spectrum"][i_height,freqindex,:],xlabel="Doppler velocity / m s-1",ylabel="spectral power / dB",xlims=xlims,ylims=ylims,logflag=0,color=colors[i_lines],linelabel="__none",linestyle='-')
        i_lines+=1
    #plot label of used frequency
    plt.text(xlims[0]+(xlims[1]-xlims[0])*0.001, ylims[1]-(ylims[1]-ylims[0])*0.07, str(pamData["frequency"][freqindex]) + 'GHz', fontsize=12)
    #add labels to distinguish the linestyle
    if not skipMC:
        axes[len(plot_vars)+raw_plots].plot(0,0,color='k',linestyle='--',label='McSnow')
    axes[len(plot_vars)+raw_plots].plot(0,0,color='k',linestyle='-',label='SB')
    #plot legend
    axes[len(plot_vars)+raw_plots].legend()
except Exception:
    print ' \n \n in except: \n \n'
    s = traceback.format_exc()
    serr = "there were errors:\n%s\n" % (s)
    sys.stderr.write(serr) 

#save figure
plt.tight_layout()
dir_distributions = '/home/mkarrer/Dokumente/plots/distributions/'
if not os.path.exists(dir_distributions + experiment): #create direktory if it does not exists
    os.makedirs(dir_distributions + experiment)
out_filestring = "/" + testcase + "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
plt.savefig(dir_distributions + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig(dir_distributions + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + dir_distributions + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince',dir_distributions + experiment + out_filestring + '.pdf'])