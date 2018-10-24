'''
create an overview panel with properties of McSnow and PAMTRA output
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
tstep = int(os.environ["tstep"])
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"]
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output
MC_dir = os.environ["MC"]
adapt_version = int(os.environ["adapt_version"]) #reading the files of the appropriate adaption version
skipMC = (os.environ["skipMC"]=="True") #allows to run the scripts also if no McSnow data is there (only 1D-SB runs) #ATTENTION: not completely implemented yet

#directory of experiments
directory = MC_dir + "/experiments/"

if not skipMC:
    #load file
    SP_file = Dataset(directory + experiment + '/mass2fr_' + str(tstep).zfill(4) + 'min_avtstep_' + str(av_tstep) + '.ncdf',mode='r')
    #create dictionary for all variables from PAMTRA
    SP = dict()
    #print pam_file.variables.keys() 

    #if necessary change name of variables
    varlist = SP_file.variables
    #read PAMTRA variables to pamData dictionary
    for var in varlist:#read files and write it with different names in Data
        SP[var] = np.squeeze(SP_file.variables[var])
    #get maximum height for plotting
    model_top = np.nanmax(SP['height'])
    #calculate values binned to here defined h-D bins
    binned_val,heightvec,d_bound_ds,d_ds,zres = __postprocess_McSnow.separate_by_height_and_diam(SP,nbins=100,diamrange=[-9,0],nheights=51,model_top=model_top)

    #calculate volume of box
    Vbox = __postprocess_McSnow.calculate_Vbox(experiment,zres)
    #divide by box volume to achieve [#]->[#/m3]
    binned_val["d_counts"] = __postprocess_McSnow.conv_num2numpm3(binned_val["d_counts"],Vbox)
    binned_val["d_counts_no_mult"] = __postprocess_McSnow.conv_num2numpm3(binned_val["d_counts"],Vbox)

    ############################
    #now: plot the binned values
    ############################
    cmap="viridis_r" #use one uniform colorbar

    #list of variables from binned_val-dictionary which should be plotted
    plot_vars = ["d_counts","av_mm"] #["d_counts","av_Frim","av_rhor","av_mm"] #"RPpSP",
    full_name = ["number density","na numb. monomer"] #["number density","rime fraction","ma rime density","na numb. monomer"] #for colorbar labelling
    var_units = ["m-3","1","kg m-3","1"] #["m-3","1","kg m-3","1"]#for colorbar labelling
    #select step of colorbar ticks for each variable
    colticksstep = [50000.,1] #[50000.,0.2,100.,1] #,20000.
    maskedcon = ['0','nan'] #['0','nan','nan','nan']#,'nan' #see if maskedcon = 'nan': within for-loop
    mincol_arr = [10,1] #[10,0.,0,1] #0, lowest number in colorbar
    maxcol_arr = [1e5,1000] #[1e5,1.0,900,1000]#,-999 heighest number in colorbar
    logcol_arr = [1,1] #[1,0,0,1] #set to one if colorbar should be in logarithmic scale
    number_of_plots = len(plot_vars)
else:
    number_of_plots = 0
    plot_vars=[]
num_pam_SB_plots = 7 #number of subplots from pamtra variables and twomoment output

figsize_height = 6.0/2.0*(number_of_plots+num_pam_SB_plots)
fig	=	plt.figure(figsize=(8.0,figsize_height))#figsize=(4, 4))

if not skipMC:
    for i,varname in enumerate(plot_vars): #let enumerate start at 0
        print '####################'
        print 'plot: ' + varname
        print '####################'


        params = {'legend.fontsize': 'x-large',
            'figure.figsize': (15, 5),
            'axes.labelsize': 'x-large', #size: Either an relative value of 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large' or an absolute font size, e.g., 12
            'axes.titlesize':'x-large',
            'xtick.labelsize':'x-large',
            'ytick.labelsize':'x-large'}
        pylab.rcParams.update(params)
        
        ax 	= 	plt.subplot2grid((len(plot_vars)+num_pam_SB_plots, 1), (i+1, 0)) #plotvars +num_pam_SB_plots -> let space for 2 PAMTRA_output plots
        #give range of colorbar and derive ticks
        mincol = mincol_arr[i];
        
        #derive step of colorbar ticks automatically for colticksstep=-999
        if colticksstep[i]==-999:
            colticksstep[i]=np.floor(np.nanmax(binned_val["av_mm"])/10)

        #calculate maxcol from array or use predefined
        if maxcol_arr[i]==-999:
            maxcol=math.ceil(np.nanmax(binned_val[varname]) / colticksstep[i]) * colticksstep[i] #maxcol is rounded to the next tickstep to have round values in the colorbar
        else:
            maxcol=maxcol_arr[i]
        
            
        #calculate colorbar ticks either linear or logarithmic
        if logcol_arr[i]==0:
            colticks = np.arange(mincol, maxcol*1.0001, colticksstep[i])
        elif logcol_arr[i]==1:
            colticks = np.logspace(np.log10(mincol),np.log10(maxcol),np.log10(maxcol)-np.log10(mincol)+1)
        if maskedcon[i] == 'nan':
            binned_val[varname] = np.ma.array(binned_val[varname],mask=np.isnan(binned_val[varname]))
        elif maskedcon[i] == '0':
            binned_val[varname] = np.ma.array(binned_val[varname],mask=(binned_val[varname]==0))
            masked = np.ma.array(binned_val[varname],mask=(binned_val[varname]==0))


        #call pcolor plot routine for height-diameter plots
        plt =  __plotting_functions.pcol_height_diam(ax,d_bound_ds,heightvec,binned_val[varname],
                                                    mincol=mincol,maxcol=maxcol,ticks=colticks,
                                                    logcol=logcol_arr[i],cmap=cmap,
                                                    collabel=full_name[i] + ' / ' + var_units[i])

'''
create plot of athmospheric variables first and add it before the histogram plots
'''
ax 	= 	plt.subplot2grid((number_of_plots+num_pam_SB_plots, 1), (0, 0))
ax2 = ax.twiny()
#read atmospheric variables
filestring_atmo = directory + experiment + "/atmo.dat"
atmo = __postprocess_McSnow.read_atmo(experiment,filestring_atmo)
#calculate rhi (relative humidity over ice)
atmo["rhi"] = __general_utilities.calc_rhi(atmo)
#from IPython.core.debugger import Tracer ; Tracer()()

#plot athmospheric variables
plt = __plotting_functions.plot_atmo(ax,ax2,atmo)
#increase i, because we added a plot before the for-loop
if skipMC:
    i=0
else:
    i+=1
    
#plot McSnow pamtra output
try:
    ##############################
    #now: plot PAMTRA output below
    ##############################

    #define file to read
    print "adaptv:",adapt_version
    if adapt_version==1:
        filename = "/adaptv1_t" + str(tstep) #+ ".nc"
    elif adapt_version==2:
        filename = "/adaptv2_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min'
    elif adapt_version==3:
        filename = "/adaptv3_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min'
        
    pam_filestring = directory + experiment + filename +".nc"

    ax = plt.subplot2grid((len(plot_vars)+num_pam_SB_plots, 1), (i+1, 0))

    #read pamtra output to pamData dictionary
    pamData = __postprocess_PAMTRA.read_pamtra(pam_filestring)

    #plot spectrogram
    plt = __plotting_functions.plot_McSnows_vt_in_spectrogram(ax,pamData,SP,experiment)
    ax = plt.subplot2grid((len(plot_vars)+num_pam_SB_plots, 1), (i+2, 0))
    plt = __plotting_functions.plot_pamtra_spectrogram(ax,pamData,freq=35.5)

    #plot reflectivities
    axrefl = plt.subplot2grid((len(plot_vars)+num_pam_SB_plots, 1), (i+3, 0))

    plt = __plotting_functions.plot_pamtra_Ze(axrefl,pamData,linestyle='--')
    
    #plot other radar-moments
    axmom = plt.subplot2grid((len(plot_vars)+num_pam_SB_plots, 1), (i+4, 0))

    plt = __plotting_functions.plot_pamtra_highermoments(axmom,pamData,linestyle='--')
except Exception:
	print ' \n \n in except: \n \n'
	s = traceback.format_exc()
   	serr = "there were errors:\n%s\n" % (s)
    	sys.stderr.write(serr) 
  
#plot SB pamtra output data
try:
    ##############################
    #now: plot PAMTRA output below
    ##############################

    #define file to read
    pam_filestring = directory + experiment + "/PAMTRA_2mom_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min.nc'
    #ax = plt.subplot2grid((len(plot_vars)+num_pam_SB_plots, 1), (i+4, 0))

    #read pamtra output to pamData dictionary
    pamData = __postprocess_PAMTRA.read_pamtra(pam_filestring)

    #plot reflectivities
    #ax = plt.subplot2grid((len(plot_vars)+num_pam_SB_plots, 1), (i+3, 0))
    plt = __plotting_functions.plot_pamtra_Ze(axrefl,pamData)
    axrefl.plot(0,0,color='k',linestyle='--',label='McSnow')
    axrefl.plot(0,0,color='k',linestyle='-',label='SB')
    axrefl.legend()
    
    #plot other radar-moments
    #axmom = plt.subplot2grid((len(plot_vars)+num_pam_SB_plots, 1), (i+4, 0))

    plt = __plotting_functions.plot_pamtra_highermoments(axmom,pamData)
    axmom.plot(0,0,color='k',linestyle='--',label='McSnow')
    axmom.plot(0,0,color='k',linestyle='-',label='SB')
    axmom.legend()
    
    ax = plt.subplot2grid((len(plot_vars)+num_pam_SB_plots, 1), (i+5, 0))
    #plot spectrogram
    plt = __plotting_functions.plot_pamtra_spectrogram(ax,pamData,freq=35.5)
    
except Exception:
    print ' \n \n in except: \n \n'
    s = traceback.format_exc()
    serr = "there were errors:\n%s\n" % (s)
    sys.stderr.write(serr) 

#load netCDF4 file
twomom_file = Dataset(directory + experiment + '/twomom_d.ncdf',mode='r')
#create dictionary for all variables from twomom_d.ncdf
twomom = dict()

#if necessary change name of variables
varlist = twomom_file.variables
#read PAMTRA variables to pamData dictionary
for var in varlist:#read files and write it with different names in Data
    twomom[var] = np.squeeze(twomom_file.variables[var])
for category in ['c','i','r','s','g','h']:
    #calculate Dmean for each category
    twomom['D_mean_' + category] =  __postprocess_SB.calc_Dmean(twomom,category)
print len(plot_vars),num_pam_SB_plots, i+6
ax = plt.subplot2grid((len(plot_vars)+num_pam_SB_plots, 1), (i+6, 0))
ax2 = ax.twiny()

i_timestep=(tstep/60)-1 #there is no output for t=0min after that there are output steps in 30 minute steps (this could vary)
plt = __plotting_functions.plot_twomom_moments(ax,ax2,twomom,i_timestep)

#save figure
plt.tight_layout()
if not os.path.exists('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment): #create direktory if it does not exists
    os.makedirs('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment)
out_filestring = filename #"/adaptv2_" + testcase + "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'])