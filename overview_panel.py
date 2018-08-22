'''
create an overview panel with properties of McSnow and PAMTRA output
'''
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

#read variables passed by shell script
tstep = os.environ["tstep"]
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"]
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output


#directory of experiments
directory = "/home/mkarrer/Dokumente/McSnow/MCSNOW/experiments/"
#experiment name (this also contains a lot of information about the run
#experiment="1d_xi100000_nz5000_lwc20_ncl0_dtc5_nrp30_rm10_rt0_vt2_h10-20_ba500"

#define filestring (including timestep)
#tstep="0300min"
filestring = directory + experiment + "/mass2fr_" + tstep + ".dat"

#read mass2fr.dat file and get SP-dictionary
#SP = __postprocess_McSnow.read_mass2frdat(experiment,filestring)
#load file
SP_file = Dataset(directory + experiment + '/mass2fr_' + tstep + '_avtstep_' + str(av_tstep) + '.ncdf',mode='r')
#create dictionary for all variables from PAMTRA
SP = dict()
#print pam_file.variables.keys() 

#if necessary change name of variables
varlist = SP_file.variables
#read PAMTRA variables to pamData dictionary
for var in varlist:#read files and write it with different names in Data
    SP[var] = np.squeeze(SP_file.variables[var])

#calculate values binned to here defined h-D bins
binned_val,heightvec,d_bound_ds,d_ds,zres = __postprocess_McSnow.seperate_by_height_and_diam(SP,nbins=100,diamrange=[-9,0],nheights=51,model_top=5000)

#calculate volume of box
Vbox = __postprocess_McSnow.calculate_Vbox(experiment,zres)
#divide by box volume to acchieve [#]->[#/m3]
binned_val["d_counts"] = __postprocess_McSnow.conv_num2numpm3(binned_val["d_counts"],Vbox)
binned_val["d_counts_no_mult"] = __postprocess_McSnow.conv_num2numpm3(binned_val["d_counts"],Vbox)

############################
#now: plot the binned values
############################
cmap="viridis_r" #use one uniform colorbar

#list of variables from binned_val-dictionary which should be plotted
plot_vars = ["d_counts","av_Frim","av_rhor","av_mm"] #"RPpSP",
full_name = ["number density","rime fraction","ma rime density","na numb. monomer"] #for colorbar labelling
var_units = ["m-3","1","kg m-3","1"]#for colorbar labelling
#select step of colorbar ticks for each variable
colticksstep = [50000.,0.2,100.,1] #,20000.
maskedcon = ['0','nan','nan','nan']#,'nan' #see if maskedcon = 'nan': within for-loop
mincol_arr = [10,0.,0,1] #0, lowest number in colorbar
maxcol_arr = [1e5,1.0,900,1000] #,-999 heighest number in colorbar
logcol_arr = [1,0,0,1] #set to one if colorbar should be in logarithmic scale

number_of_plots = len(plot_vars)
num_pamplots = 3 #number of subplots from pamtra variables

figsize_height = 6.0/2.0*(number_of_plots+num_pamplots)
fig	=	plt.figure(figsize=(8.0,figsize_height))#figsize=(4, 4))

for i,varname in enumerate(plot_vars):
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
    #plt.rc('xtick', labelsize=15) 
    #plt.rc('ytick', labelsize=15)
    #plt.rc('axes', titelsize=15)
    #plt.rc('ylabel', labelsize=15)
    #plt.rc('ytick', labelsize=15)
    
    ax 	= 	plt.subplot2grid((len(plot_vars)+num_pamplots, 1), (i, 0)) #plotvars +num_pamplots -> let space for 2 PAMTRA_output plots
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


##############################
#now: plot PAMTRA output below
##############################

#define file to read
pam_filestring = directory + experiment + "/adaptv2_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + tstep + '.nc'
ax = plt.subplot2grid((len(plot_vars)+num_pamplots, 1), (i+1, 0))

#read pamtra output to pamData dictionary
pamData = __postprocess_PAMTRA.read_pamtra(pam_filestring)

#plot reflectivities
plt = __plotting_functions.plot_pamtra_Ze(ax,pamData)

ax = plt.subplot2grid((len(plot_vars)+num_pamplots, 1), (i+2, 0))
#plot spectrogram
plt = __plotting_functions.plot_pamtra_spectrogram(ax,pamData,freq=35.5)
ax = plt.subplot2grid((len(plot_vars)+num_pamplots, 1), (i+3, 0))
plt = __plotting_functions.plot_McSnows_vt_in_spectrogram(ax,pamData,SP,experiment)

#save figure
plt.tight_layout()
if not os.path.exists('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment): #create direktory if it dies not exists
    os.makedirs('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment)
out_filestring = "/adaptv2_" + testcase + "_av_" + str(av_tstep) + "_t" + tstep
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'])