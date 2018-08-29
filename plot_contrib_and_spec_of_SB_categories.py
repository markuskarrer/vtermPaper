'''
plot the output of the PAMTRA simulations which analyze the individual categories
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
from functions import __postprocess_PAMTRA

#read variables passed by shell script
tstep = int(os.environ["tstep"])
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"]
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output

#directory of experiments
directory = "/home/mkarrer/Dokumente/McSnow/MCSNOW/experiments/"

#define figure
number_of_plots = 8
figsize_height = 6.0/2.0*(number_of_plots)
fig	=	plt.figure(figsize=(8.0,figsize_height))#figsize=(4, 4))

deactlist = ["111111","100000","010000","001000","000100","000010","000001"] #set values to zero to deactivate particle types; order: c,i,r,s,g,h
linestylelist = ["-","--",":","--",":","-.","-."]
markerlist = [" ","o","s","h","x","d","*"]
labellist = ["","cloud w.","cloud ice","rain","snow","graupel","hail"]
#plot reflectivities
ax = plt.subplot2grid((number_of_plots, 1), (0, 0))   
for i,deact in enumerate(deactlist): 
    #define file to read
    if deact=='111111':
        pam_filestring = directory + experiment + "/PAMTRA_2mom_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min.nc'
    else:
        pam_filestring = directory + experiment + "/PAMTRA_2mom_" + testcase + '_cat' + deact + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min.nc'
        ax.plot(0,-1000,color='k',linestyle=linestylelist[i],marker=markerlist[i],markerfacecolor="None",label=labellist[i]) #for the legend

    #read pamtra output to pamData dictionary
    pamData = __postprocess_PAMTRA.read_pamtra(pam_filestring)


    #plot reflectivities
    plt = __plotting_functions.plot_pamtra_Ze(ax,pamData,linestyle=linestylelist[i],marker=markerlist[i])
    
#create legend
ax.legend()

#plot spectrograms
for i,deact in enumerate(deactlist): 
    #open new panel for each category
    ax = plt.subplot2grid((number_of_plots, 1), (i+1, 0))   

    #define file to read
    if deact=='111111':
        pam_filestring = directory + experiment + "/PAMTRA_2mom_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min.nc'
    else:
        pam_filestring = directory + experiment + "/PAMTRA_2mom_" + testcase + '_cat' + deact + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min.nc'
        #ax.plot(0,-1000,color='k',linestyle=linestylelist[i],marker=markerlist[i],markerfacecolor="None",label=labellist[i]) #for the legend

    #read pamtra output to pamData dictionary
    pamData = __postprocess_PAMTRA.read_pamtra(pam_filestring)
    #plot spectrogram
    plt = __plotting_functions.plot_pamtra_spectrogram(ax,pamData,freq=35.5)
    #display category name at top left
    plt.text(-2.8,pamData["height"][-1]*0.9,labellist[i])
    
#save figure
plt.tight_layout()
if not os.path.exists('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment): #create direktory if it dies not exists
    os.makedirs('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment)
out_filestring = "/SB_PAMTRA_categories_" + testcase + "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'])