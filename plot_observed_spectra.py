'''
plot radar spectra from observations (in different frequencies) (designed for Tripex and Tripex-pol campaign); adapted from Jose's script
'''
#from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug

#import modules
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.markers as mk
import numpy as np
from netCDF4 import Dataset
import pandas as pd
import netCDF4
import glob
#self-written functions
from functions import __plotting_functions
import os 
import subprocess
import re
import __radar_postprocessing
import __postprocess_McSnow

#read time which should be analyzed from testcase string (comes from the McSnow_Pamtra_ICONinit.sh script
if "testcase" in os.environ: #this should be the case if it is called by a governing script 
    testcase = os.environ["testcase"] #"more readable" string of the experiment specifications
    date = (re.search(r'day(.*?)hour', testcase).group(1)) #this line gets the date from the testcase string
    year = date[0:4]
    month = date[4:6]
    day = date[6:8]
    hourBegin = (re.search(r'hour(.*?)min', testcase).group(1)) #this line gets the hour from the testcase string
    minBegin = (re.search(r'min(.*?)s', testcase).group(1)) #this line gets the min from the testcase string
    hourEnd = str(int(hourBegin)+1).zfill(2)  #TODO: dont have to average always one hour
    minEnd = '00' #TODO: dont have to average always full hours
    secBegin = '00'
    secEnd = '00'
    
    MC_dir = os.environ["MC"]
else: #set the date here, if you run this script directly
    year = '2018' 
    month = '11'
    day = '01' 

    hourBegin = '13'
    hourEnd = '14'

    minBegin = '00'
    minEnd =  '00'

    secBegin = '00'
    secEnd = '00'


#TODO: read campaign name from governing script (so far Tripex or Tripex-Pol)
if int(year)<2018: #thats just for lazy people
    campaign = 'tripex' #'tripex''tripex-pol'

else:
    campaign = 'tripex-pol' #'tripex''tripex-pol'


#define start and end as pandas time
start = pd.datetime(int(year), int(month), int(day),
                    int(hourBegin), int(minBegin),
                    int(secBegin))


end = pd.datetime(int(year), int(month), int(day),
                  int(hourEnd), int(minEnd), int(secEnd))

#read heightbounds to plot it consistently with the 1D-simulation
if "testcase" in os.environ:
    height_bounds = __postprocess_McSnow.read_init_vals(MC_dir) #read the heights from the init_vals.txt file

else: #no height limits if it is not run by the McSnow_Pamtra_ICONinit script
    height_bounds = [0,10000] #set some default heigh-bounds
    
##
#read the data into the obsData dictionary
##
obsData = dict() #initialize the dictionary for the observational data which gets the same nomenclatura than the pamtra Dictionary

#create a plot
num_plots=8
figsize_height = 6.0/2.0*(num_plots)
figmomwaterfall	=	plt.figure(num=1,figsize=(8.0,figsize_height))
number_of_spectrograms = 3 #3 frequencies
figsize_height = 6.0/2.0*(number_of_spectrograms)
figspectro = plt.figure(num=2,figsize=(8.0,figsize_height))

#switch to current figure
plt.figure(figmomwaterfall.number)
#define axis
axrefl = plt.subplot2grid((num_plots, 1), (0, 0))
axvDoppler = plt.subplot2grid((num_plots, 1), (1, 0))
axswidth = plt.subplot2grid((num_plots, 1), (2, 0))
axskewn = plt.subplot2grid((num_plots, 1), (3, 0))
axwaterfall = plt.subplot2grid((num_plots, 1), (4, 0), rowspan=4)
linestyle='-'

for i_freq,radarname in enumerate(['joyrad10','joyrad35','grarad94']):
    obsData["frequency"] = np.array([np.array([10.0,35.5,94.0])[i_freq]]) #this might look a bit weird, but it is working with plot_pamtra_Ze and plot_pamtra_highermoments
    #from IPython.core.debugger import Tracer ; Tracer()()
    obsData = __radar_postprocessing.read_tripex_pol_radar(obsData,year,month,day,hourBegin,hourEnd,i_freq=0,minBegin='00',minEnd='00',radar_name=radarname,onlymoments=False) 
    
    color=np.array(['b','r','g'])[i_freq]
    
    #switch to current figure
    plt.figure(figmomwaterfall.number)
    #plot reflectivities
    axrefl = __plotting_functions.plot_pamtra_Ze(axrefl,obsData,linestyle=linestyle,nolabel=True,forcedcolor=color)
    #plot other radar-moments
    axvDoppler = __plotting_functions.plot_pamtra_highermoments(axvDoppler,obsData,linestyle=linestyle,moment='vDoppler',forcedcolor=color)
    axswidth = __plotting_functions.plot_pamtra_highermoments(axswidth,obsData,linestyle=linestyle,moment='swidth',forcedcolor=color)
    axskewn = __plotting_functions.plot_pamtra_highermoments(axskewn,obsData,linestyle=linestyle,moment='skewn',forcedcolor=color)
    
    #Radar_Spectrum #Radar_Velocity #frequency #height
    axwaterfall = __plotting_functions.plot_waterfall(axwaterfall,obsData,freq=np.array([10.0,35.5,94.0])[i_freq],color=color,vel_lim=[0,3]) #,z_lim=[height_bounds[1],height_bounds[0]],linestyle=linestyle)
    
    ######
    #plot spectrograms
    ######
    #switch to current figure
    plt.figure(figspectro.number)
    #define a new axis for each frequency
    ax = plt.subplot2grid((number_of_spectrograms, 1), (i_freq, 0))
    ax = __plotting_functions.plot_pamtra_spectrogram(ax,obsData,freq=obsData["frequency"])
    ax.set_ylim([height_bounds[1],height_bounds[0]]) #this might be the range from the initialization height to the maximum of the masses
    
###########
#process and save figmomwaterfall
###########
#switch to current figure
plt.figure(figmomwaterfall.number)

#add labels to all axes and set height limits
for ax in [axrefl,axvDoppler,axswidth,axskewn]:
    ax.plot(np.nan,np.nan,color='b',linestyle='-',label='10GHz')
    ax.plot(np.nan,np.nan,color='r',linestyle='-',label='35GHz')
    ax.plot(np.nan,np.nan,color='g',linestyle='-',label='94GHz')


    ax.legend()
    ax.set_ylim([height_bounds[1],height_bounds[0]]) #this might be the range from the initialization height to the maximum of the masses
#set waterfall labels and limits
axwaterfall.plot(np.nan,np.nan,color='b',linestyle='-',label='10GHz')
axwaterfall.plot(np.nan,np.nan,color='r',linestyle='-',label='35GHz')
axwaterfall.plot(np.nan,np.nan,color='g',linestyle='-',label='94GHz')
axwaterfall.set_ylim([height_bounds[1],height_bounds[0]]) #this might be the range from the initialization height to the maximum of the masses
axwaterfall.legend()

#save figure
plt.tight_layout()
savepath = '/home/mkarrer/Dokumente/plots/observations/spectra/' + year + month + day + '/'
if not os.path.exists(savepath): #create direktory if it does not exists
    os.makedirs(savepath)
filename = 'spectra_' + hourBegin + minBegin + secBegin + '-' + hourEnd + minEnd + secEnd
plt.savefig(savepath + filename + '.pdf', dpi=400)
plt.savefig(savepath + filename + '.png', dpi=400)
print 'The pdf is at: ' + savepath +  filename + '.pdf'
subprocess.Popen(['evince',savepath + filename + '.pdf'])

###########
#process and save figspectro
###########
#switch to current figure
plt.figure(figspectro.number)

plt.tight_layout()
savepath = '/home/mkarrer/Dokumente/plots/observations/spectra/' + year + month + day + '/'

filename = 'spectrogram_' + hourBegin + minBegin + secBegin
plt.savefig(savepath + filename + '.pdf', dpi=400)
plt.savefig(savepath + filename + '.png', dpi=400)
print 'The pdf is at: ' + savepath + filename + '.pdf'
subprocess.Popen(['evince',savepath + filename + '.pdf'])