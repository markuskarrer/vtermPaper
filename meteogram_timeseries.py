'''
plot METEOGRAM output from ICON simulation as a timeseries (with the SB06 two-moment scheme)
'''

#from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug

#import modules
import netCDF4 #for reading netCDF4
import numpy as np 
import pandas as pd #necessary?
import os #for reading variables from shell script
import re #to search in string
import datetime
import csv
import matplotlib.pyplot as plt
import sys
import subprocess
import scipy.signal as scisi #import argrelextrema #to find local maxima

from functions import __plotting_functions
from functions import __general_utilities

#read variables passed by shell script
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"] #"more readable" string of the experiment specifications
MC_dir = os.environ["MC"]


#read time which should be analyzed from testcase string
date = int(re.search(r'day(.*?)hour', testcase).group(1)) #this line gets the date from the testcase string
#hour = int(re.search(r'hour(.*?)min', testcase).group(1)) #this line gets the hour from the testcase string
#minute = int(re.search(r'min(.*?)s', testcase).group(1)) #this line gets the min from the testcase string

#TODO: choose file based on experiment name
filename = "/data/inscape/icon/experiments/tripex_220km/METEOGRAM_patch002_joyce.nc"

#open netCDF4 file
nc = netCDF4.Dataset(filename)
#get formated string of date
timestrings = pd.to_datetime(netCDF4.chartostring(nc.variables['date'][:]))
unixtime = timestrings

#initialize dictionary
#twomom = dict()
#load some variables to dictionary
#twomom["heights"] = nc.variables["height_2"] #mid-level height (relevant for most variables (except: W,.. )

icon_vars = ["temp","pres","qv","qc","qnc","qr","qnr","qi","qni","qs","qns","qg","qng","qh","qnh","rho","rhw"]
derived_vars = ["rhi"]
#setup figure
number_of_plots = len(icon_vars)+len(derived_vars)
figsize_height = 6.0/2.0*(number_of_plots)
fig, axes	=	plt.subplots(nrows=number_of_plots, ncols=1, figsize=(8.0,figsize_height))

for i_plots, (new_key,icon_key) in enumerate(zip(icon_vars, #["temp","pres","qv","qc","qnc","qr","qnr","qi","qni","qs","qns","qg","qng","qh","qnh","rho","rhw"],
                                                                    ["T",   "P",   "QV","QC","QNC","QR","QNR","QI","QNI","QS","QNS","QG","QNG","QG","QNH","RHO","REL_HUM"])):
    #uncomment next 2 lines to plot only one variable for testing
    #if not i_plots==16:
    #    continue
    print "plotting", new_key, "min:",np.amin(nc.variables[icon_key][:]), "max:",np.amax(nc.variables[icon_key][:])
    axes[i_plots] == __plotting_functions.pcolor_timeseries(fig,axes[i_plots],unixtime,nc.variables["height_2"],np.transpose(nc.variables[icon_key]),varname=new_key,time_formatted=timestrings,unit=nc.variables[icon_key].unit)
##
#do additional plots with variables derived from the one in the Meteogram output
##
for j_plots, (new_key) in enumerate(derived_vars): #= ["rhi"]
    print "calculate ",new_key,"from Meteogram output"
    if new_key=="rhi":
        #calculate rhi (relative humidity over ice)
        iconData_atmo=dict(); iconData_atmo["T"]=nc.variables["T"][:,:];iconData_atmo["rh"] = nc.variables["REL_HUM"][:,:]; iconData_atmo["z"]=nc.variables["height_2"][:]; derived_var = __general_utilities.calc_rhi(iconData_atmo)
        unit = "%"
    #plot the derived variables below the one in the Meteogram output    
    axes[i_plots+j_plots+1] == __plotting_functions.pcolor_timeseries(fig,axes[i_plots+j_plots+1],unixtime,nc.variables["height_2"],np.transpose(derived_var),varname=new_key,time_formatted=timestrings,unit=unit)
        
#save figure
plt.tight_layout()
dir_timeseries = '/home/mkarrer/Dokumente/plots/Meteogram/' + str(date)
#from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug


if not os.path.exists(dir_timeseries): #create direktory if it does not exists
    os.makedirs(dir_timeseries)

out_filestring = "timeseries"
plt.savefig(dir_timeseries  + "/" + out_filestring + '.pdf', dpi=400)
plt.savefig(dir_timeseries  + "/" + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + dir_timeseries +  "/" + out_filestring + '.pdf'
subprocess.Popen(['evince',dir_timeseries + "/" + out_filestring + '.pdf'])