'''
this script checks the discrepancy of the mass (so far qi+qs) between the 3-D and the 1D-SB-scheme at the maximum of the mass in the 3D-SB scheme
'''
#from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug

#import modules
import numpy as np
import sys #for some debugging
import re #for selecting numbers from file-string
import subprocess #for using shell commands with subprocess.call()
import os #import environment variables
import csv #writing csv files
from netCDF4 import Dataset,chartostring
import pandas as pd
#self written #these files are linked here from the pythoncode/functions directory
import __postprocess_McSnow

#read variables passed by shell script
tstep = int(os.environ["tstep"])
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"]
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output
MC_dir = os.environ["MC"]


#directory of experiments
directory = MC_dir + "/experiments/"
#load netCDF4 file
twomom_file1D = Dataset(directory + experiment + '/twomom_d.ncdf',mode='r')
#create dictionary for all variables from twomom_d.ncdf
twomom1D = dict()
i_timestep_1D=(tstep/60)-1 #there is no output for t=0min after that there are output steps in 30 minute steps (this could vary)

#if necessary change name of variables
varlist = twomom_file1D.variables
#read PAMTRA variables to pamData dictionary
for var in varlist:#read files and write it with different names in Data
    twomom1D[var] = np.squeeze(twomom_file1D.variables[var][i_timestep_1D])
    
    
#####
#load the 3D-SB data
#####
#get time of initialization from governing script
date = int(re.search(r'day(.*?)hour', testcase).group(1)) #this line gets the date from the testcase string
hour = int(re.search(r'hour(.*?)min', testcase).group(1)) #this line gets the hour from the testcase string
minute = int(re.search(r'min(.*?)s', testcase).group(1)) #this line gets the min from the testcase string
#convert time to s to find nearest output timestep
analyzed_time_s = hour*3600+minute*60
#TODO: choose file based on experiment name
filename = "/data/inscape/icon/experiments/tripex_220km/METEOGRAM_patch002_joyce.nc"

#open netCDF4 file
nc = Dataset(filename)
#get formated string of date
time = pd.to_datetime(chartostring(nc.variables['date'][:]))
#get timestep in array for analysis
timestep_start = np.argmin(np.absolute(nc.variables["time"][:]-analyzed_time_s))


#initialize dictionary
twomom3D = dict()
#load some variables to dictionary
twomom3D["heights"] = nc.variables["height_2"] #mid-level height (relevant for most variables (except: W,.. )
#TODO: do not hardcode average time
average_time_s = 60*30 #120 minute average
timestep_end = np.argmin(np.absolute(nc.variables["time"][:]-analyzed_time_s-average_time_s))
#print some info in terminal
print 'analyzing time: ',time[timestep_start],' to: ',time[timestep_end]

for new_key,icon_key in zip(["temp","pres","qv","qc","qnc","qr","qnr","qi","qni","qs","qns","qg","qng","qh","qnh","rho","rhw"],
                            ["T",   "P",   "QV","QC","QNC","QR","QNR","QI","QNI","QS","QNS","QG","QNG","QG","QNH","RHO","REL_HUM"]):
    twomom3D[new_key] = np.mean(nc.variables[icon_key][timestep_start:timestep_end,:],axis=0) #ATTENTION: from here on they are temporally averaged
    if new_key in ("qc","qnc","qr","qnr","qi","qni","qs","qns","qg","qng","qh","qnh"):
        twomom3D[new_key] = twomom3D[new_key][None,:] #add singleton dimension for hydrometeors because of plotting routine
    twomom3D[new_key + '_std'] = np.std(nc.variables[icon_key][timestep_start:timestep_end,:],axis=0)
    #if new_key in ("qc","qnc","qr","qnr","qi","qni","qs","qns","qg","qng","qh","qnh"):
    #    twomom3D[new_key + '_std'] = twomom3D[new_key + '_std'][None,:] #add singleton dimension
        
#get global maximum of qi+qs in the 3D-run
i_maxq = int(np.argmax(twomom3D["qi"]+twomom3D["qs"]))

#finally check the difference at the  global maximum of qi+qs in the 3D-run
qsum1D = twomom1D["qs"][i_maxq]+twomom1D["qi"][i_maxq]
qsum3D = twomom3D["qs"][0,i_maxq]+twomom3D["qi"][0,i_maxq]
#check if qsum1D is too small (increase ssat), too large or within the tolerance interval
eps_tolerance = qsum3D/10. #10% tolerance is allowed
#write info if qsum1D is right, too high or too low in text file
f = open("/home/mkarrer/Dokumente/bash_organizing_code/tmp_qsumdir.txt", 'w')
if abs(qsum3D-qsum1D)<eps_tolerance:
    f.write('e') #equal
elif qsum1D<qsum3D:
    f.write('l') #low
elif qsum1D>qsum3D:
    f.write('h') #high
f.close()
print "qsum1D,qsum3D,eps_tolerance"
print qsum1D,qsum3D,eps_tolerance
'''
#ATTENTION: for testing searching algorithm
ssat = int(os.environ["ssat"]) #get supersaturation from governing script (this is the semi-idealized part of this setup)
qsum1D = ssat
qsum3D = 3949
print "qsum1D,qsum3D,eps_tolerance"
print qsum1D,qsum3D,eps_tolerance
raw_input("wait")

#write info if qsum1D is right, too high or too low in text file
f = open("/home/mkarrer/Dokumente/bash_organizing_code/tmp_qsumdir.txt", 'w')
#with open("/home/mkarrer/Dokumente/bash_organizing_code/tmp_qsumdir.txt","wb") as txtfile: #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in

if abs(qsum3D-qsum1D)<eps_tolerance:
    f.write('e') #equal
elif qsum1D<qsum3D:
    f.write('l') #low
elif qsum1D>qsum3D:
    f.write('h') #high
       
f.close()
'''
#from IPython.core.debugger import Tracer ; Tracer()()