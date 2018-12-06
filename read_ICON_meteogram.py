'''
read METEOGRAM output from ICON simulation (with the SB06 two-moment scheme)
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

#local functions
def add_possible_initialization_heights(ax,heights):
    '''
    add vertical lines to the axis "as" to visualize the potential heights to initialize the model
    INPUT:  ax: axes on which the line should be added
            heights: y-values (in m) where the line should be added
    '''
    for i in range(0,len(heights)):
        ax = plt.axhline(heights[i],color='grey',linestyle='--')
        
    return ax
    
#read variables passed by shell script
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"] #"more readable" string of the experiment specifications
MC_dir = os.environ["MC"]


#read time which should be analyzed from testcase string
date = (re.search(r'day(.*?)hour', testcase).group(1)) #this line gets the date from the testcase string
hour = int(re.search(r'hour(.*?)min', testcase).group(1)) #this line gets the hour from the testcase string
minute = int(re.search(r'min(.*?)s', testcase).group(1)) #this line gets the min from the testcase string
ssat = int(os.environ["ssat"]) #get supersaturation from governing script (this is the semi-idealized part of this setup)
#TODO also read in seconds
#convert time to s to find nearest output timestep
analyzed_time_s = hour*3600+minute*60

#read meteogram file
veras_tripexpol_simul = "/data/inscape/icon/experiments/juelich/testbed/testbed_"
if date=='20151124':
    filename = "/data/inscape/icon/experiments/tripex_220km/METEOGRAM_patch002_joyce.nc"
elif date.startswith('2018') or date.startswith('2019'): #TODO: this is a bit dirty in case you have other simulations in 2018/2019
    filename = veras_tripexpol_simul + date + "/METEOGRAM_patch001_" + date + "_joyce.nc"
#open netCDF4 file
nc = netCDF4.Dataset(filename)
#get formated string of date
time = pd.to_datetime(netCDF4.chartostring(nc.variables['date'][:]))
#get timestep in array for analysis
timestep_start = np.argmin(np.absolute(nc.variables["time"][:]-analyzed_time_s))


#initialize dictionary
twomom = dict()
#load some variables to dictionary
twomom["heights"] = nc.variables["height_2"] #mid-level height (relevant for most variables (except: W,.. )

#TODO: do not hardcode average time
average_time_s = 60*60 #60 minute average
timestep_end = np.argmin(np.absolute(nc.variables["time"][:]-analyzed_time_s-average_time_s))
#print some info in terminal
print 'analyzing time: ',time[timestep_start],' to: ',time[timestep_end]

for new_key,icon_key in zip(["temp","pres","qv","qc","qnc","qr","qnr","qi","qni","qs","qns","qg","qng","qh","qnh","rho","rhw"],
                            ["T",   "P",   "QV","QC","QNC","QR","QNR","QI","QNI","QS","QNS","QG","QNG","QG","QNH","RHO","REL_HUM"]):
    twomom[new_key] = np.mean(nc.variables[icon_key][timestep_start:timestep_end,:],axis=0) #ATTENTION: from here on they are temporally averaged
    if new_key in ("qc","qnc","qr","qnr","qi","qni","qs","qns","qg","qng","qh","qnh"):
        twomom[new_key] = twomom[new_key][None,:] #add singleton dimension for hydrometeors because of plotting routine
    twomom[new_key + '_std'] = np.std(nc.variables[icon_key][timestep_start:timestep_end,:],axis=0)
    if new_key in ("qc","qnc","qr","qnr","qi","qni","qs","qns","qg","qng","qh","qnh"):
        twomom[new_key + '_std'] = twomom[new_key + '_std'][None,:] #add singleton dimension

#twomom["temp"]=np.ones(twomom["temp"].shape)*273.15

###
#calculate extra variables (f.e. rhi) from other atmospheric variables
###
iconData_atmo = dict()
#calculate rhi (relative humidity over ice)
iconData_atmo["T"]=twomom["temp"];iconData_atmo["rh"] = twomom["rhw"]; iconData_atmo["z"]=twomom["heights"][:]; iconData_atmo["p"]=twomom["pres"]; twomom["rhi"] = __general_utilities.calc_rhi(iconData_atmo)
#calculate std also for rhi
twomom["rhi_std"] = np.std(nc.variables[icon_key][timestep_start:timestep_end,:],axis=0)

#ATTENTION: for testing: fix atmosphere to some value
#twomom["temp"]=np.ones(twomom["temp"].shape)*250.15
#twomom["qv"]=np.ones(twomom["qv"].shape)*0.00010
#twomom["pres"]=np.ones(twomom["pres"].shape)*101325.0
#twomom["rhi"]=np.ones(twomom["pres"].shape)*101.000
#twomom["rhi_std"]=np.ones(twomom["pres"].shape)*0.0
#print "qv_before",twomom["qv"]
#twomom["qv"] = __general_utilities.q2abs(twomom["qv"],twomom["qv"],twomom["temp"],twomom["pres"],twomom["rhw"],q_all_hydro="NaN") #qv is now in kg m-3
#print "qv_after",twomom["qv"]
#raw_input("")

        
#convert from mixing ratios [kg kg-1]/[kg-1] to absolute values [kg m-3]/[m-3]
twomom["qi_spec"] = __general_utilities.q2abs(twomom["qi"],twomom["qv"],twomom["temp"],twomom["pres"],q_all_hydro=(twomom["qc"]+twomom["qr"]+twomom["qi"]+twomom["qs"]+twomom["qg"]+twomom["qh"]))
twomom["qni_spec"] = __general_utilities.q2abs(twomom["qni"],twomom["qv"],twomom["temp"],twomom["pres"],q_all_hydro=(twomom["qc"]+twomom["qr"]+twomom["qi"]+twomom["qs"]+twomom["qg"]+twomom["qh"]))
#do you want to use the later two??
twomom["qs_spec"] = __general_utilities.q2abs(twomom["qs"],twomom["qv"],twomom["temp"],twomom["pres"],q_all_hydro=(twomom["qc"]+twomom["qr"]+twomom["qi"]+twomom["qs"]+twomom["qg"]+twomom["qh"]))
twomom["qns_spec"] = __general_utilities.q2abs(twomom["qns"],twomom["qv"],twomom["temp"],twomom["pres"],q_all_hydro=(twomom["qc"]+twomom["qr"]+twomom["qi"]+twomom["qs"]+twomom["qg"]+twomom["qh"]))
###
#find maximum in number density to initialize McSnow
###
criteria = 0 #0: choose global maximum of qni 1: choose lowest local maximum of qni 2. choose heighest level with RHi>100 3. define an arbitrary height #TODO: give this as an input by the governing script (kind of namelist)
#calculate anyway the initialization for all criterias to display them in panel but choose just one for the real initialization
num_crit=4; i_crit = np.zeros(num_crit, dtype=np.int) #set number of possible criterias to choose initialization height
qi_init = np.zeros(num_crit); qni_init = np.zeros(num_crit) ; height_init = np.zeros(num_crit) #initialize arrays which contain the values for initialization
i_crit[0] = int(np.argmax(twomom["qni"]));  #index of initialization height after criteria 0

qni_small = 1.0; i=-1
oncemore = iter([True, False])
i_loc_max_list = -1
loc_max_list = scisi.argrelextrema(twomom["qni"][0], np.greater)[0]
#from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug
while (twomom["qni"][0,loc_max_list[i_loc_max_list]]<qni_small): # or next(oncemore): #move on to next local minimum if qni>qni_small is not fullfiled; next(oncemore) achieves that the while-clause is also executed for the actual interesting minimum (where qni_small is exceeded first)
    print twomom["qni"][0,loc_max_list[i_loc_max_list]]
    i_crit[1] = loc_max_list[i_loc_max_list] #int(scisi.argrelextrema(twomom["qni"][0], np.greater)[0][-1]) #find first/next local maxima    
    #print i_loc_max_list,i_crit[1],twomom["qni"][0,loc_max_list[i_loc_max_list]],twomom["qi"][0,loc_max_list[i_loc_max_list]]
    i_loc_max_list-=1

#convert_atmovars(calc_something="calc_rhi",input_dict=iconData_atmo) #from __general_utilities
try:
    i_crit[2] = np.where(twomom["rhi"]>100)[0][0] #find heighest level with RHi> 0
except:
    print "WARNING: complete profile is subsaturated with reference to -> ice cannot find heighest level with RHi> 0 in read_ICON_meteogram.py"
    i_crit[2] = -1
    
#criteria=4: define an arbitrary height for initialization
height2init = 6000#define height for initialization here
__,i_crit[3] = __general_utilities.find_nearest(twomom["heights"],height2init)

for i in range(0,num_crit):
    #if date.startswith('2018'):
    #    qi_init[i] = twomom["qi"][0,i_crit[i]]+twomom["qs"][0,i_crit[i]]
    #    qni_init[i] = round(twomom["qni"][0,i_crit[i]]+twomom["qs"][0,i_crit[i]],0)
    #else:
    qi_init[i] = twomom["qi_spec"][0,i_crit[i]]+twomom["qs_spec"][0,i_crit[i]]
    qni_init[i] = round(twomom["qni_spec"][0,i_crit[i]]+twomom["qs_spec"][0,i_crit[i]],0)
    height_init[i] = round(twomom["heights"][i_crit[i]])

#ATTENTION: remove this (just for testing)
#qi_init[0] = 5e-4
#qni_init[0] = 1e5

#copy values from selected condition here to those who are written to file init_vals.txt later
qi_init_sel = qi_init[criteria] 
qni_init_sel = qni_init[criteria]
height_init_sel = height_init[criteria]

#get global maximum of qi+qs
i_maxq = int(np.argmax(twomom["qi"]+twomom["qs"]))
height_maxq = twomom["heights"][i_maxq]
q_maxq = twomom["heights"][i_maxq]
qn_maxq = twomom["heights"][i_maxq]

from IPython.core.debugger import Tracer ; Tracer()()
####
#in the following are the modifications (of the qv-profile) which makes the setup SEMI-idealized
####
constant_rhi_flag=True #if True define a constant rhi
if constant_rhi_flag:
    ssat_vec=np.ones_like(twomom["temp"])*ssat/100000. #ssat is in [1000%]
    ssat_vec[i_maxq:] = -5./100 #take -5% for sublimation for now

    #rhi_const=100. + ssat/1000. #TODO: this should be an INPUT from the loop in the shell script which finds the suitable rhi
    #calculate rh (relative humidity over water) from the in the previous line defined rhi (relative humidity over ice)
    #iconData_atmo["rhi"] = np.ones_like(twomom["temp"])*rhi_const
    #iconData_atmo["rhi"][i_maxq:] = np.ones_like(twomom["temp"][i_maxq:])*90.
    #twomom["rhi"]=np.ones_like(twomom["temp"])*rhi_const
    #iconData_atmo["T"]=twomom["temp"]; iconData_atmo["z"]=twomom["heights"][:]; iconData_atmo["p"]=twomom["pres"]; twomom["rh"] = __general_utilities.calc_rh(iconData_atmo)
    #calculate qv corresponding to the fix rhi and overwrite it
    #T_gt_0 = iconData_atmo["T"]>=273.15
    #T_lt_0 = iconData_atmo["T"]<273.15
    #twomom["qv_fixed_rhi"] = np.zeros_like(twomom["qv"])
    #twomom["pres"]=np.ones_like(twomom["pres"])*1e3
    #twomom["qv_fixed_rhi"][T_lt_0] = __general_utilities.rh2vap(twomom["temp"][T_lt_0],twomom["pres"][T_lt_0]/100.,iconData_atmo["rh"][T_lt_0]) #taken from PAMTRA
    #twomom["qv_fixed_rhi"][T_lt_0] = __general_utilities.rh2mixr(twomom["rh"][T_lt_0]/100.,twomom["pres"][T_lt_0],twomom["temp"][T_lt_0]) #taken from https://github.com/hendrikwout/meteo/blob/master/meteo/humidity.py
    #twomom["qv_fixed_rhi"][T_lt_0] = __general_utilities.calc_qv_from_rh(twomom["rh"][T_lt_0],twomom["pres"][T_lt_0],twomom["temp"][T_lt_0]) #taken from SB
    #twomom["qv_fixed_rhi"][T_gt_0] = twomom["qv"][T_gt_0] #treat T>0Cels differently because rhi is not defined there
    #twomom["qv_fixed_rhi"] = __general_utilities.calc_qv_from_rh(twomom["rh"],twomom["pres"],twomom["temp"]) #taken from SB
    #twomom["qv_fixed_rhi"][:i_maxq] = __general_utilities.calc_qv_from_rh(twomom["rh"][:i_maxq],twomom["pres"][:i_maxq],twomom["temp"][:i_maxq]) #taken from SB
    #twomom["qv_fixed_rhi"][i_maxq:] = __general_utilities.calc_qv_from_rh(twomom["rh"][:i_maxq],twomom["pres"][:i_maxq],twomom["temp"][:i_maxq])  #twomom["qv"][i_maxq:]*0.0 #treat heights below the maximum in mass differently because there should be sublimation
    

    #TODO: there is some bug here because rhi in the McSnow run fits not perfectly!!
    
    #from IPython.core.debugger import Tracer ; Tracer()()

#write atmospheric variables to txt file
with open(MC_dir + "/input/ecmwf_profile.txt","wb") as txtfile: #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in
    ecmwf_writer =csv.writer(txtfile, delimiter=' ', quoting=csv.QUOTE_NONE, lineterminator=os.linesep) #quoting avoids '' for formatted string; lineterminator avoids problems with system dependend lineending format https://unix.stackexchange.com/questions/309154/strings-are-missing-after-concatenating-two-or-more-variable-string-in-bash?answertab=active#tab-top
    for row in range(0,len(nc.variables["height_2"][:])):
        #ecmwf_writer.writerow([nc.variables["height_2"][row]] + [twomom["temp"][row]] + [twomom["pres"][row]] + [twomom["qv_fixed_rhi"][row]] ) #atmo_type=2 READ (unit,*,iostat=io) zz(i),tt(i),pp(i),qv(i),dummy
        #from IPython.core.debugger import Tracer ; Tracer()()
        #ecmwf_writer.writerow([nc.variables["height_2"][row]] + [270.15*np.ones_like(twomom["temp"][row])] + [twomom["pres"][row]] + [ssat_vec[row]] ) #atmo_type=5 READ(unit,*,iostat=io) pp(i), zz(i), tt(i), dummy, rh(i)
        ecmwf_writer.writerow([nc.variables["height_2"][row]] + [twomom["temp"][row]] + [twomom["pres"][row]] + [ssat_vec[row]] ) #atmo_type=2 READ (unit,*,iostat=io) zz(i),tt(i),pp(i),ssat(:)
        #print nc.variables["height_2"][row],twomom["temp"][row],twomom["pres"][row],twomom["qv"][row]
txtfile.close()
#write hydrometeor init values to txt file
with open(MC_dir + "/input/init_vals.txt","wb") as txtfile: #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in
    initvals_writer = csv.writer(txtfile, delimiter=' ', quoting=csv.QUOTE_NONE, lineterminator=os.linesep) #quoting avoids '' for formatted string; lineterminator avoids problems with system dependend lineending format https://unix.stackexchange.com/questions/309154/strings-are-missing-after-concatenating-two-or-more-variable-string-in-bash?answertab=active#tab-top
    initvals_writer.writerow(["{:010.0f}".format(height_init_sel)] + ["{:010.0f}".format(qi_init_sel*10000000)] + ["{:010.0f}".format(qni_init_sel/100)])
    initvals_writer.writerow(["{:010.0f}".format(height_init_sel)] + ["{:010.0f}".format(qi_init_sel*10000000)] + ["{:010.0f}".format(qni_init_sel/100)])

txtfile.close()
print "wrote init vals to:" + MC_dir + "/input/init_vals.txt"

number_of_plots = 4

flag_plot_icons_var = True #for debugging purposes, plot here icons variables (atmosphere + hydrometeors)
if flag_plot_icons_var:
    figsize_height = 6.0/2.0*(number_of_plots)
    fig	=	plt.figure(figsize=(8.0,figsize_height))#figsize=(4, 4))
    '''
    create plot of athmospheric variables first and add it before the histogram plots
    '''
    ax  = plt.subplot2grid((number_of_plots, 1), (0, 0))
    ax2 = ax.twiny()
    #plot_atmo in __plotting_functions.py needs special names in the dictionary
    iconData_atmo = dict()
    iconData_atmo["T"] = twomom["temp"];  iconData_atmo["rh"] = twomom["rhw"]; iconData_atmo["z"]=twomom["heights"]
    iconData_atmo["T_std"] = twomom["temp_std"]; iconData_atmo["rh_std"] = twomom["rhw_std"];
    iconData_atmo["rhi"] = twomom["rhi"]; iconData_atmo["rhi_std"] = twomom["rhi_std"];
    #plot atmospheric variables
    __plotting_functions.plot_atmo(ax,ax2,iconData_atmo)
    #add lines for the heights which can be used for initialization
    ax = add_possible_initialization_heights(ax,height_init)
    ####################################
    #plot mixing ratio + number density
    ####################################
    #number density
    mass_num_flag = 0 #0-> plot only number flux; 1-> plot only mass flux; 2-> plot both 

    ax = plt.subplot2grid((number_of_plots, 1), (1, 0))
    if mass_num_flag==2:
        ax2 = ax.twiny()
    else: #in case there is no need for a second axis, just pass the first ax twice
        ax2 = ax
    
    #plot_moments in __plotting_functions.py needs special names in the dictionary
    hei2massdens=dict() #hei2massdens stays empty here

    i_timestep = 0 #dirty workaround: the variables do have just one timesteps here, because this is choose before
    ax = __plotting_functions.plot_moments(ax,ax2,twomom,hei2massdens,i_timestep,mass_num_flag=mass_num_flag)
    #add lines for the heights which can be used for initialization
    ax = add_possible_initialization_heights(ax,height_init)
    #mass density
    mass_num_flag = 1 #0-> plot only number flux; 1-> plot only mass flux; 2-> plot both 

    ax = plt.subplot2grid((number_of_plots, 1), (2, 0))
    if mass_num_flag==2:
        ax2 = ax.twiny()
    else: #in case there is no need for a second axis, just pass the first ax twice
        ax2 = ax
        
    ax = __plotting_functions.plot_moments(ax,ax2,twomom,hei2massdens,i_timestep,mass_num_flag=mass_num_flag)
    #add lines for the heights which can be used for initialization
    ax = add_possible_initialization_heights(ax,height_init)
    #add panel with some information in text form
    ax = plt.subplot2grid((number_of_plots, 1), (3, 0))
    ax.axis("off")
    #from IPython.core.debugger import Tracer ; Tracer()() 
    plt.text(0,0,'analyzing time: {} to: {}'.format(str(time[timestep_start]),str(time[timestep_end])) +
             '\n initialize at: {:d})\n'.format(criteria) +
             '                   0) global maximum of qni: \n' + 
             '                                                 height:{:10.0f} m \n'.format(height_init[0]) + 
             '                                                 qi:       {:.2E} kg m-3 \n'.format(qi_init[0]) + 
             '                                                 qni:     {:.2E} m-3\n'.format(qni_init[0]) + 
             '                   1) lowest local maximum of qni with qni>1: \n' + 
             '                                                 height:{:6.0f}m\n'.format(height_init[1]) +
             '                                                 qi:      {:.2E}kg m-3\n'.format(qi_init[1]) + 
             '                                                 qni:     {:.2E}m-3\n'.format(qni_init[1]) +
             '                   2) heighest level with RHi>100%: \n' + 
             '                                                 height:{:6.0f}m\n'.format(height_init[2]) +
             '                                                 qi:      {:.2E}kg m-3\n'.format(qi_init[2]) + 
             '                                                 qni:     {:.2E}m-3\n'.format(qni_init[2]) +
             'the maximum of qi+qs is at :{:6.0f}m with qi+qs= {:.2E}kg m-3'.format(np.ma.getdata(twomom["heights"])[i_maxq],twomom["qi"][0,i_maxq]+twomom["qs"][0,i_maxq])
             )
    
    #save figure
    plt.tight_layout()
    if not os.path.exists('/home/mkarrer/Dokumente/plots/Meteogram/' + str(date)): # + testcase): #create direktory if it does not exists
        os.makedirs('/home/mkarrer/Dokumente/plots/Meteogram/' + str(date))

    if not os.path.exists('/home/mkarrer/Dokumente/plots/Meteogram/' + str(date) + '/' + testcase): #create direktory if it does not exists
        os.makedirs('/home/mkarrer/Dokumente/plots/Meteogram/' + str(date) + '/' + testcase)
    out_filestring = '/Meteogram_input'
    plt.savefig('/home/mkarrer/Dokumente/plots/Meteogram/' + str(date)+ '/'  + testcase + out_filestring + '.pdf', dpi=400)
    plt.savefig('/home/mkarrer/Dokumente/plots/Meteogram/' + str(date)+ '/'  + testcase + out_filestring + '.png', dpi=400)
    print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/Meteogram/' + str(date) + '/'  + testcase + out_filestring + '.pdf'
    subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/Meteogram/' + str(date) + '/'  + testcase + out_filestring + '.pdf'])

