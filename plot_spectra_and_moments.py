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
separated_by_categories = (os.environ["separated_by_categories"]=="True") #plot a line for individual categories also

#directory of experiments
directory = MC_dir + "/experiments/"

num_plots = 8 #number of subplots from pamtra variables and twomoment output
figsize_height = 6.0/2.0*(num_plots)
fig	=	plt.figure(figsize=(8.0,figsize_height))#figsize=(4, 4))


if "semiidealized" in os.environ:
    height_bounds = __postprocess_McSnow.read_init_vals(MC_dir) #read the heights from the init_vals.txt file

else: #no height limits if it is not run by the McSnow_Pamtra_ICONinit script
    height_bounds = [5000,0] #set some default heigh-bounds
    
#plot McSnow pamtra output
try:
    ##############################
    #now: plot PAMTRA output 
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

    #read pamtra output to pamData dictionary
    pamData = __postprocess_PAMTRA.read_pamtra(pam_filestring)
    try:
        axwaterfall = plt.subplot2grid((num_plots, 1), (4, 0), rowspan=4)
        #plot spectrogram
        axwaterfall = __plotting_functions.plot_waterfall(axwaterfall,pamData,freq=9.6,color='b',linestyle='--')
        axwaterfall = __plotting_functions.plot_waterfall(axwaterfall,pamData,freq=35.5,color='r',linestyle='--')
        axwaterfall = __plotting_functions.plot_waterfall(axwaterfall,pamData,freq=95.0,color='g',linestyle='--')
    except:
        print "no spectral data found in:", pam_filestring, ". Is radar_mode set to simple or moments?"
    #plot reflectivities
    axrefl = plt.subplot2grid((num_plots, 1), (0, 0))

    axrefl = __plotting_functions.plot_pamtra_Ze(axrefl,pamData,linestyle='--')
    try:
        #plot other radar-moments
        axvDoppler = plt.subplot2grid((num_plots, 1), (1, 0))
        axvDoppler = __plotting_functions.plot_pamtra_highermoments(axvDoppler,pamData,linestyle='--',moment='vDoppler')
        axswidth = plt.subplot2grid((num_plots, 1), (2, 0))
        axswidth = __plotting_functions.plot_pamtra_highermoments(axswidth,pamData,linestyle='--',moment='swidth')
        axskewn = plt.subplot2grid((num_plots, 1), (3, 0))
        axskewn = __plotting_functions.plot_pamtra_highermoments(axskewn,pamData,linestyle='--',moment='skewn')



    except Exception:
	print ' \n \n in except: \n \n'
	s = traceback.format_exc()
   	serr = "there were errors:\n%s\n" % (s)
    	sys.stderr.write(serr) 
        print "no spectral data found in:", pam_filestring, ". Is radar_mode set to simple?"
except Exception:
	print ' \n \n in except: \n \n'
	s = traceback.format_exc()
   	serr = "there were errors:\n%s\n" % (s)
    	sys.stderr.write(serr) 
  
#plot SB pamtra output data
try:
    ##############################
    #now: plot PAMTRA output 
    ##############################

    if separated_by_categories:
        category_list = ["_cat010000","_cat000100",""] #,"_cat001000","_cat000100"] #,"_cat000010","_cat000001"] #,""]
    else:
        category_list = [""]
        
    for i_cat,cat_now in enumerate(category_list): #loop over different category selections
        #define file to read
        pam_filestring = directory + experiment + "/PAMTRA_2mom_" + testcase + cat_now + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min.nc'

        #read pamtra output to pamData dictionary
        pamData = __postprocess_PAMTRA.read_pamtra(pam_filestring)
        if not "axrefl" in globals(): #in case McSnow has not been plotted above 
            axrefl = plt.subplot2grid((num_plots, 1), (0, 0))
        #plot reflectivities
        if cat_now=="": #this is the full PAMTRA run with all categories
            axrefl = __plotting_functions.plot_pamtra_Ze(axrefl,pamData)
        else: 
            axrefl = __plotting_functions.plot_pamtra_Ze(axrefl,pamData,onlyfreq=1,forcedcolor=np.array(['orange','black','red'])[i_cat],forcedlabel=np.array(['cloud ice','snow'])[i_cat])

        if cat_now=="": #this is the full PAMTRA run with all categories
            #add a label distinguishing the models
            axrefl.plot(np.nan,np.nan,color='k',linestyle='-',label='SB')
            axrefl.plot(np.nan,np.nan,color='k',linestyle='--',label='McSnow')

            axrefl.legend()
            axrefl.set_ylim([height_bounds[1],height_bounds[0]]) #this might be the range from the initialization height to the maximum of the masses

        else:
            pass

        #plot other radar-moments
        try:
            if not "axvDoppler" in globals(): #initialize axis if there was no McSnow output
                axvDoppler = plt.subplot2grid((num_plots, 1), (1, 0))
                axswidth = plt.subplot2grid((num_plots, 1), (2, 0))
                axskewn = plt.subplot2grid((num_plots, 1), (3, 0))        
                
            #plot reflectivities
            if cat_now=="": #this is the full PAMTRA run with all categories
                axvDoppler = __plotting_functions.plot_pamtra_highermoments(axvDoppler,pamData,linestyle='-',moment='vDoppler')
                axswidth = __plotting_functions.plot_pamtra_highermoments(axswidth,pamData,linestyle='-',moment='swidth')
                axskewn = __plotting_functions.plot_pamtra_highermoments(axskewn,pamData,linestyle='-',moment='skewn')            
            else: 
                axvDoppler = __plotting_functions.plot_pamtra_highermoments(axvDoppler,pamData,onlyfreq=1,forcedcolor=np.array(['orange','black','red'])[i_cat],forcedlabel=np.array(['cloud ice','snow'])[i_cat],moment='vDoppler')
                axswidth = __plotting_functions.plot_pamtra_highermoments(axswidth,pamData,onlyfreq=1,forcedcolor=np.array(['orange','black','red'])[i_cat],forcedlabel=np.array(['cloud ice','snow'])[i_cat],moment='swidth')
                axskewn = __plotting_functions.plot_pamtra_highermoments(axskewn,pamData,onlyfreq=1,forcedcolor=np.array(['orange','black','red'])[i_cat],forcedlabel=np.array(['cloud ice','snow'])[i_cat],moment='skewn')

            if cat_now=="": #this is the full PAMTRA run with all categories
                for ax in [axvDoppler,axswidth,axskewn]:         #add legend
                    ax.plot(np.nan,np.nan,color='b',linestyle='-',label='9.6GHz')
                    ax.plot(np.nan,np.nan,color='r',linestyle='-',label='35.5GHz')
                    ax.plot(np.nan,np.nan,color='g',linestyle='-',label='95.0GHz')
                    ax.plot(np.nan,np.nan,color='k',linestyle='--',label='McSnow')
                    ax.plot(np.nan,np.nan,color='k',linestyle='-',label='SB')
                    ax.legend()
                    ax.set_ylim([height_bounds[1],height_bounds[0]]) #this might be the range from the initialization height to the maximum of the masses
            else:
                for ax in [axvDoppler,axswidth,axskewn]:         #add legend
                    ax.plot(np.nan,np.nan,color=np.array(['orange','black','red'])[i_cat],linestyle='-',label=np.array(['cloud ice','snow'])[i_cat])

        except Exception:
            print ' \n \n in except: \n \n'
            s = traceback.format_exc()
            serr = "there were errors:\n%s\n" % (s)
            sys.stderr.write(serr)
            print "no spectral data found in:", pam_filestring, ". Is radar_mode set to simple?"
    try:
        if not "axwaterfall" in globals():
            axwaterfall = plt.subplot2grid((num_plots, 1), (4, 0), rowspan=4)
        #plot spectrogram
        axwaterfall = __plotting_functions.plot_waterfall(axwaterfall,pamData,freq=9.6,color='b')
        axwaterfall = __plotting_functions.plot_waterfall(axwaterfall,pamData,freq=35.5,color='r')
        axwaterfall = __plotting_functions.plot_waterfall(axwaterfall,pamData,freq=95.0,color='g')
        #add legend
        axwaterfall.plot(np.nan,np.nan,color='b',linestyle='-',label='9.6GHz')
        axwaterfall.plot(np.nan,np.nan,color='r',linestyle='-',label='35.5GHz')
        axwaterfall.plot(np.nan,np.nan,color='g',linestyle='-',label='95.0GHz')
        axwaterfall.plot(np.nan,0,color='k',linestyle='--',label='McSnow')
        axwaterfall.plot(np.nan,0,color='k',linestyle='-',label='SB')
        axwaterfall.legend()
        axwaterfall.set_ylim([height_bounds[1],height_bounds[0]]) #this might be the range from the initialization height to the maximum of the masses

    except Exception:
        #print "no spectral data found in:", pam_filestring, ". Is radar_mode set to simple or moments?"
        print ' \n \n in except: \n \n'
        s = traceback.format_exc()
        serr = "there were errors:\n%s\n" % (s)
        sys.stderr.write(serr)
except Exception:
    print ' \n \n in except: \n \n'
    s = traceback.format_exc()
    serr = "there were errors:\n%s\n" % (s)
    sys.stderr.write(serr) 


#save figure
plt.tight_layout()
if not os.path.exists('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment): #create direktory if it does not exists
    os.makedirs('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment)
out_filestring = "/spectra_" + testcase + "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'])

#plot spectrogram in new figure

#plot spectrograms
plt.clf() #clear figure



for i_scheme,scheme in enumerate(["SB","MC"]):
    if scheme=="SB":
        #define file to read
        pam_filestring = directory + experiment + "/PAMTRA_2mom_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min.nc'
        #read pamtra output to pamData dictionary
        pamData = __postprocess_PAMTRA.read_pamtra(pam_filestring)
    elif scheme=="MC":
        try:
            if adapt_version==1:
                filename = "/adaptv1_t" + str(tstep) #+ ".nc"
            elif adapt_version==2:
                filename = "/adaptv2_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min'
            elif adapt_version==3:
                filename = "/adaptv3_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min'
                
            pam_filestring = directory + experiment + filename +".nc"
            pamData = __postprocess_PAMTRA.read_pamtra(pam_filestring)

        except:
            print "no spectral data found in:", pam_filestring, ". Is radar_mode set to simple?"
            continue
    #read pamtra output to pamData dictionary
    number_of_spectrograms = 3*2 #3 frequency and two schemes (McSnow and SB)
    for i_freq,freq in enumerate([9.6,35.5,95.0]):
        ax = plt.subplot2grid((number_of_spectrograms, 1), (i_scheme*3+i_freq, 0))
        ax = __plotting_functions.plot_pamtra_spectrogram(ax,pamData,freq=freq)
        ax.set_ylim([height_bounds[1],height_bounds[0]]) #this might be the range from the initialization height to the maximum of the masses
#plt.tight_layout()
out_filestring = "/spectrogram_" + testcase + "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'])