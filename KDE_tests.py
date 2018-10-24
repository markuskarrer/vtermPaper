'''
first test to use (multivariate) KDE for PAMTRA-McSnow adaption
'''
#from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug

#import modules
import numpy as np
import sys #for some debugging
#import pyPamtra
#import re #for selecting numbers from file-string
import subprocess #for using shell commands with subprocess.call()
#import os #import environment variables
from netCDF4 import Dataset
import matplotlib.pyplot as plt

#self written #these files are linked here from the pythoncode/functions directory
import __postprocess_McSnow
#import __general_utilities

#function for formatting list of numerics to a string
def format_list2string(l): #for getting the whole list as a string
    return "["+", ".join(["%.3f" % x for x in l])+"]"


#directory of experiments
directory = "/home/mkarrer/Dokumente/McSnow/MCSNOW/experiments/"
experiment = "1d_xi100000_nz250_lwc0_ncl0_dtc5_nrp200_rm10_rt0_vt1_at0_stick1_dt0_meltt0_multt0_h10-20_ba500"
tstep = 3600
av_tstep = 5
nbins=100
#load file with superparticles (SP)
SP_file = Dataset(directory + experiment + '/mass2fr_' + str(tstep).zfill(4) + 'min_avtstep_' + str(av_tstep) + '.ncdf',mode='r')
#create dictionary for all variables from PAMTRA
SP = dict()

#select space where kde is performed
flag_space = 1 #0-> D; 1-> log(D); 2->D**4 #TODO: this implementation does not work yet


#read PAMTRA variables to pamData dictionary
for var in SP_file.variables:#read files and write it with different names in Data
    SP[var] = np.squeeze(SP_file.variables[var])
#subsample SP to reduce computational time in testing (only specific heights)
heightvec = np.array([2500.,2600])

i=0 #TODO: expand this to all heights (for loop and bigger heightvec)
condition_in_height = np.logical_and(heightvec[i]<=SP["height"],SP["height"]<heightvec[i+1])
SP_in_heightbin = dict()
for var in SP_file.variables:
    SP_in_heightbin[var] = SP[var][condition_in_height] #TODO: expand this to all heights
#constant for bandwidth from Shima 2009
sigma0         = 0.62     #! Shima 2009, Sec. 5.1.4

#varnames = ["m_tot","Frim","height","d_rime","vt","xi",    "rhor","a_rime","mr_crit","diam",    "proj_A",   "mm",         "m_rime",   "m_wat"]
#get necessary parameter of m-D and A-D relationship
mth,unr_alf,unr_bet,rhoi,rhol = __postprocess_McSnow.return_parameter_mD_AD_rel()[0:4]
#determine number of SP
number_ofSP = SP_in_heightbin['m_tot'].shape[0]; print "number_ofSP",number_ofSP

#transform to logspace and set target array
if flag_space==0:
	x =  SP_in_heightbin["diam"]
	x_grid = np.exp(np.arange(-12,-2,nbins)) #np.linspace(1e-6,1e-2,500)
elif flag_space==1:
	x = SP_in_heightbin["diam"] 
	x_grid = np.logspace(-12,-2,nbins) #np.logspace(-8,-2,100) #this must be equidistant (in log space) to calculate x_grid_logdiff!!
	x_log = np.log(x)
	x_grid_log = np.log(x_grid)
	x_grid_logdiff=x_grid_log[1]-x_grid_log[0]
elif flag_space==2:
	x = SP_in_heightbin["diam"]**4
	x_grid = np.arange(-12,-2,nbins)**4


# sigma for kernel estimate, sigma = sigma0/N_s**(1/5), see Shima Sec 5.1.4
sigma0_lin = 1e-8 #np.exp(sigma0)**2
sigmai = (sigma0 / number_ofSP**0.2) #ATTENTION:  this is defined **-1 in McSnow's mo_output.f90 
sigmai_lin = (sigma0_lin / number_ofSP**0.2)

#####
#calculate pdfs with different kde methods
#####
#perform kde with searched bandwidth
guessed_bandwidths = np.linspace(sigmai*0.01,sigmai*2,5) #list with possible bandwidths
sel_bandwidth,pdf = __postprocess_McSnow.kde_bandwidth_estimated(x_log,x_grid_log,guessed_bandwidths);pdf=pdf*x_grid_logdiff
#perform kde with referenc rule bandwidth
pdf_fixed_bandwith = __postprocess_McSnow.kde_statsmodels_m(x_log, x_grid_log, bandwidth=sigmai);pdf_fixed_bandwith=pdf_fixed_bandwith*x_grid_logdiff
#perform kde self-written routine
pdf_fixed_bandwith_self_written = __postprocess_McSnow.kernel_estimate(x,x_grid,sigmai);pdf_fixed_bandwith_self_written=pdf_fixed_bandwith_self_written*x_grid_logdiff
#perform kde self-written routine (take multiplicity into account)
#SP_in_heightbin["xi"][:]=0
pdf_fixed_bandwith_self_written_weighted = __postprocess_McSnow.kernel_estimate(x,x_grid,sigmai,weight=SP_in_heightbin["xi"]);pdf_fixed_bandwith_self_written_weighted=pdf_fixed_bandwith_self_written_weighted*x_grid_logdiff
#perform kde self-written routine in D space
D_min=1e-5; D_max=1e-3
x_grid_rsqeqi = np.linspace(D_min**(1./2.),D_max**(1./2.),nbins*10) #this must be equidistant (in D**2 space) to calculate x_grid_lindiff!!
x_grid_D2diff = x_grid_rsqeqi[1]-x_grid_rsqeqi[0]
x_grid_rsqeqi = x_grid_rsqeqi**2
pdf_fixed_bandwith_self_written_D2KDE = __postprocess_McSnow.kernel_estimate(x,x_grid_rsqeqi,sigmai_lin,space='D2',weight=SP_in_heightbin["xi"]);pdf_fixed_bandwith_self_written_D2KDE=pdf_fixed_bandwith_self_written_D2KDE*x_grid_D2diff
#from IPython.core.debugger import Tracer ; Tracer()()

#get bin width for xgrid
del_x_grid = x_grid[1::]-x_grid[:-1]
if flag_space==0:
    x_grid_lin = x_grid
elif flag_space==1:
    #transform back to lin space
    x_grid_lin = x_grid #np.exp(x_grid)
elif flag_space==2:
    x_grid_lin = x_grid**(1./4.)

flag_plotting = 1 #plot the calculated pdfs
flag_log = 2 #0-> semilogx 1-> semilogy 2-> loglog 3->linlin
if flag_plotting==1:
    #plot for visualization
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    #plotting the pdfs
    strFormat = len(guessed_bandwidths) * '{:3.2f} ' #format for the list
    #plot pdf for automatically selected bandwidth
    #ATTENTION: replace /1.0 with /del_x_grid  to normalize it
    #print x_grid_lin[:-1], pdf[:-1]/1.0, "pdf_searched_bandwidth {:3.3f}ln(m) \n choosen from \n {}".format(sel_bandwidth,format_list2string(guessed_bandwidths))

    #plot the different KDE-results in one plot
    ax.plot(x_grid_lin[:-1], pdf[:-1]/1.0, color='blue', alpha=0.3, lw=3,label="pdf_searched_bandwidth {:3.3f}ln(m) \n choosen from \n {}".format(sel_bandwidth,format_list2string(guessed_bandwidths)))
    #plot pdf with reference rule bandwidth
    ax.plot(x_grid_lin[:-1], pdf_fixed_bandwith[:-1]/1.0,color='red', alpha=1.0, lw=1,linestyle='-',label="pdf_reference_rule_bandwidth {:3.3f}ln(m)".format(sigmai))
    #plot pdf with self-writted kde
    ax.plot(x_grid_lin[:-1], pdf_fixed_bandwith_self_written[:-1]/1.0,color='green', alpha=1.0, lw=1,linestyle='--',label="pdf_reference_rule_bandwidth (sw) {:3.3f}ln(m)".format(sigmai))
    #plot pdf with self-writted kde including weights
    ax.plot(x_grid_lin[:-1], pdf_fixed_bandwith_self_written_weighted[:-1]/1.0,color='orange', alpha=1.0, lw=1,linestyle='--',label="pdf_reference_rule_bandwidth \n (sw-w) {:3.3f}ln(m)".format(sigmai))
    #plot pdf with self-writted kde (in lin space) including weights
    #from IPython.core.debugger import Tracer ; Tracer()()
    ax.plot(x_grid_rsqeqi[:-1], pdf_fixed_bandwith_self_written_D2KDE[:-1]/1.0,color='purple', alpha=1.0, lw=1,linestyle='--',label="pdf_reference_rule_bandwidth \n (sw-w-lin) {:.2e}m**2".format(sigmai_lin))

    #change the scaling strategy if not both axes should be linear
    if flag_log==0:
        ax.set_xscale("log")
    elif flag_log==1:
        ax.set_yscale("log")
    elif flag_log==2:
        ax.set_xscale("log")
        ax.set_yscale("log")
    elif flag_log==3:
        pass
    #define limits
    ymin = 1e-4
    ymax = 1e4
    if flag_log==0:
        print "np.nanmax(pdf)",np.nanmax(pdf)
        ax.set_ylim([0,2*np.nanmax(pdf)])
    if flag_log==1:
        ax.set_xlim([0,x_grid_lin[max(np.argwhere(pdf> ymin)[-1],np.argwhere(pdf_fixed_bandwith> ymin)[-1])]]) #adaptively select xlim depending on ymin
        ax.set_ylim([ymin,ymax])
    elif flag_log==2:
        ax.set_xlim([10**-6,10**-2]) #x_grid_lin[max(np.argwhere(pdf> ymin)[-1],np.argwhere(pdf_fixed_bandwith> ymin)[-1])]]) #adaptively select xlim depending on ymin    
        ax.set_ylim([ymin,ymax])
    elif flag_log==3:
        ax.set_xlim([0,x_grid_lin[max(np.argwhere(pdf> ymin)[-1],np.argwhere(pdf_fixed_bandwith> ymin)[-1])]]) #adaptively select xlim depending on ymin
        ax.set_ylim([0,2*np.nanmax(pdf)])

    #get and rescale xticklabels
    xticklabels = ax.get_xticks()
    ax.set_xticklabels(xticklabels*1000) #here the xticklabels are converted from m to mm
    #label
    ax.set_xlabel("diameter D / mm")
    ax.set_ylabel("f(ln D)")
    plt.legend()
    #save figure
    plt.tight_layout()
    plt.savefig('/home/mkarrer/Dokumente/plots/kde_test/kde_test_flaglog' + str(flag_log) + '.pdf', dpi=400)
    subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/kde_test/kde_test_flaglog' + str(flag_log) + '.pdf'])
