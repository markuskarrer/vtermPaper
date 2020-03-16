# coding: utf-8
#import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.pylab as pylab
from matplotlib.ticker import FormatStrFormatter
import os
import subprocess
import random
import csv
#import other self-defined functions
import __postprocess_McSnow
import __postprocess_SB
import __fallspeed_relations
import __tools_for_processing_Jagg
import __plotting_functions
import __setup_mDAD_frommodels
import generate_2Dhist_of_N_D_Nmono_from_MC_and_Jagg 
from IPython.core.debugger import Tracer ; debug=Tracer()
from matplotlib import rc


'''
this code reads in properties from Jussis aggregate model
and fits different functions to each property (m-D(monomer dependent and binary),A-D(monomer dependent and binary),v-D for all fallspeed models) and displays them together
'''


def fit_powerlaw(xdata,ydata):
    from scipy import optimize
    ##########
    # Fitting the data -- Least Squares Method
    ##########

    # Power-law fitting is best done by first converting
    # to a linear equation and then fitting to a straight line.
    # Note that the `logyerr` term here is ignoring a constant prefactor.
    #
    #  y = a * x^b
    #  log(y) = log(a) + b*log(x)
    #
    '''
    ##########
    # Generate data points with noise
    ##########
    num_points = 20

    # Note: all positive, non-zero data
    xdata = np.linspace(1.1, 10.1, num_points)
    powerlaw = lambda x, amp, index: amp * (x**index)
    ydata = powerlaw(xdata, 10.0, -2.0)     # simulated perfect data
    yerr = 0.01 * ydata                      # simulated errors (10%)

    ydata += np.random.randn(num_points) * yerr       # simulated noisy data
    '''


    logx = np.log10(xdata)
    logy = np.log10(ydata)
    logyerr = 1.*np.ones_like(ydata)

    # define our (line) fitting function
    fitfunc = lambda p, x: p[0] + p[1] * x
    errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

    pinit = [0.9, 1.0]
    out = optimize.leastsq(errfunc, pinit,
                           args=(logx, logy, logyerr), full_output=1)

    pfinal = out[0]
    covar = out[1]
    print "pfinal",pfinal
    print "covar",covar

    index = pfinal[1]
    amp = 10.0**pfinal[0]

    indexErr = np.sqrt( covar[1][1] )
    ampErr = np.sqrt( covar[0][0] ) * amp
    print "amp,index",amp,index 
    return amp,index

tumbling=False #so far only the use of non-tumbling projected area is implemented here
show_Nmono_fits = True #show fit lines for some selected monomer numbers
v_small_notuse = 0.0  #ATTENTION: particle smaller than .5ms-1 are not used #ATTENTION: applied to B?hm only


def read_and_plot(fig,axes,particle_types,called_by_main=False,plot_diagnostic=""):
    '''
    this performs all the processing and plotting
    
    ARGUMENTS:
    fig: figure handle
    axes: axes handles
    particle type: monomer type whcih is processed and plotted
    called_by_main: if this script has been executed directly save the figures here
    plot_diagnostic: what should be plotted?
    '''
    #flatten axes in case we have more than one column
    try:
        axes = axes.flat
    except:
        print axes," cannot be flattened any more"

    for i_particle_type,particle_type in enumerate(particle_types):    



        
        #read the properties of the individual particles from the files
        particle_dic,N_mono_list = __tools_for_processing_Jagg.read_particle_prop_files(
        #prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation/",
        #prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation_local/samesizeparam/",
        prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation_local/diffsizeparam/",
                                                                                D_small_notuse=1e-4, #ATTENTION: particle smaller than D_small_notuse are not considered (e.g. because of resolution issues)
                                                                                N_small_notuse=1, #ATTENTION: THIS HAS NO EFFECT!!  monomer numbers smaller than N_small_notuse are not considered
                                                                                #grid_res_array = [5e-6,10e-6], #array of grid resolutions (defines the files which are read in)
                                                                                grid_res_array = [40e-6], #array of grid resolutions (defines the files which are read in)

                                                                                particle_type = particle_type, #define the habit
                                                                                test_with_small_sample = False
                                                                                )
        #show the dictionary (e.g. to check that it's not empty because the path was wrong)
        print particle_dic;
        #particle_dic["diam"] = particle_dic["Dmax_xy_aligned"]
        print particle_dic.keys()
        #particle_dic["diam"] = particle_dic["rad_gyr"]*2.0
        
        if plot_diagnostic=="RG_vs_Dmax":
            number_of_plots = 2 #8
            #optimize the appearance of the plot (figure size, fonts)
            [fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots,legend_fontsize='small')
            axes.set_xscale("log")
            axes.set_yscale("log")
            axes.scatter(particle_dic["diam"],particle_dic["rad_gyr"])
            diam = np.linspace(min(particle_dic["diam"]),max(particle_dic["diam"]))
            a,b = fit_powerlaw(particle_dic["diam"],particle_dic["rad_gyr"])
            axes.plot(diam,a*diam**b,color='k',linestyle='--')
            axes.text(
                0.1, 0.9,'{:.3f}'.format(a) + "$D_{max}$^" + '{:.3f}'.format(b),
                horizontalalignment='center',
verticalalignment='center',
                transform = axes.transAxes)
            axes.set_xlabel(r"$D_{max}$ [m]")
            axes.set_ylabel(r"$RG$ [m]")
        if plot_diagnostic=="Dmax_Dcirc":
            number_of_plots = 2 #8
            #optimize the appearance of the plot (figure size, fonts)
            [fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots,legend_fontsize='small')

            diam = np.linspace(min(particle_dic["diam"]),max(particle_dic["diam"]))
            #fit
            a,b = fit_powerlaw(particle_dic["diam"],particle_dic["D_circum_proj_area"])
            
            axes[0].set_xscale("log")
            axes[0].set_yscale("log")
            axes[0].scatter(particle_dic["diam"],particle_dic["D_circum_proj_area"],s=1)
            axes[1].set_xlabel(r"$D_{max}$ [m]")
            axes[1].set_ylabel(r"$D_{proj}/D_{max}$")
            axes[0].plot(diam,a*diam**b,color='k',linestyle='--')
            axes[0].plot(diam,diam,color='r',linestyle='--')
            axes[0].text(
                0.1, 0.9,'{:.3f}'.format(a) + "$D_{max}$^" + '{:.3f}'.format(b),
                horizontalalignment='center',
                verticalalignment='center',
                transform = axes[1].transAxes)
            axes[0].set_xlabel(r"$D_{max}$ [m]")
            axes[0].set_ylabel(r"$D_{proj}$")


            #relative
            axes[1].set_xscale("log")
            axes[1].set_yscale("linear")
            axes[1].scatter(particle_dic["diam"],(particle_dic["D_circum_proj_area"])/particle_dic["diam"],s=1)
            axes[1].set_xlabel(r"$D_{max}$ [m]")
            axes[1].set_ylabel(r"$D_{proj}/D_{max}$")
        #save figure
        plt.tight_layout()
        out_filestring = "otherdiag_" + particle_type            
        dir_save = "/home/mkarrer/Dokumente/plots/Jagg/"
        plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
        plt.savefig(dir_save + out_filestring + '.png', dpi=100)
        print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
        subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])


if __name__ == '__main__':
    #from sys import argv 
   
    particle_types=["needle"] #,"dendrite","column","needle","rosette","mixcoldend1","mixcolumndend"]

    #optimize the appearance of the plot (figure size, fonts)
    fig=0 #will be overwritten in read_and_plot
    axes=0 #will be overwritten in read_and_plot
    
    #MAIN PART
    read_and_plot(fig,axes,particle_types,called_by_main=True,plot_diagnostic="RG_vs_Dmax")
    #read_and_plot(fig,axes,particle_types,called_by_main=True,plot_diagnostic="Dmax_Dcirc")
    
           
