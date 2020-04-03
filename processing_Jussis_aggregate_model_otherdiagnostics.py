# coding: utf-8
#import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.pylab as pylab
from matplotlib.ticker import FormatStrFormatter
import os
import sys
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


def read_and_plot(fig,axes,particle_types,called_by_main=False,plot_diagnostic="",N_small_notuse=1):
    '''
    this performs all the processing and plotting
    
    ARGUMENTS:
    fig: figure handle
    axes: axes handles
    particle type: monomer type whcih is processed and plotted
    called_by_main: if this script has been executed directly save the figures here
    plot_diagnostic: what should be plotted?
    N_small_notuse: monomer numbers smaller than N_small_notuse are not considered
    '''
    #flatten axes in case we have more than one column
    try:
        axes = axes.flat
    except:
        print axes," cannot be flattened any more"
    if "monoonly" in plot_diagnostic:
        monoonly = True
    else:
        monoonly = False
        
    for i_particle_type,particle_type in enumerate(particle_types):    



        
        #read the properties of the individual particles from the files
        particle_dic,N_mono_list = __tools_for_processing_Jagg.read_particle_prop_files(
        #prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation/",
        #prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation_local/",
        prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation_local_rot/",
        #prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation_local/samesizeparam/",
        #prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation_local/diffsizeparam/",
                                                                                D_small_notuse=1e-4, #ATTENTION: particle smaller than D_small_notuse are not considered (e.g. because of resolution issues)
                                                                                N_small_notuse=N_small_notuse, #ATTENTION: THIS HAS NO EFFECT!!  monomer numbers smaller than N_small_notuse are not considered
                                                                                #grid_res_array = [5e-6,10e-6], #array of grid resolutions (defines the files which are read in)
                                                                                grid_res_array = [10e-6], #array of grid resolutions (defines the files which are read in)

                                                                                particle_type = particle_type, #define the habit
                                                                                test_with_small_sample = False,
         use_only_monomers = monoonly, #get only Nmono=1 particle properties back
                                                                                )
        #show the dictionary (e.g. to check that it's not empty because the path was wrong)
        #print particle_dic;
        #particle_dic["diam"] = particle_dic["Dmax_xy_aligned"]
        #print particle_dic.keys()
        #particle_dic["diam"] = particle_dic["rad_gyr"]*2.0
        
        if plot_diagnostic=="RG_vs_Dmax":
            number_of_plots = 1 #8
            #optimize the appearance of the plot (figure size, fonts)
            [fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots,legend_fontsize='small')
            axes.set_xscale("log")
            axes.set_yscale("log")
            axes.scatter(particle_dic["diam"],particle_dic["rad_gyr"])
            
            #set up diameter-array for analyzing the fit
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

            #set up diameter-array for analyzing the fit
            diam = np.linspace(min(particle_dic["diam"]),max(particle_dic["diam"]))
            #fit
            a,b = fit_powerlaw(particle_dic["diam"],particle_dic["D_circum_proj_area"])
            #print "D_circum_proj_Area=",a,"*diam**",b
            
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
            axes[0].set_ylabel(r"$D_{proj}$ [m]")


            #relative
            axes[1].set_xscale("log")
            axes[1].set_yscale("linear")
            axes[1].scatter(particle_dic["diam"],(particle_dic["D_circum_proj_area"])/particle_dic["diam"],s=1)
            axes[1].set_xlabel(r"$D_{max}$ [m]")
            axes[1].set_ylabel(r"$D_{proj}/D_{max}$")
        if "projarea_side" in plot_diagnostic:

            number_of_plots = 1
            #optimize the appearance of the plot (figure size, fonts)
            if "comp2obs" in plot_diagnostic:
                [fig,axis] = __plotting_functions.proper_font_and_fig_size(number_of_plots,aspect_ratio=1./8.,legend_fontsize='small')
            else:
                [fig,axis] = __plotting_functions.proper_font_and_fig_size(number_of_plots,legend_fontsize='small')

            #show the data
            axis.set_xscale("log")
            #axis.set_yscale("log")
            i=-3
            
            axis.scatter(particle_dic["diam"],particle_dic["area"]/(np.pi/4.*particle_dic["diam"]**2),s=1,label="(x-y)")
            if not ("onlytop" in plot_diagnostic) and not("siderandom" in plot_diagnostic):
                axis.scatter(particle_dic["Dmax_yz_aligned"],particle_dic["area_aligned_side_yz"]/(np.pi/4.*particle_dic["Dmax_yz_aligned"]**2),s=1,label="(y-z)")
                axis.scatter(particle_dic["Dmax_xz_aligned"],particle_dic["area_aligned_side_xz"]/(np.pi/4.*particle_dic["Dmax_xz_aligned"]**2),s=1,label="(x-z)")
            elif "siderandom" in plot_diagnostic:
                axis.scatter(particle_dic["Dmax_random_aligned"],particle_dic["area_aligned_side_random"]/(np.pi/4.*particle_dic["Dmax_random_aligned"]**2),s=1,label="random side")

            #set up diameter-array for analyzing the fit
            diam = np.linspace(min(particle_dic["diam"]),max(particle_dic["diam"]),100)
            diam_all = np.linspace(1e-4,4e-2,100)
            diam1 = np.linspace(min(particle_dic["Dmax_yz_aligned"]),max(particle_dic["Dmax_yz_aligned"]),100)
            diam2 = np.linspace(min(particle_dic["Dmax_xz_aligned"]),max(particle_dic["Dmax_xz_aligned"]),100)
            diam_rand = np.linspace(min(particle_dic["Dmax_random_aligned"]),max(particle_dic["Dmax_random_aligned"]),100)

            #fit the area ratio (x-y and y-z)
            a,b = fit_powerlaw(particle_dic["diam"],particle_dic["area"]/(np.pi/4.*particle_dic["diam"]**2))
            a_side1,b_side1 = fit_powerlaw(particle_dic["Dmax_yz_aligned"],particle_dic["area_aligned_side_yz"]/(np.pi/4.*particle_dic["Dmax_yz_aligned"]**2))
            a_side2,b_side2 = fit_powerlaw(particle_dic["Dmax_xz_aligned"],particle_dic["area_aligned_side_xz"]/(np.pi/4.*particle_dic["Dmax_xz_aligned"]**2))
            a_side_rand,b_side_rand = fit_powerlaw(particle_dic["Dmax_random_aligned"],particle_dic["area_aligned_side_random"]/(np.pi/4.*particle_dic["Dmax_random_aligned"]**2))

            #display the coefficients
            #print "aratio(x-y)",a,"*D(3D)**",b
            #print "aratio(y-z)",a_side1,"*D(y-z)**",b_side1
            #print "aratio(x-z)",a_side2,"*D(x-z)**",b_side2
            
            #show the fit lines
            axis.semilogx(diam,a*diam**b,label="fit (x-y) a={0:4.3f},b={1:4.2f}".format(a,b))
            if not ("onlytop" in plot_diagnostic) and not ("siderandom" in plot_diagnostic):
                axis.semilogx(diam1,a_side1*diam1**b_side1,label="fit (y-z) a={0:4.3f},b={1:4.2f}".format(a_side1,b_side1),linestyle='--')
                axis.semilogx(diam2,a_side2*diam2**b_side2,label="fit (x-z) a={0:4.3f},b={1:4.2f}".format(a_side2,b_side2),linestyle='--')
            elif "siderandom" in plot_diagnostic:
                diam_median_rand = np.logspace(-4,-2,20)
                arat_median = np.ones_like(diam_median_rand)*np.nan
                for i,diam_now in enumerate(diam_median_rand[:-1]):
                    arat_now = particle_dic["area_aligned_side_random"]/(np.pi/4.*particle_dic["Dmax_random_aligned"]**2) 
                    arat_now = arat_now[np.logical_and(particle_dic["Dmax_random_aligned"]>diam_median_rand[i],particle_dic["Dmax_random_aligned"]<diam_median_rand[i+1])]
                    arat_median[i] = np.median(arat_now)
                axis.semilogx(diam_median_rand,arat_median,label="median (side random)",linestyle='-.',color="black")
                axis.semilogx(diam_rand,a_side_rand*diam_rand**b_side_rand,label="fit (side random) a={0:4.3f},b={1:4.2f}".format(a_side_rand,b_side_rand),linestyle='--')
            
            #overlay the area ratio derived from the A-D fit
            a_A = 0.066; b_A = 1.79
            a_ratio_from_paper_fit = a_A*diam_all**b_A/(np.pi/4.*diam_all**2)
            if not ("monoonly" in plot_diagnostic) and ("onlytop" in plot_diagnostic):
                axis.semilogx(diam_all,a_ratio_from_paper_fit,linestyle="--",label="from A-D fit",color="r")
            
            if ("monoonly" in plot_diagnostic):
                axis.text(0.5,1.0, "Nmono=" + str(N_small_notuse),
                         horizontalalignment='right',
                         verticalalignment='top',
                         transform = axis.transAxes,
                         fontsize=18)
            else:
                axis.text(0.5,1.0, "Nmono>=" + str(N_small_notuse),
                         horizontalalignment='right',
                         verticalalignment='top',
                         transform = axis.transAxes,
                         fontsize=18)
            #legend
            axis.legend()

            #labels
            axis.set_xlabel(r"$D_{max}$ [m]")
            axis.set_ylabel(r"area ratio")
            
            #limits
            if "monoonly" in plot_diagnostic:
                axis.set_xlim([1e-4,1e-3])
            if "comp2obs" in plot_diagnostic:
                if "comp2obsKor" in plot_diagnostic:
                    axis.set_xscale("linear")
                    axis.set_xlim([0,0.001])
                else:
                    axis.set_xlim([1e-4,1e-2])
            else:
                axis.set_xlim([1e-4,4e-2])
            axis.set_ylim([0.0,1.0])
            
            #grid
            axis.grid(b=True,which="both")
        
        if plot_diagnostic=="filled_vs_unfilled_area":
            number_of_plots = 1 #8
            #optimize the appearance of the plot (figure size, fonts)
            [fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots,legend_fontsize='small')
            axes.set_xscale("log")
            axes.set_yscale("linear")

            ratio = particle_dic["area_aligned_filled"]/particle_dic["area"]
            ratio[ratio<0] = 0.9
            axes.scatter(particle_dic["diam"],ratio)
            
            axes.set_xlabel(r"$D_{max}$ [m]")
            axes.set_ylabel(r"area(filled)/area(unfilled)")

        #save figure
        plt.tight_layout()
        if "monoonly" in plot_diagnostic:
            out_filestring = "otherdiag_" + particle_type + plot_diagnostic + "_Nmonoeq" + str(N_small_notuse)
        elif "projarea_side" in plot_diagnostic:
            out_filestring = "otherdiag_" + particle_type + plot_diagnostic + "_Nmonogt" + str(N_small_notuse)
        else:
            out_filestring = "otherdiag_" + particle_type + plot_diagnostic
            
        dir_save = "/home/mkarrer/Dokumente/plots/Jagg/"
        plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
        plt.savefig(dir_save + out_filestring + '.png', dpi=100)
        print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
        if "monoonly" in plot_diagnostic:
            sys.exit()
        #subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])


if __name__ == '__main__':
   
    particle_types=["mixcolumndend"] #,"dendrite","column","needle","rosette","mixcoldend1","mixcolumndend"]

    #optimize the appearance of the plot (figure size, fonts)
    fig=0 #will be overwritten in read_and_plot
    axes=0 #will be overwritten in read_and_plot
    
    #MAIN PART
    #read_and_plot(fig,axes,particle_types,called_by_main=True,plot_diagnostic="RG_vs_Dmax")
    #read_and_plot(fig,axes,particle_types,called_by_main=True,plot_diagnostic="Dmax_Dcirc")
    for N in [1,2,10,20,50,100]:
        #read_and_plot(fig,axes,particle_types,called_by_main=True,plot_diagnostic="projarea_side_onlytop",N_small_notuse=N)
        #read_and_plot(fig,axes,particle_types,called_by_main=True,plot_diagnostic="projarea_side_monoonly",N_small_notuse=N)
        #read_and_plot(fig,axes,particle_types,called_by_main=True,plot_diagnostic="projarea_side",N_small_notuse=N)
        #read_and_plot(fig,axes,particle_types,called_by_main=True,plot_diagnostic="projarea_siderandom_comp2obs",N_small_notuse=N)
        pass #read_and_plot(fig,axes,particle_types,called_by_main=True,plot_diagnostic="projarea_siderandom_comp2obsKor",N_small_notuse=N)
     
    read_and_plot(fig,axes,particle_types,called_by_main=True,plot_diagnostic="filled_vs_unfilled_area")
