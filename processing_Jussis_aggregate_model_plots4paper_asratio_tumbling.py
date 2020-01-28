# coding: utf-8
#import packages
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.pylab as pylab
import csv #to read the txt files
import os
import glob #to get filename to list
import subprocess
import itertools#to read only certain lines of the txt files
import re
from matplotlib.colors import LogNorm
import matplotlib.pylab as pylab
a=0
import argparse; parser = argparse.ArgumentParser() #get optional arguments from shell
#import other self-defined functions
import __postprocess_McSnow
import __postprocess_SB
import __fallspeed_relations
import __tools_for_processing_Jagg
import __plotting_functions
import __setup_mDAD_frommodels
import generate_2Dhist_of_N_D_Nmono_from_MC_and_Jagg
from matplotlib import rc

from IPython.core.debugger import Tracer ; debug=Tracer()
'''
this code reads in properties from Jussis aggregate model, calculate + plots the following fallspeed
and fits different functions to each property (aspect ratio, fall speed for different tumbling settings)
'''

#from IPython.core.debugger import Tracer ; Tracer()()
#switch on/of fitting of terminal velocity (removes the fit lines)
fitvel=1 #0-> no fitting of velocity at all; 1-> fit only monomer and all aggregates; 2-> fit everything (including seperated for different monomer numbers)
tumbling=False #True: take the rotated projected area; False: take the area calculated from the aligned particle
take_mono_prop_theor=True #take theoretical values (as assumed in Jagg) instead of fits to N_mono=1

#define where the txt files with properties of the aggregates are stored
#prop_file_folder = "/data/optimice/Jussis_aggregates_bugfixedrotation/tumbling_asratio/"
prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation/"

#define which particles (with which resolution to read in)
grid_res_array = [5e-6,10e-6] #resolution for big medium and smaller particles in mum [10e-6,5e-5]


#calculate the arrays with masses and areas corresponding to the common fit_dic["diam"] and based on the (piecewise) power-law fits
fit_dic = dict() #initialize fit dictionary
#set up array of diameters for bin edges, the KDE-grid and the fitting
low_diam_log=-4; 
#highest diameter for fit
high_diam_log=-1.3979400086720375 #np.log10(3e-2)=-1.5228787452803376 #np.log10(4e-2)=-1.3979400086720375

diam = np.logspace(low_diam_log,high_diam_log,50) #set up array of diameters (for bin edges and KDE-grid)
fit_dic["diam"] = diam


#loop over different particle habits
particle_types = ["plate"] #"needle","dendrite","plate","mixdendneedle","column"] #,"column","needle","dendrite","rosette","bullet"] #"needle","dendrite","plate","column","rosette","bullet","plate_dendrite",...]] # ,"bullet"rosette
for particle_type_comb in particle_types:

    print "########################"
    print "#####" + particle_type_comb + "########"
    print "########################"
    particle_dic = __tools_for_processing_Jagg.init_particle_dict() #initialize dictionary which contains information about the type of the pristine crystals (dentrites, needles, plates, ...) the underlying distribution of the
    #allow also combination of habits
    if particle_type_comb not in ("needle","dendrite","plate","column","rosette","bullet"):
        #particle_type_comb should be a combination of particles in particle type (e.g. "plate_dendrite")
        particle_type_list = re.split(r'_',particle_type_comb) #split the combination to read in all habits
    else: #just one habit is analyzed
        particle_type_list = [particle_type_comb]
    #########
    #START: read the properties of the aggregates into the particle_dic dictionary
    #########

    for particle_type in particle_type_list: #this is a single particle type (no combination allowed after that)

            
            
        #read the properties of the individual particles from the files
        particle_dic,N_mono_list = __tools_for_processing_Jagg.read_particle_prop_files(prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation/",
                                                                                D_small_notuse=1e-4, #ATTENTION: particle smaller than D_small_notuse are not considered (e.g. because of resolution issues)
                                                                                N_small_notuse=1, #ATTENTION:  monomer numbers smaller than N_small_notuse are not considered
                                                                                grid_res_array = [5e-6,10e-6], #[1e-6,5e-6,10e-6], #array of grid resolutions (defines the files which are read in)
                                                                                particle_type = particle_type, #define the habit
                                                                                test_with_small_sample = False
                                                                                )
    #show the dictionary (e.g. to check that it's not empty because the path was wrong)
    print particle_dic

    #########
    #END: read the properties of the aggregates into the particle_dic dictionary
    #########     

   
    ######
    ###plotting
    ###
    #initialize the plot
    ###
    number_of_plots = 2 #28
    
    #optimize the appearance of the plot (figure size, fonts)
    [fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots,aspect_ratio=0.25)
    params = {'legend.fontsize': 'medium',
                'axes.labelsize': 'medium', #size: Either a relative value of 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large' or an absolute font size, e.g., 12
        'axes.titlesize':'medium',
        'xtick.labelsize':'medium',
        'ytick.labelsize':'medium'}
    pylab.rcParams.update(params)
    i_ax=0
    #define size of colorbar labels
    colbarlabelsize=11
    # define the colormap and the discrete boundaries (defined by the considered monomer numbers) and set up an array of the same colors (for line-plotting)
    cmap = plt.cm.hot #brg
    #modify lowest color in colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmaplist = cmaplist[30:256-30] #crop colormap to avoid too dark/light scatterer #out of 256 colors

    #colors for the binary relevant
    colormonomer = colors.ColorConverter.to_rgb("b") #(1.0,0.0,0.0,1.0) #make Nmono=1 blue
    colorallagg =  colors.ColorConverter.to_rgb("g") #(1.0,0.0,0.0,1.0) #make Nmono=1 blue
    
    #reflect the selected color also in the colormap
    cmaplist[0] = colormonomer
    cmap = colors.LinearSegmentedColormap.from_list('mcm',cmaplist, cmap.N)
    bounds = np.append(N_mono_list,[9999])
    norm = colors.BoundaryNorm(bounds, cmap.N)
    usecolors = cmap(np.linspace(0,1,N_mono_list.shape[0])) #get colors from colormap to make consistent line colors for the monomer number
    
    cmaplist[0] = (0.0,0.0,1.0,1.0)
    cmap = colors.LinearSegmentedColormap.from_list('mcm',cmaplist, cmap.N)

    bounds = np.append(N_mono_list,[9999])
    norm = colors.BoundaryNorm(bounds, cmap.N)
    usecolors = cmap(np.linspace(0,1,N_mono_list.shape[0])) #get colors from colormap to make consistent line colors for the KDE-lines #pylab.cm.Greens(np.linspace(0,1,N_mono_list.shape[0]))
    
    ##define gray colormap
    # define the colormap and the discrete boundaries (defined by the considered monomer numbers) and set up an array of the same colors (for line-plotting)
    cmap2 = plt.cm.gray #brg
    #modify lowest color in colormap
    cmaplist2 = [cmap2(i) for i in range(cmap2.N)]
    cmaplist2 = cmaplist2[60:256-30] #crop colormap to avoid too dark/light scatterer #out of 256 colors

    #colors for the binary relevant
    #colormonomer = colors.ColorConverter.to_rgb("b") #(1.0,0.0,0.0,1.0) #make Nmono=1 blue
    #colorallagg =  colors.ColorConverter.to_rgb("g") #(1.0,0.0,0.0,1.0) #make Nmono=1 blue
    
    #span the colormaps
    cmap2 = colors.LinearSegmentedColormap.from_list('mcm',cmaplist2, cmap2.N)
    bounds2 = np.append(N_mono_list,[9999])
    norm2 = colors.BoundaryNorm(bounds2, cmap2.N)
    #####
    #START: plots in m-D/A-D space
    #####

    #initialize an array which contains all fit results
    fitresult_arr_dir_powerlaw_string = np.zeros((5,N_mono_list.shape[0],2)) #fit results (saved as floats) #dimensions: 0-> [mass,array]; 1-> monomer number; 2-> parameters (a,b) in fit
    
    #fit the m-D and A-D properties
    for i_prop,(prop,prop_unit,prop_short) in enumerate(zip(["mass","area","area_partal10std","area_partal20std","area_partal40std","area_partal60std"],["kg","m2","m2","m2","m2","m2"],["m","A","A","A","A","A"])): #,"HWJussi","KCJussi"]):
        print "fitting and plotting: ",prop
        
        plot_allfits_bool=False #plot fits and display fit coefficients in text for all monomer number
        plot_binary_fits= False #plot fits and display fit coefficients in text for monomer and allagg
        
        #fit m-D and A-D relations
        axes[i_ax],particle_dic,fitline_allagg = __tools_for_processing_Jagg.fitting_wrapper(axes[i_ax],particle_dic,prop,N_mono_list,usecolors,fit_dic,diam,function="powerlaw",weight="None",linewidthinvalid=0.9,linewidthvalid=1.7,plot_fits=plot_allfits_bool,plot_binary_fits=plot_binary_fits)
        
        if take_mono_prop_theor: #overwrite fit for monomers with theoretical value
            if "mass" in prop:
                fit_dic[prop + '_coeff_Nmono_1'] = __tools_for_processing_Jagg.calc_mD_AD_coeffs(particle_type)[0:2]
            if "area"==prop: #dont overwrite non-aligned areas
                fit_dic[prop + '_coeff_Nmono_1'] = __tools_for_processing_Jagg.calc_mD_AD_coeffs(particle_type)[2:4]
            fit_dic[prop + "_Nmono_1"] = fit_dic[prop + "_coeff_Nmono_1"][0]*fit_dic["diam"]**fit_dic[prop + "_coeff_Nmono_1"][1] #recalculate arrays of masses and areas 
            
    for i_prop,prop in enumerate(["vterm_B92","vterm_HW10","vterm_KC05","vterm_mitch_heym","vterm_bohm89"]):

        #calculate vterm based on the m-D/A-D fits
        for partal_add in ["","_partal10std","_partal20std","_partal40std","_partal60std"]:
            #calculate the terminal velocity based on the m-D and A-D fits and plot this line # use the colormap and the bounds which where defined before
            if prop=="vterm_B92":
                velocity_model = "bohm"
            else:
                velocity_model = prop[6:] #get the name of the terminal velocity model from the property name (prop) which starts with "vterm"
            #change the scale (according to xscale_vec and yscale_vec)
            axes[i_ax].set_xscale("log") #axes[i_ax].set_xscale(xscale_vec[i_prop]) 
            #axes[i_ax].set_yscale("linear") #axes[i_ax].set_yscale(yscale_vec[i_prop])
            for i,N_mono_now in enumerate(["1","allagg"]): #N_mono_list):
                fit_dic[prop + partal_add + "_fitted_via_mD_AD_Nmono" + str(N_mono_now)] = __fallspeed_relations.calc_vterm(velocity_model,fit_dic["mass" + "_Nmono_" + str(N_mono_now)],fit_dic["diam"],fit_dic["area" + partal_add + "_Nmono_" + str(N_mono_now)])

    ###
    #START: plot v-D lines according to m/A-D fits for different terminal velocity models and tumbling settings
    ###
    #i_ax+=1
    i_ax_mono=i_ax
    i_ax_agg=i_ax+1
    for i_align,partal_add in enumerate(["","_partal20std","_partal40std"]): #,"_partal60std"]): #,"_partal10std","_partal20std","_partal40std","_partal60std"]):
        for i_prop,prop in enumerate(["vterm_B92","vterm_HW10","vterm_KC05"]): #,"vterm_bohm89"]): #,"vterm_mitch_heym"]): #"vterm_B92","vterm_HW10","vterm_KC05"]): #,"vterm_mitch_heym"]): #,"HWJussi","KCJussi"]):
            velocity_model = prop[6:] #get the name of the terminal velocity model from the property name (prop) which starts with "vterm"
            
            #lines for monomers
            axes[i_ax_mono].semilogx(fit_dic["diam"],fit_dic["vterm_" + velocity_model + partal_add + "_fitted_via_mD_AD_Nmono1"],marker='None',color=np.array(["k","r","y","b"])[i_prop],label=np.array([velocity_model,"__None","__None","__None","__None"])[i_align],linestyle=np.array(["-","--","-.",":"])[i_align])
            
            #all agregates
            axes[i_ax_agg].semilogx(fit_dic["diam"],fit_dic["vterm_" + velocity_model + partal_add + "_fitted_via_mD_AD_Nmonoallagg"],marker='None',color=np.array(["k","r","y","b"])[i_prop],label=np.array([velocity_model,"__None","__None","__None","__None"])[i_align],linestyle=np.array(["-","--","-.",":"])[i_align])
            
        for i_ax_now in [i_ax_mono,i_ax_agg]:
            #add a legend for the tumbling
            if i_align==0:  #no tumbling
                axes[i_ax_now].plot(np.nan,np.nan,c='k',linestyle=np.array(["-","--","-.",":"])[i_align],label="no tumbling")
            else:
                axes[i_ax_now].plot(np.nan,np.nan,c='k',linestyle=np.array(["-","--","-.",":"])[i_align],label="tumbling std: "+ partal_add[7:9] + "$^\circ$")
    show_rel_diff=True
    if show_rel_diff:
        for Nmonotype in ["Nmono1","Nmonoallagg"]:
            print "#######compare tumbling"
            for i_prop,prop in enumerate(["vterm_B92","vterm_HW10","vterm_KC05"]):
                partal_add = "_partal40std"
                print "\n"
                print Nmonotype,prop,partal_add,
                print "absol vdiff",np.max(fit_dic[prop + partal_add + "_fitted_via_mD_AD_"+ Nmonotype]-fit_dic[prop + "_fitted_via_mD_AD_" + Nmonotype])
                print "rel vdiff",np.max((fit_dic[prop + partal_add + "_fitted_via_mD_AD_"+ Nmonotype]-fit_dic[prop + "_fitted_via_mD_AD_" + Nmonotype])/fit_dic[prop + "_fitted_via_mD_AD_" + Nmonotype])
                print "\n"
            print "#######compare hydro models"
            for i_prop,prop in enumerate(["vterm_HW10","vterm_KC05"]):
                print "\n"
                print Nmonotype,prop,partal_add + "-vterm_B92"
                dconsidered=1e-2 #diameter bigger than this are not considered
                print "absol vdiff",np.min(fit_dic[prop + "_fitted_via_mD_AD_"+ Nmonotype][fit_dic["diam"]<dconsidered]-fit_dic["vterm_B92_fitted_via_mD_AD_" + Nmonotype][fit_dic["diam"]<dconsidered]),np.max(fit_dic[prop + "_fitted_via_mD_AD_"+ Nmonotype][fit_dic["diam"]<dconsidered]-fit_dic["vterm_B92_fitted_via_mD_AD_" + Nmonotype][fit_dic["diam"]<dconsidered]) 
                print "rel vdiff",np.min((fit_dic[prop + "_fitted_via_mD_AD_"+ Nmonotype][fit_dic["diam"]<dconsidered]-fit_dic["vterm_B92_fitted_via_mD_AD_" + Nmonotype][fit_dic["diam"]<dconsidered])/fit_dic["vterm_B92_fitted_via_mD_AD_" + Nmonotype][fit_dic["diam"]<dconsidered]),np.max((fit_dic[prop + "_fitted_via_mD_AD_"+ Nmonotype][fit_dic["diam"]<dconsidered]-fit_dic["vterm_B92_fitted_via_mD_AD_" + Nmonotype][fit_dic["diam"]<dconsidered])/fit_dic["vterm_B92_fitted_via_mD_AD_" + Nmonotype][fit_dic["diam"]<dconsidered]) 
        #debug()
        #raw_input("seen it? then press a key")

    for i_ax_now in [i_ax_mono,i_ax_agg]:
        
        
        #make labels
        axes[i_ax_now].set_xlabel("$D_{max}$ [m]")
        axes[i_ax_now].set_ylabel("$v_{term}$ [m/s]" ) #TODO: plot also the untis of these properties

        #change the axis
        axes[i_ax_now].set_xlim([10**low_diam_log,10**high_diam_log]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
        axes[i_ax_now].set_ylim([0,2.5]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
        axes[i_ax_now].grid(which="major")
    #show legend
    axes[i_ax_mono].legend(loc="lower right") 
    #add titles
    axes[i_ax_mono].text(1.1e-4,2.3,"a) $N_{mono}=1$",fontsize=14)
    axes[i_ax_agg].text(1.1e-4,2.3,"b) $N_{mono}>1$",fontsize=14)
    ###
    #END: plot v-D lines according to m/A-D fits for different terminal velocity models and tumbling settings
    ###
    
    #save the plot (and open it)
    plt.tight_layout()
    dir_save = '/home/mkarrer/Dokumente/plots/4paper/'
    if not os.path.exists(dir_save): #create direktory if it does not exists
        os.makedirs(dir_save)

    out_filestring = "tumbling"
    plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
    plt.savefig(dir_save + out_filestring + '.png', dpi=400)
    print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
    subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])

    
