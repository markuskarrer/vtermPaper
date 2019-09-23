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

#from IPython.core.debugger import Tracer ; Tracer()()

'''
this code reads in properties from Jussis aggregate model, calculate + plots the following fallspeed
and fits different functions to each property (aspect ratio, fall speed for different tumbling settings)
'''

show_lines = dict()
#define which lines to show (ATTENTION: default, can be overwritten by optional arguments below):
show_lines["SB_mD"] = False #show the m-D relation of the SB scheme
show_lines["vel_fits"] = True #show the velocity based on the fit
show_lines["powerlaw_fits"] = False
show_lines["Atlas_fit"] = False #show the Atlas-type fit for cloud ice & snow
#show the fits as numbers
show_lines["fittext_vel"] = False #show text for the velocity fits
#from previous model settings
show_lines["SB_powerlaw"] = False #show also the old fit of Axels power-law
show_lines["SB_Atlas"] = False #show also the old fit of Axels Atlas-type
#from other models
show_lines["P3"] = False
show_lines["MC"] = False
show_lines["del_vterm"] = False #plot the velocity difference to all (ATTENTION: check implementation!!) velocity parameterizations 
show_lines["no_simparticle"] = False #if True, dont show the simulated particles (just compare parameterizations)
scale_vel_d = 0.001 #scale the velocity by some factor so it gets visible (TODO: add a real meaning to that value and a second y-axis)
show_lines["highest_diam_for_fit"] = -999 #see high_diam_log= ... below
prefer_vterm_model = "vterm_bohm" #due to the amount of plots, most plots are only done for one specific velocity model
#overwrite with optional arguments to this function
parser.add_argument("--show_SB_mD", help="")
parser.add_argument("--show_vel_fits", help="show the velocity fits (power-law/Atlas)")
parser.add_argument("--show_powerlaw_fits", help="")
parser.add_argument("--show_Atlas_fit", help="")
parser.add_argument("--show_fittext_vel", help="")
parser.add_argument("--show_SB_powerlaw", help="")
parser.add_argument("--show_SB_Atlas", help="")
parser.add_argument("--show_P3", help="")
parser.add_argument("--show_MC", help="")
parser.add_argument("--del_vterm", help="")
parser.add_argument("--show_no_simparticles", help="dont show the simulated particles (just compare parameterizations)")
parser.add_argument("--highest_diam_for_fit", help="highest diameter used for fitting in m")

args = parser.parse_args()
if args.show_vel_fits:
    show_lines["SB_mD"] = (args.show_SB_mD=="True")
    show_lines["vel_fits"] = (args.show_vel_fits=="True")
    show_lines["powerlaw_fits"] = (args.show_powerlaw_fits=="True")
    show_lines["Atlas_fit"] = (args.show_Atlas_fit=="True")
    show_lines["fittext_vel"] = (args.show_fittext_vel=="True")
    show_lines["SB_powerlaw"] = (args.show_SB_powerlaw=="True")
    show_lines["SB_Atlas"] = (args.show_SB_Atlas=="True")
    show_lines["P3"] = (args.show_P3=="True")
    show_lines["MC"] = (args.show_MC=="True")
    show_lines["del_vterm"] = (args.del_vterm=="True")
    show_lines["no_simparticle"] = (args.show_no_simparticles=="True")
    show_lines["highest_diam_for_fit"] = args.highest_diam_for_fit


#add displayed lines to string to distinguish them in the saved files
add_displayed_lines2string = '' 
for key in show_lines.keys():
    if show_lines[key]:
        add_displayed_lines2string+= '_' + key
#print add_displayed_lines2string
#from IPython.core.debugger import Tracer ; Tracer()()
#switch on/of fitting of terminal velocity (removes the fit lines)
fitvel=1 #0-> no fitting of velocity at all; 1-> fit only monomer and all aggregates; 2-> fit everything (including seperated for different monomer numbers)
tumbling=False #True: take the rotated projected area; False: take the area calculated from the aligned particle
take_mono_prop_theor=True #take theoretical values (as assumed in Jagg) instead of fits to N_mono=1
as_ratio_as_input2bohm=False #show the difference between considering the aspect ratio in B?hm and neglecting it

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
if show_lines["highest_diam_for_fit"]==-999:
    high_diam_log=-1.3979400086720375 #np.log10(3e-2)=-1.5228787452803376 #np.log10(4e-2)=-1.3979400086720375
else:
    high_diam_log=np.log10(float(show_lines["highest_diam_for_fit"])) #
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

    #START: get the 2D weighting function (derived from the comparison of the aggregate model output and McSnow runs)
    ######
    #calculate the histogram
    diam_edges,Nmono_edges,H_Jagg = generate_2Dhist_of_N_D_Nmono_from_MC_and_Jagg.calc_histogram(particle_dic["diam"],particle_dic["N_monomer"],low_diam_log=low_diam_log, high_diam_log=high_diam_log,nbins=40)
    #recalculate the weighting factor to be uniform
    H = np.nan*np.ones_like(H_Jagg)
    for i_diam in range(H_Jagg.shape[0]): #first dimension is the diameter
        ###calculate the weighting factor
        #we normalize the weighting factor at each diameter to one
        sum_valid_Jagg = np.nansum(H_Jagg[i_diam,:])
        if sum_valid_Jagg>0: #check if there is any Jagg data
            H[i_diam,:] = np.divide(1,H_Jagg[i_diam,:], out=np.zeros_like(H_Jagg[i_diam,:]), where=H_Jagg[i_diam,:]!=0)/sum(np.divide(1,H_Jagg[i_diam,:], out=np.zeros_like(H_Jagg[i_diam,:]), where=H_Jagg[i_diam,:]!=0))#sum_valid_Jagg
            #the next line might help to understand how the weighting is performed: 1. weight each diameter equally 2. homogenize each monomer bin for each diameter
            #print diam[i_diam],H_Jagg[i_diam,:],np.divide(1,H_Jagg[i_diam,:], out=np.zeros_like(H_Jagg[i_diam,:]), where=H_Jagg[i_diam,:]!=0),sum(np.divide(1,H_Jagg[i_diam,:], out=np.zeros_like(H_Jagg[i_diam,:]), where=H_Jagg[i_diam,:]!=0)); raw_input()
    weighting_factor=H

    
    #initialize the weighting attribute
    particle_dic["weight"] = np.zeros_like(particle_dic["diam"])
    print particle_dic
    #set the weight for each simulated aggregate
    for i_agg in range(0,particle_dic["diam"].shape[0]):#loop over all aggregates to give each aggregate the right weighting factor
        #TODO: this is probably very inefficiently programmed with nested loops and ifs, but maybe this is not important here
        for i_diam_now,diam_now in enumerate(diam_edges[:-1]):
            for i_N_mono_now,N_mono_now in enumerate(Nmono_edges[:-1]):
                if diam_edges[i_diam_now]<=particle_dic["diam"][i_agg]<diam_edges[i_diam_now+1]: #check if we are in the right diameter bin
                    if Nmono_edges[i_N_mono_now]<=particle_dic["N_monomer"][i_agg]<Nmono_edges[i_N_mono_now+1]: #check if we are in the right monomer number bin
                        particle_dic["weight"][i_agg] = weighting_factor[i_diam_now,i_N_mono_now] ##from IPython.core.debugger import Tracer ; Tracer()()    
    ######
    #END: calculate a weigting of the simulated particle to get an accurate fit for the whole diameter range and give equally importance to all monomer numbers
    ######
    ###plotting
    ###
    #initialize the plot
    ###
    number_of_plots = 8 #28
    
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
    #self-defined discrete colormap
    # define the colormap and the discrete boundaries (defined by the considered monomer numbers) and set up an array of the same colors (for line-plotting)
    '''
        cmapN = plt.cm.hot #brg
        #modify lowest color in colormap
        cmaplistN = [cmapN(i) for i in range(cmapN.N)]
        cmaplistN = cmaplistN[30:256-30] #crop colormap to avoid too dark/light scatterer #out of 256 colors
    '''
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
            #from IPython.core.debugger import Tracer ; Tracer()()
    #calculate the terminal velocitys of the simulated particles individually for each available velocity model
    for velocity_model in ["HW10","KC05","bohm","bohmar"]: #,"mitch_heym"]:
        for partal_add in ["","_partal10std","_partal20std","_partal40std","_partal60std"]:
            particle_dic["vterm_" + velocity_model + partal_add] = __fallspeed_relations.calc_vterm(velocity_model,particle_dic["mass"],particle_dic["diam"],particle_dic["area" + partal_add],particle_dic["as_ratio"])
            
    print particle_dic #overview of the dictionary, might be helpful sometimes

    
   
    ####
    #START: plot aspect_ratio
    ###
    if not show_lines["no_simparticle"]:
        #change the scale
        axes[i_ax].set_xscale("log") 
        axes[i_ax].set_yscale("linear") 
        #plot the scatter plot of the individual particle # use the colormap and the bounds which where defined before
        im = axes[i_ax].scatter(particle_dic["diam"],particle_dic["as_ratio"],s=1,c=particle_dic["N_monomer"],rasterized=True,norm=norm,cmap=cmap,marker='o')
    #make labels
    axes[i_ax].set_xlabel("diameter D / m")
    axes[i_ax].set_ylabel("aspect ratio as / 1" ) #TODO: plot also the untis of these properties

    #change the axis
    axes[i_ax].set_xlim([10**low_diam_log,10**high_diam_log]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
    axes[i_ax].set_ylim([0,1]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
    axes[i_ax].grid(which="major")

    #add colorbar
    cbar = fig.colorbar(im,ax=axes[i_ax])
    cbar.set_label("monomer number")
    #shift ticks to the middle of the color
    tick_locs = (bounds[:-1:3]) + 0.5*(bounds[1::3]-bounds[:-1:3])
    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels(N_mono_list[::3])

    #shrink colorbar lables
    cbar.ax.tick_params(labelsize=colbarlabelsize)
    ####
    #END: plot aspect_ratio
    ###
    
    #loop over different terminal velocity models in which 1. the velocity of each individual particle is plotted 2. the velocity according to the m-D and A-D fits 3. different possible fits (power-law, Atlas type) are calculated (and plotted)
    if as_ratio_as_input2bohm:
        vterm_array=["vterm_bohm","vterm_bohmar"]
    else:
        vterm_array=["vterm_bohm","vterm_HW10","vterm_KC05"]
    for i_prop,prop in enumerate(vterm_array): #,"vterm_HW10","vterm_KC05"]): #,"vterm_mitch_heym"]): #,"HWJussi","KCJussi"]):
        print "fitting and plotting: ",prop
        i_ax=0
        #####
        #1. the velocity of each individual particle is plotted 
        #####
        for partal_add in ["","_partal10std","_partal20std","_partal40std","_partal60std"]:
            if True: #prop==prefer_vterm_model:
                i_ax+=1
                #change the scale (according to xscale_vec and yscale_vec)
                axes[i_ax].set_xscale("log") #axes[i_ax].set_xscale(xscale_vec[i_prop]) 
                axes[i_ax].set_yscale("linear") #axes[i_ax].set_yscale(yscale_vec[i_prop])
                if prop=="vterm_bohmar":
                    #plot the scatter plot of the individual particle # use the colormap and the bounds which where defined before
                    im = axes[i_ax].scatter(particle_dic["diam"],particle_dic[prop + partal_add],s=1,c=particle_dic["N_monomer"],rasterized=True,norm=norm,cmap=cmap,marker='o') #,norm=colors.LogNorm(),cmap="gist_ncar")+
                    im = axes[i_ax].scatter(particle_dic["diam"][particle_dic["N_monomer"]==1],particle_dic[prop][particle_dic["N_monomer"]==1],s=1,c=particle_dic["N_monomer"][particle_dic["N_monomer"]==1],rasterized=True,norm=norm,cmap=cmap,marker='o') #,norm=colors.LogNorm(),cmap="gist_ncar")+
                elif prop=="vterm_bohm":
                    #plot the scatter plot of the individual particle # use the colormap and the bounds which where defined before
                    im = axes[i_ax].scatter(particle_dic["diam"],particle_dic[prop + partal_add],s=1,c=particle_dic["N_monomer"],rasterized=True,norm=norm2,cmap=cmap,marker='o') #,norm=colors.LogNorm(),cmap="gist_ncar")+
                    im = axes[i_ax].scatter(particle_dic["diam"][particle_dic["N_monomer"]==1],particle_dic[prop][particle_dic["N_monomer"]==1],s=1,c=particle_dic["N_monomer"][particle_dic["N_monomer"]==1],rasterized=True,norm=norm2,cmap=cmap,marker='o') #,norm=colors.LogNorm(),cmap="gist_ncar")+
                    
            if not as_ratio_as_input2bohm:
                
                ####
                #2. the velocity according to the m-D and A-D fits
                ####
                #calculate the terminal velocity based on the m-D and A-D fits and plot this line # use the colormap and the bounds which where defined before
                velocity_model = prop[6:] #get the name of the terminal velocity model from the property name (prop) which starts with "vterm"
                #change the scale (according to xscale_vec and yscale_vec)
                axes[i_ax].set_xscale("log") #axes[i_ax].set_xscale(xscale_vec[i_prop]) 
                axes[i_ax].set_yscale("linear") #axes[i_ax].set_yscale(yscale_vec[i_prop])
                for i,N_mono_now in enumerate([1]): #N_mono_list):
                    fit_dic["vterm_" + velocity_model + partal_add + "_fitted_via_mD_AD_Nmono" + str(N_mono_now)] = __fallspeed_relations.calc_vterm(velocity_model,fit_dic["mass" + "_Nmono_" + str(N_mono_now)],fit_dic["diam"],fit_dic["area" + partal_add + "_Nmono_" + str(N_mono_now)])
                    if prop==prefer_vterm_model:
                        if show_lines["vel_fits"]:
                            if N_mono_now<2 or fitvel>=2:
                                fitline_Nmonospecific, = axes[i_ax].plot(fit_dic["diam"],fit_dic["vterm_" + velocity_model + partal_add  + "_fitted_via_mD_AD_Nmono" + str(N_mono_now)],marker='None',c=usecolors[i],label="Nmono=1") #,norm=colors.LogNorm(),cmap="gist_ncar")
                                print partal_add,fit_dic["vterm_" + velocity_model + partal_add  + "_fitted_via_mD_AD_Nmono" + str(N_mono_now)]; raw_input()
                
                #now for all aggregates
                fit_dic["vterm_" + velocity_model + partal_add + "_fitted_via_mD_AD_Nmonoallagg"] = __fallspeed_relations.calc_vterm(velocity_model,fit_dic["mass" + "_Nmono_allagg"],fit_dic["diam"],fit_dic["area" + partal_add + "_Nmono_allagg"])
                if prop==prefer_vterm_model:
                    if show_lines["vel_fits"]:
                        if fitvel>=1:
                            if prop==prefer_vterm_model:
                                fitline_Nmonoallagg, = axes[i_ax].plot(fit_dic["diam"],fit_dic["vterm_" + velocity_model + partal_add + "_fitted_via_mD_AD_Nmonoallagg"],marker='None',c="g",label="Nmono>1")
                            #fitline_Nmonoallagg_dense_enough, = axes[i_ax].plot(fit_dic["diam"][fit_dic["diam_dens_enough_allagg"]],fit_dic["vterm_" + velocity_model + "_fitted_via_mD_AD_Nmonoallagg"][fit_dic["diam_dens_enough_allagg"]],marker='None',c="black",label="Nmono>1",linewidth=linewidthvalid)
                #####
                #3. different fits (power-law, Atlas type) are calculated (and plotted)
                #####
                if fitvel>=1:
                    #fit an Atlas type to the "mean v-D" relations (derived from the A-D and m-D fits)
                    for i,(str_monomer_agg,corr_cat) in enumerate(zip(["1","allagg"],["Nmono=1","Nmono>1"])): 
                        ##Atlas-type
                        [fitresult,covar] = __tools_for_processing_Jagg.fit_data(fit_dic["diam"],fit_dic[prop +  "_fitted_via_mD_AD_Nmono" + str_monomer_agg],func="Atlas") #, weight=fit_dic["dens" + "_Nmono_" + str_monomer_agg])
                        #from IPython.core.debugger import Tracer ; Tracer()()


                        [fitresult_Deq,covar_Deq] = __tools_for_processing_Jagg.fit_data(((6.*fit_dic["mass_coeff_Nmono_" + str_monomer_agg][0]*fit_dic["diam"]**fit_dic["mass_coeff_Nmono_" + str_monomer_agg][1]/(np.pi*1000.))**(1./3.)),fit_dic[prop +  "_fitted_via_mD_AD_Nmono" + str_monomer_agg],func="Atlas") #, weight=fit_dic["dens" + "_Nmono_" + str_monomer_agg])
                        
                        fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg] = fitresult #copy the Atlas-fit (as a function of D_max) coefficients in the fit_dic dictionary
                        fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_Deq"] = fitresult_Deq #copy the Atlas-fit (as a function of D_eq) coefficients in the fit_dic dictionary

                        #recalculate v(D)= A-B*exp(-gam*D) from fit
                        fit_dic[prop + "_Nmono_" + str_monomer_agg] = fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg][0]-fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg][1]*np.exp(-fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg][2]*fit_dic["diam"]) #calculate arrays of masses and areas based on the m-D fits
                        
                        #power-law
                        [fitresult,covar] = __tools_for_processing_Jagg.fit_data(fit_dic["diam"],fit_dic[prop +  "_fitted_via_mD_AD_Nmono" + str_monomer_agg],func="powerlaw") #, weight=fit_dic["dens" + "_Nmono_" + str_monomer_agg])
                        fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_powerlaw"] = fitresult #copy the power-law coefficients in the fit_dic dictionary

                        #recalculate powerlaw from fit
                        fit_dic[prop + "_Nmono_" + str_monomer_agg + "_powerlaw"] = fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_powerlaw"][0]*fit_dic["diam"]**fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_powerlaw"][1] #calculate arrays of masses and areas based on the m-D fits
                        if prop==prefer_vterm_model:
                            if show_lines["Atlas_fit"]:
                                Atlasfit_Nmonoallagg, = axes[i_ax].semilogx(fit_dic["diam"],fit_dic[prop + "_Nmono_" + str_monomer_agg],marker='None',c=np.array(["b","g"])[i],linestyle="-",label="Atlas-type " + corr_cat,linewidth=linewidthvalid)
                                if show_lines["del_vterm"]:
                                    Atlasfit_Nmonoallagg_delv, = axes[i_ax].semilogx(fit_dic["diam"][:-1],scale_vel_d*1./np.diff(fit_dic["diam"])*np.diff(fit_dic[prop + "_Nmono_" + str_monomer_agg]),marker='x',c=np.array(["b","g"])[i],linestyle="-",label="Atlas-type " + corr_cat + " del v",linewidth=linewidthvalid)
                            if show_lines["powerlaw_fits"]:
                                powerlaw_handle, = axes[i_ax].semilogx(fit_dic["diam"], fit_dic[prop + "_Nmono_" + str_monomer_agg + "_powerlaw"],marker='None',c=np.array(["b","g"])[i],linestyle="-.",label="powerlaw " + corr_cat,linewidth=linewidthvalid)
                if prop==prefer_vterm_model:
                    if fitvel>=1:
                        if show_lines["fittext_vel"]:
                            #add labels for the fit result to the plot
                            axes[i_ax] = __tools_for_processing_Jagg.add_fitresult_text(axes[i_ax],fit_dic,prop,[1],function="Atlas")
                            axes[i_ax] = __tools_for_processing_Jagg.add_fitresult_text(axes[i_ax],fit_dic,prop,[1],function="powerlaw_fallspeed")

                
            if len(vterm_array)==(i_prop+1):
                if partal_add=='':
                    axes[i_ax].text(0.999,0.00,"no tumbling",horizontalalignment='right',verticalalignment='bottom',transform = axes[i_ax].transAxes)

                else:
                    #add a description of the tumbling settings
                    axes[i_ax].text(0.999,0.00,"tumbling std: "+ partal_add[7:9] + "$^\circ$",horizontalalignment='right',verticalalignment='bottom',transform = axes[i_ax].transAxes)

                #show legend
                axes[i_ax].legend() 
                #expl1, = axes[i_ax].semilogx(np.nan,np.nan,label="model assumptions:",linestyle=""); expl2, = axes[i_ax].semilogx(np.nan,np.nan,label="fit to simulated particles:",linestyle="");axes[i_ax].legend(handles=[expl2,Atlasfit_Nmono1,Atlasfit_Nmonoallagg,expl1,SB_ci_handle,SB_snow_handle,MC_handle,P3_handle]) #uncomment this line for a special order
                
                #make labels
                axes[i_ax].set_xlabel("diameter D / m")
                axes[i_ax].set_ylabel(prop + " v / m s-1" ) #TODO: plot also the untis of these properties

                #change the axis
                axes[i_ax].set_xlim([10**low_diam_log,10**high_diam_log]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                axes[i_ax].set_ylim([0,2.5]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                axes[i_ax].grid(which="major")

                #add colorbar
                cbar = fig.colorbar(im,ax=axes[i_ax])
                cbar.set_label("monomer number")
                #shift ticks to the middle of the color
                tick_locs = (bounds[:-1:3]) + 0.5*(bounds[1::3]-bounds[:-1:3])
                cbar.set_ticks(tick_locs)
                cbar.set_ticklabels(N_mono_list[::3])

                #shrink colorbar lables
                cbar.ax.tick_params(labelsize=colbarlabelsize)
    if not as_ratio_as_input2bohm:
        
        ###
        #START: plot v-D lines according to m/A-D fits for different terminal velocity models and tumbling settings
        ###
        i_ax+=1
        i_ax_mono=i_ax
        i_ax_agg=i_ax+1
        for i_align,partal_add in enumerate(["","_partal20std","_partal40std"]): #,"_partal10std","_partal20std","_partal40std","_partal60std"]):
            for i_prop,prop in enumerate(["vterm_bohm","vterm_HW10","vterm_KC05"]): #"vterm_bohm","vterm_HW10","vterm_KC05"]): #,"vterm_mitch_heym"]): #,"HWJussi","KCJussi"]):
                velocity_model = prop[6:] #get the name of the terminal velocity model from the property name (prop) which starts with "vterm"
                
                #lines for monomers
                axes[i_ax_mono].semilogx(fit_dic["diam"],fit_dic["vterm_" + velocity_model + partal_add + "_fitted_via_mD_AD_Nmono1"],marker='None',color=np.array(["k","r","y"])[i_prop],label=np.array([velocity_model + " monomer","__None","__None","__None","__None"])[i_align],linestyle=np.array(["-","--","-.",":"])[i_align])

                #all agregates
                axes[i_ax_agg].semilogx(fit_dic["diam"],fit_dic["vterm_" + velocity_model + partal_add + "_fitted_via_mD_AD_Nmonoallagg"],marker='None',color=np.array(["k","r","y"])[i_prop],label=np.array([velocity_model + " all agg","__None","__None","__None","__None"])[i_align],linestyle=np.array(["-","--","-.",":"])[i_align])
                

            for i_ax_now in [i_ax_mono,i_ax_agg]:
                #add a legend for the tumbling
                if i_align==0:  #no tumbling
                    axes[i_ax_now].plot(np.nan,np.nan,c='k',linestyle=np.array(["-","--","-.",":"])[i_align],label="no tumbling")
                else:
                    axes[i_ax_now].plot(np.nan,np.nan,c='k',linestyle=np.array(["-","--","-.",":"])[i_align],label="tumbling std: "+ partal_add[7:9] + "deg")
        for i_ax_now in [i_ax_mono,i_ax_agg]:
            #show legend
            axes[i_ax_now].legend() 
            
            #make labels
            axes[i_ax_now].set_xlabel("diameter D / m")
            axes[i_ax_now].set_ylabel("terminal velocity v / m s-1" ) #TODO: plot also the untis of these properties

            #change the axis
            axes[i_ax_now].set_xlim([10**low_diam_log,10**high_diam_log]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
            axes[i_ax_now].set_ylim([0,2.5]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
            axes[i_ax_now].grid(which="major")
        
        ###
        #END: plot v-D lines according to m/A-D fits for different terminal velocity models and tumbling settings
        ###
    
    #save the plot (and open it)
    plt.tight_layout()
    dir_save = '/home/mkarrer/Dokumente/plots/Jagg/'
    if not os.path.exists(dir_save): #create direktory if it does not exists
        os.makedirs(dir_save)

    if as_ratio_as_input2bohm:
        add_asratio_in_hydro = "_asratio_in_bohm"
    else:
        add_asratio_in_hydro = ""
    out_filestring = "Jagg_asratio_tumbling" + particle_type_comb +  "_compscat2models" +  add_displayed_lines2string + add_asratio_in_hydro
    plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
    plt.savefig(dir_save + out_filestring + '.png', dpi=400)
    print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
    subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
    #plt.clf()
    #plt.close()
    

    
    ax_description_array=["asratio","vterm_0deg","vterm_10deg","vterm_20deg","vterm_40deg","vterm_60deg","vterm_comb_mono","vterm_comb_agg"]
    for save_axis,ax_description in enumerate(ax_description_array):
        if "vterm_comb" in ax_description:
            print "saving subfigure",save_axis,ax_description
            try: #clear label of axis before and recover current one
                #pass
                axes[save_axis].set_xlabel("diameter D / m")
                axes[save_axis-1].set_xlabel("")
            except:
                pass
            extent = axes[save_axis].get_window_extent().transformed(fig.dpi_scale_trans.inverted())

            fig.savefig('/home/mkarrer/Dokumente/plots/tmp.pdf',bbox_inches=extent.expanded(1.6, 1.3),dpi=400) #(1.6, 1.4) are width and height expanded around center

            subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/tmp.pdf /home/mkarrer/Dokumente/plots/4paper/' + out_filestring + '_' + ax_description + '.pdf',shell=True)


    # Save just the portion _inside_ the third axis's boundaries
    save_axis=3
    for i_axis in range(0,len(axes)):#clear all other axes
        if i_axis==save_axis:
            continue
        axes[i_axis].set_xlabel("")
        axes[i_axis].axes.get_xaxis().set_ticklabels([])
    extent = axes[save_axis].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    axes[save_axis].set_ylabel("terminal velocity v / m s-1")
    fig.savefig(dir_save + out_filestring + 'spec_ax_figure.png', bbox_inches=extent.expanded(1.5, 1.4),dpi=400)
    fig.savefig(dir_save + out_filestring + 'spec_ax_figure.pdf', bbox_inches=extent.expanded(1.5, 1.4),dpi=400)