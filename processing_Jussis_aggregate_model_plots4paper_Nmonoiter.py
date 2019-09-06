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
#import other self-defined functions
import __postprocess_McSnow
import __postprocess_SB
import __fallspeed_relations
import __tools_for_processing_Jagg
import __plotting_functions
import __setup_mDAD_frommodels
import generate_2Dhist_of_N_D_Nmono_from_MC_and_Jagg 
#from IPython.core.debugger import Tracer ; Tracer()()
from matplotlib import rc


'''
this code reads in properties from Jussis aggregate model
and fits different functions to each property (m-D(monomer dependent and binary),A-D(monomer dependent and binary),v-D for all fallspeed models) and displays them together
'''
tumbling=False #so far only the use of non-tumbling projected area is implemented here
only_spec_Nmono=True #plot only specific monomer numbers

particle_types = ["needle"] #["plate","dendrite","mixdendneedle","needle"] #["needle","column","plate","dendrite","bullet","rosette"] # ,"bullet"rosette
for particle_type in particle_types:

    #read the properties of the individual particles from the files
    particle_dic,N_mono_list = __tools_for_processing_Jagg.read_particle_prop_files(prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_Nmonoiter/",
                                                                            D_small_notuse=1e-4, #ATTENTION: particle smaller than D_small_notuse are not considered (e.g. because of resolution issues)
                                                                            N_small_notuse=1, #ATTENTION:  monomer numbers smaller than N_small_notuse are not considered
                                                                            grid_res_array = [10e-6], #array of grid resolutions (defines the files which are read in)
                                                                            particle_type = particle_type, #define the habit
                                                                            test_with_small_sample = True
                                                                            )
    #crop N_mono_list to not show to large N_mono/N_initial rations
    N_mono_list[(N_mono_list<=100)] # & (N_mono_list>=50)]
    particle_dic["N_ratio"] = particle_dic["N_initial"]/particle_dic["N_monomer"] #calculate the ratio between initial population and N_mono
    #from IPython.core.debugger import Tracer ; Tracer()()

    if only_spec_Nmono:
        #particle_dic["diam"][((particle_dic["N_ratio"])>=1) && ((particle_dic["N_ratio"])<=1)] = np.nan
        particle_dic["diam"][(particle_dic["N_monomer"]!=2)] = np.nan
        #get rid of other than N_initial/N_mono=[1,100]
    #show the dictionary (e.g. to check that it's not empty because the path was wrong)
    print particle_dic
    #randomize order in dictionary (for an unbiased view in scatter-plot)
    N_particles = particle_dic["area"].shape[0]
    particle_indices_array = range(0,N_particles)
    random.shuffle(particle_indices_array) #shuffle the indice array (random order)
    for key in particle_dic.keys():
        particle_dic[key] = particle_dic[key][particle_indices_array] #shuffle the arrays
    #get m-D from the assumptions in the aggregate model (this is used at various places)
    a,b,c,d = __tools_for_processing_Jagg.calc_mD_AD_coeffs(particle_type)

    if tumbling:
        particle_dic["area"] = particle_dic["area_partal40std"] #use the tumbled area instead
    
    #calculate the terminal velocitys of the simulated particles individually for each available velocity model
    for velocity_model in ["HW10","KC05","bohm"]: #,"mitch_heym"]:
            particle_dic["vterm_" + velocity_model] = __fallspeed_relations.calc_vterm(velocity_model,particle_dic["mass"],particle_dic["diam"],particle_dic["area"])
    
    #####################
    ######START: fitting
    #####################
    #initialize an array which contains all fit results
    fitresult_arr_dir_powerlaw_string = np.zeros((5,N_mono_list.shape[0],2)) #fit results (saved as floats) #dimensions: 0-> [mass,array]; 1-> monomer number; 2-> parameters (a,b) in fit

    fit_dic = dict() #initialize fit dictionary
    #set up array of diameters for displaying the fits
    #low_diam_log_detailed= -4; high_diam_log_detailed=-1.3979400086720375
    low_diam_log_detailed= np.log10(5e-4); high_diam_log_detailed=np.log10(5e-3)

    diam = np.logspace(low_diam_log_detailed,high_diam_log_detailed,500) #set up array of diameters (for bin edges and KDE-grid)
    diam_log_detailed = np.log10(diam)
    diam_center = (diam[:-1]+diam[1:])/2 #center of the diameter bins
    diam_log = np.log10(diam_center)
    diam_log_center = 10**((np.log10(diam[:-1])+np.log10(diam[1:]))/2) #logarithmic center of the diameter bins
    #define at which monomer numbers the fit should be evaluated
    Nmono_fit_eval = np.array([1,2,5,10,100]) 
    Nmono_log = np.log10(Nmono_fit_eval)
    
    #pass already the diameter to the fit-array
    fit_dic["diam"] = diam
    '''
    ###
    #START: monomer number dependent fitting 
    ###

    #fit the m-D and A-D properties
    print "monomer number dependent fitting"
    for i_prop,(prop,prop_unit) in enumerate(zip(["mass","area"],["kg","m2"])): #,"HWJussi","KCJussi"]):

        #fit the modified power-law m(D,N_mono)= a(N_mono)*D**(b(N_mono)* a_0*D**b0
        ab_func="powerlaw_rational1" # "powerlaw_polynom1"powerlaw_rational1" "powerlaw_rational2" and "powerlaw_rational3" need some more advanced initial guesses #select the functional form of a and b 
        if prop=="mass" or (not tumbling): #for mass and area without tumbling we can constrain more parameters

            if ab_func=="powerlaw_rational1":
                fitresult_mod_powerlaw = __tools_for_processing_Jagg.fit_2D_rational_function_transformwrapper(particle_dic["diam"],particle_dic["N_monomer"],particle_dic[prop],func='powerlaw_rational1_fixp0_fixr0',method='lm',weight='None',habit=particle_type,prop=prop)
            else: #look for other functions in process_Jussis_aggregate_model_mod_powerlaw
                print "func: ", ab_func, " not (carefully) implemented in process_Jussis_aggregate_model_plots4paper (and its modules)"; sys.exit()
        elif prop=="area" and tumbling:
            #fitresult_mod_powerlaw = __tools_for_processing_Jagg.fit_2D_rational_function_transformwrapper(particle_dic["diam"],particle_dic["N_monomer"],particle_dic[prop],func='rational2_2_pade_fixq00',method='lm',weight='None',habit=particle_type,prop=prop)
            if ab_func=="powerlaw_rational1":
                fitresult_mod_powerlaw = __tools_for_processing_Jagg.fit_2D_rational_function_transformwrapper(particle_dic["diam"],particle_dic["N_monomer"],particle_dic[prop],func='powerlaw_rational1',method='lm',weight='None',habit=particle_type,prop=prop)
            else: #TODO: include different fitting methods
                print "func: ", ab_func, " not (carefully) implemented in process_Jussis_aggregate_model_plots4paper (and its modules)"; sys.exit()
                
        #save the fit results in a dictionary
        fit_dic[prop +"_coeff_mod_powerlaw"] = fitresult_mod_powerlaw
        
        ###
        #END: monomer number dependent fitting 
        ###
        
        
        ###
        #START: get the fitted properties from the fit coefficients
        ###

        #spanned grid from D-N
        D_N_grid =  np.meshgrid(diam_log,Nmono_log)
        D_N_grid_detailed = np.meshgrid(diam_log_detailed,Nmono_log)
        #apply the fit coefficients to the above spanned D,N field
        if ab_func=="powerlaw_rational1":
            if prop=="mass" or (not tumbling):
                prop_fitted = __tools_for_processing_Jagg.powerlaw_rational1_fixp0_fixr0(D_N_grid,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3])#this is the logscaled property (mass/area) divided by the monomer property
                prop_fitted_detailed = __tools_for_processing_Jagg.powerlaw_rational1_fixp0_fixr0(D_N_grid_detailed,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3])#this is the logscaled property (mass/area) divided by the monomer property
            elif prop=="area" and tumbling:
                prop_fitted = __tools_for_processing_Jagg.powerlaw_rational1(D_N_grid,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3],fit_dic[prop +"_coeff_mod_powerlaw"][4],fit_dic[prop +"_coeff_mod_powerlaw"][5])#this is the logscaled property (mass/area) divided by the monomer property
                prop_fitted_detailed = __tools_for_processing_Jagg.powerlaw_rational1(D_N_grid_detailed,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3],fit_dic[prop +"_coeff_mod_powerlaw"][4],fit_dic[prop +"_coeff_mod_powerlaw"][5])#this is the logscaled property (mass/area) divided by the monomer property
            else: #TODO: include different fitting methods
                    print "func: ", ab_func, " not (carefully) implemented in process_Jussis_aggregate_model_mod_powerlaw (and its modules)"; sys.exit()


        #transform the normed property back   
        if prop=="mass":
            mass_fit = 10**(prop_fitted)*a*diam_center**b
            mass_fit_detailed = 10**(prop_fitted_detailed)*a*diam**b#for plotting the velocity more bins are nicer
        if prop=="area":
            area_fit = 10**(prop_fitted)*c*diam_center**d
            area_fit_detailed = 10**(prop_fitted_detailed)*c*diam**d#for plotting the velocity more bins are nicer

        ###
        #END: get the fitted properties from the fit coefficients
        ###
        
        ###
        #START: get fits for all aggregates
        ###
        
        ##
        # make a fit for all aggregates in N_mono_list
        ##
        #fit to m-D and A-D coefficients for all particles with N_monomer>1
        [fit_dic[prop + "_coeff_Nmono_allagg"],covar] = __tools_for_processing_Jagg.fit_data(particle_dic["diam"],particle_dic[prop],func='powerlaw',weight="None")
        #calculate the values corresponding to the spanned diameter from the fit-coefficients
        fit_dic[prop + "_Nmono_allagg"] = fit_dic[prop + "_coeff_Nmono_allagg"][0]*fit_dic["diam"]**fit_dic[prop + "_coeff_Nmono_allagg"][1] #calculate arrays of masses and areas based on the m-D fits
        
        ###
        #END: get fits for all aggregates
        ###
        
    #####################
    ######END: fitting
    #####################
    '''
    #####################
    ######START: plotting figure1
    #####################
    
    ##general settings for the plotting
    number_of_plots = 7

    #optimize the appearance of the plot (figure size, fonts)
    [fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots,legend_fontsize='medium')
    i_ax=-1
    linewidth=0.7
    linewidth_white=1.0
    ##define colormaps
    #self-defined discrete colormap (for the monomer number N)
    # define the colormap and the discrete boundaries (defined by the considered monomer numbers) and set up an array of the same colors (for line-plotting)
    cmapN = plt.cm.viridis #brg
    #modify lowest color in colormap
    cmaplistN = [cmapN(i) for i in range(cmapN.N)]
    #cmaplistN = cmaplistN[30:256-30] #crop colormap to avoid too dark/light scatterer #out of 256 colors

    #colors for the binary relevant
    colormonomer = colors.ColorConverter.to_rgb("r") #(1.0,0.0,0.0,1.0) #make Nmono=1 blue
    colorallagg =  colors.ColorConverter.to_rgb("g") #(1.0,0.0,0.0,1.0) #make Nmono=1 blue
    
    #reflect the selected color also in the colormap
    cmaplistN[0] = colormonomer
    cmapN = colors.LinearSegmentedColormap.from_list('mcm',cmaplistN, cmapN.N)
    boundsN = np.append(N_mono_list,[9999])
    normN = colors.BoundaryNorm(boundsN, cmapN.N)
    usecolorsN = cmapN(np.linspace(0,1,N_mono_list.shape[0])) #get colors from colormap to make consistent line colors for the monomer number

    ###
    #START: plot the prop(m or D) vs diameter as scatter and the 2D-fit
    ###
    
    #do the same plots for mass and area     
    for i_prop,(prop,prop_unit,prop_short) in enumerate(zip(["mass","area"],["kg","m2"],["m","A"])):
        i_ax+=1 #start at first axis
        print "fitting ",prop
        #plot the individual particle
        axes[i_ax].set_yscale('log')
        axes[i_ax].set_xscale('log')
        im = axes[i_ax].scatter(particle_dic["diam"],particle_dic[prop],s=1,c=particle_dic["N_ratio"],rasterized=True,norm=normN,cmap=cmapN) #all used particles        
        '''
        #overlay the fit
        for i_Nmono,N_mono_now in enumerate(Nmono_fit_eval): #range(0,N_mono_list.shape[0],5):
            if prop=="mass":
                fitlines = axes[i_ax].plot(diam_center,mass_fit[i_Nmono,:],c="white",linewidth=linewidth_white,linestyle="-")
                fitlines = axes[i_ax].plot(diam_center,mass_fit[i_Nmono,:],c=usecolorsN[N_mono_now==N_mono_list][0],linewidth=linewidth,linestyle="-")

            if prop=="area":
                fitlines = axes[i_ax].plot(diam_center,area_fit[i_Nmono,:],c="white",linewidth=linewidth_white,linestyle="-")
                fitlines = axes[i_ax].plot(diam_center,area_fit[i_Nmono,:],c=usecolorsN[N_mono_now==N_mono_list][0],linewidth=linewidth,linestyle="-")

        #also plot all aggregates in this plot
        fitlines = axes[i_ax].plot(np.nan,np.nan,c="b",linewidth=2.*linewidth,label="monomers")
        #fitlines = axes[i_ax].plot(fit_dic["diam"],fit_dic[prop + "_Nmono_allagg"],c="w",linewidth=2.*linewidth_white,label="__None")
        fitlines = axes[i_ax].plot(fit_dic["diam"],fit_dic[prop + "_Nmono_allagg"],c="g",linewidth=2.*linewidth,label="all aggregates",linestyle='-.')
        '''
        #make labels
        axes[i_ax].set_xlabel("diameter D / m")
        axes[i_ax].set_ylabel((" / ").join((prop+ ' ' + prop_short,prop_unit))) 

        #change the axis
        axes[i_ax].set_xlim([1e-4,4e-2])
        axes[i_ax].set_ylim([0,np.array([1e-4,1e-3])[i_prop]]) #define the upper limit of the displayed axis
        axes[i_ax].grid(which="minor")
        #add colorbar
        cbar = fig.colorbar(im,ax=axes[i_ax])
        cbar.set_label(r"$N_{initial}/N_{mono}$")
        #shift ticks to the middle of the color
        tick_locs = (boundsN[:-1:3]) + 0.5*(boundsN[1::3]-boundsN[:-1:3])
        cbar.set_ticks(tick_locs)
        cbar.set_ticklabels(N_mono_list[::3])
        '''
        #add fit result text
        if prop=="mass" or (not tumbling):
            
            if ab_func=="powerlaw_rational1":
                
                if prop=="mass":
                    fit_string_coeff1 = "a= 10^(({:.2f}".format(fit_dic[prop +"_coeff_mod_powerlaw"][0])+"*log(N))/(1{:+.2f}".format(fit_dic[prop +"_coeff_mod_powerlaw"][1]) + "*log(N)))*{:.4f}".format(a)
                    fit_string_coeff2 = "b= ({:.4f}".format(fit_dic[prop +"_coeff_mod_powerlaw"][2])+"*log(N))/(1+{:.2f}".format(fit_dic[prop +"_coeff_mod_powerlaw"][3]) + "*log(N))+{:.2f}".format(b)
                elif prop =="area":
                    fit_string_coeff1 = "c= 10^(({:.2f}".format(fit_dic[prop +"_coeff_mod_powerlaw"][0])+"*log(N))/(1{:+.2f}".format(fit_dic[prop +"_coeff_mod_powerlaw"][1]) + "*log(N)))*{:.4f}".format(c)
                    fit_string_coeff2 = "d= ({:.4f}".format(fit_dic[prop +"_coeff_mod_powerlaw"][2])+"*log(N))/(1{:+.2f}".format(fit_dic[prop +"_coeff_mod_powerlaw"][3]) + "*log(N))+{:.2f}".format(d)
                
                fitted_samples = ""
                for Nmono in Nmono_fit_eval: 
                    if prop=="mass":
                        a_fitted = 10**(fit_dic[prop +"_coeff_mod_powerlaw"][0]*np.log10(Nmono)/(1+fit_dic[prop +"_coeff_mod_powerlaw"][1]*np.log10(Nmono)))*a
                        b_fitted = (fit_dic[prop +"_coeff_mod_powerlaw"][2]*np.log10(Nmono)/(1+fit_dic[prop +"_coeff_mod_powerlaw"][3]*np.log10(Nmono)))+b
                        fitted_samples+= "$N_{mono}$: " + str(Nmono) + " a: {:.4f}".format(a_fitted) + " b: {:.2f}".format(b_fitted) + "\n"
                    elif prop=="area":
                        c_fitted = 10**(fit_dic[prop +"_coeff_mod_powerlaw"][0]*np.log10(Nmono)/(1+fit_dic[prop +"_coeff_mod_powerlaw"][1]*np.log10(Nmono)))*c
                        d_fitted = (fit_dic[prop +"_coeff_mod_powerlaw"][2]*np.log10(Nmono)/(1+fit_dic[prop +"_coeff_mod_powerlaw"][3]*np.log10(Nmono)))+d
                        fitted_samples+= "$N_{mono}$: " + str(Nmono) + " c: {:.4f}".format(c_fitted) + " d: {:.2f}".format(d_fitted) + "\n"
                        
                #add also rows for all aggregates
                if prop=="mass":
                    fitted_samples+= "all agg: " "a: {:.4f}".format(fit_dic[prop + "_coeff_Nmono_allagg"][0]) + " b: {:.2f}".format(fit_dic[prop + "_coeff_Nmono_allagg"][1]) + "\n" 
                if prop=="area":
                    fitted_samples+= "all agg: " "c: {:.4f}".format(fit_dic[prop + "_coeff_Nmono_allagg"][0]) + " d: {:.2f}".format(fit_dic[prop + "_coeff_Nmono_allagg"][1]) + "\n" 
 
            axes[i_ax].text(0.02,0.98,fit_string_coeff1 + "\n" +fit_string_coeff2 + "\n" + fitted_samples,horizontalalignment='left',verticalalignment='top',transform=axes[i_ax].transAxes,fontsize=6)

        elif prop=="area" and tumbling:
            
            if ab_func=="powerlaw_rational1":
                
                fit_string_coeff1 = "c= 10^(({:.2f}".format(fit_dic[prop +"_coeff_mod_powerlaw"][0])+"{:+.2f}".format(fit_dic[prop +"_coeff_mod_powerlaw"][1])+"*log(N))/(1+{:.4f}".format(fit_dic[prop +"_coeff_mod_powerlaw"][2]) + "*log(N)))*{:.4f}".format(c)
                fit_string_coeff2 = "d= ({:.2f}".format(fit_dic[prop +"_coeff_mod_powerlaw"][3])+"{:+.2f}".format(fit_dic[prop +"_coeff_mod_powerlaw"][4])+"*log(N))/(1+{:.2f}".format(fit_dic[prop +"_coeff_mod_powerlaw"][5]) + "*log(N)))+{:.4f}".format(d)
                
                fitted_samples = ""
                for Nmono in Nmono_fit_eval: 
                    c_fitted = 10**((fit_dic[prop +"_coeff_mod_powerlaw"][0]+fit_dic[prop +"_coeff_mod_powerlaw"][1]*np.log10(Nmono))/(1+fit_dic[prop +"_coeff_mod_powerlaw"][2]*np.log10(Nmono)))*c
                    d_fitted = (fit_dic[prop +"_coeff_mod_powerlaw"][3]+fit_dic[prop +"_coeff_mod_powerlaw"][4]*np.log10(Nmono))/(1+fit_dic[prop +"_coeff_mod_powerlaw"][5]*np.log10(Nmono))+d
                    fitted_samples+= "$N_{mono}$: " + str(Nmono) + " c: {:.4f}".format(c_fitted) + " d: {:.2f}".format(d_fitted) + "\n"
                        
                #add also rows for all aggregates
                if prop=="mass":
                    fitted_samples+= "all agg: " "a: {:.4f}".format(fit_dic[prop + "_coeff_Nmono_allagg"][0]) + " b: {:.2f}".format(fit_dic[prop + "_coeff_Nmono_allagg"][1]) + "\n" 
                if prop=="area":
                    fitted_samples+= "all agg: " "c: {:.4f}".format(fit_dic[prop + "_coeff_Nmono_allagg"][0]) + " d: {:.2f}".format(fit_dic[prop + "_coeff_Nmono_allagg"][1]) + "\n" 
 
            axes[i_ax].text(0.02,0.98,fit_string_coeff1 + "\n" +fit_string_coeff2 + "\n" + fitted_samples,horizontalalignment='left',verticalalignment='top',transform=axes[i_ax].transAxes,fontsize=6)
        '''
        ###
        #END: plot the prop(m or D) vs diameter as scatter and the 2D-fit
        ###
        
        ###
        #START: plot the prop(m or D) vs diameter as scatter and the 2D-fit (with transformed values D->log10(D) N->log10(N) prop->prop/prop_monomer)
        ###
        #mass vs diameter
        i_ax+=1

        #plot the individual particle normalized by the properties of the monomer
        match_Nmono = (particle_dic["N_monomer"]<=100)
        if prop=="mass":
            #im = axes[i_ax].scatter(np.log10(particle_dic["diam"][match_Nmono]),np.log10(particle_dic[prop][match_Nmono]/(a*particle_dic["diam"][match_Nmono]**b)),s=1,c=particle_dic["N_ratio"][match_Nmono],rasterized=True,norm=normN,cmap=cmapN) #old with log10(..) as label
            im = axes[i_ax].scatter(particle_dic["diam"][match_Nmono],particle_dic[prop][match_Nmono]/(a*particle_dic["diam"][match_Nmono]**b),s=1,c=particle_dic["N_ratio"][match_Nmono],rasterized=True,norm=normN,cmap=cmapN) 

        elif prop=="area":
            #im = axes[i_ax].scatter(np.log10(particle_dic["diam"][match_Nmono]),np.log10(particle_dic[prop][match_Nmono]/(c*particle_dic["diam"][match_Nmono]**d)),s=1,c=particle_dic["N_ratio"][match_Nmono],rasterized=True,norm=normN,cmap=cmapN) #old with log10(..) as label
            im = axes[i_ax].scatter(particle_dic["diam"][match_Nmono],particle_dic[prop][match_Nmono]/(c*particle_dic["diam"][match_Nmono]**d),s=1,c=particle_dic["N_ratio"][match_Nmono],rasterized=True,norm=normN,cmap=cmapN)
        '''
        #overlay the scatterer with the fit at some selected monomer numbers
        for i_Nmono,N_mono_now in enumerate(Nmono_fit_eval):
            #fitlines = axes[i_ax].plot(np.log10(diam_center),(prop_fitted[i_Nmono,:]),c="white",linewidth=linewidth_white,linestyle="-")
            fitlines = axes[i_ax].semilogx(diam_center,10**prop_fitted[i_Nmono,:],c="white",linewidth=linewidth_white,linestyle="-") #old with log10(..) as label
            #fitlines = axes[i_ax].plot(np.log10(diam_center),(prop_fitted[i_Nmono,:]),c=usecolorsN[N_mono_now==N_mono_list][0],linewidth=linewidth,linestyle="-") #old with log10(..) as label
            fitlines = axes[i_ax].semilogx(diam_center,10**prop_fitted[i_Nmono,:],c=usecolorsN[N_mono_now==N_mono_list][0],linewidth=linewidth,linestyle="-")
            
        #also plot all aggregates in this plot
        if prop=="mass":
            #fitlines = axes[i_ax].plot(np.log10(fit_dic["diam"]),np.log10(fit_dic[prop + "_Nmono_allagg"]/(a*fit_dic["diam"]**b)),c="g",linewidth=2.*linewidth,label="all aggregates")#old with log10(..) as label
            fitlines = axes[i_ax].plot(np.nan,np.nan,c="b",linewidth=linewidth,label="monomers")
            #fitlines_white = axes[i_ax].semilogx(fit_dic["diam"],fit_dic[prop + "_Nmono_allagg"]/(a*fit_dic["diam"]**b),c="w",linewidth=2.*linewidth_white,label="__None")
            fitlines = axes[i_ax].semilogx(fit_dic["diam"],fit_dic[prop + "_Nmono_allagg"]/(a*fit_dic["diam"]**b),c="g",linewidth=2.*linewidth,label="all aggregates",linestyle='-.')

        elif prop=="area":
            #fitlines = axes[i_ax].plot(np.log10(fit_dic["diam"]),np.log10(fit_dic[prop + "_Nmono_allagg"]/(c*fit_dic["diam"]**d)),c="g",linewidth=2.*linewidth,label="all aggregates")#old with log10(..) as label
            fitlines_white = axes[i_ax].plot(fit_dic["diam"],fit_dic[prop + "_Nmono_allagg"]/(c*fit_dic["diam"]**d),c="w",linewidth=2.*linewidth,label="__None")
            fitlines = axes[i_ax].plot(np.nan,np.nan,c="b",linewidth=linewidth,label="monomers")
            fitlines = axes[i_ax].plot(fit_dic["diam"],fit_dic[prop + "_Nmono_allagg"]/(c*fit_dic["diam"]**d),c="g",linewidth=2.*linewidth,label="all aggregates",linestyle='-.')
        '''
        #set yticks and yticklabels
        #axes[i_ax].set_yticks([0.1,0.2,0.3,0.4,0.5,0.7,1.0,2.0])
        #axes[i_ax].set_yticklabels([0.1,0.2,0.3,0.4,0.5,0.7,1.0,2.0])

        axes[i_ax].set_xscale("log")
        #make labels
        axes[i_ax].set_xlabel("diameter D / m")
        axes[i_ax].set_ylabel(prop + "/ \n" + prop + " monomer / 1")
        #change the axis
        axes[i_ax].set_xlim([1e-4,4e-2])
        #axes[i_ax].set_ylim([0,np.array([1e-4,1e-3])[i_prop]]) #define the upper limit of the displayed axis
        axes[i_ax].grid(which="minor")
        #add colorbar
        cbar = fig.colorbar(im,ax=axes[i_ax])
        cbar.set_label(r"$N_{initial}/N_{mono}$")
        #shift ticks to the middle of the color
        tick_locs = (boundsN[:-1:3]) + 0.5*(boundsN[1::3]-boundsN[:-1:3])
        cbar.set_ticks(tick_locs)
        cbar.set_ticklabels(N_mono_list[::3])
        ###
        #END: plot the prop(m or D) vs diameter as scatter and the 2D-fit (with transformed values D->log10(D) N->log10(N) prop->prop/prop_monomer)
        ###
        
    ####
    #plot vterm and its fits
    ####
    for i_prop,prop in enumerate(["vterm_bohm","vterm_HW10","vterm_KC05"]): #,"vterm_mitch_heym"]): #,"HWJussi","KCJussi"]):
        print "fitting and plotting: ",prop
        
        #define current axis
        i_ax+=1         
        #change the scale (according to xscale_vec and yscale_vec)
        axes[i_ax].set_xscale("log") #axes[i_ax].set_xscale(xscale_vec[i_prop]) 
        axes[i_ax].set_yscale("linear") #axes[i_ax].set_yscale(yscale_vec[i_prop])
        
        #####
        #1. the velocity of each individual particle is plotted 
        #####
        #plot the scatter plot of the individual particle # use the colormap and the bounds which where defined before
        im = axes[i_ax].scatter(particle_dic["diam"],particle_dic[prop],s=1,c=particle_dic["N_ratio"],rasterized=True,norm=normN,cmap=cmapN,marker='o')
        
        #show legend
        axes[i_ax].legend() 
        
        #make labels
        axes[i_ax].set_xlabel("diameter D / m")
        axes[i_ax].set_ylabel(prop + " / m s-1" ) #TODO: plot also the untis of these properties

        #change the axis
        axes[i_ax].set_xlim([10**low_diam_log_detailed,10**high_diam_log_detailed]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
        axes[i_ax].set_ylim([0,2.5]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
        axes[i_ax].grid(which="major")

        #add colorbar
        cbar = fig.colorbar(im,ax=axes[i_ax])
        cbar.set_label(r"$N_{initial}$/$N_{mono}$")
        #shift ticks to the middle of the color
        tick_locs = (boundsN[:-1:3]) + 0.5*(boundsN[1::3]-boundsN[:-1:3])
        cbar.set_ticks(tick_locs)
        cbar.set_ticklabels(N_mono_list[::3])
        '''
        #plot v-D corresponding to some selected monomer numbers
        for i_Nmono,N_mono_now in enumerate(Nmono_fit_eval):
            #if N_mono_now in [1,2,5,10,100,1000]:
            #calculate vterm from the m/A-D relations and the corresponding fall speeds
            vterm_fit = __fallspeed_relations.calc_vterm(prop[6:],mass_fit_detailed[i_Nmono,:],diam,area_fit_detailed[i_Nmono,:])
            if i_Nmono==0 and i_prop==0:
                print "area_fitted",area_fit_detailed[i_Nmono,:]
            #plot the lines (with a white surrounding)
            fitlines_white = axes[i_ax].semilogx(diam,vterm_fit,c="white",linestyle="-",label="__None",linewidth=linewidth_white)
            fitlines = axes[i_ax].semilogx(diam,vterm_fit,c=usecolorsN[N_mono_now==N_mono_list][0],linestyle="-",label="__None",linewidth=linewidth)
        
        #add a lines for the aggregates
        vterm_fit_allagg = __fallspeed_relations.calc_vterm(prop[6:],fit_dic["mass_Nmono_allagg"],fit_dic["diam"],fit_dic["area_Nmono_allagg"])

        fitlines_monomers = axes[i_ax].semilogx(np.nan,np.nan,c="b",linestyle="-",label="monomers",linewidth=2.*linewidth)
        #fitlines_allagg = axes[i_ax].semilogx(fit_dic["diam"],vterm_fit_allagg,c="w",linestyle="-",label="__None",linewidth=2.*linewidth_white)
        fitlines_allagg = axes[i_ax].semilogx(fit_dic["diam"],vterm_fit_allagg,c="g",label="all aggregates",linewidth=2.*linewidth,linestyle='-.')
        '''
    '''
    #add legends
    for ax in axes: #,7,8,9,10]:
        ax.legend(loc="lower right")
    '''

    #save the plot (and open it)
    plt.tight_layout()
    dir_save = '/home/mkarrer/Dokumente/plots/Jagg/'
    if not os.path.exists(dir_save): #create direktory if it does not exists
        os.makedirs(dir_save)
    if tumbling:
        tumbling_add="_wtumbling"
    else: 
        tumbling_add="_wotumbling"


    out_filestring = "mAVvsDfits_" + particle_type + '_' + tumbling_add
    
    plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
    plt.savefig(dir_save + out_filestring + '.png', dpi=100)
    print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
    subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
    #plt.clf()
    #plt.close()
    '''
    # Save just the portion _inside_ the second axis's boundaries
    save_axis=14
    for i_axis in range(0,len(axes)):#clear all other axes
        if i_axis==save_axis:
            continue
        axes[i_axis].set_xlabel("")
        axes[i_axis].axes.get_xaxis().set_ticklabels([])
    extent = axes[save_axis].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    axes[save_axis].set_ylabel("terminal velocity / m s-1")
    fig.savefig(dir_save + out_filestring + 'spec_ax_figure.png', bbox_inches=extent.expanded(1.5, 1.6),dpi=400)
    fig.savefig(dir_save + out_filestring + 'spec_ax_figure.pdf', bbox_inches=extent.expanded(1.5, 1.6),dpi=400)
    '''
    #####################
    ######END: plotting figure1
    #####################
    #from IPython.core.debugger import Tracer ; Tracer()()
        #####################
    ######START: plotting figure2: ratios m/A**0.25 m/A**0.5 m/A 
    #####################
    
    ##general settings for the plotting
    number_of_plots = 3

    #optimize the appearance of the plot (figure size, fonts)
    [fig2,axes2] = __plotting_functions.proper_font_and_fig_size(number_of_plots,legend_fontsize='medium')
    
    ###add the ratios to the dictionary
    #individual particle
    particle_dic["mABohm"] = particle_dic["mass"]/(a*particle_dic["diam"][match_Nmono]**b)/(particle_dic["area"]/(particle_dic["diam"][match_Nmono]**2*c*particle_dic["diam"][match_Nmono]**d))**0.25
    particle_dic["mAHW10"] = particle_dic["mass"]/(a*particle_dic["diam"][match_Nmono]**b)/(particle_dic["area"]/(particle_dic["diam"][match_Nmono]**2*c*particle_dic["diam"][match_Nmono]**d))**0.5
    particle_dic["mAKC"] = particle_dic["mass"]/(a*particle_dic["diam"][match_Nmono]**b)/(particle_dic["area"]/(particle_dic["diam"][match_Nmono]**2*c*particle_dic["diam"][match_Nmono]**d))
    #fits
    particle_dic["mABohm_fit"] = mass_fit/(a*diam_center**b)/(area_fit/(diam_center**2*c*diam_center**d))**0.25
    particle_dic["mAHW10_fit"] = mass_fit/(a*diam_center**b)/(area_fit/(diam_center**2*c*diam_center**d))**0.5
    particle_dic["mAKC_fit"] = mass_fit/(a*diam_center**b)/(area_fit/(diam_center**2*c*diam_center**d))
    ###
    #START: plot the prop(m or D) vs diameter as scatter and the 2D-fit
    ###
    i_ax=-1
    #do the same plots for mass and area     
    for i_prop,(prop,prop_unit,prop_short) in enumerate(zip(["mABohm","mAHW10","mAKC"],[r"$ 1$",r"$1$",r"$1$"],[r"$\frac{(m/m_{mono})}{(A_r/A_{mono})^{0.25}}$",r"$\frac{(m/m_{mono})}{(A_r/A_{mono})^{0.5}}$",r"$\frac{(m/m_{mono})}{(A_r/A_{mono})}$"])):#,r"$(m/m_{mono})/(A/A_{mono})^{0.5}$","$(m/m_{mono})/(A/A_{mono})"])):
        i_ax+=1 #start at first axis
        print "plotting ",prop
        #plot the individual particle
        axes2[i_ax].set_yscale('log')
        axes2[i_ax].set_xscale('log')
        im = axes2[i_ax].scatter(particle_dic["diam"],particle_dic[prop],s=1,c=particle_dic["N_ratio"],rasterized=True,norm=normN,cmap=cmapN) #all used particles        
        
        #overlay the fit
        for i_Nmono,N_mono_now in enumerate(Nmono_fit_eval): #range(0,N_mono_list.shape[0],5):
            fitlines = axes2[i_ax].plot(diam_center,particle_dic[prop + "_fit"][i_Nmono,:],c="white",linewidth=linewidth_white,linestyle="-")
            fitlines = axes2[i_ax].plot(diam_center,particle_dic[prop + "_fit"][i_Nmono,:],c=usecolorsN[N_mono_now==N_mono_list][0],linewidth=linewidth,linestyle="-")


        #also plot all aggregates in this plot
        #fitlines = axes2[i_ax].plot(np.nan,np.nan,c="b",linewidth=2.*linewidth,label="monomers")
        #fitlines = axes2[i_ax].plot(fit_dic["diam"],fit_dic[prop + "_Nmono_allagg"],c="w",linewidth=2.*linewidth_white,label="__None")
        #fitlines = axes2[i_ax].plot(fit_dic["diam"],fit_dic[prop + "_Nmono_allagg"],c="g",linewidth=2.*linewidth,label="all aggregates",linestyle='-.')

        #make labels
        axes2[i_ax].set_xlabel("diameter D / m")
        axes2[i_ax].set_ylabel((" / ").join((prop_short,prop_unit))) 

        #change the axis
        #axes2[i_ax].set_xlim([1e-4,4e-2])
        #axes2[i_ax].set_ylim([0,np.array([1e-4,1e-3])[i_prop]]) #define the upper limit of the displayed axis
        axes2[i_ax].grid(which="minor")
        #add colorbar
        cbar = fig.colorbar(im,ax=axes2[i_ax])
        cbar.set_label(r"$N_{initial}/N_{mono}$")
        #shift ticks to the middle of the color
        tick_locs = (boundsN[:-1:3]) + 0.5*(boundsN[1::3]-boundsN[:-1:3])
        cbar.set_ticks(tick_locs)
        cbar.set_ticklabels(N_mono_list[::3])
    
    #add legends
    for ax in axes2: #,7,8,9,10]:
        ax.legend(loc="lower right")
   

    #save the plot (and open it)
    plt.tight_layout()
    dir_save = '/home/mkarrer/Dokumente/plots/Jagg/'
    if not os.path.exists(dir_save): #create direktory if it does not exists
        os.makedirs(dir_save)
    if tumbling:
        tumbling_add="_wtumbling"
    else: 
        tumbling_add="_wotumbling"


    out_filestring = "mA_ratios_vsDfits_" + particle_type + '_' + tumbling_add
    
    plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
    plt.savefig(dir_save + out_filestring + '.png', dpi=100)
    print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
    subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])