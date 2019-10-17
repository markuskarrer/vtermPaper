# coding: utf-8
#import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.pylab as pylab
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
#from IPython.core.debugger import Tracer ; Tracer()()
from matplotlib import rc

import compare_fallspeed_parameterizations

'''
this code reads in properties from Jussis aggregate model
and fits different functions to each property (m-D(monomer dependent and binary),A-D(monomer dependent and binary),v-D for all fallspeed models) and displays them together
'''

def read_and_plot(fig,axes,particle_types,powerlawAtlas="Atlas",hydro_model
    ="all"):
    '''
    read and plot the aggregate data along different fits to these data and some microphysics schemes assumptions
    ARGUMENTS:
    fig: figure handle
    axes: axes handles
    particle type: monomer type whcih is processed and plotted
    powerlawAtlas: which fits should be overplotted [powerlaw,Atlas,schemes]
    '''

    tumbling=False 
    take_mono_prop_theor = True #so far only True is implemented

    rhow=1000. #density of water

    #show which types of fit
    show_lines=dict()
    #powerlawAtlas="schemes" #"Atlas","powerlaw" "schemes" #this is defined in __main__
    if powerlawAtlas=="powerlaw":
        show_lines["powerlaw_fits"] = True
        show_lines["highest_diam_for_fit"] = 0.01 #-999 for all -2 for increasing v range only
        show_lines["Atlas_fit"] = False #show the Atlas-type fit for cloud ice & snow
        #from previous model settings
        show_lines["SB_powerlaw"] = False #show also the old fit of Axels power-law    
        show_lines["selected_schemes"] = False

    elif powerlawAtlas=="Atlas":
        show_lines["powerlaw_fits"] = False
        show_lines["highest_diam_for_fit"] = -999 #-999 for all -2 for increasing v range only
        show_lines["Atlas_fit"] = True #show the Atlas-type fit for cloud ice & snow
        #from previous model settings
        show_lines["SB_powerlaw"] = False #show also the old fit of Axels power-law
        show_lines["selected_schemes"] = False

    elif powerlawAtlas=="schemes":
        show_lines["powerlaw_fits"] = False
        show_lines["highest_diam_for_fit"] = -999 #-999 for all -2 for increasing v range only
        show_lines["Atlas_fit"] = False #show the Atlas-type fit for cloud ice & snow
        #from previous model settings
        show_lines["SB_powerlaw"] = False #show also the old fit of Axels power-law
        scheme_dic  = compare_fallspeed_parameterizations.get_model_params()
        show_lines["selected_schemes"] =True


    #add displayed lines to string to distinguish them in the saved files
    add_displayed_lines2string = '' 
    for key in show_lines.keys():
        if show_lines[key]:
            add_displayed_lines2string+= '_' + key
    ##open file for writing fit-parameters
    with open("/home/mkarrer/Dokumente/plots/fit_results_vterm.txt","w") as txtfile: #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in
        fit_param_writer = csv.writer(txtfile, delimiter='&', quoting=csv.QUOTE_NONE, lineterminator=os.linesep, escapechar=" ") #quoting avoids '' for formatted string; lineterminator avoids problems with system dependend lineending format https://unix.stackexchange.com/questions/309154/strings-are-missing-after-concatenating-two-or-more-variable-string-in-bash?answertab=active#tab-top
        fit_param_writer.writerow(["particle_type","am","bm","aA","bA","vterm_Atlas_Dmax_A","vterm_Atlas_Dmax_B","vterm_Atlas_Dmax_C","vterm_Atlas_Deq_A","vterm_Atlas_Deq_B","vterm_Atlas_Deq_C","vterm_pow_Dmax_a","verm_pow_Dmax_b","vterm_pow_m_a","verm_pow_m_b"])
        #particle_types = ["plate"] #["plate","dendrite","mixdendneedle","needle","column","mixcolumndend"] #"mixofdentneed""needle","column","plate","dendrite"] #["needle","column","plate","dendrite","bullet","rosette"] # ,"bullet"rosette
        for particle_type in particle_types:

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
            show_largest_Dmax = False #output largest D_max
            if show_largest_Dmax:
                print "max(Dmax[Nmono<=10])",np.max(particle_dic["diam"][particle_dic["N_monomer"]<=10])
                print "max(Dmax[Nmono<=100])",np.max(particle_dic["diam"][particle_dic["N_monomer"]<=100])
                print "max(Dmax)",np.max(particle_dic["diam"])
                raw_input()
                continue
            #randomize order in dictionary (for an unbiased view in scatter-plot)
            N_particles = particle_dic["area"].shape[0]
            particle_indices_array = range(0,N_particles)
            random.shuffle(particle_indices_array) #shuffle the indice array (random order)
            for key in particle_dic.keys():
                if any(particle_dic[key]): #an empty dictionary will return false
                    particle_dic[key] = particle_dic[key][particle_indices_array] #shuffle the arrays

            #from IPython.core.debugger import Tracer ; Tracer()()
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
            low_diam_log_detailed= -4; high_diam_log_detailed_plotrange=-1.3979400086720375
            if show_lines["highest_diam_for_fit"]==-999:
                high_diam_log_detailed=-1.3979400086720375 #np.log10(3e-2)=-1.5228787452803376 #np.log10(4e-2)=-1.3979400086720375
            else:
                #from IPython.core.debugger import Tracer ; Tracer()()
                high_diam_log_detailed=np.log10(show_lines["highest_diam_for_fit"])

            diam = np.logspace(low_diam_log_detailed,high_diam_log_detailed,500) #set up array of diameters (only the range to which the plot is adjusted)
            diam_full = np.logspace(low_diam_log_detailed,high_diam_log_detailed_plotrange,500) #set up array of diameters 

            diam_log_detailed = np.log10(diam)
            diam_center = (diam[:-1]+diam[1:])/2 #center of the diameter bins
            diam_log = np.log10(diam_center)
            diam_log_center = 10**((np.log10(diam[:-1])+np.log10(diam[1:]))/2) #logarithmic center of the diameter bins
            #define at which monomer numbers the fit should be evaluated
            Nmono_fit_eval = np.array([1,2,5,10,100,1000]) 
            Nmono_log = np.log10(Nmono_fit_eval)

            
            #pass already the diameter to the fit-array
            fit_dic["diam"] = diam
            fit_dic["diam_full"] = diam_full
            fit_dic["diam_wide"] = np.logspace(-5,-1,1000) #for plotting equivalent diameter we need a wider diameter range
                
            ###
            #START: monomer number dependent fitting 
            ###

            #fit the m-D and A-D properties
            print "binary m-D and A-D fitting"
            for i_prop,(prop,prop_unit) in enumerate(zip(["mass","area"],["kg","m2"])): #,"HWJussi","KCJussi"]):
            
                
                ###
                #START: get fits for all aggregates
                ###
                
                ##
                # make a fit for all aggregates in N_mono_list
                ##
                #fit to m-D and A-D coefficients for all particles with N_monomer>1
                [fit_dic[prop + "_coeff_Nmono_allagg"],covar] = __tools_for_processing_Jagg.fit_data(particle_dic["diam"],particle_dic[prop],func='powerlaw',weight="None")
                
                if take_mono_prop_theor: #overwrite fit for monomers with theoretical value
                    if prop=="mass":
                        fit_dic['mass_coeff_Nmono_1'] = __tools_for_processing_Jagg.calc_mD_AD_coeffs(particle_type)[0:2]
                    if prop=="area":
                        fit_dic['area_coeff_Nmono_1'] = __tools_for_processing_Jagg.calc_mD_AD_coeffs(particle_type)[2:4]
                    fit_dic[prop + "_Nmono_1"] = fit_dic[prop + "_coeff_Nmono_1"][0]*fit_dic["diam"]**fit_dic[prop + "_coeff_Nmono_1"][1] #recalculate arrays of masses and areas 
                #calculate the values corresponding to the spanned diameter from the fit-coefficients
                fit_dic[prop + "_Nmono_allagg"] = fit_dic[prop + "_coeff_Nmono_allagg"][0]*fit_dic["diam"]**fit_dic[prop + "_coeff_Nmono_allagg"][1] #calculate arrays of masses and areas based on the m-D fits
                
                ###
                #END: get fits for all aggregates
                ###
                
            ###
            #START: calculate power-law and Atlas-fit based on the exact v-D line
            ###
            
            #first calculate the exact v-D
            


            
            #loop over different terminal velocity models in which different possible fits (power-law, Atlas type) are calculated
            for i_prop,prop in enumerate(["vterm_bohm","vterm_HW10","vterm_KC05"]): #,"vterm_mitch_heym"]): #,"HWJussi","KCJussi"]):
                print "fitting and plotting: ",prop
                #first calculate the exact v-D
                for i_mono,(str_monomer_agg,corr_cat) in enumerate(zip(["1","allagg"],["cloud ice","snow"])):
                    fit_dic[prop +  "_fitted_via_mD_AD_Nmono" + str_monomer_agg] = __fallspeed_relations.calc_vterm(prop[6:],fit_dic["mass_Nmono_" + str_monomer_agg],fit_dic["diam"],fit_dic["area_Nmono_" + str_monomer_agg])

                #fit an Atlas type to the "mean v-D" relations (derived from the A-D and m-D fits)
                for i,(str_monomer_agg,corr_cat) in enumerate(zip(["1","allagg"],["$N_{mono}$=1","$N_{mono}$>1"])): 
                    
                    ##Atlas-type
                    [fitresult,covar] = __tools_for_processing_Jagg.fit_data(fit_dic["diam"],fit_dic[prop +  "_fitted_via_mD_AD_Nmono" + str_monomer_agg],func="Atlas") #, weight=fit_dic["dens" + "_Nmono_" + str_monomer_agg])
                    #from IPython.core.debugger import Tracer ; Tracer()()


                    [fitresult_Deq,covar_Deq] = __tools_for_processing_Jagg.fit_data(((6.*fit_dic["mass_coeff_Nmono_" + str_monomer_agg][0]*fit_dic["diam"]**fit_dic["mass_coeff_Nmono_" + str_monomer_agg][1]/(np.pi*rhow))**(1./3.)),fit_dic[prop +  "_fitted_via_mD_AD_Nmono" + str_monomer_agg],func="Atlas") #, weight=fit_dic["dens" + "_Nmono_" + str_monomer_agg])
                    
                    fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg] = fitresult #copy the Atlas-fit (as a function of D_max) coefficients in the fit_dic dictionary
                    fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_Deq"] = fitresult_Deq #copy the Atlas-fit (as a function of D_eq) coefficients in the fit_dic dictionary

                    #recalculate v(D)= A-B*exp(-gam*D) from fit (Dmax)
                    fit_dic[prop + "_Nmono_" + str_monomer_agg] = fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg][0]-fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg][1]*np.exp(-fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg][2]*fit_dic["diam_full"]) #calculate arrays of masses and areas based on the m-D fits
                    #recalculate v(D)= A-B*exp(-gam*D) from fit (Dmax)
                    fit_dic[prop + "_Nmono_" + str_monomer_agg + "_Deq"] = fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_Deq"][0]-fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_Deq"][1]*np.exp(-fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_Deq"][2]*fit_dic["diam_wide"]) #calculate arrays of masses and areas based on the m-D fits
                    
                    #power-law
                    [fitresult,covar] = __tools_for_processing_Jagg.fit_data(fit_dic["diam"],fit_dic[prop +  "_fitted_via_mD_AD_Nmono" + str_monomer_agg],func="powerlaw") #, weight=fit_dic["dens" + "_Nmono_" + str_monomer_agg])

                    fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_powerlaw"] = fitresult #copy the power-law coefficients in the fit_dic dictionary

                    #recalculate powerlaw from fit
                    fit_dic[prop + "_Nmono_" + str_monomer_agg + "_powerlaw"] = fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_powerlaw"][0]*fit_dic["diam_full"]**fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_powerlaw"][1] #calculate arrays of masses and areas based on the m-D fits
                
                        
            ###
            #END: calculate power-law and Atlas-fit based on the exact v-D line
            ###
                
            #####################
            ######END: fitting
            #####################
            
            #####################
            ######START: plotting
            #####################
            

            i_ax=-1
            linewidth=1.0
            linewidth_white=1.3
            ##define colormaps
            #self-defined discrete colormap (for the monomer number N)
            # define the colormap and the discrete boundaries (defined by the considered monomer numbers) and set up an array of the same colors (for line-plotting)
            cmapN = plt.cm.gray #brg
            #modify lowest color in colormap
            cmaplistN = [cmapN(i) for i in range(cmapN.N)]
            cmaplistN = cmaplistN[60:256-30] #crop colormap to avoid too dark/light scatterer #out of 256 colors

            #colors for the binary relevant
            colormonomer = colors.ColorConverter.to_rgb("b") #(1.0,0.0,0.0,1.0) #make Nmono=1 blue
            colorallagg =  colors.ColorConverter.to_rgb("g") #(1.0,0.0,0.0,1.0) #make Nmono=1 blue
            
            #span the colormaps
            cmapN = colors.LinearSegmentedColormap.from_list('mcm',cmaplistN, cmapN.N)
            boundsN = np.append(N_mono_list,[9999])
            normN = colors.BoundaryNorm(boundsN, cmapN.N)
            usecolorsN = cmapN(np.linspace(0,1,N_mono_list.shape[0])) #get colors from colormap to make consistent line colors for the monomer number

            ###
            #START: plot the prop(m or D) vs diameter as scatter and the 2D-fit
            ###
            
            
            ####
            #plot vterm and its fits
            ####
            for i_prop,prop in enumerate(["vterm_bohm","vterm_HW10","vterm_KC05"]): #,"vterm_mitch_heym"]): #,"HWJussi","KCJussi"]):
                print "fitting and plotting: ",prop
                if hydro_model!=all and hydro_model!=prop:
                    continue
                #define current axis
                i_ax+=1         
                #change the scale (according to xscale_vec and yscale_vec)
                axes[i_ax].set_xscale("log") #axes[i_ax].set_xscale(xscale_vec[i_prop]) 
                axes[i_ax].set_yscale("linear") #axes[i_ax].set_yscale(yscale_vec[i_prop])
                
                #####
                #1. the velocity of each individual particle is plotted 
                #####
                #plot the scatter plot of the individual particle # use the colormap and the bounds which where defined before
                im = axes[i_ax].scatter(particle_dic["diam"],particle_dic[prop],s=1,c=particle_dic["N_monomer"],rasterized=True,norm=normN,cmap=cmapN,marker='o')
                
                #show legend
                axes[i_ax].legend() 
                
                #make labels
                axes[i_ax].set_xlabel("diameter $D_{max}$ / m")
                axes[i_ax].set_ylabel(prop + " / m s-1" ) #TODO: plot also the untis of these properties

                #change the axis
                axes[i_ax].set_xlim([10**low_diam_log_detailed,10**high_diam_log_detailed_plotrange]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                axes[i_ax].set_ylim([0,2.5]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                axes[i_ax].grid(which="major")

                #add colorbar
                cbar = fig.colorbar(im,ax=axes[i_ax])
                cbar.set_label(r"$N_{mono}$")
                #shift ticks to the middle of the color
                tick_locs = (boundsN[:-1:3]) + 0.5*(boundsN[1::3]-boundsN[:-1:3])
                cbar.set_ticks(tick_locs)
                cbar.set_ticklabels(N_mono_list[::3])
                if not show_lines["selected_schemes"]:

                    #add a line for the monomers
                    fitlines_allagg = axes[i_ax].semilogx(fit_dic["diam"],fit_dic[prop +  "_fitted_via_mD_AD_Nmono1"],c="b",linestyle="-",label="$N_{mono}$=1",linewidth=linewidth)

                    #add a line for the aggregates
                    fitlines_allagg = axes[i_ax].semilogx(fit_dic["diam"],fit_dic[prop +  "_fitted_via_mD_AD_Nmonoallagg"],c="g",linestyle="-",label="$N_{mono}$>1",linewidth=linewidth)

                ###add fitted lines
                ##monomers
                if show_lines["Atlas_fit"]:
                    #Atlas
                    #Atlasfit_monomers, = axes[i_ax].semilogx(fit_dic["diam_full"],fit_dic[prop + "_Nmono_1"],c="deepskyblue",linestyle="--",label="Atlas-type Nmono=1 (Dmax)",linewidth=linewidth)
                    Atlasfit_monomersDeq = axes[i_ax].semilogx((np.pi*rhow*fit_dic["diam_wide"]**3/(6*fit_dic["mass_coeff_Nmono_1"][0]))**(1./fit_dic["mass_coeff_Nmono_1"][1]),fit_dic[prop +  "_Nmono_1_Deq"],c="b",linestyle="--",label="$N_{mono}$=1 fit",linewidth=linewidth)
                    

                if show_lines["powerlaw_fits"]:
                    #power-law
                    powerlaw_monomers, = axes[i_ax].semilogx(fit_dic["diam_full"], fit_dic[prop + "_Nmono_1_powerlaw"],marker='None',c="b",linestyle="-.",label="$N_{mono}$=1 fit",linewidth=linewidth)
                    

                
                ##all aggregates
                if show_lines["Atlas_fit"]:
                    #Atlas
                    #Atlasfit_Nmonoallagg, = axes[i_ax].semilogx(fit_dic["diam_full"],fit_dic[prop + "_Nmono_allagg"],c="limegreen",linestyle="--",label="Atlas-type " + corr_cat + " (Dmax)" ,linewidth=linewidth)
                    Atlasfit_NmonoallaggDeq = axes[i_ax].semilogx((np.pi*rhow*fit_dic["diam_wide"]**3/(6*fit_dic["mass_coeff_Nmono_allagg"][0]))**(1./fit_dic["mass_coeff_Nmono_allagg"][1]),fit_dic[prop +  "_Nmono_allagg" + "_Deq"],c="g",linestyle="--",label="$N_{mono}$>1 fit",linewidth=linewidth)

                if show_lines["powerlaw_fits"]:
                    #power-law
                    powerlaw_handle, = axes[i_ax].semilogx(fit_dic["diam_full"], fit_dic[prop + "_Nmono_allagg" + "_powerlaw"],marker='None',c="g",linestyle="-.",label="$N_{mono}>1$ fit",linewidth=linewidth)
                if show_lines["selected_schemes"]:
                    #SB
                    SB_ci_handle,   = axes[i_ax].plot(scheme_dic["diam"],scheme_dic["mDADvD_dict_SBcloudice"]["v(D_array)"],color='blue', label="SB cloud ice",linestyle=':',linewidth=linewidth)
                    SB_snow_handle, = axes[i_ax].plot(scheme_dic["diam"],scheme_dic["mDADvD_dict_SBsnow"]["v(D_array)"]    ,color='green',label="SB snow",     linestyle=':',linewidth=linewidth)
                    
                    #P3
                    D_min_P3=20e-6     #see do jj = 1,1000 ...d1 = real(jj)*20.*1.e-6 - 10.e-6 in create_p3_lookupTable_1.f90
                    D_max_P3=20.*1000.*1e-6       #see do jj = 1,1000 ...d1 = real(jj)*20.*1.e-6 - 10.e-6 in create_p3_lookupTable_1.f90
                    scheme_dic["mDADvD_dict_P3"]["v_mitch_heym(D_array)"] = np.where(scheme_dic["diam"]<D_min_P3,np.nan,scheme_dic["mDADvD_dict_P3"]["v_mitch_heym(D_array)"])
                    scheme_dic["mDADvD_dict_P3"]["v_mitch_heym(D_array)"] = np.where(scheme_dic["diam"]>D_max_P3,np.nan,scheme_dic["mDADvD_dict_P3"]["v_mitch_heym(D_array)"])
                    
                    P3_handle, = axes[i_ax].plot(scheme_dic["diam"],scheme_dic["mDADvD_dict_P3"]["v_mitch_heym(D_array)"],color='orange',label="P3 unrimed",linestyle='-',linewidth=linewidth) #mitch heym is hardcoded because there are no other options yet
                    
                    #Morrison
                    Morr_ci_handle,   = axes[i_ax].plot(scheme_dic["diam"],scheme_dic["morr2mom_cloudice"]["v(D_array)"],color='blue', label="Morr cloud ice",linestyle='--',linewidth=linewidth)
                    Morr_snow_handle, = axes[i_ax].plot(scheme_dic["diam"],scheme_dic["morr2mom_snow"]["v(D_array)"]    ,color='green',label="Morr snow",     linestyle='--',linewidth=linewidth)
                ####
                #START: get model parameters
                #for comparing with models 
                ####

                #SB
                mDADvD_dict_SBcloudice = __setup_mDAD_frommodels.get_model_mDADs(model="SBcloudice",verbose=False)
                mDADvD_dict_SBsnow = __setup_mDAD_frommodels.get_model_mDADs(model="SBsnow",verbose=False)
                #get the m,A,v-array corresponding to diam and save them in separate dictionaries
                for dict_now in (mDADvD_dict_SBcloudice,mDADvD_dict_SBsnow):
                    dict_now = __setup_mDAD_frommodels.calc_area_mass_vterm_arrays(fit_dic["diam_full"],dict_now)

                ####
                #END: get model parameters
                #for comparing with models 
                ####
                
                #plot the 
                if show_lines["SB_powerlaw"]: #the SB dont have to specify the area
                    SB_cloudicehandle, = axes[i_ax].semilogx(fit_dic["diam_full"],mDADvD_dict_SBcloudice["v(D_array)"],color='b',label="SB cloud ice",linestyle=':')
                    SB_snowhandle, = axes[i_ax].semilogx(fit_dic["diam_full"],mDADvD_dict_SBsnow["v(D_array)"],color='g',label="SB snow",linestyle=':')
                #add legends
                for ax in axes: #,7,8,9,10]:
                    ax.legend(ncol=1,bbox_to_anchor=[0.02, 0.9], loc='upper left')
        

            #save the plot (and open it)
            plt.tight_layout()
            dir_save = '/home/mkarrer/Dokumente/plots/Jagg/'
            if not os.path.exists(dir_save): #create direktory if it does not exists
                os.makedirs(dir_save)
            if tumbling:
                tumbling_add="_wtumbling"
            else: 
                tumbling_add="_wotumbling"


            out_filestring = "vvsDfits_" + particle_type + '_' + tumbling_add + '_' + add_displayed_lines2string
            
            plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
            plt.savefig(dir_save + out_filestring + '.png', dpi=100)
            print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
            subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
            
            for i_prop,prop in enumerate(["vterm_bohm"]): #,"vterm_HW10","vterm_KC05"]): #,"vterm_mitch_heym"]): #,"HWJussi","KCJussi"]):

                ###
                #calculate the relations which can be put into the SB-scheme and display it here
                ###
                print "####################"
                print "#coefficients for SB derived from the here investigated particles (hydro_model:" + prop + ")"
                print "####################"
                #cloud ice
                print "ice: m(D):a= ",fit_dic["mass_coeff_Nmono_1"][0]," b= ",fit_dic["mass_coeff_Nmono_1"][1] ," v(D):a= ",fit_dic[prop + "_coeff_Nmono_1_powerlaw"][0]," b= ",fit_dic[prop + "_coeff_Nmono_1_powerlaw"][1]
                print "ice: Atlas-type (D_max): ",fit_dic[prop + "_coeff_Nmono_1"]
                print "ice: A(D):a= ",fit_dic["area_coeff_Nmono_1"][0]," b= ",fit_dic["area_coeff_Nmono_1"][1]
                
                #calculate the D(m) and v(m) coefficients
                fit_dic["mass_coeff_Nmono_1_N(m)"] = __postprocess_SB.convert_ND_to_Nm_from_coeff(fit_dic["mass_coeff_Nmono_1"][0],fit_dic["mass_coeff_Nmono_1"][1],fit_dic[prop + "_coeff_Nmono_1_powerlaw"][0],fit_dic[prop + "_coeff_Nmono_1_powerlaw"][1])
                print "ice: D(m):a= ", fit_dic["mass_coeff_Nmono_1_N(m)"][0]," b= ",fit_dic["mass_coeff_Nmono_1_N(m)"][1], " v(m):a= ",fit_dic["mass_coeff_Nmono_1_N(m)"][2]," b= ",fit_dic["mass_coeff_Nmono_1_N(m)"][3]
                print "ice: Atlas-type (D_eq): ",fit_dic[prop + "_coeff_Nmono_1_Deq"]

                #snow
                print "snow: m(D):a= ", fit_dic["mass_coeff_Nmono_allagg"][0]," b= ",fit_dic["mass_coeff_Nmono_allagg"][1], " v(D):a= ",fit_dic[prop + "_coeff_Nmono_allagg_powerlaw"][0]," b= ",fit_dic[prop + "_coeff_Nmono_allagg_powerlaw"][1]
                print "snow: Atlas-type (D_max): ",fit_dic[prop + "_coeff_Nmono_allagg"]
                print "snow: A(D):a= ",fit_dic["area_coeff_Nmono_allagg"][0]," b= ",fit_dic["area_coeff_Nmono_allagg"][1]

                #calculate the D(m) and v(m) coefficients
                fit_dic["mass_coeff_Nmono_allagg_N(m)"] = __postprocess_SB.convert_ND_to_Nm_from_coeff(fit_dic["mass_coeff_Nmono_allagg"][0],fit_dic["mass_coeff_Nmono_allagg"][1],fit_dic[prop + "_coeff_Nmono_allagg_powerlaw"][0],fit_dic[prop + "_coeff_Nmono_allagg_powerlaw"][1])
                print "snow: D(m):a= ", fit_dic["mass_coeff_Nmono_allagg_N(m)"][0]," b= ",fit_dic["mass_coeff_Nmono_allagg_N(m)"][1], " v(m):a= ",fit_dic["mass_coeff_Nmono_allagg_N(m)"][2]," b= ",fit_dic["mass_coeff_Nmono_allagg_N(m)"][3]
                print "snow: Atlas-type (D_eq): ",fit_dic[prop + "_coeff_Nmono_allagg_Deq"]
                ##write fit-parameter to txt-file
                #convert arrays first to strings
                if prop=="vterm_bohm":

                    for coeff_array in ["mass_coeff_Nmono_allagg","area_coeff_Nmono_allagg","mass_coeff_Nmono_1","area_coeff_Nmono_1",prop + "_coeff_Nmono_1",prop + "_coeff_Nmono_1_Deq",prop + "_coeff_Nmono_1_powerlaw","mass_coeff_Nmono_allagg","area_coeff_Nmono_allagg",prop + "_coeff_Nmono_allagg",prop + "_coeff_Nmono_allagg_Deq",prop + "_coeff_Nmono_allagg_powerlaw"]:
                        print coeff_array
                        fit_dic[coeff_array + "_str_prec"] = ["{:8.6f}".format(tmp_str) for tmp_str in fit_dic[coeff_array][:]]
                        fit_dic[coeff_array + "_str"] = ["{:4.3f}".format(tmp_str) for tmp_str in fit_dic[coeff_array][:]]
                    for i,(str_monomer_agg,fullstring) in enumerate(zip(["1","allagg"],["monomers","aggregates"])): 
                        
                        #more digits
                        #["particle_type","am","bm","aA","bA","vterm_Atlas_Dmax_A","vterm_Atlas_Dmax_B","vterm_Atlas_Dmax_C","vterm_Atlas_Deq_A","vterm_Atlas_Deq_B","vterm_Atlas_Deq_C","vterm_pow_Dmax_a","verm_pow_Dmax_b","vterm_pow_m_a","verm_pow_m_b"]
                        fit_param_writer.writerow("")
                        fit_param_writer.writerow([fullstring])
                        fit_param_writer.writerow([particle_type + "_prec",
                                                fit_dic["mass_coeff_Nmono_" + str_monomer_agg + "_str_prec"][0],fit_dic["mass_coeff_Nmono_" + str_monomer_agg + "_str_prec"][1], #m(D)
                                                fit_dic["area_coeff_Nmono_" + str_monomer_agg + "_str_prec"][0],fit_dic["area_coeff_Nmono_" + str_monomer_agg + "_str_prec"][1], #A(D)
                                                fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg +  "_str_prec"][0],fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg+ "_str_prec"][1],fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg+ "_str_prec"][2], #Atlas Dmax
                                                fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_Deq" + "_str_prec"][0],fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_Deq"+ "_str_prec"][1],fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_Deq" + "_str_prec"][2], #Atlas Deq
                                                fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_powerlaw"+ "_str_prec"][0],fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_powerlaw"][1]]) #power law Dmax
                        #less digits for Table
                        fit_param_writer.writerow([particle_type,
                                                fit_dic["mass_coeff_Nmono_" + str_monomer_agg + "_str"][0],fit_dic["mass_coeff_Nmono_" + str_monomer_agg + "_str"][1], #m(D)
                                                fit_dic["area_coeff_Nmono_" + str_monomer_agg + "_str"][0],fit_dic["area_coeff_Nmono_" + str_monomer_agg + "_str"][1], #A(D)
                                                fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_str"][0],fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_str"][1],fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_str"][2], #Atlas Dmax
                                                    fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_Deq" + "_str"][0],fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_Deq" + "_str"][1],fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_Deq" + "_str"][2], #Atlas Deq
                                                #Atlas Deq
                                                fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_powerlaw" + "_str"][0],fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_powerlaw" + "_str"][1]]) #power law
                    '''
                    #aggregates
                    #more digits
                    fit_param_writer.writerow("")
                    fit_param_writer.writerow(["aggregates"])
                    fit_param_writer.writerow([particle_type + "_prec",fit_dic["mass_coeff_Nmono_allagg" + "_str_prec"][0],fit_dic["mass_coeff_Nmono_allagg" + "_str_prec"][1],fit_dic["area_coeff_Nmono_allagg" + "_str_prec"][0],fit_dic["area_coeff_Nmono_allagg" + "_str_prec"][1],fit_dic[prop + "_coeff_Nmono_allagg" + "_str_prec"][0],fit_dic[prop + "_coeff_Nmono_allagg"+ "_str_prec"][1],fit_dic[prop + "_coeff_Nmono_allagg"+ "_str_prec"][2],fit_dic[prop + "_coeff_Nmono_allagg_powerlaw"+ "_str_prec"][0],fit_dic[prop + "_coeff_Nmono_allagg_powerlaw"][1]])
                    #less digits for Table
                    fit_param_writer.writerow([particle_type,fit_dic["mass_coeff_Nmono_allagg" + "_str"][0],fit_dic["mass_coeff_Nmono_allagg" + "_str"][1],fit_dic["area_coeff_Nmono_allagg" + "_str"][0],fit_dic["area_coeff_Nmono_allagg" + "_str"][1],fit_dic[prop + "_coeff_Nmono_allagg" + "_str"][0],fit_dic[prop + "_coeff_Nmono_allagg" + "_str"][1],fit_dic[prop + "_coeff_Nmono_allagg" + "_str"][2],fit_dic[prop + "_coeff_Nmono_allagg_powerlaw" + "_str"][0],fit_dic[prop + "_coeff_Nmono_allagg_powerlaw" + "_str"][1]])
                    #fit_param_writer.writerow([particle_type,"{:4.2f}".format(a),'&'.join(fit_dic["mass_coeff_mod_powerlaw" + "_str"][0:2]),"{:4.2f}".format(b),'&'.join(fit_dic["mass_coeff_mod_powerlaw" + "_str"][2:]),"{:4.2f}".format(c),'&'.join(fit_dic["area_coeff_mod_powerlaw" + "_str"][0:2]),"{:4.2f}".format(d),'&'.join(fit_dic["area_coeff_mod_powerlaw" + "_str"][2:]),'&'.join(fit_dic["mass_coeff_Nmono_allagg" + "_str"]),'&'.join(fit_dic["area_coeff_Nmono_allagg" + "_str"])])
                    '''
                #if i_prop==0:
                #    raw_input("wait to check or copy coefficients")
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
            ######END: plotting
            #####################
            #from IPython.core.debugger import Tracer ; Tracer()()
            if particle_type=="plate":
                ax_description_array=["bohm","HW10","KC05"]
                for save_axis,ax_description in enumerate(ax_description_array):
                    if ax_description=="bohm":
                        print "saving subfigure",save_axis,ax_description
                        try: #clear label of axis before and recover current one
                            axes[save_axis].set_xlabel("diameter $D_{max}$ / m")
                            axes[save_axis-1].set_xlabel("")
                        except:
                            pass
                        extent = axes[save_axis].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
                        #from IPython.core.debugger import Tracer ; Tracer()()

                        fig.savefig('/home/mkarrer/Dokumente/plots/tmp.pdf',bbox_inches=extent.expanded(1.6, 1.45),dpi=400) #(1.6, 1.4) are width and height expanded around center

                        subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/tmp.pdf /home/mkarrer/Dokumente/plots/4paper/' + out_filestring + '_' + ax_description + '.pdf',shell=True)
    return axes
                        
if __name__ == '__main__':
    
    ##general settings for the plotting
    number_of_plots = 3

    #optimize the appearance of the plot (figure size, fonts)
    [fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots,legend_fontsize='medium')
    fig,axes = read_and_plot(fig,axes,["plate"],powerlawAtlas="powerlaw",hydro_model
    ="all")