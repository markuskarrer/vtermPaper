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
tumbling=False #so far only the use of non-tumbling projected area is implemented here
show_Nmono_fits = True #show fit lines for some selected monomer numbers
v_small_notuse = 0.0  #ATTENTION: particle smaller than .5ms-1 are not used #ATTENTION: applied to B?hm only


def read_and_plot(fig,axes,particle_types,MAE=None,txtfile="/home/mkarrer/Dokumente/plots/fit_results.txt",write_fit_params=False,plot_vars=None,called_by_main=False,use_Deq=False,forKamil=False,skip_plotting=False):
    '''
    this performs all the processing and plotting
    
    ARGUMENTS:
    fig: figure handle
    axes: axes handles
    particle type: monomer type whcih is processed and plotted
    MAE: dictionary to save mean absolute error results
    txtfile: here the coefficients are writteAn to
    write_fit_params: should coefficients be written
    plot_vars: which variables should be processes and plotted (fm comes always with mass; fA comes always with area)? (if "None" all are plotted)
    called_by_main: if this script has been executed directly save the figures here
    use_Deq: use mass equivalent mass instead of Dmax
    forKamil: some different format for Kamils closure paper
    skip_plotting: dont plot but only calculate the coefficients
    '''
    if MAE==None:
        MAE= dict()
    #flatten axes in case we have more than one column
    try:
        axes = axes.flat
    except:
        print axes," cannot be flattened any more"
    if use_Deq:
        show_Nmono_fits=False
    else: 
        show_Nmono_fits=True        
    with open(txtfile,"w") as txtfile: #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in

        for i_particle_type,particle_type in enumerate(particle_types):
            if called_by_main:
                ##general settings for the plotting
                number_of_plots = 9 #8

                #optimize the appearance of the plot (figure size, fonts)
                [fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots,legend_fontsize='small')
    
            if write_fit_params:
                fit_param_writer = csv.writer(txtfile, delimiter='&', quoting=csv.QUOTE_NONE, lineterminator=os.linesep, escapechar=" ") #quoting avoids '' for formatted string; lineterminator avoids problems with system dependend lineending format https://unix.stackexchange.com/questions/309154/strings-are-missing-after-concatenating-two-or-more-variable-string-in-bash?answertab=active#tab-top
                fit_param_writer.writerow(["particle_type","am1","am2","bm1","bm2","aA1","aA2","bA1","bA2","am_allagg","bm_allagg","aA_allagg","bA_allagg"])
            rimed=False 
            if forKamil:
                if ("simultaneous" in particle_type) or ("subsequent" in particle_type) or ("rimeonly" in particle_type): #rimed particles
                    prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_ssrga/Jussis_aggregates_rimed/"
                    rimed=True
                else:
                    prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_ssrga/" + "Jussis_aggregates_" + particle_type + "/"

                grid_res_array = [10e-6]
            else:
                prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation/"
                grid_res_array = [5e-6,10e-6]
            #read the properties of the individual particles from the files
            if rimed:
                particle_dic,N_mono_list = __tools_for_processing_Jagg.read_particle_prop_files_rimed(   prop_file_folder, particle_type)
            
            else:
                particle_dic,N_mono_list = __tools_for_processing_Jagg.read_particle_prop_files(
            prop_file_folder = prop_file_folder,#"/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation/",
            #prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation_local/samesizeparam/",
            #prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation_local/diffsizeparam/",
                                                                                    D_small_notuse=1e-4, #ATTENTION: particle smaller than D_small_notuse are not considered (e.g. because of resolution issues)
                                                                                    N_small_notuse=1, #ATTENTION: THIS HAS NO EFFECT!!  monomer numbers smaller than N_small_notuse are not considered
                                                                                    grid_res_array = grid_res_array, #[5e-6,10e-6], #array of grid resolutions (defines the files which are read in)
                                                                                    #grid_res_array = [40e-6], #array of grid resolutions (defines the files which are read in)

                                                                                    particle_type = particle_type, #define the habit
                                                                                    test_with_small_sample = False,
                                                                                    forKamil = forKamil
                                                                                    )
            #show the dictionary (e.g. to check that it's not empty because the path was wrong)
            print particle_dic;
            #particle_dic["diam"] = particle_dic["Dmax_xy_aligned"]
            print particle_dic.keys()
            #particle_dic["diam"] = particle_dic["rad_gyr"]*2.0
            #randomize order in dictionary (for an unbiased view in scatter-plot)
            N_particles = particle_dic["area"].shape[0]
            particle_indices_array = range(0,N_particles)
            random.shuffle(particle_indices_array) #shuffle the indice array (random order)
            for key in particle_dic.keys():
                if len(particle_dic[key])>0 and not forKamil: 
                    particle_dic[key] = particle_dic[key][particle_indices_array] #shuffle the arrays
            #get m-D from the assumptions in the aggregate model (this is used at various places)
            a,b,c,d = __tools_for_processing_Jagg.calc_mD_AD_coeffs(particle_type)

            if tumbling:
                particle_dic["area"] = particle_dic["area_partal40std"] #use the tumbled area instead
            #calculate the terminal velocitys of the simulated particles individually for each available velocity model
            for velocity_model in ["HW10","KC05","bohm"]: #,"mitch_heym"]:
                    particle_dic["vterm_" + velocity_model] = __fallspeed_relations.calc_vterm(velocity_model,particle_dic["mass"],particle_dic["diam"],particle_dic["area"]) #+np.random.randn(particle_dic["mass"].shape[0])*0.4 #test noise on data
            #calculate the area ratio
            particle_dic["Aratio"] = particle_dic["area"]/(np.pi/4*particle_dic["diam"]**2)
            if use_Deq:
                rho_i = 917.6 #bulk ice density
                particle_dic["diam"] = (6.*particle_dic["mass"]/(np.pi*rho_i))**(1./3.)    

            #if : #TODO: do this with a flag for all hydro-models?
            #for i_particle,vterm_particle_now in enumerate(particle_dic["vterm_" + velocity_model]):
            #    if vterm_particle_now
            #####################
            ######START: fitting
            #####################
            #initialize an array which contains all fit results
            fitresult_arr_dir_powerlaw_string = np.zeros((5,N_mono_list.shape[0],2)) #fit results (saved as floats) #dimensions: 0-> [mass,array]; 1-> monomer number; 2-> parameters (a,b) in fit

            fit_dic = dict() #initialize fit dictionary
            #set up array of diameters for displaying the fits
            low_diam_log_detailed= -4; high_diam_log_detailed=-1.3979400086720375
            
            diam = np.logspace(low_diam_log_detailed,high_diam_log_detailed,500) #set up array of diameters (for bin edges and KDE-grid)
            diam_log_detailed = np.log10(diam)
            diam_center = (diam[:-1]+diam[1:])/2 #center of the diameter bins
            diam_log = np.log10(diam_center)
            diam_log_center = 10**((np.log10(diam[:-1])+np.log10(diam[1:]))/2) #logarithmic center of the diameter bins
            #define at which monomer numbers the fit should be evaluated
            Nmono_fit_eval = np.array([1,2,5,10,100,1000]) 
            Nmono_log = np.log10(Nmono_fit_eval)
            
            #pass already the diameter to the fit-array
            fit_dic["diam"] = diam
                
            ###
            #START: monomer number dependent fitting 
            ###

            #fit the m-D and A-D properties
            print "monomer number dependent fitting"
            for i_prop,(prop,prop_unit) in enumerate(zip(["mass","area"],["kg","m2"])): #,"HWJussi","KCJussi"]):

                #fit the modified power-law m(D,N_mono)= a(N_mono)*D**(b(N_mono)* a_0*D**b0
                ab_func="powerlaw_rational1" # "powerlaw_polynom1"powerlaw_rational1" "powerlaw_rational2" and "powerlaw_rational3" need some more advanced initial guesses #select the functional form of a and b 
                if prop=="mass" or (not tumbling): #for mass and area without tumbling we can constrain more parameters
                    if particle_type=="dendrite" and prop=="area":
                        #if prop=="mass":
                        #    guess=[-0.1,0.0,-0.0,0.0]
                        guess=[-0.05,0.0,-0.0,0.0]
                        ftol=1e-2
                    elif particle_type=="rosette" and prop=="area":
                        #if prop=="mass":
                        #    guess=[-0.1,0.0,-0.0,0.0]
                        guess=[0.00,0.0,0.0,0.0]
                        ftol=1e-2
                    elif particle_type=="column" and prop=="mass":
                        #if prop=="mass":
                        #    guess=[-0.1,0.0,-0.0,0.0]
                        guess=[0.00,0.0,0.0,0.0]
                        ftol=1e-2
                    elif particle_type=="needle" and prop=="mass":
                        #if prop=="mass":
                        #    guess=[-0.1,0.0,-0.0,0.0]
                        guess=[0.00,0.0,0.0,0.0]
                        ftol=1e-2
                    else:
                        guess=[0.0,0.0,0.0,0.0]
                        ftol=1e-08
                    if ab_func=="powerlaw_rational1":
                        fitresult_mod_powerlaw = __tools_for_processing_Jagg.fit_2D_rational_function_transformwrapper(particle_dic["diam"][particle_dic["vterm_bohm"]>v_small_notuse],particle_dic["N_monomer"][particle_dic["vterm_bohm"]>v_small_notuse],particle_dic[prop][particle_dic["vterm_bohm"]>v_small_notuse],func='powerlaw_rational1_fixp0_fixr0',method='lm',weight='None',habit=particle_type,prop=prop,guess=guess,ftol=ftol)
                    else: #look for other functions in process_Jussis_aggregate_model_mod_powerlaw
                        print "func: ", ab_func, " not (carefully) implemented in process_Jussis_aggregate_model_plots4paper (and its modules)"; sys.exit()
                elif prop=="area" and tumbling:
                    #fitresult_mod_powerlaw = __tools_for_processing_Jagg.fit_2D_rational_function_transformwrapper(particle_dic["diam"],particle_dic["N_monomer"],particle_dic[prop],func='rational2_2_pade_fixq00',method='lm',weight='None',habit=particle_type,prop=prop)
                    if ab_func=="powerlaw_rational1":
                        fitresult_mod_powerlaw = __tools_for_processing_Jagg.fit_2D_rational_function_transformwrapper(particle_dic["diam"][particle_dic["vterm_bohm"]>v_small_notuse],particle_dic["N_monomer"][particle_dic["vterm_bohm"]>v_small_notuse],particle_dic[prop][particle_dic["vterm_bohm"]>v_small_notuse],func='powerlaw_rational1',method='lm',weight='None',habit=particle_type,prop=prop)
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
                    if prop=="mass":
                        mass_fitted = __tools_for_processing_Jagg.powerlaw_rational1_fixp0_fixr0(D_N_grid,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3])#this is the logscaled property (mass/area) divided by the monomer property
                        mass_fitted_detailed = __tools_for_processing_Jagg.powerlaw_rational1_fixp0_fixr0(D_N_grid_detailed,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3])#this is the logscaled property (mass/area) divided by the monomer property
                        
                    elif prop=="area" and (not tumbling):
                        #ATTENTION: fits are overwritten here
                        #if particle_type=="dendrite":
                        #    fit_dic[prop +"_coeff_mod_powerlaw"] = [-0.8,0.0,-0.1,0.0]
                        area_fitted = __tools_for_processing_Jagg.powerlaw_rational1_fixp0_fixr0(D_N_grid,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3])#this is the logscaled property (mass/area) divided by the monomer property
                        area_fitted_detailed = __tools_for_processing_Jagg.powerlaw_rational1_fixp0_fixr0(D_N_grid_detailed,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3])#this is the logscaled property (mass/area) divided by the monomer property
                    elif prop=="area" and tumbling:
                        area_fitted = __tools_for_processing_Jagg.powerlaw_rational1(D_N_grid,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3],fit_dic[prop +"_coeff_mod_powerlaw"][4],fit_dic[prop +"_coeff_mod_powerlaw"][5])#this is the logscaled property (mass/area) divided by the monomer property
                        area_fitted_detailed = __tools_for_processing_Jagg.powerlaw_rational1(D_N_grid_detailed,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3],fit_dic[prop +"_coeff_mod_powerlaw"][4],fit_dic[prop +"_coeff_mod_powerlaw"][5])#this is the logscaled property (mass/area) divided by the monomer property
                    else: #TODO: include different fitting methods
                            print "func: ", ab_func, " not (carefully) implemented in process_Jussis_aggregate_model_mod_powerlaw (and its modules)"; sys.exit()

                if prop=="mass":
                    mass_fit = 10**(mass_fitted)*a*diam_center**b
                    mass_fit_detailed = 10**(mass_fitted_detailed)*a*diam**b#for plotting the velocity more bins are nicer
                    
                if prop=="area":
                    area_fit = 10**(area_fitted)*c*diam_center**d
                    area_fit_detailed = 10**(area_fitted_detailed)*c*diam**d#for plotting the velocity more bins are nicer

                
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
                [fit_dic[prop + "_coeff_Nmono_allagg"],covar] = __tools_for_processing_Jagg.fit_data(particle_dic["diam"][particle_dic["vterm_bohm"]>v_small_notuse],particle_dic[prop][particle_dic["vterm_bohm"]>v_small_notuse],func='powerlaw',weight="None")
                print "\n",prop,fit_dic[prop + "_coeff_Nmono_allagg"] 
                #if prop=="area":
                #    fit_dic["area_coeff_Nmono_allagg"][0] = np.pi/4 #test: use opaque particles as suggested by Westbrook, PhD thesis
                #    fit_dic["area_coeff_Nmono_allagg"][1] = 2.0 #test: use opaque particles as suggested by Westbrook, PhD thesis

                #calculate the values corresponding to the spanned diameter from the fit-coefficients
                fit_dic[prop + "_Nmono_allagg"] = fit_dic[prop + "_coeff_Nmono_allagg"][0]*fit_dic["diam"]**fit_dic[prop + "_coeff_Nmono_allagg"][1] #calculate arrays of masses and areas based on the m-D fits
                
                ###
                #END: get fits for all aggregates
                ###

                ###
                #START: calculate the mean absolute error                
                ###
                #evaluate the m(D,Nmono) fit for the whole particle_dic to get the "predicted" property
                fit_eval_on_particle_dic = __tools_for_processing_Jagg.powerlaw_rational1_fixp0_fixr0([[np.log10(particle_dic["diam"])],[np.log10(particle_dic["N_monomer"])]],fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3])[0] #this is the logscaled property (mass/area) divided by the monomer property
                particle_dic[prop + "_predicted_monodep"] = 10**(fit_eval_on_particle_dic)
                if prop=="mass":
                    particle_dic[prop + "_predicted_binary"] = fit_dic[prop + "_coeff_Nmono_allagg"][0]*particle_dic["diam"]**fit_dic[prop + "_coeff_Nmono_allagg"][1]/(a*particle_dic["diam"]**b)
                    actual = particle_dic[prop][particle_dic["N_monomer"]>1]/(a*particle_dic["diam"][particle_dic["N_monomer"]>1]**b)
                elif prop=="area":
                    particle_dic[prop + "_predicted_binary"] = fit_dic[prop + "_coeff_Nmono_allagg"][0]*particle_dic["diam"]**fit_dic[prop + "_coeff_Nmono_allagg"][1]/(c*particle_dic["diam"]**d)
                    actual = particle_dic[prop][particle_dic["N_monomer"]>1]/(c*particle_dic["diam"][particle_dic["N_monomer"]>1]**d)

                predicted_monodep = particle_dic[prop + "_predicted_monodep"][particle_dic["N_monomer"]>1]
                predicted_binary = particle_dic[prop + "_predicted_binary"][particle_dic["N_monomer"]>1]
                MAE[particle_type +prop + "_monodep"] = np.mean(np.abs(actual-predicted_monodep))
                MAE[particle_type +prop + "_binary"] = np.mean(np.abs(actual-predicted_binary))
                #print prop,MAE; debug()

                ###
                #END: calculate the mean absolute error                
                ###

            if write_fit_params:
                ##write fit-parameter to txt-file
                #convert arrays first to strings
                
                for coeff_array in ["mass_coeff_Nmono_allagg","mass_coeff_mod_powerlaw","area_coeff_Nmono_allagg","area_coeff_mod_powerlaw"]:
                    fit_dic[coeff_array + "_str_prec"] = ["{:8.6f}".format(tmp_str) for tmp_str in fit_dic[coeff_array][:]]
                    fit_dic[coeff_array + "_str"] = ["{:4.3f}".format(tmp_str) for tmp_str in fit_dic[coeff_array][:]]
                    

                #fit_param_writer.writerow([particle_type + "_prec "," {:10.6f} ".format(a),' &'.join(fit_dic["mass_coeff_mod_powerlaw" + "_str_prec"][0:2])," {:10.6f} ".format(b),' &'.join(fit_dic["mass_coeff_mod_powerlaw" + "_str_prec"][2:])," {:10.6f} ".format(c),' &'.join(fit_dic["area_coeff_mod_powerlaw" + "_str_prec"][0:2])," {:10.6f} ".format(d),' &'.join(fit_dic["area_coeff_mod_powerlaw" + "_str_prec"][2:]),' &'.join(fit_dic["mass_coeff_Nmono_allagg" + "_str_prec"]),' &'.join(fit_dic["area_coeff_Nmono_allagg" + "_str_prec"])])
                fit_param_writer.writerow([particle_type + "_prec","{:8.6f}".format(a),fit_dic["mass_coeff_mod_powerlaw" + "_str_prec"][0],fit_dic["mass_coeff_mod_powerlaw" + "_str_prec"][1],"{:8.6f}".format(b),fit_dic["mass_coeff_mod_powerlaw" + "_str_prec"][2],fit_dic["mass_coeff_mod_powerlaw" + "_str_prec"][3],"{:8.6f}".format(c),fit_dic["area_coeff_mod_powerlaw" + "_str_prec"][0],fit_dic["area_coeff_mod_powerlaw" + "_str_prec"][1],"{:8.6f}".format(d),fit_dic["area_coeff_mod_powerlaw" + "_str_prec"][2],fit_dic["area_coeff_mod_powerlaw" + "_str_prec"][3],fit_dic["mass_coeff_Nmono_allagg" + "_str_prec"][0],fit_dic["mass_coeff_Nmono_allagg" + "_str_prec"][1],fit_dic["area_coeff_Nmono_allagg" + "_str_prec"][0],fit_dic["area_coeff_Nmono_allagg" + "_str_prec"][1]])
                fit_param_writer.writerow([particle_type,'&'.join(fit_dic["mass_coeff_mod_powerlaw" + "_str"][0:2]),'&'.join(fit_dic["mass_coeff_mod_powerlaw" + "_str"][2:]),'&'.join(fit_dic["area_coeff_mod_powerlaw" + "_str"][0:2]),'&'.join(fit_dic["area_coeff_mod_powerlaw" + "_str"][2:]),'&'.join(fit_dic["mass_coeff_Nmono_allagg" + "_str"]),'&'.join(fit_dic["area_coeff_Nmono_allagg" + "_str"])])
                #fit_param_writer.writerow([particle_type,'&'.join(fit_dic["area_coeff_Nmono_allagg" + "_str"]),'&'.join(fit_dic["area_coeff_mod_powerlaw" + "_str"])])
                print ["particle_type","am1","am2","bm1","bm2","aA1","aA2","bA1","bA2","am_allagg","bm_allagg","aA_allagg","bA_allagg"]
                print [particle_type,'&'.join(fit_dic["mass_coeff_mod_powerlaw" + "_str"][0:2]),'&'.join(fit_dic["mass_coeff_mod_powerlaw" + "_str"][2:]),'&'.join(fit_dic["area_coeff_mod_powerlaw" + "_str"][0:2]),'&'.join(fit_dic["area_coeff_mod_powerlaw" + "_str"][2:]),'&'.join(fit_dic["mass_coeff_Nmono_allagg" + "_str"]),'&'.join(fit_dic["area_coeff_Nmono_allagg" + "_str"])]
            #####################
            ######END: fitting
            #####################
             
            #####################
            ######START: plotting figure1
            #####################
            if skip_plotting:
                continue

            i_ax=-1
            linewidth=1.0
            linewidth_white=1.3
            ##define colormaps
            #self-defined discrete colormap (for the monomer number N)
            # define the colormap and the discrete boundaries (defined by the considered monomer numbers) and set up an array of the same colors (for line-plotting)
            cmapN = plt.cm.hot #brg
            #modify lowest color in colormap
            cmaplistN = [cmapN(i) for i in range(cmapN.N)]
            cmaplistN = cmaplistN[30:256-30] #crop colormap to avoid too dark/light scatterer #out of 256 colors

            #colors for the binary relevant
            colormonomer = colors.ColorConverter.to_rgb("b") #(1.0,0.0,0.0,1.0) #make Nmono=1 blue
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
            for i_prop,(prop,prop_unit,prop_short) in enumerate(zip(["mass","area","Aratio"],["kg","m2","1"],["m","A","Ar"])):
                #skip variables which are not requested by plot_vars
                #if plot_vars!=None and not (prop in plot_vars) or prop=="Aratio": 
                if plot_vars!=None and not (prop in plot_vars): 
                    continue
                
                
                i_ax+=1 #start at first axis
                print "fitting ",prop
                #plot the individual particle
                axes[i_ax].set_yscale('log')
                axes[i_ax].set_xscale('log')
                im = axes[i_ax].scatter(particle_dic["diam"][particle_dic["vterm_bohm"]>v_small_notuse],particle_dic[prop][particle_dic["vterm_bohm"]>v_small_notuse],s=1,c=particle_dic["N_monomer"][particle_dic["vterm_bohm"]>v_small_notuse],rasterized=True,norm=normN,cmap=cmapN) #all used particles        
                fitlines = axes[i_ax].semilogx(np.nan,np.nan,linestyle="",label='__None',linewidth=linewidth) #prop_short + "-D power \nlaw fits",linewidth=linewidth)

                if show_Nmono_fits and (not prop=="Aratio"):
                    #overlay the fit
                    for i_Nmono,N_mono_now in enumerate(Nmono_fit_eval): #range(0,N_mono_list.shape[0],5):
                        if prop=="mass":
                            fitlines = axes[i_ax].plot(diam_center,mass_fit[i_Nmono,:],c="white",linewidth=linewidth_white,linestyle="-")
                            fitlines = axes[i_ax].plot(diam_center,mass_fit[i_Nmono,:],c=usecolorsN[N_mono_now==N_mono_list][0],linewidth=linewidth,linestyle="-",label="$N_{mono}=$" + "{}".format(N_mono_now))

                        if prop=="area":
                            #fitlines = axes[i_ax].plot(diam_center,area_fit[i_Nmono,:],c="white",linewidth=linewidth_white,linestyle="-")
                            fitlines = axes[i_ax].plot(diam_center,area_fit[i_Nmono,:],c=usecolorsN[N_mono_now==N_mono_list][0],linewidth=linewidth,linestyle="-",label="$N_{mono}=$" + "{}".format(N_mono_now))
                        if prop=="Aratio":
                            #fitlines = axes[i_ax].plot(diam_center,area_fit[i_Nmono,:],c="white",linewidth=linewidth_white,linestyle="-")
                            fitlines = axes[i_ax].plot(diam_center,area_fit[i_Nmono,:]/(np.pi/4*diam_center**2),c=usecolorsN[N_mono_now==N_mono_list][0],linewidth=linewidth,linestyle="-",label="$N_{mono}=$" + "{}".format(N_mono_now))
                if not prop=="Aratio":
                    #also plot all aggregates in this plot
                    #fitlines = axes[i_ax].plot(np.nan,np.nan,c="b",linewidth=2.*linewidth,label="monomers")
                    #fitlines = axes[i_ax].plot(fit_dic["diam"],fit_dic[prop + "_Nmono_allagg"],c="w",linewidth=2.*linewidth_white,label="__None")
                    fitlines = axes[i_ax].plot(fit_dic["diam"],fit_dic[prop + "_Nmono_allagg"],c="g",linewidth=2.*linewidth,label=r"$N_{mono}$>1",linestyle='--')

                #make labels
                if use_Deq:
                    axes[i_ax].set_xlabel("$D_{eq}$ [m]")
                else:
                    axes[i_ax].set_xlabel("$D_{max}$ [m]")
                if not prop=="Aratio":
                    axes[i_ax].set_ylabel(prop_short + " [" +prop_unit+ "]") 
                else:
                    axes[i_ax].set_ylabel(' arearatio '+ ' ' + prop_short + " [" +prop_unit+ "]") 
                #change the axis
                if use_Deq:
                    axes[i_ax].set_xlim([1e-5,1e-3])
                else:
                    axes[i_ax].set_xlim([1e-4,4e-2])

                axes[i_ax].set_ylim([0,np.array([1e-4,1e-3,1e-0])[i_prop]]) #define the upper limit of the displayed axis
                axes[i_ax].grid(which="both")
                #add colorbar (not when plotting only selected variables)
                if plot_vars==None:
                    cbar = fig.colorbar(im,ax=axes[i_ax])
                    cbar.set_label("$N_{mono}$")
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
                #match_Nmono = (particle_dic["N_monomer"]<=1000)
                if prop=="mass":
                    #im = axes[i_ax].scatter(np.log10(particle_dic["diam"]),np.log10(particle_dic[prop]/(a*particle_dic["diam"]**b)),s=1,c=particle_dic["N_monomer"],rasterized=True,norm=normN,cmap=cmapN) #old with log10(..) as label
                    im = axes[i_ax].scatter(particle_dic["diam"][particle_dic["vterm_bohm"]>v_small_notuse],particle_dic[prop][particle_dic["vterm_bohm"]>v_small_notuse]/(a*particle_dic["diam"][particle_dic["vterm_bohm"]>v_small_notuse]**b),s=1,c=particle_dic["N_monomer"][particle_dic["vterm_bohm"]>v_small_notuse],rasterized=True,norm=normN,cmap=cmapN) 
                    #im = axes[i_ax].scatter(particle_dic["diam"][particle_dic["vterm_bohm"]>v_small_notuse],particle_dic[prop + "_predicted_binary"][particle_dic["vterm_bohm"]>v_small_notuse],s=1,c=particle_dic["N_monomer"][particle_dic["vterm_bohm"]>v_small_notuse],rasterized=True,norm=normN,cmap=cmapN) #to test the MAE implementation

                elif prop=="area":
                    #im = axes[i_ax].scatter(np.log10(particle_dic["diam"]),np.log10(particle_dic[prop]/(c*particle_dic["diam"]**d)),s=1,c=particle_dic["N_monomer"],rasterized=True,norm=normN,cmap=cmapN) #old with log10(..) as label
                    im = axes[i_ax].scatter(particle_dic["diam"][particle_dic["vterm_bohm"]>v_small_notuse],particle_dic[prop][particle_dic["vterm_bohm"]>v_small_notuse]/(c*particle_dic["diam"][particle_dic["vterm_bohm"]>v_small_notuse]**d),s=1,c=particle_dic["N_monomer"][particle_dic["vterm_bohm"]>v_small_notuse],rasterized=True,norm=normN,cmap=cmapN)
                elif prop=="Aratio":
                    #im = axes[i_ax].scatter(np.log10(particle_dic["diam"]),np.log10(particle_dic[prop]/(c*particle_dic["diam"]**d)),s=1,c=particle_dic["N_monomer"],rasterized=True,norm=normN,cmap=cmapN) #old with log10(..) as label
                    im = axes[i_ax].scatter(particle_dic["diam"][particle_dic["vterm_bohm"]>v_small_notuse],particle_dic[prop][particle_dic["vterm_bohm"]>v_small_notuse]/(c*particle_dic["diam"][particle_dic["vterm_bohm"]>v_small_notuse]**d/(np.pi/4*particle_dic["diam"]**2)),s=1,c=particle_dic["N_monomer"][particle_dic["vterm_bohm"]>v_small_notuse],rasterized=True,norm=normN,cmap=cmapN)
                if show_Nmono_fits:
                    fitlines = axes[i_ax].semilogx(np.nan,np.nan,linestyle="",label='__None',linewidth=linewidth) #prop_short + "-D power \nlaw fits",linewidth=linewidth)

                    #overlay the scatterer with the fit at some selected monomer numbers
                    for i_Nmono,N_mono_now in enumerate(Nmono_fit_eval):
                        if prop=="mass":
                            prop_fitted=mass_fitted
                        elif prop=="area":
                            prop_fitted=area_fitted
                        #fitlines = axes[i_ax].plot(np.log10(diam_center),(prop_fitted[i_Nmono,:]),c="white",linewidth=linewidth_white,linestyle="-")
                        fitlines = axes[i_ax].semilogx(diam_center,10**prop_fitted[i_Nmono,:],c="white",linewidth=linewidth_white,linestyle="-") #old with log10(..) as label
                        #fitlines = axes[i_ax].plot(np.log10(diam_center),(prop_fitted[i_Nmono,:]),c=usecolorsN[N_mono_now==N_mono_list][0],linewidth=linewidth,linestyle="-") #old with log10(..) as label
                        fitlines = axes[i_ax].semilogx(diam_center,10**prop_fitted[i_Nmono,:],c=usecolorsN[N_mono_now==N_mono_list][0],linewidth=linewidth,linestyle="-",label="$N_{mono}=$" + "{}".format(N_mono_now))
                    
                #also plot all aggregates in this plot
                if prop=="mass":
                    #fitlines = axes[i_ax].plot(np.log10(fit_dic["diam"]),np.log10(fit_dic[prop + "_Nmono_allagg"]/(a*fit_dic["diam"]**b)),c="g",linewidth=2.*linewidth,label="all aggregates")#old with log10(..) as label
                    #fitlines = axes[i_ax].plot(np.nan,np.nan,c="b",linewidth=linewidth,label="monomers")
                    #fitlines_white = axes[i_ax].semilogx(fit_dic["diam"],fit_dic[prop + "_Nmono_allagg"]/(a*fit_dic["diam"]**b),c="w",linewidth=2.*linewidth_white,label="__None")
                    fitlines = axes[i_ax].semilogx(fit_dic["diam"],fit_dic[prop + "_Nmono_allagg"]/(a*fit_dic["diam"]**b),c="g",linewidth=2.*linewidth,label=r"$N_{mono}$>1",linestyle='--')

                elif prop=="area":
                    #fitlines = axes[i_ax].plot(np.log10(fit_dic["diam"]),np.log10(fit_dic[prop + "_Nmono_allagg"]/(c*fit_dic["diam"]**d)),c="g",linewidth=2.*linewidth,label="all aggregates")#old with log10(..) as label
                    #fitlines_white = axes[i_ax].plot(fit_dic["diam"],fit_dic[prop + "_Nmono_allagg"]/(c*fit_dic["diam"]**d),c="w",linewidth=2.*linewidth,label="__None")
                    #fitlines = axes[i_ax].plot(np.nan,np.nan,c="b",linewidth=linewidth,label="monomers")
                    fitlines = axes[i_ax].plot(fit_dic["diam"],fit_dic[prop + "_Nmono_allagg"]/(c*fit_dic["diam"]**d),c="g",linewidth=2.*linewidth,label="$N_{mono}$>1",linestyle='--')
                
                #set yticks and yticklabels
                #axes[i_ax].set_yticks([0.1,0.2,0.3,0.4,0.5,0.7,1.0,2.0])
                #axes[i_ax].set_yticklabels([0.1,0.2,0.3,0.4,0.5,0.7,1.0,2.0])
                #set xscale to log
                axes[i_ax].set_xscale("log")

                #make labels
                axes[i_ax].set_xlabel("$D_{max}$ [m]")
                if prop=='area':
                    axes[i_ax].set_ylabel("$f_A$")
                elif prop=='mass':
                    axes[i_ax].set_ylabel("$f_m$")

                #change the axis
                if use_Deq:
                    axes[i_ax].set_xlim([1e-5,1e-3])
                else:
                    axes[i_ax].set_xlim([1e-4,4e-2])
                #axes[i_ax].set_ylim([0,np.array([1e-4,1e-3])[i_prop]]) #define the upper limit of the displayed axis
                axes[i_ax].grid(which="both")
                
                #add colorbar (not when plotting only selected variables)
                tick_locs = (boundsN[:-1:3]) + 0.5*(boundsN[1::3]-boundsN[:-1:3])
                if plot_vars==None:
                    cbar = fig.colorbar(im,ax=axes[i_ax])
                    cbar.set_label("$N_{mono}$")
                    #shift ticks to the middle of the color
                    cbar.set_ticks(tick_locs)
                    cbar.set_ticklabels(N_mono_list[::3])
                ###
                #END: plot the prop(m or D) vs diameter as scatter and the 2D-fit (with transformed values D->log10(D) N->log10(N) prop->prop/prop_monomer)
                ###
                
            ####
            #plot vterm and its fits
            ####
            for i_prop,(prop,vmodel_short) in enumerate(zip(["vterm_bohm","vterm_HW10","vterm_KC05"],["B92","HW10","KC05"])): #,"vterm_mitch_heym"]): #,"HWJussi","KCJussi"]):
                print "fitting and plotting: ",prop
                #skip variables which are not requested by plot_vars
                if plot_vars!=None and not (prop in plot_vars): 
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
                im = axes[i_ax].scatter(particle_dic["diam"][particle_dic["vterm_bohm"]>v_small_notuse],particle_dic[prop][particle_dic["vterm_bohm"]>v_small_notuse],s=1,c=particle_dic["N_monomer"][particle_dic["vterm_bohm"]>v_small_notuse],rasterized=True,norm=normN,cmap=cmapN,marker='o')
                
                #show legend
                axes[i_ax].legend() 
                
                #make labels
                if use_Deq:
                    axes[i_ax].set_xlabel("$D_{eq}$ [m]")
                else:
                    axes[i_ax].set_xlabel("$D_{max}$ [m]")
                axes[i_ax].set_ylabel(prop + " [m s-1]" ) #TODO: plot also the untis of these properties

                #change the axis
                #axes[i_ax].set_xlim([10**low_diam_log_detailed,10**high_diam_log_detailed]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                axes[i_ax].set_ylim([0,2.5]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                axes[i_ax].grid(which="both")

                #add colorbar
                cbar = fig.colorbar(im,ax=axes[i_ax])
                cbar.set_label("$N_{mono}$")
                #shift ticks to the middle of the color
                tick_locs = (boundsN[:-1:3]) + 0.5*(boundsN[1::3]-boundsN[:-1:3])
                cbar.set_ticks(tick_locs)
                cbar.set_ticklabels(N_mono_list[::3])
                if show_Nmono_fits:
                    fitlines = axes[i_ax].semilogx(np.nan,np.nan,linestyle="",label="$v_{term}$ using m-D\n& A-D fits and "+ vmodel_short,linewidth=linewidth)

                    #plot v-D corresponding to some selected monomer numbers
                    for i_Nmono,N_mono_now in enumerate(Nmono_fit_eval):
                        #if N_mono_now in [1,2,5,10,100,1000]:
                        #calculate vterm from the m/A-D relations and the corresponding fall speeds
                        vterm_fit = __fallspeed_relations.calc_vterm(prop[6:],mass_fit_detailed[i_Nmono,:],diam,area_fit_detailed[i_Nmono,:])
                        #if i_Nmono==0 and i_prop==0:
                        #    print "area_fitted",area_fit_detailed[i_Nmono,:]
                        #plot the lines (with a white surrounding)
                        fitlines_white = axes[i_ax].semilogx(diam,vterm_fit,c="white",linestyle="-",label="__None",linewidth=linewidth_white)
                        fitlines = axes[i_ax].semilogx(diam,vterm_fit,c=usecolorsN[N_mono_now==N_mono_list][0],linestyle="-",label="$N_{mono}=$" + "{}".format(N_mono_now),linewidth=linewidth)
                
                #add a lines for the aggregates
                vterm_fit_allagg = __fallspeed_relations.calc_vterm(prop[6:],fit_dic["mass_Nmono_allagg"],fit_dic["diam"],fit_dic["area_Nmono_allagg"])

                #fitlines_monomers = axes[i_ax].semilogx(np.nan,np.nan,c="b",linestyle="-",label="monomers",linewidth=2.*linewidth)
                #fitlines_allagg = axes[i_ax].semilogx(fit_dic["diam"],vterm_fit_allagg,c="w",linestyle="-",label="__None",linewidth=2.*linewidth_white)
                if not use_Deq:
                    fitlines_allagg = axes[i_ax].semilogx(fit_dic["diam"],vterm_fit_allagg,c="g",label="$N_{mono}$>1",linewidth=2.*linewidth,linestyle='--')

                
            #add legends
            for i_ax in [0,6,7,8]: #,7,8,9,10]:
                #skip variables which are not requested by plot_vars
                if plot_vars!=None and not ("vterm" in plot_vars) and i_ax>0: 
                    continue
                if i_ax=="0":
                    axes[i_ax].legend(ncol=2,loc="upper left")
                else:
                    axes[i_ax].legend(ncol=2,loc="lower right")
            if called_by_main:
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
                #subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
                
                # Save just the portion _inside_ the .. axis's boundaries
                save_axis=2
                
                for i_ax in range(0,len(axes)):#clear all other axes
                    if i_ax==save_axis:
                        continue
                    axes[i_ax].set_xlabel("")
                    axes[i_ax].axes.get_xaxis().set_ticklabels([])
                ''' 
                ax_description_array=["mD","mDnorm","AD","ADnorm","AratioD","AratioDnorm","vterm_bohm","vterm_HW10","vterm_KC05"]
                for save_axis,ax_description in enumerate(ax_description_array):
                    print "saving subfigure",save_axis,ax_description
                    try: #clear label of axis before and recover current one
                        axes[save_axis].set_xlabel("$D_{max}$ [m]")
                        axes[save_axis-1].set_xlabel("")
                    except:
                        pass
                    if "vterm" in ax_description:
                        axes[save_axis].set_ylabel("$v_{term}$ [$ms^{-1}$]")

                    extent = axes[save_axis].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
                    #from IPython.core.debugger import Tracer ; Tracer()()

                    fig.savefig('/home/mkarrer/Dokumente/plots/tmp.pdf',bbox_inches=extent.expanded(1.6, 1.45),dpi=400) #(1.6, 1.4) are width and height expanded around center
                '''
                #subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/tmp.pdf /home/mkarrer/Dokumente/plots/4paper/' + out_filestring + '_' + ax_description + '.pdf',shell=True)
                fig.clf()
            else:
                if i_particle_type==len(particle_types)-1:
                    return fig,axes,im,tick_locs,N_mono_list #im is for the colorbar
    return MAE

if __name__ == '__main__':
    
   
    particle_types=["plate","needle","dendrite","column","mixcoldend1","mixcolumndend"]

    txtfile="/home/mkarrer/Dokumente/plots/fit_results.txt"
    #optimize the appearance of the plot (figure size, fonts)
    fig=0 #will be overwritten in read_and_plot
    axes=0 #will be overwritten in read_and_plot
    
    MAE = dict() #dictionary for some mean absolute error measures

    #MAIN PART
    MAE = read_and_plot(fig,axes,particle_types,MAE=MAE,txtfile=txtfile,write_fit_params=True,called_by_main=True)
    print "MAE (mean absolute error)"
    print "mass_monodep mass_binary area_monodep area_binary"
    for particle_type in particle_types: 
        print particle_type," & {:8.3f} & {:8.3f} & {:8.3f} & {:8.3f} ".format(MAE[particle_type + "mass" + "_monodep"],MAE[particle_type + "mass" + "_binary"],MAE[particle_type + "area" + "_monodep"],MAE[particle_type + "area" + "_binary"])
