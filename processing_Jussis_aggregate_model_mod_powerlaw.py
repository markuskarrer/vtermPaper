# coding: utf-8
#import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.pylab as pylab
import csv #to read the txt files
import os
import sys
import glob #to get filename to list
import subprocess
import itertools#to read only certain lines of the txt files
from matplotlib.colors import LogNorm
#import other self-defined functions
import __postprocess_McSnow
import __postprocess_SB
import __fallspeed_relations
import __tools_for_processing_Jagg
import __plotting_functions
import generate_2Dhist_of_N_D_Nmono_from_MC_and_Jagg 
#from IPython.core.debugger import Tracer ; Tracer()()

'''
this code reads in properties from Jussis aggregate model
and fits different functions to each property (m-D,A-D,v-D for all fallspeed models), displays them as a function of monomer and makes an appropriate multivariate fit
'''


#define where the txt files with properties of the aggregates are stored
#prop_file_folder = "/data/optimice/Jussis_aggregates/"
prop_file_folder = "/data/optimice/Jussis_aggregates/fromHPC/"

grid_res = 10e-6
if grid_res==40e-6:
    print "processing particle with resolution: ",grid_res
    sensrun_folder = '/'
else:
    print "processing particle with resolution: ",grid_res
    sensrun_folder = 'res_'+ str(int(grid_res*1e6)) + 'mum/'# if '/' this will be in /data/optimice/Jussis_aggregates

    #sensrun_folder = '/' #res_10mum/'# if '/' this will be in /data/optimice/Jussis_aggregates


#read the properties of the aggregates into the particle_dic dictionary
particle_types = ["needle","column","plate","dendrite"] #["needle","column","plate","dendrite","bullet","rosette"] # ,"bullet"rosette
for particle_type in particle_types:
    #get m-D from the assumptions in the aggregate model (this is used at various places)
    a,b,c,d = __tools_for_processing_Jagg.calc_mD_AD_coeffs(particle_type)
    ###
    #plot the m-D, A-D and v-D relations
    ###
    number_of_plots = 14 #28

    #optimize the appearance of the plot (figure size, fonts)
    [fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots)
    
    print "########################"
    print "#####" + particle_type + "########"
    print "########################"
    particle_dic = __tools_for_processing_Jagg.init_particle_dict() #initialize dictionary which contains information about the type of the pristine crystals (dentrites, needles, plates, ...) the underlying distribution of the   
    if not os.path.exists(prop_file_folder + sensrun_folder):
        print prop_file_folder + sensrun_folder, " doesnt exist (exit in process_Jussis_aggregate_model_mod_powerlaw.py"
        sys.exit()
    for filename in glob.glob(prop_file_folder + sensrun_folder +  particle_type + '*properties.txt'): #loop over all property files (with different mean sizes)
        print "reading: " + filename
        #if float(filename[39:44])<9.: continue #"2." in filename: continue #print "dont take this" #continue
        #print filename[39:44],float(filename[39:44])<5
        with open(filename,"rb") as txtfile: #TODO: this is just one file! #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in
            prop_reader = csv.reader(filter(lambda row: row[0]!='#',txtfile), delimiter=' ', quoting=csv.QUOTE_NONNUMERIC, lineterminator=os.linesep) #row[0]!='#': skip the header; quoting avoids '' for formatted string; lineterminator avoids problems with system dependend lineending format https://unix.stackexchange.com/questions/309154/strings-are-missing-after-concatenating-two-or-more-variable-string-in-bash?answertab=active#tab-top

            #define the monomer numbers to read in
            take_all_Nmono_bool = True
            if take_all_Nmono_bool:
                N_mono_list = np.append(np.append(np.array(range(1,10)),np.array(range(10,100,10))),np.array(range(100,1001,100)))
            else: #select certain monomer numbers
                N_mono_list = np.array([1,2,3,4,5,6,7,8,9,10]) #10,50]) #,2,3,5,10,20,50,100,200,500,1000]) #-1 is the index shift
            #for N_mono_now in N_mono_list: # magnitude in [1,10,100]: #itertools.islice (which is used for reading only certain rows in the txt files) just works with range() like input, so this is the workaround to be able to select all kind of monomer numbers

            #read the property file in to the particle_dic dictionary
            for i_row,row_content in enumerate(prop_reader):
                if row_content[0] in N_mono_list: #i_row==0 or i_row==4 or i_row==9:
                    particle_dic["particle_type"] = np.append(particle_dic["particle_type"],particle_type)
                    particle_dic["N_monomer"] =     np.append(particle_dic["N_monomer"],row_content[0]) #python indices start with 0
                    particle_dic["mass"] =          np.append(particle_dic["mass"],     row_content[1])
                    particle_dic["area"] =          np.append(particle_dic["area"],     row_content[2])
                    particle_dic["diam"] =          np.append(particle_dic["diam"],     row_content[3])
                    
    print particle_dic #overview of the dictionary, might be helpful sometimes

    #set up array of diameters for displaying the fit
    diam = np.logspace(-4,-1.5,50) #set up array of diameters (for bin edges and KDE-grid)


    fit_dic = dict() #initialize fit dictionary

    #initialize an array which contains all fit results
    fitresult_arr_dir_powerlaw_string = np.zeros((5,N_mono_list.shape[0],2)) #fit results (saved as floats) #dimensions: 0-> [mass,array]; 1-> monomer number; 2-> parameters (a,b) in fit
    
    #fit the m-D and A-D properties
    for i_prop,(prop,prop_unit) in enumerate(zip(["mass","area"],["kg","m2"])): #,"HWJussi","KCJussi"]):
        print "fitting and plotting: ",prop
        
        ###
        #START: fitting 
        ###
        #fit m-D and A-D relations
        dummy_axes=0;dummy_colors=0 #we dont have any axes here to plot on
        particle_dic["weight_by_MC"] = np.ones_like(particle_dic["diam"]) #WORKAROUND: we dont need the weighting here
        
        #fit a power-law for each monomer number individually
        dummy_axes,particle_dic,fitline_allagg = __tools_for_processing_Jagg.fitting_wrapper(dummy_axes,particle_dic,prop,N_mono_list,dummy_colors,fit_dic,diam,function="powerlaw",linewidthinvalid=0.05,linewidthvalid=0.7,plot_fits=False)


        #fit the modified power-law m(D,N_mono)= a(N_mono)*D**(b(N_mono)* a_0*D**b0
        ab_func="powerlaw_rational1" # "powerlaw_polynom1"powerlaw_rational1" "powerlaw_rational2" and "powerlaw_rational3" need some more advanced initial guesses #select the functional form of a and b 
        if prop=="mass":
            #fitresult_mod_powerlaw = _
            #__tools_for_processing_Jagg.fit_2D_rational_function_transformwrapper(particle_dic["diam"],particle_dic["N_monomer"],particle_dic[prop],func='rational1_1_pade_fixp00_and_fixq00',method='lm',weight='None',habit=particle_type,prop=prop)
            if ab_func=="powerlaw_polynom1":
                fitresult_mod_powerlaw = __tools_for_processing_Jagg.fit_2D_rational_function_transformwrapper(particle_dic["diam"],particle_dic["N_monomer"],particle_dic[prop],func='powerlaw_polynom1_powerlaw_fixp0_fixr0',method='lm',weight='None',habit=particle_type,prop=prop)
            elif ab_func=="powerlaw_rational1":
                fitresult_mod_powerlaw = __tools_for_processing_Jagg.fit_2D_rational_function_transformwrapper(particle_dic["diam"],particle_dic["N_monomer"],particle_dic[prop],func='powerlaw_rational1_fixp0_fixr0',method='lm',weight='None',habit=particle_type,prop=prop)
            elif ab_func=="powerlaw_rational2":
                fitresult_mod_powerlaw = __tools_for_processing_Jagg.fit_2D_rational_function_transformwrapper(particle_dic["diam"],particle_dic["N_monomer"],particle_dic[prop],func='powerlaw_rational2_fixp0_fixr0',method='lm',weight='None',habit=particle_type,prop=prop)
            elif ab_func=="powerlaw_rational3":
                fitresult_mod_powerlaw = __tools_for_processing_Jagg.fit_2D_rational_function_transformwrapper(particle_dic["diam"],particle_dic["N_monomer"],particle_dic[prop],func='powerlaw_rational3_fixp0_fixr0',method='lm',weight='None',habit=particle_type,prop=prop)
            else: #TODO: include different fitting methods
                print "func: ", ab_func, " not (carefully) implemented in process_Jussis_aggregate_model_mod_powerlaw (and its modules)"; sys.exit()
        elif prop=="area":
            #fitresult_mod_powerlaw = __tools_for_processing_Jagg.fit_2D_rational_function_transformwrapper(particle_dic["diam"],particle_dic["N_monomer"],particle_dic[prop],func='rational2_2_pade_fixq00',method='lm',weight='None',habit=particle_type,prop=prop)
            if ab_func=="powerlaw_polynom1":
                fitresult_mod_powerlaw = __tools_for_processing_Jagg.fit_2D_rational_function_transformwrapper(particle_dic["diam"],particle_dic["N_monomer"],particle_dic[prop],func='powerlaw_polynom1_powerlaw',method='lm',weight='None',habit=particle_type,prop=prop)
            elif ab_func=="powerlaw_rational1":
                fitresult_mod_powerlaw = __tools_for_processing_Jagg.fit_2D_rational_function_transformwrapper(particle_dic["diam"],particle_dic["N_monomer"],particle_dic[prop],func='powerlaw_rational1',method='lm',weight='None',habit=particle_type,prop=prop)
            elif ab_func=="powerlaw_rational2":
                fitresult_mod_powerlaw = __tools_for_processing_Jagg.fit_2D_rational_function_transformwrapper(particle_dic["diam"],particle_dic["N_monomer"],particle_dic[prop],func='powerlaw_rational2',method='lm',weight='None',habit=particle_type,prop=prop)
            elif ab_func=="powerlaw_rational3":
                fitresult_mod_powerlaw = __tools_for_processing_Jagg.fit_2D_rational_function_transformwrapper(particle_dic["diam"],particle_dic["N_monomer"],particle_dic[prop],func='powerlaw_rational3',method='lm',weight='None',habit=particle_type,prop=prop)
            else: #TODO: include different fitting methods
                print "func: ", ab_func, " not (carefully) implemented in process_Jussis_aggregate_model_mod_powerlaw (and its modules)"; sys.exit()
                
        #save the fit results in a dictionary
        fit_dic[prop +"_coeff_mod_powerlaw"] = fitresult_mod_powerlaw
        
        ###
        #END: fitting 
        ###
        
    #do the same plots for mass and area     
    i_ax=-1
    for i_prop,(prop,prop_unit,prop_short) in enumerate(zip(["mass","area"],["kg","m2"],["m","A"])):
        print "fitting ",prop
        
        #####
        # START: plot a,b as a function of the monomer number (fit conducted for the subset of particle with the identical monomer number)
        #####
        
        #extract the fitted parameters from the dictionary to arrays
        fit_dic[prop +"_coeff_Nmono_merged"] = np.ones((len(N_mono_list),2)) # a parameter in m=a*D**b

        #merge the a-b fits for individual N_mono in one array
        for i_N_mono,N_mono in enumerate(N_mono_list):
            fit_dic[prop +"_coeff_Nmono_merged"][i_N_mono] = fit_dic[prop +"_coeff_Nmono_" + str(N_mono)]
        
        #plot the fit coefficients as a function of the monomer number
        i_ax+=1
        axes[i_ax].semilogx(N_mono_list,fit_dic[prop +"_coeff_Nmono_merged"][:,0],label=particle_type)
        #make labels
        axes[i_ax].set_xlabel("number of monomers / 1")
        if prop=="mass":
            axes[i_ax].set_ylabel("a coefficient in m=a*D**b")
        if prop=="area":
            axes[i_ax].set_ylabel("c coefficient in A=c*D**d")
        
        #now b in m=a*D**b
        i_ax+=1
        axes[i_ax].semilogx(N_mono_list,fit_dic[prop +"_coeff_Nmono_merged"][:,1],label=particle_type)
        #make labels
        axes[i_ax].set_xlabel("number of monomers / 1")
        if prop=="mass":
            axes[i_ax].set_ylabel("b coefficient in m=a*D**b")
        if prop=="area":
            axes[i_ax].set_ylabel("d coefficient in A=c*D**d")
        
        
        #####
        # END: plot a,b as a function of the monomer number (fit conducted for the subset of particle with the identical monomer number)
        #####
        
        ###
        #START: plot the property (normalized against sphere or the monomer) at selected diameters
        ###
        i_ax+=1
        
        #define the diameter array for the binning of the diameter (used for the calculation of the mean)
        low_diam_log=-4; high_diam_log=np.log10(0.032);nbins=6
        diam_edges = np.logspace(low_diam_log,high_diam_log,nbins) #set up array of diameters
        
        #loop over different diameter bins
        diam_selected1_handle_normmono = [None]*diam_edges[:-1].shape[0]
        diam_center = (diam_edges[:-1]+diam_edges[1:])/2
        diam_log_center = 10**((np.log10(diam_edges[:-1])+np.log10(diam_edges[1:]))/2)
        for i_diam,diam_lower in enumerate(diam_edges[:-1]):
            #get the center of the diameter bin
            #diam_log_center =#10**((np.log10(diam_edges[i_diam])+np.log10(diam_edges[i_diam+1]))/2)
            
            #calculate the properties of the monomer (for a large range of diameter)
            if prop=="mass":
                prop_monomer = a*diam_log_center[i_diam]**b
            elif prop=="area":
                prop_monomer = c*diam_log_center[i_diam]**d

            #initialize array containing the mean property (mass/area) at different diameters
            prop_agg_at_diam = np.ones_like(N_mono_list)*np.nan
            debug_piecewise_fitting=False #give some additional output later
            
            #perform a piecewise linear fitting
            for i_N_mono,N_mono in enumerate(N_mono_list):
                if debug_piecewise_fitting:
                    print i_N_mono,N_mono
                #select by the monomer number and the diameter bin (TODO: log-uniform diameter bins)
                fullfill_bin = np.logical_and(particle_dic["N_monomer"]==N_mono,np.logical_and(diam_edges[i_diam]<particle_dic["diam"],particle_dic["diam"]<diam_edges[i_diam+1]))
                if np.sum(fullfill_bin)>10: #dont do this for a too small sample
                    #make a piecewise fit (in the limits of the above defined diam_edges)
                    [fitresult_mod_powerlaw,covar] = __tools_for_processing_Jagg.fit_data(particle_dic["diam"][fullfill_bin],particle_dic[prop][fullfill_bin],func="powerlaw",method="leastsq",habit=particle_type,prop=prop)
                    prop_agg_at_diam[i_N_mono] = fitresult_mod_powerlaw[0]*diam_log_center[i_diam]**fitresult_mod_powerlaw[1]
                    if debug_piecewise_fitting:
                        print "result",np.nanmean(particle_dic[prop][fullfill_bin]),prop_agg_at_diam[i_N_mono] ,particle_dic[prop][fullfill_bin]
                        print "fitresult",fitresult_mod_powerlaw
                        print "ratio",prop_agg_at_diam[i_N_mono]/prop_monomer
            

            #normalize with the sphere or the monomer properties, respectively
            if prop=="mass":
                prop_at_diam_normed_by_sphere = prop_agg_at_diam/(1./6.*np.pi*diam_log_center[i_diam]**3) #->leads to effective density
                prop_at_diam_normed_by_monomer = prop_agg_at_diam/prop_monomer
            if prop=="area":
                prop_at_diam_normed_by_sphere = prop_agg_at_diam/(np.pi*diam_log_center[i_diam]**2)
                prop_at_diam_normed_by_monomer = prop_agg_at_diam/prop_monomer 
            
            ##plot the property normalized by a sphere
            diam_selected1_handle_normsphere = axes[i_ax].loglog(N_mono_list,prop_at_diam_normed_by_sphere,label="{:.2f}mm-{:.2f}mm".format(diam_edges[i_diam]*1e3,diam_edges[i_diam+1]*1e3))
            #make labels
            axes[i_ax].set_xlabel("number of monomers / 1")
            if prop=="mass":
                axes[i_ax].set_ylabel("eff. density / kg m-3")
            if prop=="area":
                axes[i_ax].set_title("proj. area ratio Ae/Aspher")
                
            ##same as above but normalized to the monomer properties
            diam_selected1_handle_normmono[i_diam] = axes[i_ax+1].plot(N_mono_list,prop_at_diam_normed_by_monomer,label="{:.2f}mm-{:.2f}mm".format(diam_edges[i_diam]*1e3,diam_edges[i_diam+1]*1e3))
        
        

        ### add the modified powerlaw plot to the above plotted data from the particles
        
        ##set up the D and N_mono array for analyzing the fit
        #D
        #diam_fit_eval = diam_center #link this to the selected diameter bins above #np.logspace(-4,-1.5,20)
        diam_log = np.log10(diam_center)
        #N_mono
        Nmono_fit_eval = N_mono_list #evaluate the 
        Nmono_log = np.log10(Nmono_fit_eval)
        
        #self-defined discrete colormap (for the monomer number N)
        # define the colormap and the discrete boundaries (defined by the considered monomer numbers) and set up an array of the same colors (for line-plotting)
        cmapN = plt.cm.viridis #brg
        #modify lowest color in colormap
        cmaplistN = [cmapN(i) for i in range(cmapN.N)]
        #cmaplist = cmaplist[80:] #crop colormap to avoid to light scatterer
        cmaplistN[0] = (1.0,0.0,0.0,1.0)
        cmapN = colors.LinearSegmentedColormap.from_list('mcm',cmaplistN, cmapN.N)
        boundsN = np.append(N_mono_list,[9999])
        normN = colors.BoundaryNorm(boundsN, cmapN.N)
        usecolorsN = cmapN(np.linspace(0,1,N_mono_list.shape[0])) #get colors from colormap to make consistent line colors for the monomer number
        
        #self-defined discrete colormap (for the diameter D)
        # define the colormap and the discrete boundaries (defined by the considered monomer numbers) and set up an array of the same colors (for line-plotting)
        cmapD = plt.cm.viridis #brg
        #modify lowest color in colormap
        #cmaplistD = [cmapD(i) for i in range(cmapD.N)]
        #cmapD = colors.LinearSegmentedColormap.from_list('mcm',cmaplistD, cmapD.N)
        boundsD = np.linspace(-4,-1.5,20) #this is log10(D)!!
        normD = colors.BoundaryNorm(boundsD, cmapD.N)
        usecolorsD = cmapD(np.linspace(0,1,boundsD.shape[0])) #get colors from colormap to make consistent line colors for the monomer number
        #from IPython.core.debugger import Tracer ; Tracer()()

        #spanned grid from D-N
        D_N_grid =  np.meshgrid(diam_log,Nmono_log)
        
        #apply the fit coefficients to the above spanned D,N field
        if ab_func=="powerlaw_polynom1":
            if prop=="mass":
                prop_fitted = __tools_for_processing_Jagg.powerlaw_polynom1_powerlaw_fixp0_fixr0(D_N_grid,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1])#this is the logscaled property (mass/area) divided by the monomer property
            if prop=="area":
                prop_fitted = __tools_for_processing_Jagg.powerlaw_polynom1_powerlaw(D_N_grid,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3])#this is the logscaled property (mass/area) divided by the monomer property
        elif ab_func=="powerlaw_rational1":
            if prop=="mass":
                prop_fitted = __tools_for_processing_Jagg.powerlaw_rational1_fixp0_fixr0(D_N_grid,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3])#this is the logscaled property (mass/area) divided by the monomer property
            if prop=="area":
                prop_fitted = __tools_for_processing_Jagg.powerlaw_rational1(D_N_grid,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3],fit_dic[prop +"_coeff_mod_powerlaw"][4],fit_dic[prop +"_coeff_mod_powerlaw"][5])#this is the logscaled property (mass/area) divided by the monomer property     
        elif ab_func=="powerlaw_rational2":
            if prop=="mass":
                prop_fitted = __tools_for_processing_Jagg.powerlaw_rational2_fixp0_fixr0(D_N_grid,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3],fit_dic[prop +"_coeff_mod_powerlaw"][4],fit_dic[prop +"_coeff_mod_powerlaw"][5],fit_dic[prop +"_coeff_mod_powerlaw"][6],fit_dic[prop +"_coeff_mod_powerlaw"][7])#this is the logscaled property (mass/area) divided by the monomer property
            if prop=="area":
                prop_fitted = __tools_for_processing_Jagg.powerlaw_rational2(D_N_grid,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3],fit_dic[prop +"_coeff_mod_powerlaw"][4],fit_dic[prop +"_coeff_mod_powerlaw"][5],fit_dic[prop +"_coeff_mod_powerlaw"][6],fit_dic[prop +"_coeff_mod_powerlaw"][7],fit_dic[prop +"_coeff_mod_powerlaw"][8],fit_dic[prop +"_coeff_mod_powerlaw"][9])#this is the logscaled property (mass/area) divided by the monomer property      
        elif ab_func=="powerlaw_rational3":
            if prop=="mass":
                prop_fitted = __tools_for_processing_Jagg.powerlaw_rational2_fixp0_fixr0(D_N_grid,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3],fit_dic[prop +"_coeff_mod_powerlaw"][4],fit_dic[prop +"_coeff_mod_powerlaw"][5],fit_dic[prop +"_coeff_mod_powerlaw"][6],fit_dic[prop +"_coeff_mod_powerlaw"][7])#this is the logscaled property (mass/area) divided by the monomer property
            if prop=="area":
                prop_fitted = __tools_for_processing_Jagg.powerlaw_rational2(D_N_grid,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3],fit_dic[prop +"_coeff_mod_powerlaw"][4],fit_dic[prop +"_coeff_mod_powerlaw"][5],fit_dic[prop +"_coeff_mod_powerlaw"][6],fit_dic[prop +"_coeff_mod_powerlaw"][7],fit_dic[prop +"_coeff_mod_powerlaw"][8],fit_dic[prop +"_coeff_mod_powerlaw"][9])#this is the logscaled property (mass/area) divided by the monomer property      
        else: #TODO: include different fitting methods
                print "func: ", ab_func, " not (carefully) implemented in process_Jussis_aggregate_model_mod_powerlaw (and its modules)"; sys.exit()
        #prop_fitted = __tools_for_processing_Jagg.polynom1_powerlaw_fixp0_fixr0(D_N_grid,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1]) #,fit_dic[prop +"_coeff_mod_powerlaw"][11])
        #prop_fitted = __tools_for_processing_Jagg.rational1_1_pade_fixp00_and_fixq00(D_N_grid,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3],fit_dic[prop +"_coeff_mod_powerlaw"][4],fit_dic[prop +"_coeff_mod_powerlaw"][5])
        
        #prop_fitted = __tools_for_processing_Jagg.rational2_2_pade_fixp00_and_fixq00(D_N_grid,fit_dic[prop +"_coeff_mod_powerlaw"][0],fit_dic[prop +"_coeff_mod_powerlaw"][1],fit_dic[prop +"_coeff_mod_powerlaw"][2],fit_dic[prop +"_coeff_mod_powerlaw"][3],fit_dic[prop +"_coeff_mod_powerlaw"][4],fit_dic[prop +"_coeff_mod_powerlaw"][5],fit_dic[prop +"_coeff_mod_powerlaw"][6],fit_dic[prop +"_coeff_mod_powerlaw"][8],fit_dic[prop +"_coeff_mod_powerlaw"][9]) #,fit_dic[prop +"_coeff_mod_powerlaw"][10])
    
        #transform the normed property back   
        if prop=="mass":
            mass_fit = 10**(prop_fitted)*a*diam_center**b
        if prop=="area":
            area_fit = 10**(prop_fitted)*c*diam_center**d
        
            
        for i_diam,diam_center_now in enumerate(diam_center):
            #from IPython.core.debugger import Tracer ; Tracer()()
            axes[i_ax+1].plot(Nmono_fit_eval,10**prop_fitted[:,i_diam],label="fit@{:.2f}mm".format(diam_center_now*1.e3),linestyle="--",color=diam_selected1_handle_normmono[i_diam][0].get_color())

        #make labels
        axes[i_ax+1].set_xlabel("number of monomers / 1")
        axes[i_ax+1].set_ylabel(prop + " / " + prop + " monomer / 1")

        ###
        #END: plot the property (normalized against sphere or the monomer) at selected diameters
        ###
            
        ###
        #plot the modified powerlaw fit (additional dependency on the monomer number)
        ###

        ###
        #START: plot the prop(m or D) vs diameter as scatter and the 2D-fit
        ###
        i_ax+=1
        diam_array = np.logspace(-4,-1,1000)
        
        #change the scale (according to xscale_vec and yscale_vec)
        axes[i_ax].set_xscale("log") #axes[i_ax].set_xscale(xscale_vec[i_prop]) 
        axes[i_ax].set_yscale("log") #axes[i_ax].set_yscale(yscale_vec[i_prop])
        


        i_ax+=1
        #plot the individual particle normalized by the properties of the monomer
        axes[i_ax].set_yscale('log')
        axes[i_ax].set_xscale('log')
        im = axes[i_ax].scatter(particle_dic["diam"],particle_dic[prop],s=1,c=particle_dic["N_monomer"],rasterized=True,norm=normN,cmap=cmapN) #all used particles        
        axes[i_ax].legend(loc="lower right")
        
        #overlay the fit
        for i_Nmono in range(0,Nmono_fit_eval.shape[0],5):
            #from IPython.core.debugger import Tracer ; Tracer()()
            if prop=="mass":
                fitlines = axes[i_ax].plot(diam_center,mass_fit[i_Nmono,:],c=usecolorsN[i_Nmono],linestyle="--")
            if prop=="area":
                fitlines = axes[i_ax].plot(diam_center,area_fit[i_Nmono,:],c=usecolorsN[i_Nmono],linestyle="--")
        #make labels
        axes[i_ax].set_xlabel("diameter / m")
        axes[i_ax].set_ylabel((" / ").join((prop,prop_unit))) 

        #change the axis
        axes[i_ax].set_ylim([0,np.array([1e-4,1e-3])[i_prop]]) #define the upper limit of the displayed axis
        axes[i_ax].grid(which="minor")
        #add colorbar
        cbar = fig.colorbar(im,ax=axes[i_ax])
        cbar.set_label("monomer number")
        #shift ticks to the middle of the color
        tick_locs = (boundsN[:-1:3]) + 0.5*(boundsN[1::3]-boundsN[:-1:3])
        cbar.set_ticks(tick_locs)
        cbar.set_ticklabels(N_mono_list[::3])
        
        
        ###
        #END: plot the prop(m or D) vs diameter as scatter and the 2D-fit
        ###
        
        ###
        #START: plot the prop(m or D) vs diameter as scatter and the 2D-fit (with transformed values D->log10(D) N->log10(N) prop->prop/prop_monomer)
        ###
        #mass vs diameter
        i_ax+=1
        #plot the individual particle normalized by the properties of the monomer
        #axes[i_ax].set_yscale('log')
        #axes[i_ax].set_xscale('log')
        #im = axes[i_ax].scatter(particle_dic["diam"],particle_dic[prop]/(a*particle_dic["diam"]**b),s=1,c=particle_dic["N_monomer"],rasterized=True,norm=normN,cmap=cmapN) #all used particles #TODO also for area
        match_Nmono = (particle_dic["N_monomer"]<=1000)
        if prop=="mass":
            im = axes[i_ax].scatter(np.log10(particle_dic["diam"][match_Nmono]),np.log10(particle_dic[prop][match_Nmono]/(a*particle_dic["diam"][match_Nmono]**b)),s=1,c=particle_dic["N_monomer"][match_Nmono],rasterized=True,norm=normN,cmap=cmapN) #all used particles #TODO also for area
        elif prop=="area":
            im = axes[i_ax].scatter(np.log10(particle_dic["diam"][match_Nmono]),np.log10(particle_dic[prop][match_Nmono]/(c*particle_dic["diam"][match_Nmono]**d)),s=1,c=particle_dic["N_monomer"][match_Nmono],rasterized=True,norm=normN,cmap=cmapN) #all used particles #TODO also for area
        #overlay lines for the fit
        
        #transform the normed property back
        
        #if prop=="mass":
        #    mass_rat_initial = 10**initial_guess*a*diam_fit_eval**b
        #if prop=="area":
        #    area_rat_fit = 10**initial_guess*c*diam_fit_eval**d
        
        #overlay the fit
        for i_Nmono in range(0,Nmono_fit_eval.shape[0],5):
            fitlines = axes[i_ax].plot(np.log10(diam_center),(prop_fitted[i_Nmono,:]),c=usecolorsN[i_Nmono],linestyle="--")
            
        '''
        for i_Nmono in range(0,Nmono_fit_eval.shape[0]):
            if not (i_Nmono % 5)==0:
                continue
            fitlines = axes[i_ax].plot(np.log10(diam_center),initial_guess[i_Nmono,:],linestyle='--',color=usecolorsN[i_Nmono])
        '''
        axes[i_ax].legend(loc="lower right")
        
        #make labels
        axes[i_ax].set_xlabel("log10(diameter)")
        axes[i_ax].set_ylabel("log10(" +prop + "/ " + prop + " monomer)")

        #change the axis
        #axes[i_ax].set_ylim([0,np.array([1e-4,1e-3])[i_prop]]) #define the upper limit of the displayed axis
        axes[i_ax].grid(which="minor")
        #add colorbar
        cbar = fig.colorbar(im,ax=axes[i_ax])
        cbar.set_label("monomer number")
        #shift ticks to the middle of the color
        tick_locs = (boundsN[:-1:3]) + 0.5*(boundsN[1::3]-boundsN[:-1:3])
        cbar.set_ticks(tick_locs)
        cbar.set_ticklabels(N_mono_list[::3])
        ###
        #START: plot the prop(m or D) vs N_monomer as scatter and the 2D-fit (with transformed values D->log10(D) N->log10(N) prop->prop/prop_monomer)
        ###
        i_ax+=1
        #plot the individual particle normalized by the properties of the monomer
        #axes[i_ax].set_yscale('log')
        axes[i_ax].set_xscale('log')
        if prop=="mass":
            im = axes[i_ax].scatter(particle_dic["N_monomer"],np.log10(particle_dic[prop]/(a*particle_dic["diam"]**b)),s=1,norm=normD,c=np.log10(particle_dic["diam"]),rasterized=True) 
        #,norm=normN,cmap=cmapN) #all used particles 
        if prop=="area":
            im = axes[i_ax].scatter(particle_dic["N_monomer"],np.log10(particle_dic[prop]/(c*particle_dic["diam"]**d)),s=1,c=np.log10(particle_dic["diam"]),rasterized=True) 
        #,norm=normN,cmap=cmapN) #all used particles
        
        #overlay the fit
        for i_diam in range(0,diam_center.shape[0]):
            #get the i-th element of the colorbar in which diam_center lies
            i_color=0
            while(np.log10(diam_center[i_diam])>boundsD[i_color+1]):
                i_color+=1
                
            if prop=="mass":
                fitlines = axes[i_ax].plot(Nmono_fit_eval,(prop_fitted[:,i_diam]),c=usecolorsD[i_color],linestyle="--")
            if prop=="area":
                fitlines = axes[i_ax].plot(Nmono_fit_eval,(prop_fitted[:,i_diam]),c=usecolorsD[i_color],linestyle="--")

        #make labels
        axes[i_ax].set_xlabel("monomer number / m")
        axes[i_ax].set_ylabel("log10(" +prop + "/ " + prop + " monomer)")

        #change the axis
        #axes[i_ax].set_ylim([0,np.array([1e-4,1e-3])[i_prop]]) #define the upper limit of the displayed axis
        #axes[i_ax].grid(which="minor")
        #add colorbar
        cbar = fig.colorbar(im,ax=axes[i_ax])
        cbar.set_label("log10(diameter)")
        #shift ticks to the middle of the color
        #tick_locs = (boundsN[:-1:3]) + 0.5*(boundsN[1::3]-boundsN[:-1:3])
        #cbar.set_ticks(tick_locs)
        cbar.set_ticklabels(10**cbar.get_ticks())
        '''
        for i_diam in range(0,diam_center.shape[0],5):
            fitlines = axes[i_ax].plot(Nmono_fit_eval,initial_guess[:,i_diam],linestyle='--')
        '''
        ###
        #END: plot the prop(m or D) vs N_monomer as scatter and the 2D-fit (with transformed values D->log10(D) N->log10(N) prop->prop/prop_monomer)
        ###
    for i_ax in range(0,12):
        axes[i_ax].legend(bbox_to_anchor=(1.02, 1.02))
        #ax.legend(bbox_to_anchor=(1.1, 1.05))
    #save the plot (and open it)
    plt.tight_layout()
    dir_save = '/home/mkarrer/Dokumente/plots/Jagg/'
    if not os.path.exists(dir_save): #create direktory if it does not exists
        os.makedirs(dir_save)
    out_filestring = "mod_powerlaw_fit_Jagg_mD_AD_vD_" + particle_type + "_gridres" + str(grid_res)
    plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
    plt.savefig(dir_save + out_filestring + '.png', dpi=100)
    print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
    subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
    plt.clf()
    plt.close()
