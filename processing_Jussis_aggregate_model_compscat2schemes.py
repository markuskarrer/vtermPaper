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
import argparse; parser = argparse.ArgumentParser() #get optional arguments from shell
#import other self-defined functions
import __postprocess_McSnow
import __postprocess_SB
import __fallspeed_relations
import __tools_for_processing_Jagg
import __plotting_functions
import __setup_mDAD_frommodels
import generate_2Dhist_of_N_D_Nmono_from_MC_and_Jagg
#from IPython.core.debugger import Tracer ; Tracer()()

'''
this code reads in properties from Jussis aggregate model, calculate + plots the following fallspeed
and fits different functions to each property (m-D,A-D,v-D for all fallspeed models)
'''

show_lines = dict()
#define which lines to show (ATTENTION: default, can be overwritten by optional arguments below):
show_lines["SB_mD"] = True #show the m-D relation of the SB scheme
show_lines["vel_fits"] = True #show the velocity fits (power-law/Atlas)
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
show_lines["highest_diam_for_fit"] = -999
scale_vel_d = 0.001 #scale the velocity by some factor so it gets visible (TODO: add a real meaning to that value and a second y-axis)
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
if args.highest_diam_for_fit:
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

#define where the txt files with properties of the aggregates are stored
#prop_file_folder = "/data/optimice/Jussis_aggregates/tumbling_asratio_complete/"
prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation"
#prop_file_folder = "/data/optimice/Jussis_aggregates_bugfixedrotation/tumbling_asratio/"

#define which particles (with which resolution to read in)
grid_res_array = [10e-6] #[1e-6,5e-6,10e-6] #resolution for big medium and smaller particles in mum [10e-6,5e-5]

####
#START: get model parameters
#for comparing with models 
####

#McSnow
mDADvD_dict_MC = __setup_mDAD_frommodels.get_model_mDADs(model="MC")
#SB
mDADvD_dict_SBcloudice = __setup_mDAD_frommodels.get_model_mDADs(model="SBcloudice")
mDADvD_dict_SBsnow = __setup_mDAD_frommodels.get_model_mDADs(model="SBsnow")
mDADvD_dict_SBcloudiceAtlas = __setup_mDAD_frommodels.get_model_mDADs(model="SBcloudice_Atlas")
mDADvD_dict_SBsnowAtlas = __setup_mDAD_frommodels.get_model_mDADs(model="SBsnow_Atlas")
#P3
mDADvD_dict_P3 = __setup_mDAD_frommodels.get_model_mDADs(model="P3")
#Morrison scheme (from WRF 4.0)
morr2mom_cloudice  = __setup_mDAD_frommodels.get_model_mDADs(model="morr2mom_cloudice")
morr2mom_snow  = __setup_mDAD_frommodels.get_model_mDADs(model="morr2mom_snow")

#calculate the arrays with masses and areas corresponding to the common fit_dic["diam"] and based on the (piecewise) power-law fits
fit_dic = dict() #initialize fit dictionary
#set up array of diameters for bin edges, the KDE-grid and the fitting
low_diam_log=-4; 
#highest diameter for fit
print show_lines["highest_diam_for_fit"]
if show_lines["highest_diam_for_fit"]==-999:
    high_diam_log=-1.3979400086720375 #np.log10(3e-2)=-1.5228787452803376 #np.log10(4e-2)=-1.3979400086720375
else:
    high_diam_log=np.log10(float(show_lines["highest_diam_for_fit"])) #
diam = np.logspace(low_diam_log,high_diam_log,50) #set up array of diameters (for bin edges and KDE-grid)
fit_dic["diam"] = diam
#get the m,A,v-array corresponding to diam and save them in separate dictionaries
for dict_now in (mDADvD_dict_MC,mDADvD_dict_P3,mDADvD_dict_SBcloudice,mDADvD_dict_SBsnow,mDADvD_dict_SBcloudiceAtlas,mDADvD_dict_SBsnowAtlas,morr2mom_cloudice,morr2mom_snow):
    dict_now = __setup_mDAD_frommodels.calc_area_mass_vterm_arrays(fit_dic["diam"],dict_now)

####
#END: get model parameters
#for comparing with models 
####

#loop over different particle habits
particle_types = ["bullet"] #"plate","needle","dendrite","column","rosette","bullet","plate_dendrite",...]] # ,"bullet"rosette
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
        for grid_res in grid_res_array: #loop over different grid-resolutions in order to have small monomer numbers modelled woth higher resolution
            
            print "processing particle with resolution: ",grid_res
            sensrun_folder = 'res_'+ str(int(grid_res*1e6)) + 'mum/'
            for filename in glob.glob(prop_file_folder + sensrun_folder +  particle_type + '*properties.txt'): #loop over all property files (with different mean sizes)
                #read size parameter from filename in order to disregard some size parameter
                m=re.search(prop_file_folder + sensrun_folder + particle_type + '_(.+?)_properties.txt', filename) #search for patter
                size_param_now = float(m.group(1))
                if size_param_now>1000: #ignore all sizeparameter above ... or below ...
                    continue
                #verbose reading of files
                print "reading: " + filename

                with open(filename,"rb") as txtfile: #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in
                    prop_reader = csv.reader(filter(lambda row: row[0]!='#',txtfile), delimiter=' ', quoting=csv.QUOTE_NONNUMERIC, lineterminator=os.linesep) #row[0]!='#': skip the header; quoting avoids '' for formatted string; lineterminator avoids problems with system dependend lineending format https://unix.stackexchange.com/questions/309154/strings-are-missing-after-concatenating-two-or-more-variable-string-in-bash?answertab=active#tab-top

                    #define the monomer numbers to read in
                    take_all_Nmono_bool = True
                    if grid_res==1e-6 or len(grid_res_array)==1: #particles < 10mum are modelled only for small monomer numbers #len(grid_res_array)==1 is for some quick test
                        N_mono_list_tmp = np.array(range(1,10))
                    elif grid_res==5e-6: #particles < 10mum are modelled only for small monomer numbers
                        N_mono_list_tmp = np.array(range(10,101,10))
                    elif grid_res==10e-6:#take only the large monomer numbers from the low-resolution runs
                        N_mono_list_tmp = np.array(range(100,1001,100))
                    
                    
                    for i_row,row_content in enumerate(prop_reader): #TODO: the row is not any more equal to the monomer number #read the row with N_mono_now (the currently considered monomer number)
                        if row_content[0] in N_mono_list_tmp: #i_row==0 or i_row==4 or i_row==9:
                            particle_dic["particle_type"] = np.append(particle_dic["particle_type"],particle_type)
                            particle_dic["N_monomer"] = np.append(particle_dic["N_monomer"],row_content[0]) #python indices start with 0
                            particle_dic["mass"] = np.append(particle_dic["mass"],row_content[1])
                            if tumbling:
                                particle_dic["area"] = np.append(particle_dic["area"],row_content[2]) #,row_content[2]) -[2] for tumbling [4] for disabled tumbling
                            elif not tumbling:
                                particle_dic["area"] = np.append(particle_dic["area"],row_content[4]) #,row_content[2]) -[2] for tumbling [4] for disabled tumbling
                            particle_dic["diam"] = np.append(particle_dic["diam"],row_content[3])
            if len(grid_res_array)==1: #some quick tests
                N_mono_list = np.array(range(1,10))
            else:
                N_mono_list = np.append(np.append(np.array(range(1,10)),np.array(range(10,100,10))),np.array(range(100,1001,100))) #the full monomer number list
            
    #########
    #END: read the properties of the aggregates into the particle_dic dictionary
    #########                

    #calculate the terminal velocitys of the simulated particles individually for each available velocity model
    for velocity_model in ["HW10","KC05","bohm"]: #,"mitch_heym"]:
            particle_dic["vterm_" + velocity_model] = __fallspeed_relations.calc_vterm(velocity_model,particle_dic["mass"],particle_dic["diam"],particle_dic["area"])
    print particle_dic #overview of the dictionary, might be helpful sometimes

    ######
    #calculate a weigting of the simulated particle to get an accurate fit for the whole diameter range and give equally importance to all monomer numbers
    ######
    #START: get the 2D weighting function (derived from the comparison of the aggregate model output and McSnow runs)
    ######
    weight_uniform_or_MC = 0 #0-> uniform 1-> MC
    if weight_uniform_or_MC==0:
        #diam_edges,Nmono_edges,H_MC,H_Jagg,H = generate_2Dhist_of_N_D_Nmono_from_MC_and_Jagg.N_D_Dmono_from_MC_and_Jagg(particle_type=particle_type)
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
    elif weight_uniform_or_MC==1:
        diam_edges,Nmono_edges,H_MC,H_Jagg,H = generate_2Dhist_of_N_D_Nmono_from_MC_and_Jagg.N_D_Dmono_from_MC_and_Jagg(particle_type=particle_type)
    weighting_factor=H
    
    #initialize the weighting attribute
    particle_dic["weight"] = np.zeros_like(particle_dic["diam"])
    
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
    
    ###
    #initialize the plot
    ###
    number_of_plots = 6 #28
    
    #optimize the appearance of the plot (figure size, fonts)
    [fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots,aspect_ratio=0.25)
    params = {'legend.fontsize': 'medium',
                'axes.labelsize': 'medium', #size: Either a relative value of 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large' or an absolute font size, e.g., 12
        'axes.titlesize':'medium',
        'xtick.labelsize':'medium',
        'ytick.labelsize':'medium'}
    pylab.rcParams.update(params)
    #define size of colorbar labels
    colbarlabelsize=11
    
    ###
    #START: subplot 1: histogram and KDE of the data
    ###
    
    #self-defined discrete colormap
    # define the colormap and the discrete boundaries (defined by the considered monomer numbers) and set up an array of the same colors (for line-plotting)
    cmap = plt.cm.BuGn #brg
    #modify lowest color in colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmaplist = cmaplist[80:] #crop colormap to avoid to light scatterer
    #cmap = cmap(np.linspace(0.2, 0.8, 100))

    cmaplist[0] = (0.0,0.0,1.0,1.0)
    cmap = colors.LinearSegmentedColormap.from_list('mcm',cmaplist, cmap.N)
    
    bounds = np.append(N_mono_list,[9999])
    norm = colors.BoundaryNorm(bounds, cmap.N)
    usecolors = cmap(np.linspace(0,1,N_mono_list.shape[0])) #get colors from colormap to make consistent line colors for the KDE-lines #pylab.cm.Greens(np.linspace(0,1,N_mono_list.shape[0]))

    #overlay the histograms corresponding to specific monomer numbers
    i_ax=0
    ax_dens = axes[i_ax].twinx() #second y-axis for density (derived by KDE)


    #separating sparse from valid data
    dens_lim = 1e-3 #0.005 #value of separating sparse from valid data
    linewidthinvalid=0.9;linewidthvalid=1.7 #define with which linewidth valid and invalid data should be plotted
    #calculate the histogram and KDE for the monomer numbers in N_mono_list
    for i,N_mono_now in enumerate(N_mono_list): 
        #calculate and plot histogram
        #print particle_dic["N_monomer"] #==N_mono_now
        hist = axes[i_ax].hist(particle_dic["diam"][particle_dic["N_monomer"]==N_mono_now],bins=diam,alpha=0.2,normed=False,color=usecolors[i])
        axes[i_ax].set_xscale("log") #set x-axis to logarithmic
        #estimate the PDF with KDE 
        sigmai = 0.2 #0.62 #TODO: is there a way to better constrain this value
        dens = __postprocess_McSnow.kernel_estimate(particle_dic["diam"][particle_dic["N_monomer"]==N_mono_now],diam,sigmai,weight="None",space='loge')   #TODO: is space='loge' appropriate
        fit_dic["dens" + "_Nmono_" + str(N_mono_now)] = dens
        
        #separate diameter region with reasonably high enough density (available particle) from roughly dense /no particle
        fit_dic["diam_dens_enough_" + str(N_mono_now)] = fit_dic["dens" + "_Nmono_" + str(N_mono_now)]>dens_lim
        fit_dic["diam_dens_low_" + str(N_mono_now)] = fit_dic["dens" + "_Nmono_" + str(N_mono_now)]<dens_lim
        
        #plot the KDE-result
        dens_handle = ax_dens.plot(diam,dens,color=usecolors[i])
    #calculate also the density for all aggregates with Nmono>1
    dens = __postprocess_McSnow.kernel_estimate(particle_dic["diam"][particle_dic["N_monomer"]>1],diam,sigmai,weight="None",space='loge')   #TODO: is space='loge' appropriate
    fit_dic["dens" + "_Nmono_allagg"] = dens
    #separate diameter region with reasonably high enough density (available particle) from roughly dense /no particle
    fit_dic["diam_dens_enough_allagg"] = fit_dic["dens" + "_Nmono_allagg"]>dens_lim
    fit_dic["diam_dens_low_allagg"] = fit_dic["dens" + "_Nmono_allagg"]<dens_lim
    #plot the KDE-result also for the Nmono>1 case
    dens_handle_Nmonoallagg, = ax_dens.plot(diam,dens,color="black",label="Nmono>1")
    ax_dens.tick_params(axis="y",direction="in", pad=-27) #direction="in" is for the ticks; pad defines where the label is located (horizontally)
    #add a legend for the fit lines which are not explained by the colorbar
    axes[i_ax].legend(handles=[dens_handle_Nmonoallagg],loc="upper left", prop={'size': 8})
    
    #make labels for the histogram and density plot
    axes[i_ax].set_xlabel("diameter D / m")
    axes[i_ax].set_ylabel("absolute counts")
    ax_dens.set_ylabel("PDF")

    ###
    #END: subplot 1: histogram and KDE of the data
    ###
    
    #####
    #START: plots in m-D/A-D space
    #####

    #initialize an array which contains all fit results
    fitresult_arr_dir_powerlaw_string = np.zeros((5,N_mono_list.shape[0],2)) #fit results (saved as floats) #dimensions: 0-> [mass,array]; 1-> monomer number; 2-> parameters (a,b) in fit
    
    #fit the m-D and A-D properties
    for i_prop,(prop,prop_unit,prop_short) in enumerate(zip(["mass","area"],["kg","m2"],["m","A"])): #,"HWJussi","KCJussi"]):
        print "fitting and plotting: ",prop
        
        i_ax = i_prop+1 #+1 because we plotted already the histogram
        #change the scale (according to xscale_vec and yscale_vec)
        axes[i_ax].set_xscale("log") #axes[i_ax].set_xscale(xscale_vec[i_prop]) 
        axes[i_ax].set_yscale("log") #axes[i_ax].set_yscale(yscale_vec[i_prop])
        
        #plot the scatter plot (of mass and area for single particle) # use the colormap and the bounds which where defined before
        im = axes[i_ax].scatter(particle_dic["diam"][particle_dic["N_monomer"]>0],particle_dic[prop][particle_dic["N_monomer"]>0],s=1,c=particle_dic["N_monomer"][particle_dic["N_monomer"]>0],rasterized=True,norm=norm,cmap=cmap) #all used particles

        #fit m-D and A-D relations
        #axes[i_ax],particle_dic,fitline_allagg = __tools_for_processing_Jagg.fitting_wrapper(axes[i_ax],particle_dic,prop,N_mono_list,usecolors,fit_dic,diam,function="powerlaw",linewidthinvalid=0.05,linewidthvalid=0.7)
        plot_allfits_bool=False #plot fits and display fit coefficients in text for all monomer number
        plot_binary_fits= True #plot fits and display fit coefficients in text for monomer and allagg
        #axes[i_ax],particle_dic,fitline_allagg = __tools_for_processing_Jagg.fitting_wrapper(axes[i_ax],particle_dic,prop,N_mono_list,usecolors,fit_dic,diam,function="powerlaw",weight=particle_dic["weight"][N_mono_now],linewidthinvalid=0.9,linewidthvalid=1.7,plot_fits=plot_allfits_bool,plot_binary_fits=plot_binary_fits)
        axes[i_ax],particle_dic,fitline_allagg =__tools_for_processing_Jagg.fitting_wrapper(axes[i_ax],particle_dic,prop,N_mono_list,usecolors,fit_dic,diam,function="powerlaw",weight="None",linewidthinvalid=0.9,linewidthvalid=1.7,plot_fits=plot_allfits_bool,plot_binary_fits=plot_binary_fits)
        
        if take_mono_prop_theor: #overwrite fit for monomers with theoretical value
            if prop=="mass":
                fit_dic['mass_coeff_Nmono_1'] = __tools_for_processing_Jagg.calc_mD_AD_coeffs(particle_type)[0:2]
            if prop=="area":
                fit_dic['area_coeff_Nmono_1'] = __tools_for_processing_Jagg.calc_mD_AD_coeffs(particle_type)[2:4]
            fit_dic[prop + "_Nmono_1"] = fit_dic[prop + "_coeff_Nmono_1"][0]*fit_dic["diam"]**fit_dic[prop + "_coeff_Nmono_1"][1] #recalculate arrays of masses and areas 
            #from IPython.core.debugger import Tracer ; Tracer()()
            
        ####
        #compare to models
        ####
       
        #plot the properties as assumed in the schemes in the current plot
        if prop=="mass" and show_lines["SB_mD"]: #the SB dont have to specify the area
            axes[i_ax].semilogx(fit_dic["diam"],mDADvD_dict_SBcloudice[prop_short + "(D_array)"],color='blue',label="SB cloud ice",linestyle='--')
            axes[i_ax].semilogx(fit_dic["diam"],mDADvD_dict_SBsnow[prop_short + "(D_array)"],color='y',label="SB snow",linestyle='--')
        if show_lines["MC"]:
            axes[i_ax].semilogx(fit_dic["diam"],mDADvD_dict_MC[prop_short + "(D_array)"],color='red',label="McSnow",linestyle='--')
        if show_lines["P3"]:
            axes[i_ax].semilogx(fit_dic["diam"],mDADvD_dict_P3["v_mitch_heym(D_array)"],color='orange',label="P3unrimed",linestyle='-.') #mitch heym is hardcoded because there are no other options yet

        #show legend
        axes[i_ax].legend(loc="lower right")
        #make labels
        axes[i_ax].set_xlabel("diameter D / m")
        axes[i_ax].set_ylabel((" / ").join((prop,prop_unit))) 
        #change the axis
        axes[i_ax].set_xlim([10**(low_diam_log),10**(high_diam_log)]) #define the upper limit of the displayed axis
        axes[i_ax].set_ylim([0,np.array([1e-4,2e-2])[i_prop]]) #define the upper limit of the displayed axis
        axes[i_ax].grid(which="minor")
        #add colorbar
        cbar = fig.colorbar(im,ax=axes[i_ax])
        cbar.set_label("monomer number")
        #shift ticks to the middle of the color
        tick_locs = (bounds[:-1:3]) + 0.5*(bounds[1::3]-bounds[:-1:3])
        cbar.set_ticks(tick_locs)
        cbar.set_ticklabels(N_mono_list[::3])
        if not plot_allfits_bool:
            #add labels for the fit result to the plot
            axes[i_ax] = __tools_for_processing_Jagg.add_fitresult_text(axes[i_ax],fit_dic,prop,N_mono_list,function="powerlaw",hide_Nmono_fits=take_all_Nmono_bool)
        
        #shrink colorbar lables
        cbar.ax.tick_params(labelsize=colbarlabelsize) 
    #####
    #END: plots in m-D/A-D space
    #####



    #loop over different terminal velocity models in which 1. the velocity of each individual particle is plotted 2. the velocity according to the m-D and A-D fits 3. different possible fits (power-law, Atlas type) are calculated (and plotted)
    for i_prop,prop in enumerate(["vterm_bohm","vterm_HW10","vterm_KC05"]): #,"vterm_mitch_heym"]): #,"HWJussi","KCJussi"]):
        print "fitting and plotting: ",prop
        
        #define current axis
        i_ax = i_prop+1+2 #+1 because we plotted already the histogram +2 because we plotted already mass and area
        
        #change the scale (according to xscale_vec and yscale_vec)
        axes[i_ax].set_xscale("log") #axes[i_ax].set_xscale(xscale_vec[i_prop]) 
        axes[i_ax].set_yscale("linear") #axes[i_ax].set_yscale(yscale_vec[i_prop])
        
        #####
        #1. the velocity of each individual particle is plotted 
        #####
        if not show_lines["no_simparticle"]:
            #plot the scatter plot of the individual particle # use the colormap and the bounds which where defined before
            im = axes[i_ax].scatter(particle_dic["diam"][particle_dic["N_monomer"]>0],particle_dic[prop][particle_dic["N_monomer"]>0],s=1,c=particle_dic["N_monomer"][particle_dic["N_monomer"]>0],rasterized=True,norm=norm,cmap=cmap,marker='o') #,norm=colors.LogNorm(),cmap="gist_ncar")+
            #im = axes[i_ax].scatter(particle_dic["diam"][particle_dic["N_monomer"]==1],particle_dic[prop][particle_dic["N_monomer"]==1],s=1,c=particle_dic["N_monomer"][particle_dic["N_monomer"]==1],rasterized=True,norm=norm,cmap=cmap,marker='o') #,norm=colors.LogNorm(),cmap="gist_ncar")+

        ####
        #2. the velocity according to the m-D and A-D fits
        ####
        #calculate the terminal velocity based on the m-D and A-D fits and plot this line # use the colormap and the bounds which where defined before
        velocity_model = prop[6:] #get the name of the terminal velocity model from the property name (prop) which starts with "vterm"
        for i,N_mono_now in enumerate(N_mono_list):
            fit_dic["vterm_" + velocity_model + "_fitted_via_mD_AD_Nmono" + str(N_mono_now)] = __fallspeed_relations.calc_vterm(velocity_model,fit_dic["mass" + "_Nmono_" + str(N_mono_now)],fit_dic["diam"],fit_dic["area" + "_Nmono_" + str(N_mono_now)])
            if show_lines["vel_fits"]:
                if N_mono_now<2 or fitvel>=2:
                    fitline_Nmonospecific, = axes[i_ax].semilogx(fit_dic["diam"],fit_dic["vterm_" + velocity_model + "_fitted_via_mD_AD_Nmono" + str(N_mono_now)],marker='None',c=usecolors[i],linewidth=linewidthinvalid,label="Nmono=1") #,norm=colors.LogNorm(),cmap="gist_ncar")

        #now for all aggregates
        fit_dic["vterm_" + velocity_model + "_fitted_via_mD_AD_Nmonoallagg"] = __fallspeed_relations.calc_vterm(velocity_model,fit_dic["mass" + "_Nmono_allagg"],fit_dic["diam"],fit_dic["area" + "_Nmono_allagg"])
        if show_lines["vel_fits"]:
            if fitvel>=1:
                fitline_Nmonoallagg, = axes[i_ax].semilogx(fit_dic["diam"],fit_dic["vterm_" + velocity_model + "_fitted_via_mD_AD_Nmonoallagg"],marker='None',c='y',label="Nmono>1",linewidth=linewidthinvalid)
                #fitline_Nmonoallagg_dense_enough, = axes[i_ax].plot(fit_dic["diam"][fit_dic["diam_dens_enough_allagg"]],fit_dic["vterm_" + velocity_model + "_fitted_via_mD_AD_Nmonoallagg"][fit_dic["diam_dens_enough_allagg"]],marker='None',c="black",label="Nmono>1",linewidth=linewidthvalid)

        #####
        #3. different fits (power-law, Atlas type) are calculated (and plotted)
        #####
        if fitvel>=1:
            #fit an Atlas type to the "mean v-D" relations (derived from the A-D and m-D fits)
            for i,(str_monomer_agg,corr_cat) in enumerate(zip(["1","allagg"],["cloud ice","snow"])): 
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
                
                if show_lines["Atlas_fit"]:
                    Atlasfit_Nmonoallagg, = axes[i_ax].semilogx(fit_dic["diam"],fit_dic[prop + "_Nmono_" + str_monomer_agg],marker='None',c=np.array(["b",'y'])[i],linestyle="--",label="Atlas-type " + corr_cat,linewidth=linewidthvalid)
                    if show_lines["del_vterm"]:
                        Atlasfit_Nmonoallagg_delv, = axes[i_ax].semilogx(fit_dic["diam"][:-1],scale_vel_d*1./np.diff(fit_dic["diam"])*np.diff(fit_dic[prop + "_Nmono_" + str_monomer_agg]),marker='x',c=np.array(["b","g"])[i],linestyle="--",label="Atlas-type " + corr_cat + " del v",linewidth=linewidthvalid)
                if show_lines["powerlaw_fits"]:
                    powerlaw_handle, = axes[i_ax].semilogx(fit_dic["diam"], fit_dic[prop + "_Nmono_" + str_monomer_agg + "_powerlaw"],marker='None',c=np.array(["b",'y'])[i],linestyle="-.",label="powerlaw " + corr_cat,linewidth=linewidthvalid)
        if fitvel>=1:
            if show_lines["fittext_vel"]:
                #add labels for the fit result to the plot
                axes[i_ax] = __tools_for_processing_Jagg.add_fitresult_text(axes[i_ax],fit_dic,prop,[1],function="Atlas")
                axes[i_ax] = __tools_for_processing_Jagg.add_fitresult_text(axes[i_ax],fit_dic,prop,[1],function="powerlaw_fallspeed")

        


        ####
        #compare to models
        ####
       
        #plot the fall speed in the current plot
        if show_lines["SB_powerlaw"]:
            SB_ci_handle,   = axes[i_ax].semilogx(fit_dic["diam"],mDADvD_dict_SBcloudice["v(D_array)"],color='blue', label="SB cloud ice",linestyle=':') #,marker='x')
            SB_snow_handle, = axes[i_ax].semilogx(fit_dic["diam"],mDADvD_dict_SBsnow["v(D_array)"]    ,color='y',label="SB snow",     linestyle=':') #,marker='x')
            if show_lines["del_vterm"]:
                SB_snow_handle_diff, = axes[i_ax].semilogx(fit_dic["diam"][:-1],scale_vel_d*1./np.diff(fit_dic["diam"])*np.diff(mDADvD_dict_SBsnow["v(D_array)"]),color='y',label="SB snow del v",marker='x',     linestyle='-.') #optionally add the gradients of the terminal velocity
        if show_lines["MC"]:
            MC_handle, = axes[i_ax].semilogx(fit_dic["diam"],mDADvD_dict_MC["v_" + velocity_model + "(D_array)"],color='red',label="McSnow(" + velocity_model + ")",linestyle='--')
        if show_lines["P3"]:
            P3_handle, = axes[i_ax].semilogx(fit_dic["diam"],mDADvD_dict_P3["v_mitch_heym(D_array)"],color='orange',label="P3unrimed",linestyle='-.') #mitch heym is hardcoded because there are no other options yet
        if show_lines["SB_Atlas"]:
            Axel_iceAtlas_handle, = axes[i_ax].semilogx(mDADvD_dict_SBcloudiceAtlas["D_max_from_moltenD"],mDADvD_dict_SBcloudiceAtlas["v(D_array)"],color='blue',label="SB cloud ice Atlas",linestyle='--',marker='x')
            Axel_snowAtlas_handle, = axes[i_ax].semilogx(mDADvD_dict_SBsnowAtlas["D_max_from_moltenD"],mDADvD_dict_SBsnowAtlas["v(D_array)"],color='y',label="SB snow Atlas",linestyle='--',marker='x')
            if show_lines["del_vterm"]:
                Axel_snowAtlas_handle, = axes[i_ax].semilogx(mDADvD_dict_SBsnowAtlas["D_max_from_moltenD"][:-1],scale_vel_d*1./np.diff(mDADvD_dict_SBsnowAtlas["D_max_from_moltenD"])*np.diff(mDADvD_dict_SBsnowAtlas["v(D_array)"]),color='y',label="SB snow Atlas del v",linestyle='--',marker="x")
                

        #show legend
        axes[i_ax].legend() 
        #expl1, = axes[i_ax].semilogx(np.nan,np.nan,label="model assumptions:",linestyle=""); expl2, = axes[i_ax].semilogx(np.nan,np.nan,label="fit to simulated particles:",linestyle="");axes[i_ax].legend(handles=[expl2,Atlasfit_Nmono1,Atlasfit_Nmonoallagg,expl1,SB_ci_handle,SB_snow_handle,MC_handle,P3_handle]) #uncomment this line for a special order
        
        #make labels
        axes[i_ax].set_xlabel("diameter D / m")
        axes[i_ax].set_ylabel(prop + " / m s-1" ) #TODO: plot also the untis of these properties

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
        
        ###
        #calculate the relations which can be put into the SB-scheme and display it here
        ###
        print "####################"
        print "#coefficients for SB derived from the here investigated particles (hydro_model:" + prop + ")"
        print "####################"
        #cloud ice
        print "ice: m(D):a= ",fit_dic["mass_coeff_Nmono_1"][0]," b= ",fit_dic["mass_coeff_Nmono_1"][1] ," v(D):a= ",fit_dic[prop + "_coeff_Nmono_1_powerlaw"][0]," b= ",fit_dic[prop + "_coeff_Nmono_1_powerlaw"][1]
        print "ice: Atlas-type (D_max): ",fit_dic[prop + "_coeff_Nmono_1"]

        #calculate the D(m) and v(m) coefficients
        fit_dic["mass_coeff_Nmono_1_N(m)"] = __postprocess_SB.convert_ND_to_Nm_from_coeff(fit_dic["mass_coeff_Nmono_1"][0],fit_dic["mass_coeff_Nmono_1"][1],fit_dic[prop + "_coeff_Nmono_1_powerlaw"][0],fit_dic[prop + "_coeff_Nmono_1_powerlaw"][1])
        print "ice: D(m):a= ", fit_dic["mass_coeff_Nmono_1_N(m)"][0]," b= ",fit_dic["mass_coeff_Nmono_1_N(m)"][1], " v(m):a= ",fit_dic["mass_coeff_Nmono_1_N(m)"][2]," b= ",fit_dic["mass_coeff_Nmono_1_N(m)"][3]
        print "ice: Atlas-type (D_eq): ",fit_dic[prop + "_coeff_Nmono_1_Deq"]

        #snow
        print "snow: m(D):a= ", fit_dic["mass_coeff_Nmono_allagg"][0]," b= ",fit_dic["mass_coeff_Nmono_allagg"][1], " v(D):a= ",fit_dic[prop + "_coeff_Nmono_allagg_powerlaw"][0]," b= ",fit_dic[prop + "_coeff_Nmono_allagg_powerlaw"][1]
        print "snow: Atlas-type (D_max): ",fit_dic[prop + "_coeff_Nmono_allagg"]
        #calculate the D(m) and v(m) coefficients
        fit_dic["mass_coeff_Nmono_allagg_N(m)"] = __postprocess_SB.convert_ND_to_Nm_from_coeff(fit_dic["mass_coeff_Nmono_allagg"][0],fit_dic["mass_coeff_Nmono_allagg"][1],fit_dic[prop + "_coeff_Nmono_allagg_powerlaw"][0],fit_dic[prop + "_coeff_Nmono_allagg_powerlaw"][1])
        print "snow: D(m):a= ", fit_dic["mass_coeff_Nmono_allagg_N(m)"][0]," b= ",fit_dic["mass_coeff_Nmono_allagg_N(m)"][1], " v(m):a= ",fit_dic["mass_coeff_Nmono_allagg_N(m)"][2]," b= ",fit_dic["mass_coeff_Nmono_allagg_N(m)"][3]
        print "snow: Atlas-type (D_eq): ",fit_dic[prop + "_coeff_Nmono_allagg_Deq"]

        
    #add a colorbar also for the histogram (which is "stolen" from the scatter plots)
    for ax in [axes[0],ax_dens]:
        cbar = fig.colorbar(im,ax=ax)
        cbar.set_label("monomer number")
        #shift ticks to the middle of the color
        tick_locs = (bounds[:-1]) + 0.5*(bounds[1:]-bounds[:-1])
        cbar.set_ticks(tick_locs)
        cbar.set_ticklabels(N_mono_list)
    
        #shrink colorbar lables
        cbar.ax.tick_params(labelsize=colbarlabelsize)
    
    #save the plot (and open it)
    plt.tight_layout()
    dir_save = '/home/mkarrer/Dokumente/plots/Jagg/'
    if not os.path.exists(dir_save): #create direktory if it does not exists
        os.makedirs(dir_save)
    if tumbling:
        tumbling_add="_wtumbling"
    else: 
        tumbling_add="_wotumbling"
    out_filestring = "Jagg_mD_AD_vD_" + particle_type_comb + "_gridres" + str(grid_res) + "_compscat2models" + tumbling_add + add_displayed_lines2string
    plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
    plt.savefig(dir_save + out_filestring + '.png', dpi=400)
    print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
    subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
    #plt.clf()
    #plt.close()
    
    # Save just the portion _inside_ the second axis's boundaries
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
