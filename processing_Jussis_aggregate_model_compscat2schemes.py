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

#define which lines to show:

show_vel_fits = True #show the velocity fits (power-law/Atlas)
show_powerlaw_fits = False
show_Atlas_fit = False #show the Atlas-type fit for cloud ice & snow
#show the fits as numbers
show_fittext = False #show text for the fit
#from previous model settings
show_Axels_powerlaw = False #show also the old fit of Axels power-law
show_Axels_atlas = False #show also the old fit of Axels Atlas-type
#from other models
show_P3 = False
show_MC = False

#switch on/of fitting of terminal velocity (removes the fit lines)
fitvel=1 #0-> no fitting of velocity at all; 1-> fit only monomer and all aggregates; 2-> fit everything (including seperated for different monomer numbers)


#define where the txt files with properties of the aggregates are stored
prop_file_folder = "/data/optimice/Jussis_aggregates/fromHPC/"
#define which particles (with which resolution to read in)
grid_res = 10e-6
if grid_res==40e-6:
    print "processing particle with resolution: ",grid_res
    sensrun_folder = '/'
else:
    print "processing particle with resolution: ",grid_res
    sensrun_folder = 'res_'+ str(int(grid_res*1e6)) + 'mum/'

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

#calculate the arrays with masses and areas corresponding to the common fit_dic["diam"] and based on the (piecewise) power-law fits
fit_dic = dict() #initialize fit dictionary
#set up array of diameters (for bin edges and KDE-grid) and the fitting
diam = np.logspace(-4,np.log10(3e-2),50) #set up array of diameters (for bin edges and KDE-grid)
fit_dic["diam"] = diam
#
for dict_now in (mDADvD_dict_MC,mDADvD_dict_P3,mDADvD_dict_SBcloudice,mDADvD_dict_SBsnow,mDADvD_dict_SBcloudiceAtlas,mDADvD_dict_SBsnowAtlas):
    dict_now = __setup_mDAD_frommodels.calc_area_mass_vterm_arrays(fit_dic["diam"],dict_now)

####
#END: get model parameters
#for comparing with models 
####

#read the properties of the aggregates into the particle_dic dictionary
particle_types = ["plate_dendrite_column"] #"needle","dendrite","plate","column","rosette","bullet"] # ,"bullet"rosette
for particle_type_comb in particle_types:
    print "########################"
    print "#####" + particle_type_comb + "########"
    print "########################"
    particle_dic = __tools_for_processing_Jagg.init_particle_dict() #initialize dictionary which contains information about the type of the pristine crystals (dentrites, needles, plates, ...) the underlying distribution of the
    #allow also combination of habits
    if particle_type_comb not in ("needle","dendrite","plate","column","rosette","bullet"):
        #particle_type_comb should be a combination of particles in particle type (e.g. "plate_dendrite")
        particle_type_list = re.split(r'_',particle_type_comb) #split the combination to read in all habits
    for particle_type in particle_type_list: #this is a single particle type (no combination allowed after that)
        for filename in glob.glob(prop_file_folder + sensrun_folder +  particle_type + '*properties.txt'): #loop over all property files (with different mean sizes)
            #read size parameter from filename in order to disregard some size parameter
            m=re.search(prop_file_folder + sensrun_folder + particle_type + '_(.+?)_properties.txt', filename) #search for patter
            size_param_now = float(m.group(1))
            if size_param_now>1000: #ignore all sizeparameter above ... or below ...
                continue
            
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
                    N_mono_list = np.array([1,2,3,5,7,10,20,30,50,70,100,200,300,400,500,700,1000]) #10,50]) #,2,3,5,10,20,50,100,200,500,1000]) #-1 is the index shift

                for i_row,row_content in enumerate(prop_reader): #TODO: the row is not any more equal to the monomer number #read the row with N_mono_now (the currently considered monomer number)
                    if row_content[0] in N_mono_list: #i_row==0 or i_row==4 or i_row==9:
                        particle_dic["particle_type"] = np.append(particle_dic["particle_type"],particle_type)
                        particle_dic["N_monomer"] = np.append(particle_dic["N_monomer"],row_content[0]) #python indices start with 0
                        particle_dic["mass"] = np.append(particle_dic["mass"],row_content[1])
                        particle_dic["area"] = np.append(particle_dic["area"],row_content[2])
                        particle_dic["diam"] = np.append(particle_dic["diam"],row_content[3])
                    

    #calculate the terminal velocitys for each available velocity model
    for velocity_model in ["HW10","KC05","bohm"]: #,"mitch_heym"]:
            particle_dic["vterm_" + velocity_model] = __fallspeed_relations.calc_vterm(velocity_model,particle_dic["mass"],particle_dic["diam"],particle_dic["area"])
    print particle_dic #overview of the dictionary, might be helpful sometimes


    #get the 2D weighting function (derived from the comparison of the aggregate model output and McSnow runs)
    weight_uniform_or_MC = 0 #0-> uniform 1-> MC
    if weight_uniform_or_MC==0:
        diam_edges,Nmono_edges,H_MC,H_Jagg,H = generate_2Dhist_of_N_D_Nmono_from_MC_and_Jagg.N_D_Dmono_from_MC_and_Jagg(particle_type=particle_type)
        #recalculate the weighting factor to be uniform
        for i_diam in range(H_MC.shape[0]): #first dimension is the diameter
            ###calculate the weighting factor
            #we normalize the weighting factor at each diameter to one
            sum_valid_Jagg = np.nansum(H_Jagg[i_diam,:])
            if sum_valid_Jagg>0: #check if there is any Jagg data
                H[i_diam,:] = np.divide(1,H_Jagg[i_diam,:], out=np.zeros_like(H_Jagg[i_diam,:]), where=H_Jagg[i_diam,:]!=0)/sum_valid_Jagg
    elif weight_uniform_or_MC==1:
        diam_edges,Nmono_edges,H_MC,H_Jagg,H = generate_2Dhist_of_N_D_Nmono_from_MC_and_Jagg.N_D_Dmono_from_MC_and_Jagg(particle_type=particle_type)
    weighting_factor=H

    #initialize the weighting attribute
    particle_dic["weight_by_MC"] = np.zeros_like(particle_dic["diam"])
    
    #set the weight for each simulated aggregate
    for i_agg in range(0,particle_dic["diam"].shape[0]):#loop over all aggregates to give each aggregate the right weighting factor
        
        #TODO: this is probably very inefficiently programmed with nested loops and ifs, but maybe this is not important here
        for i_diam_now,diam_now in enumerate(diam_edges[:-1]):
            for i_N_mono_now,N_mono_now in enumerate(Nmono_edges[:-1]):
                if diam_edges[i_diam_now]<=particle_dic["diam"][i_agg]<diam_edges[i_diam_now+1]: #check if we are in the right diameter bin
                    if Nmono_edges[i_N_mono_now]<=particle_dic["N_monomer"][i_agg]<Nmono_edges[i_N_mono_now+1]: #check if we are in the right monomer number bin
                        particle_dic["weight_by_MC"][i_agg] = weighting_factor[i_diam_now,i_N_mono_now] ##from IPython.core.debugger import Tracer ; Tracer()()    
    
    ###
    #plot the m-D, A-D and v-D relations
    ###
    number_of_plots = 6 #28
    
    #optimize the appearance of the plot (figure size, fonts)
    [fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots)
    params = {'legend.fontsize': 'medium',
                'axes.labelsize': 'medium', #size: Either a relative value of 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large' or an absolute font size, e.g., 12
        'axes.titlesize':'medium',
        'xtick.labelsize':'medium',
        'ytick.labelsize':'medium'}
    pylab.rcParams.update(params)
    #define size of colorbar labels
    colbarlabelsize=11
    
    ###
    #subplot 1: histogram and KDE of the data
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
    axes[i_ax].set_xlabel("diameter / m")
    axes[i_ax].set_ylabel("absolute counts")
    ax_dens.set_ylabel("PDF")

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
        im = axes[i_ax].scatter(particle_dic["diam"],particle_dic[prop],s=1,c=particle_dic["N_monomer"],rasterized=True,norm=norm,cmap=cmap) #all used particles

        #fit m-D and A-D relations
        #axes[i_ax],particle_dic,fitline_allagg = __tools_for_processing_Jagg.fitting_wrapper(axes[i_ax],particle_dic,prop,N_mono_list,usecolors,fit_dic,diam,function="powerlaw",linewidthinvalid=0.05,linewidthvalid=0.7)
        plot_allfits_bool=False #plot fits and display fit coefficients in text for all monomer number
        plot_binary_fits= True #plot fits and display fit coefficients in text for monomer and allagg
        axes[i_ax],particle_dic,fitline_allagg = __tools_for_processing_Jagg.fitting_wrapper(axes[i_ax],particle_dic,prop,N_mono_list,usecolors,fit_dic,diam,function="powerlaw",weight='None',linewidthinvalid=0.9,linewidthvalid=1.7,plot_fits=plot_allfits_bool,plot_binary_fits=plot_binary_fits)
        
       
        ####
        #compare to models
        ####
       
        #plot the properties as assumed in the schemes in the current plot
        if prop=="mass": #the SB dont have to specify the area
            axes[i_ax].semilogx(fit_dic["diam"],mDADvD_dict_SBcloudice[prop_short + "(D_array)"],color='blue',label="SB cloud ice",linestyle='--')
            axes[i_ax].semilogx(fit_dic["diam"],mDADvD_dict_SBsnow[prop_short + "(D_array)"],color='green',label="SB snow",linestyle='--')
        if show_MC:
            axes[i_ax].semilogx(fit_dic["diam"],mDADvD_dict_MC[prop_short + "(D_array)"],color='red',label="McSnow",linestyle='--')
        if show_P3:
            axes[i_ax].semilogx(fit_dic["diam"],mDADvD_dict_P3["v_mitch_heym(D_array)"],color='orange',label="P3unrimed",linestyle='-.') #mitch heym is hardcoded because there are no other options yet



        #show legend
        axes[i_ax].legend(loc="lower right")
        
        #make labels
        axes[i_ax].set_xlabel("diameter / m")
        axes[i_ax].set_ylabel((" / ").join((prop,prop_unit))) 

        #change the axis
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
            #add a legend for the fit lines which are not explained by the colorbar
            #axes[i_ax].legend(handles=[fitline_allagg],loc="lower right", prop={'size': 8})
        
        #shrink colorbar lables
        cbar.ax.tick_params(labelsize=colbarlabelsize) 

    #loop over different terminal velocity models
    for i_prop,prop in enumerate(["vterm_bohm","vterm_HW10","vterm_KC05"]): #,"vterm_mitch_heym"]): #,"HWJussi","KCJussi"]):
        print "fitting and plotting: ",prop
        
        
        i_ax = i_prop+1+2 #+1 because we plotted already the histogram +2 because we plotted already mass and area
        #change the scale (according to xscale_vec and yscale_vec)
        axes[i_ax].set_xscale("log") #axes[i_ax].set_xscale(xscale_vec[i_prop]) 
        axes[i_ax].set_yscale("linear") #axes[i_ax].set_yscale(yscale_vec[i_prop])
        
        #plot the scatter plot of the individual particle # use the colormap and the bounds which where defined before
        im = axes[i_ax].scatter(particle_dic["diam"],particle_dic[prop],s=1,c=particle_dic["N_monomer"],rasterized=True,norm=norm,cmap=cmap,marker='o') #,norm=colors.LogNorm(),cmap="gist_ncar")+
        im = axes[i_ax].scatter(particle_dic["diam"][particle_dic["N_monomer"]==1],particle_dic[prop][particle_dic["N_monomer"]==1],s=1,c=particle_dic["N_monomer"][particle_dic["N_monomer"]==1],rasterized=True,norm=norm,cmap=cmap,marker='o') #,norm=colors.LogNorm(),cmap="gist_ncar")+

        #calculate the terminal velocity based on the m-D and A-D fits and plot this line # use the colormap and the bounds which where defined before
        velocity_model = prop[6:] #get the name of the terminal velocity model from the property name (prop) which starts with "vterm"
        for i,N_mono_now in enumerate(N_mono_list):
            fit_dic["vterm_" + velocity_model + "_fitted_via_mD_AD_Nmono" + str(N_mono_now)] = __fallspeed_relations.calc_vterm(velocity_model,fit_dic["mass" + "_Nmono_" + str(N_mono_now)],fit_dic["diam"],fit_dic["area" + "_Nmono_" + str(N_mono_now)])
            if show_vel_fits:
                if N_mono_now<2 or fitvel>=2:
                    fitline_Nmonospecific, = axes[i_ax].plot(fit_dic["diam"],fit_dic["vterm_" + velocity_model + "_fitted_via_mD_AD_Nmono" + str(N_mono_now)],marker='None',c=usecolors[i],linewidth=linewidthinvalid,label="Nmono=1") #,norm=colors.LogNorm(),cmap="gist_ncar")
                    #fitline_Nmonospecific_dense_enough, = axes[i_ax].plot(fit_dic["diam"][fit_dic["diam_dens_enough_" + str(N_mono_now)]],fit_dic["vterm_" + velocity_model + "_fitted_via_mD_AD_Nmono" + str(N_mono_now)][fit_dic["diam_dens_enough_" + str(N_mono_now)]],marker='None',c=usecolors[i],linewidth=linewidthvalid) 

        #now for all aggregates
        fit_dic["vterm_" + velocity_model + "_fitted_via_mD_AD_Nmonoallagg"] = __fallspeed_relations.calc_vterm(velocity_model,fit_dic["mass" + "_Nmono_allagg"],fit_dic["diam"],fit_dic["area" + "_Nmono_allagg"])
        if show_vel_fits:
            if fitvel>=1:
                fitline_Nmonoallagg, = axes[i_ax].plot(fit_dic["diam"],fit_dic["vterm_" + velocity_model + "_fitted_via_mD_AD_Nmonoallagg"],marker='None',c="g",label="Nmono>1",linewidth=linewidthinvalid)
                #fitline_Nmonoallagg_dense_enough, = axes[i_ax].plot(fit_dic["diam"][fit_dic["diam_dens_enough_allagg"]],fit_dic["vterm_" + velocity_model + "_fitted_via_mD_AD_Nmonoallagg"][fit_dic["diam_dens_enough_allagg"]],marker='None',c="black",label="Nmono>1",linewidth=linewidthvalid)

        if fitvel>=1:
            #fit an Atlas type to the "mean v-D" relations (derived from the A-D and m-D fits)
            for i,str_monomer_agg in enumerate(["1","allagg"]): # make a power-law fit for the pristine crystals and the aggregate
                [fitresult,covar] = __tools_for_processing_Jagg.fit_data(fit_dic["diam"],fit_dic[prop +  "_fitted_via_mD_AD_Nmono" + str_monomer_agg],func="Atlas") #, weight=fit_dic["dens" + "_Nmono_" + str_monomer_agg])
                fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg] = fitresult #copy the Atlas-fit coefficients in the fit_dic dictionary
                #fit additionally a powerlaw
                [fitresult,covar] = __tools_for_processing_Jagg.fit_data(fit_dic["diam"],fit_dic[prop +  "_fitted_via_mD_AD_Nmono" + str_monomer_agg],func="powerlaw") #, weight=fit_dic["dens" + "_Nmono_" + str_monomer_agg])
                fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_powerlaw"] = fitresult #copy the Atlas-fit coefficients in the fit_dic dictionary

                #recalculate A-B*exp(-gam*D) from fit
                fit_dic[prop + "_Nmono_" + str_monomer_agg] = fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg][0]-fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg][1]*np.exp(-fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg][2]*fit_dic["diam"]) #calculate arrays of masses and areas based on the m-D fits
                #recalculate powerlaw from fit
                fit_dic[prop + "_Nmono_" + str_monomer_agg + "_powerlaw"] = fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_powerlaw"][0]*fit_dic["diam"]**fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_powerlaw"][1] #calculate arrays of masses and areas based on the m-D fits
                
                
                if show_Atlas_fit:
                    if str_monomer_agg=="1": #no label
                        Atlasfit_Nmono1, = axes[i_ax].plot(fit_dic["diam"],fit_dic[prop + "_Nmono_" + str_monomer_agg],marker='None',c=np.array(["b","g"])[i],linestyle="-",label="Atlas-type (cloud ice)",linewidth=linewidthvalid)
                        #from IPython.core.debugger import Tracer ; Tracer()()
                        #Atlasfit_Nmono1_dense_enough, = axes[i_ax].plot(fit_dic["diam"][fit_dic["diam_dens_enough_" + str(1)]],fit_dic[prop + "_Nmono_" + str_monomer_agg][fit_dic["diam_dens_enough_" + str(1)]],marker='None',c=np.array(["b","black"])[i],linestyle="--",label="Nmono=1(weighted Atlas fit)",linewidth=linewidthvalid)  
                if show_Atlas_fit:
                    if str_monomer_agg=="allagg":
                        Atlasfit_Nmonoallagg, = axes[i_ax].plot(fit_dic["diam"],fit_dic[prop + "_Nmono_" + str_monomer_agg],marker='None',c=np.array(["b","g"])[i],linestyle="-",label="Atlas-type (snow)",linewidth=linewidthvalid)
                        #Atlasfit_Nmonoallagg_dense_enough, = axes[i_ax].plot(fit_dic["diam"][fit_dic["diam_dens_enough_allagg"]],fit_dic[prop + "_Nmono_" + str_monomer_agg][fit_dic["diam_dens_enough_allagg"]],marker='None',c=np.array(["b","black"])[i],linestyle="--",label="Nmono>1(weighted Atlas fit)",linewidth=linewidthvalid)
                if show_powerlaw_fits:
                    powerlaw_handle, = axes[i_ax].plot(fit_dic["diam"], fit_dic[prop + "_Nmono_" + str_monomer_agg + "_powerlaw"],marker='None',c=np.array(["b","g"])[i],linestyle="-.",label="powerlaw",linewidth=linewidthvalid)
        if fitvel>=1:
            if show_fittext:
                #add labels for the fit result to the plot
                axes[i_ax] = __tools_for_processing_Jagg.add_fitresult_text(axes[i_ax],fit_dic,prop,[1],function="Atlas")
                axes[i_ax] = __tools_for_processing_Jagg.add_fitresult_text(axes[i_ax],fit_dic,prop,[1],function="powerlaw_fallspeed")

        


        ####
        #compare to models
        ####
       
        #plot the fall speed in the current plot
        if show_Axels_powerlaw:
            SB_ci_handle, = axes[i_ax].semilogx(fit_dic["diam"],mDADvD_dict_SBcloudice["v(D_array)"],color='blue',label="SB cloud ice",linestyle='--')
            SB_snow_handle, = axes[i_ax].semilogx(fit_dic["diam"],mDADvD_dict_SBsnow["v(D_array)"],color='green',label="SB snow",linestyle='--')
        if show_MC:
            MC_handle, = axes[i_ax].semilogx(fit_dic["diam"],mDADvD_dict_MC["v_" + velocity_model + "(D_array)"],color='red',label="McSnow(" + velocity_model + ")",linestyle='--')
        if show_P3:
            P3_handle, = axes[i_ax].semilogx(fit_dic["diam"],mDADvD_dict_P3["v_mitch_heym(D_array)"],color='orange',label="P3unrimed",linestyle='-.') #mitch heym is hardcoded because there are no other options yet
        if show_Axels_atlas:
            Axel_iceAtlas_handle, = axes[i_ax].semilogx(mDADvD_dict_SBcloudiceAtlas["D_max_from_moltenD"],mDADvD_dict_SBcloudiceAtlas["v(D_array)"],color='blue',label="SB cloud ice Atlas",linestyle=':')
            Axel_snowAtlas_handle, = axes[i_ax].semilogx(mDADvD_dict_SBsnowAtlas["D_max_from_moltenD"],mDADvD_dict_SBsnowAtlas["v(D_array)"],color='green',label="SB snow Atlas",linestyle=':')
        #show legend
        axes[i_ax].legend() 
        #axes[i_ax].legend(handles=[SB_ci_handle,SB_snow_handle,Atlasfit_Nmono1,Atlasfit_Nmonoallagg,MC_handle,P3_handle]) #uncomment this line for a special order
        #make labels
        axes[i_ax].set_xlabel("diameter / m")
        axes[i_ax].set_ylabel(prop + " / m s-1" ) #TODO: plot also the untis of these properties

        #change the axis
        axes[i_ax].set_xlim([1e-4,2e-2]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
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
    out_filestring = "Jagg_mD_AD_vD_" + particle_type + "_gridres" + str(grid_res) + "_compscat2models"
    plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
    plt.savefig(dir_save + out_filestring + '.png', dpi=400)
    print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
    subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
    #plt.clf()
    #plt.close()
    
    # Save just the portion _inside_ the second axis's boundaries
    extent = axes[3].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    axes[3].set_ylabel("terminal velocity / m s-1")
    fig.savefig(dir_save + 'spec_ax_figure.png', bbox_inches=extent.expanded(1.4, 1.5),dpi=400)
