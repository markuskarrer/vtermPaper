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
from matplotlib.colors import LogNorm
#import other self-defined functions
import __postprocess_McSnow
import __postprocess_SB
import __fallspeed_relations
import __tools_for_processing_Jagg
import __plotting_functions
#from IPython.core.debugger import Tracer ; Tracer()()

'''
this code reads in properties from Jussis aggregate model, calculate + plots the following fallspeed
and fits different functions to each property (m-D,A-D,v-D for all fallspeed models)
'''


#define where the txt files with properties of the aggregates are stored
prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation"

grid_res = 10e-6
if grid_res==40e-6:
    print "processing particle with resolution: ",grid_res
    sensrun_folder = '/'
else:
    print "processing particle with resolution: ",grid_res
    sensrun_folder = 'res_'+ str(int(grid_res*1e6)) + 'mum/'# if '/' this will be in /data/optimice/Jussis_aggregates

    #sensrun_folder = '/' #res_10mum/'# if '/' this will be in /data/optimice/Jussis_aggregates


#read the properties of the aggregates into the particle_dic dictionary
particle_types = ["plate"] #,"dendrite","plate","column","rosette"] # ,"bullet"rosette
for particle_type in particle_types:
    print "########################"
    print "#####" + particle_type + "########"
    print "########################"
    particle_dic = __tools_for_processing_Jagg.init_particle_dict() #initialize dictionary which contains information about the type of the pristine crystals (dentrites, needles, plates, ...) the underlying distribution of the   
    for i_file,filename in enumerate(glob.glob(prop_file_folder + sensrun_folder +  particle_type + '*properties.txt')): #loop over all property files (with different mean sizes)
        print "reading: " + filename
        #if float(filename[39:44])<9.: continue #"2." in filename: continue #print "dont take this" #continue
        #print filename[39:44],float(filename[39:44])<5
        with open(filename,"rb") as txtfile: #TODO: this is just one file! #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in
            prop_reader = csv.reader(filter(lambda row: row[0]!='#',txtfile), delimiter=' ', quoting=csv.QUOTE_NONNUMERIC, lineterminator=os.linesep) #row[0]!='#': skip the header; quoting avoids '' for formatted string; lineterminator avoids problems with system dependend lineending format https://unix.stackexchange.com/questions/309154/strings-are-missing-after-concatenating-two-or-more-variable-string-in-bash?answertab=active#tab-top

            #define the monomer numbers to read in
            take_all_Nmono_bool = False
            if take_all_Nmono_bool:
                N_mono_list = np.append(np.array(range(1,10)),np.array(range(10,100,10)))
            else: #select certain monomer numbers
                N_mono_list = np.array([1,2,3,5,10,20,50,100]) #-1 is the index shift
            #for N_mono_now in N_mono_list: # magnitude in [1,10,100]: #itertools.islice (which is used for reading only certain rows in the txt files) just works with range() like input, so this is the workaround to be able to select all kind of monomer numbers
            #    #for row_content in itertools.islice(prop_reader,N_mono_now-1,N_mono_now): #prop_reader: read one line from the property file and add this to the particle_dic dictionary
            #    #print N_mono_now
            for i_row,row_content in enumerate(prop_reader): #TODO: the row is not any more equal to the monomer number #read the row with N_mono_now (the currently considered monomer number)
                if row_content[0] in N_mono_list: #i_row==0 or i_row==4 or i_row==9:
                    particle_dic["particle_type"] = np.append(particle_dic["particle_type"],particle_type)
                    particle_dic["N_monomer"] = np.append(particle_dic["N_monomer"],row_content[0]) #python indices start with 0
                    particle_dic["mass"] = np.append(particle_dic["mass"],row_content[1])
                    particle_dic["area"] = np.append(particle_dic["area"],row_content[2])
                    particle_dic["diam"] = np.append(particle_dic["diam"],row_content[3])
    print particle_dic #overview of the dictionary, might be helpful sometimes

    #calculate the fall speeds for each available velocity model
    for velocity_model in ["HW10","KC05","bohm","mitch_heym"]:
            #try: #sometimes there is an empty array
            particle_dic["vterm_" + velocity_model] = __fallspeed_relations.calc_vterm(velocity_model,particle_dic["mass"],particle_dic["diam"],particle_dic["area"])
            #except Exception as e:
            #    print e
            #    pass

    #print particle_dic #overview of the dictionary, might be helpful sometimes

    ###
    #plot the m-D, A-D and v-D relations
    ###
    number_of_plots = 7 #28
    
    #optimize the appearance of the plot (figure size, fonts)
    [fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots)

    ###
    #subplot 1: histogram and KDE of the data
    ###
    
    #self-defined discrete colormap
    # define the colormap and the discrete boundaries (defined by the considered monomer numbers) and set up an array of the same colors (for line-plotting)
    cmap = plt.cm.brg
    bounds = np.append(N_mono_list,[9999])
    norm = colors.BoundaryNorm(bounds, cmap.N)
    usecolors = pylab.cm.brg(np.linspace(0,1,N_mono_list.shape[0]))

    #overlay the histograms corresponding to specific monomer numbers
    i_ax=0
    ax_dens = axes[i_ax].twinx() #second y-axis for density (derived by KDE)

    #set up array of diameters (for bin edges and KDE-grid) and the fitting
    diam = np.logspace(-4,-2,50) #set up array of diameters (for bin edges and KDE-grid)


    fit_dic = dict() #initialize fit dictionary
    #separating sparse from valid data
    dens_lim = 1e-3 #0.005 #value of separating sparse from valid data
    linewidthinvalid=0.5;linewidthvalid=0.7 #define with which linewidth valid and invalid data should be plotted
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
    for i_prop,(prop,prop_unit) in enumerate(zip(["mass","area"],["kg","m2"])): #,"HWJussi","KCJussi"]):
        print "fitting and plotting: ",prop
        
        i_ax = i_prop+1 #+1 because we plotted already the histogram
        #change the scale (according to xscale_vec and yscale_vec)
        axes[i_ax].set_xscale("log") #axes[i_ax].set_xscale(xscale_vec[i_prop]) 
        axes[i_ax].set_yscale("log") #axes[i_ax].set_yscale(yscale_vec[i_prop])
        
        #plot the scatter plot (of mass and area for single particle) # use the colormap and the bounds which where defined before
        im = axes[i_ax].scatter(particle_dic["diam"],particle_dic[prop],s=1,c=particle_dic["N_monomer"],rasterized=True,norm=norm,cmap=cmap)
        
        #fit m-D and A-D relations
        axes[i_ax],particle_dic,fitline_allagg = __tools_for_processing_Jagg.fitting_wrapper(axes[i_ax],particle_dic,prop,N_mono_list,usecolors,fit_dic,diam,function="powerlaw",linewidthinvalid=0.05,linewidthvalid=0.7)
        
        #make labels
        axes[i_ax].set_xlabel("diameter / m")
        axes[i_ax].set_ylabel((" / ").join((prop,prop_unit))) #TODO: plot also the untis of these properties

        #change the axis
        axes[i_ax].set_ylim([0,np.array([1e-4,1e-3])[i_prop]]) #define the upper limit of the displayed axis
        axes[i_ax].grid(which="minor")
        #add colorbar
        cbar = fig.colorbar(im,ax=axes[i_ax])
        cbar.set_label("monomer number")
        #shift ticks to the middle of the color
        tick_locs = (bounds[:-1]) + 0.5*(bounds[1:]-bounds[:-1])
        cbar.set_ticks(tick_locs)
        cbar.set_ticklabels(N_mono_list)
        
        #add labels for the fit result to the plot
        axes[i_ax] = __tools_for_processing_Jagg.add_fitresult_text(axes[i_ax],fit_dic,prop,N_mono_list,function="powerlaw",hide_Nmono_fits=take_all_Nmono_bool)
        #add a legend for the fit lines which are not explained by the colorbar
        axes[i_ax].legend(handles=[fitline_allagg],loc="lower right", prop={'size': 8})
    #loop over different fall speed models
    for i_prop,prop in enumerate(["vterm_bohm","vterm_HW10","vterm_KC05","vterm_mitch_heym"]): #,"HWJussi","KCJussi"]):
        print "fitting and plotting: ",prop
        
        
        i_ax = i_prop+1+2 #+1 because we plotted already the histogram +2 because we plotted already mass and area
        #change the scale (according to xscale_vec and yscale_vec)
        axes[i_ax].set_xscale("log") #axes[i_ax].set_xscale(xscale_vec[i_prop]) 
        axes[i_ax].set_yscale("linear") #axes[i_ax].set_yscale(yscale_vec[i_prop])
        
        #plot the scatter plot of the individual particle # use the colormap and the bounds which where defined before
        im = axes[i_ax].scatter(particle_dic["diam"],particle_dic[prop],s=1,c=particle_dic["N_monomer"],rasterized=True,norm=norm,cmap=cmap) #,norm=colors.LogNorm(),cmap="gist_ncar")
        
        #calculate the terminal velocity based on the m-D and A-D fits and plot this line # use the colormap and the bounds which where defined before
        velocity_model = prop[6:] #get the name of the fall speed model from the property name (prop) which starts with "vterm"
        for i,N_mono_now in enumerate(N_mono_list):
            fit_dic["vterm_" + velocity_model + "_fitted_via_mD_AD_Nmono" + str(N_mono_now)] = __fallspeed_relations.calc_vterm(velocity_model,fit_dic["mass" + "_Nmono_" + str(N_mono_now)],fit_dic["diam"],fit_dic["area" + "_Nmono_" + str(N_mono_now)])
            fitline_Nmonospecific, = axes[i_ax].plot(fit_dic["diam"],fit_dic["vterm_" + velocity_model + "_fitted_via_mD_AD_Nmono" + str(N_mono_now)],marker='None',c=usecolors[i],linewidth=linewidthinvalid) #,norm=colors.LogNorm(),cmap="gist_ncar")
            fitline_Nmonospecific_dense_enough, = axes[i_ax].plot(fit_dic["diam"][fit_dic["diam_dens_enough_" + str(N_mono_now)]],fit_dic["vterm_" + velocity_model + "_fitted_via_mD_AD_Nmono" + str(N_mono_now)][fit_dic["diam_dens_enough_" + str(N_mono_now)]],marker='None',c=usecolors[i],linewidth=linewidthvalid) 
        #now for all aggregates
        fit_dic["vterm_" + velocity_model + "_fitted_via_mD_AD_Nmonoallagg"] = __fallspeed_relations.calc_vterm(velocity_model,fit_dic["mass" + "_Nmono_allagg"],fit_dic["diam"],fit_dic["area" + "_Nmono_allagg"])
        fitline_Nmonoallagg, = axes[i_ax].plot(fit_dic["diam"],fit_dic["vterm_" + velocity_model + "_fitted_via_mD_AD_Nmonoallagg"],marker='None',c="black",label="Nmono>1",linewidth=linewidthinvalid)
        fitline_Nmonoallagg_dense_enough, = axes[i_ax].plot(fit_dic["diam"][fit_dic["diam_dens_enough_allagg"]],fit_dic["vterm_" + velocity_model + "_fitted_via_mD_AD_Nmonoallagg"][fit_dic["diam_dens_enough_allagg"]],marker='None',c="black",label="Nmono>1",linewidth=linewidthvalid)

        #fit an Atlas type to the "mean v-D" relations (derived from the A-D and m-D fits)
        for i,str_monomer_agg in enumerate(["1","allagg"]): # make a power-law fit for the pristine crystals and the aggregate
            [fitresult,covar] = __tools_for_processing_Jagg.fit_data(fit_dic["diam"],fit_dic[prop +  "_fitted_via_mD_AD_Nmono" + str_monomer_agg],func="Atlas", weight=fit_dic["dens" + "_Nmono_" + str_monomer_agg])
            fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg] = fitresult #copy the Atlas-fit coefficients in the fit_dic dictionary
            #recalculate A-B*exp(-gam*D) from fit
            fit_dic[prop + "_Nmono_" + str_monomer_agg] = fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg][0]-fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg][1]*np.exp(-fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg][2]*fit_dic["diam"]) #calculate arrays of masses and areas based on the m-D fits
            
            if str_monomer_agg=="1": #no label
                Atlasfit_Nmono1, = axes[i_ax].plot(fit_dic["diam"],fit_dic[prop + "_Nmono_" + str_monomer_agg],marker='None',c=np.array(["b","black"])[i],linestyle="--",label="Nmono=1(weighted Atlas fit)",linewidth=linewidthinvalid)
                #from IPython.core.debugger import Tracer ; Tracer()()
                Atlasfit_Nmono1_dense_enough, = axes[i_ax].plot(fit_dic["diam"][fit_dic["diam_dens_enough_" + str(1)]],fit_dic[prop + "_Nmono_" + str_monomer_agg][fit_dic["diam_dens_enough_" + str(1)]],marker='None',c=np.array(["b","black"])[i],linestyle="--",label="Nmono=1(weighted Atlas fit)",linewidth=linewidthvalid)                
            elif str_monomer_agg=="allagg":
                Atlasfit_Nmonoallagg, = axes[i_ax].plot(fit_dic["diam"],fit_dic[prop + "_Nmono_" + str_monomer_agg],marker='None',c=np.array(["b","black"])[i],linestyle="--",label="Nmono>1(Atlas fit)",linewidth=linewidthinvalid)
                Atlasfit_Nmonoallagg_dense_enough, = axes[i_ax].plot(fit_dic["diam"][fit_dic["diam_dens_enough_allagg"]],fit_dic[prop + "_Nmono_" + str_monomer_agg][fit_dic["diam_dens_enough_allagg"]],marker='None',c=np.array(["b","black"])[i],linestyle="--",label="Nmono>1(weighted Atlas fit)",linewidth=linewidthvalid)

        #add labels for the fit result to the plot
        axes[i_ax] = __tools_for_processing_Jagg.add_fitresult_text(axes[i_ax],fit_dic,prop,[1],function="Atlas")
        #add a legend for the fit lines which are not explained by the colorbar
        axes[i_ax].legend(handles=[fitline_Nmonoallagg_dense_enough,Atlasfit_Nmonoallagg_dense_enough],loc="lower right", prop={'size': 8})
        
        #make labels
        axes[i_ax].set_xlabel("diameter / m")
        axes[i_ax].set_ylabel(prop + " / m s-1" ) #TODO: plot also the untis of these properties

        #change the axis
        axes[i_ax].set_ylim([0,2.0]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
        axes[i_ax].grid(which="major")

        #add colorbar
        cbar = fig.colorbar(im,ax=axes[i_ax])
        cbar.set_label("monomer number")
        #shift ticks to the middle of the color
        tick_locs = (bounds[:-1]) + 0.5*(bounds[1:]-bounds[:-1])
        cbar.set_ticks(tick_locs)
        cbar.set_ticklabels(N_mono_list)

    #add a colorbar also for the histogram (which is "stolen" from the scatter plots)
    for ax in [axes[0],ax_dens]:
        cbar = fig.colorbar(im,ax=ax)
        cbar.set_label("monomer number")
        #shift ticks to the middle of the color
        tick_locs = (bounds[:-1]) + 0.5*(bounds[1:]-bounds[:-1])
        cbar.set_ticks(tick_locs)
        cbar.set_ticklabels(N_mono_list)
    
    
    #save the plot (and open it)
    plt.tight_layout()
    dir_save = '/home/mkarrer/Dokumente/plots/Jagg/'
    if not os.path.exists(dir_save): #create direktory if it does not exists
        os.makedirs(dir_save)
    out_filestring = "Jagg_mD_AD_vD_" + particle_type + "_gridres" + str(grid_res)
    plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
    plt.savefig(dir_save + out_filestring + '.png', dpi=100)
    print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
    subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
    plt.clf()
    plt.close()
