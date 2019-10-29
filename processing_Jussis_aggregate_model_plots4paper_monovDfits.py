# coding: utf-8
#import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.pylab as pylab
import os
import subprocess
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


'''
this code reads in properties from Jussis aggregate model
and fits different functions to each property (m-D(monomer dependent and binary),A-D(monomer dependent and binary),v-D for all fallspeed models) and displays only v-D lines
'''
def read_and_plot(axes,hydro_model
    ="all",forw_mod=False,tumbling=False,called_by_main=False):
    '''
    read and plot the aggregate data and plot vterm (1. based on m/A-D power laws: forw_mod=False; 2. as seen by PIP: forw_mode=True)
    ARGUMENTS:
    axes: axes handles
    particle type: monomer type whcih is processed and plotted
    forw_mod: use this script as PIP forw model?
    '''

    def rchop(thestring, ending): #remove ending from a string
        if thestring.endswith(ending):
            return thestring[:-len(ending)]
        return thestring

    take_mono_prop_theor = True #so far only True is implemented
    skip_mono = True #dont use the monomer property at all #if you set this to true check the properties in __tools_for_processing_Jagg....

    if forw_mod==True:
        no_monomers_in_plot = True #skip plotting lines for monomers
        no_agg_in_plot = False #skip aggregates in plot

    else:
        no_monomers_in_plot = False
        no_agg_in_plot = True #skip aggregates in plot
        

    #show which types of fit
    show_lines=dict()
    show_lines["powerlaw_fits"] = False
    show_lines["highest_diam_for_fit"] = -999 #-999 for all -2 for increasing v range only
    show_lines["Atlas_fit"] = False #show the Atlas-type fit for cloud ice & snow
    #from previous model settings
    show_lines["SB_powerlaw"] = False #show also the old fit of Axels power-law
    v_small_notuse = 0.5  #ATTENTION: particle smaller than .5ms-1 are not used #ATTENTION: applied to B?hm only #ATTENTION: applied to percentiles only
    vnoisestdval = 0.2
    vnoisestdval2 = 0.4



    #add displayed lines to string to distinguish them in the saved files
    add_displayed_lines2string = '' 
    for key in show_lines.keys():
        if show_lines[key]:
            add_displayed_lines2string+= '_' + key
    particle_types = ["plate","dendrite","column","needle","rosette"]
    particle_typelabel = particle_types #["plate","dendrite","needle","column","rosette"]#,"Seifert14"] #,"mixcolumndend"] #"mixdendneedle""needle","column","plate","dendrite"] #["needle","column","plate","dendrite","bullet","rosette"] # ,"bullet"rosette
        
    for i_particle_type,particle_type in enumerate(particle_types):
        #get m-D from the assumptions in the aggregate model (this is used at various places)
        a,b,c,d = __tools_for_processing_Jagg.calc_mD_AD_coeffs(particle_type,skip_mono=skip_mono)

            
        #####################
        ######START: fitting
        #####################
        #initialize an array which contains all fit results
        fit_dic = dict() #initialize fit dictionary
        #set up array of diameters for displaying the fits
        low_diam_log_detailed= -4; high_diam_log_detailed_plotrange=-1.3979400086720375
        if show_lines["highest_diam_for_fit"]==-999:
            high_diam_log_detailed=-1.3979400086720375 #np.log10(3e-2)=-1.5228787452803376 #np.log10(4e-2)=-1.3979400086720375
        else:
            high_diam_log_detailed=show_lines["highest_diam_for_fit"]
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
        #calculate monomer properties   
        fit_dic["mass_Nmono_1"]=a*fit_dic["diam"]**b
        fit_dic["area_Nmono_1"]=c*fit_dic["diam"]**d
        ###
        #START: calculate power-law and Atlas-fit based on the exact v-D line
        ###
        
        #first calculate the exact v-D
        
        
        #loop over different terminal velocity models in which different possible fits (power-law, Atlas type) are calculated
        for i_prop,prop in enumerate(["vterm_bohm","vterm_HW10","vterm_KC05","vterm_mitch_heym"]): #,"vterm_mitch_heym"]): #,"HWJussi","KCJussi"]):
            print "fitting and plotting: ",prop
            #first calculate the exact v-D
            for i_mono,(str_monomer_agg,corr_cat) in enumerate(zip(["1"],["cloud ice"])):
                fit_dic[prop +  "_fitted_via_mD_AD_Nmono" + str_monomer_agg] = __fallspeed_relations.calc_vterm(prop[6:],fit_dic["mass_Nmono_" + str_monomer_agg],fit_dic["diam"],fit_dic["area_Nmono_" + str_monomer_agg])
                
            '''
            #calculate velocity also for percentiles of m and A
            for percentile in [25,50,75]:
                fit_dic[prop + "_" + str(percentile) + "_via_mD_AD_Nmono" + str_monomer_agg] = __fallspeed_relations.calc_vterm(prop[6:],fit_dic["mass_" + str(percentile) + "perc_Nmono_allagg"],diam_percentile,fit_dic["area_" + str(percentile) + "perc_Nmono_allagg"]) 
            '''
            #fit an Atlas type to the "mean v-D" relations (derived from the A-D and m-D fits)
            for i_mono,(str_monomer_agg,corr_cat) in enumerate(zip(["1"],["cloud ice"])): 
                if (no_monomers_in_plot or particle_type=="Seifert14") and i_mono==0:
                    continue
                ##Atlas-type
                [fitresult,covar] = __tools_for_processing_Jagg.fit_data(fit_dic["diam"],fit_dic[prop +  "_fitted_via_mD_AD_Nmono" + str_monomer_agg],func="Atlas") #, weight=fit_dic["dens" + "_Nmono_" + str_monomer_agg])
                #from IPython.core.debugger import Tracer ; Tracer()()


                [fitresult_Deq,covar_Deq] = __tools_for_processing_Jagg.fit_data(((6.*fit_dic["mass_Nmono_" + str_monomer_agg][0]/(np.pi*1000.))**(1./3.)),fit_dic[prop +  "_fitted_via_mD_AD_Nmono" + str_monomer_agg],func="Atlas") #, weight=fit_dic["dens" + "_Nmono_" + str_monomer_agg])
                
                fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg] = fitresult #copy the Atlas-fit (as a function of D_max) coefficients in the fit_dic dictionary
                fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg + "_Deq"] = fitresult_Deq #copy the Atlas-fit (as a function of D_eq) coefficients in the fit_dic dictionary

                #recalculate v(D)= A-B*exp(-gam*D) from fit
                fit_dic[prop + "_Nmono_" + str_monomer_agg] = fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg][0]-fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg][1]*np.exp(-fit_dic[prop + "_coeff_Nmono_" + str_monomer_agg][2]*fit_dic["diam_full"]) #calculate arrays of masses and areas based on the m-D fits
                
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
        linewidth_white=1.0
        markersize=2.0
        #         plate,dendrite,column,needle,rosette,mix1,mix2
        linestyles=["-","--",      ":",  "-.",    ":"  ,"-","-"]
        markers=   ["", ""  ,      "" ,   ""  ,    "o" ,"o","x"]
        ###
        #START: plot the prop(m or D) vs diameter as scatter and the 2D-fit
        ###
        
        
        ####
        #plot vterm (and its fits)
        ####
        for i_prop,prop in enumerate(["vterm_bohm","vterm_HW10","vterm_KC05","vterm_mitch_heym"]): #,"vterm_mitch_heym"]): #,"HWJussi","KCJussi"]):
            if hydro_model!="all" and hydro_model!=prop:
                continue
            print "fitting and plotting: ",prop
            
            #define current axis
            i_ax+=1         
            #change the scale (according to xscale_vec and yscale_vec)
            axes[i_ax].set_xscale("log") #axes[i_ax].set_xscale(xscale_vec[i_prop]) 
            axes[i_ax].set_yscale("linear") #axes[i_ax].set_yscale(yscale_vec[i_prop])
            
            #add a line for the monomers
            #fitlines_monomers = axes[i_ax].semilogx(fit_dic["diam"],fit_dic[prop +  "_fitted_via_mD_AD_Nmono1"],c="b",linestyle=["-","--","-.",":"][i_particle_type],linewidth=linewidth,marker='o',markevery=20)
            fitlines_monomers = axes[i_ax].semilogx(fit_dic["diam"],fit_dic[prop +  "_fitted_via_mD_AD_Nmono1"],c="b",linestyle=linestyles[i_particle_type],marker=markers[i_particle_type],markersize=markersize,markevery=50,linewidth=linewidth)

            
            #add a line for the monomers
            #fitlines_monomer = axes[i_ax].semilogx(fit_dic["diam"],fit_dic[prop +  "_fitted_via_mD_AD_Nmono1"],c=["b","b","b","b","b","b"][i_particle_type],linestyle=["-","--","--",":","-.","--","-","-","-."][i_particle_type],marker=["","","o","","","x","d","P",""][i_particle_type],markevery=50,markersize=markersize,linewidth=linewidth)
                ##add fitted lines
                #monomers
            if show_lines["Atlas_fit"]:
                #Atlas
                Atlasfit_monomers, = axes[i_ax].semilogx(fit_dic["diam_full"],fit_dic[prop + "_Nmono_1"],marker='None',markersize=markersize,c="b",linestyle="--",label="Atlas-type cl. ice "  + particle_type,linewidth=2*linewidth)
            if show_lines["powerlaw_fits"]:
                #power-law
                powerlaw_monomers, = axes[i_ax].semilogx(fit_dic["diam_full"], fit_dic[prop + "_Nmono_1_powerlaw"],marker='None',c="b",linestyle="-.",label="powerlaw cl. ice "  + particle_type,linewidth=2*linewidth)
                
            #make labels
            axes[i_ax].set_xlabel("$D_{max}$ [m]")
            axes[i_ax].set_ylabel("$v_{term}$ [m/s]" ) #TODO: plot also the untis of these properties

            #change the axis
            axes[i_ax].set_xlim([10**low_diam_log_detailed,10**high_diam_log_detailed_plotrange]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
            axes[i_ax].set_ylim([0,2.5]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
            axes[i_ax].grid(b=True,which="both",axis="both")

    for i_ax,ax_now in enumerate(axes):
        for i_particle_type,particle_type in enumerate(particle_types):
            particle_type_label=particle_type

            ax_now.plot(np.nan,np.nan,color="b",linestyle=linestyles[i_particle_type],marker=markers[i_particle_type],linewidth=linewidth,label=particle_typelabel[i_particle_type],markevery=50, markersize=markersize)
                
        #show legend
        if called_by_main: 
            axes[i_ax].legend(ncol=1,loc='upper left') #, bbox_to_anchor=(1.0, 0.5)) #ncol=1,loc='center left', bbox_to_anchor=(1.0, 0.5))


    return axes






if __name__ == '__main__':

    ##general settings for the plotting
    number_of_plots = 4

    #optimize the appearance of the plot (figure size, fonts)
    aspect_ratio=1./4.
    legend_fontsize='x-large' #'medium'
    #increase font sizes
    params = {'legend.fontsize': legend_fontsize,
        'axes.labelsize': 'xx-large', #size: Either a relative value of 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large' or an absolute font size, e.g., 12
        'axes.titlesize':'xx-large',
        'xtick.labelsize':'xx-large',
        'ytick.labelsize':'xx-large'}
    pylab.rcParams.update(params)
    #define figure
    figsize_height = 1./aspect_ratio*number_of_plots
    fig,axes = plt.subplots(nrows=number_of_plots, ncols=1, figsize=(12.5,figsize_height))
    
    #CALL THE MAIN PART
    tumbling=False; forw_mod=False
    axes = read_and_plot(axes,hydro_model="all",forw_mod=forw_mod,tumbling=tumbling,called_by_main=False)

    #save the plot (and open it)
    plt.tight_layout()

    dir_save = '/home/mkarrer/Dokumente/plots/Jagg/'
    if not os.path.exists(dir_save): #create direktory if it does not exists
        os.makedirs(dir_save)
    if tumbling:
        tumbling_add="_wtumbling"
    else: 
        tumbling_add="_wotumbling"
        if forw_mod:
            sideproj_add="_side_proj"
        else: 
            sideproj_add=""
        #all_particle_types = "".join(particle_types)
        out_filestring = "vaggvsDfits_" + '_' + tumbling_add + '_' +  sideproj_add

        plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
        plt.savefig(dir_save + out_filestring + '.png', dpi=100)
        print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
        subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])


        #save the axis individually
        ax_description_array=["vterm_bohm","vterm_HW10","vterm_KC05"]
        for save_axis,ax_description in enumerate(ax_description_array):
            if forw_mod==True: #in this case we need only bohm
                continue
            print "saving subfigure",save_axis,ax_description
            try: #clear label of axis before and recover current one
                if forw_mod:
                    axes[save_axis].set_xlabel("diameter $D_{max,side}$ / m")            
                else:
                    axes[save_axis].set_xlabel("diameter $D_{max}$ / m")
                axes[save_axis-1].set_xlabel("")
            except:
                pass
            extent = axes[save_axis].get_window_extent().transformed(fig.dpi_scale_trans.inverted())

            fig.savefig('/home/mkarrer/Dokumente/plots/tmp.pdf',bbox_inches=extent.expanded(2.3, 1.45),dpi=400) #,bbox_inches=extent.expanded(2.3, 1.3),dpi=400)bbox_inches=extent.expanded(1.6, 1.45),dpi=400) are width and height expanded around center

            subprocess.call('cp ' + '/home/mkarrer/Dokumente/plots/tmp.pdf /home/mkarrer/Dokumente/plots/4paper/' + out_filestring + '_' + ax_description + '.pdf',shell=True)
