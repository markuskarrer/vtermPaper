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
from IPython.core.debugger import Tracer ; debug = Tracer()
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
        no_monomers_in_plot = True
        no_agg_in_plot = False #skip aggregates in plot
        

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


    ####
    #START: read in Hyytiala data from Annakaisa
    ####

    hyytiala_data = "/home/mkarrer/Dokumente/insitu_vterm/Hyytiala_MedianV_perc25_75.csv"
    vD_hyyt_dic = dict() #dictionary which contains arrays of the observed variables from the Hyytiala site
    vD_hyyt_dic["Dmax"] = []        #maximum dimension [mm]
    vD_hyyt_dic["v_25perc"] = []    #25% percenticle of the velocity
    vD_hyyt_dic["v_75perc"] = []    #75% percenticle of the velocity
    vD_hyyt_dic["Dmax_wNAN"] = []   #maximum dimension, when including NaN velocity
    vD_hyyt_dic["v"] = []           #median velocity
    with open(hyytiala_data,"rb") as txtfile: #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in
        prop_reader = csv.reader(filter(lambda row: row[0]!='D',txtfile), delimiter=',', quoting=csv.QUOTE_NONNUMERIC, lineterminator=os.linesep) #row[0]!='#': skip the header; quoting avoids '' for formatted string; lineterminator avoids problems with system dependend lineending format https://unix.stackexchange.com/questions/309154/strings-are-missing-after-concatenating-two-or-more-variable-string-in-bash?answertab=active#tab-top
        #for i_row,row_content in enumerate(prop_reader):
        for row_content in prop_reader:
        
            #replace missing values and matlab "nan"'s with python "NaN"s
            for i_entry,entry in enumerate(row_content):
                if entry in ("","nan"):
                    row_content[i_entry] = "NaN"

            vD_hyyt_dic["Dmax"] = np.append(vD_hyyt_dic["Dmax"],float(row_content[0]))
            vD_hyyt_dic["v_25perc"] = np.append(vD_hyyt_dic["v_25perc"],float(row_content[1]))
            vD_hyyt_dic["v_75perc"] = np.append(vD_hyyt_dic["v_75perc"],float(row_content[2]))
            vD_hyyt_dic["Dmax_wNAN"] = np.append(vD_hyyt_dic["Dmax_wNAN"],float(row_content[3]))
            vD_hyyt_dic["v"] = np.append(vD_hyyt_dic["v"],float(row_content[4]))

        
    ####
    #END: read in Hyytiala data from Annakaisa
    ####

    ####
    #START: read in CARE data from Annakaisa
    ####

    hyytiala_data = "/home/mkarrer/Dokumente/insitu_vterm/CARE_2014_2017_V_D_median_25_75_larger05.csv"
    vD_CARE_dic = dict() #dictionary which contains arrays of the observed variables from the Hyytiala site
    vD_CARE_dic["Dmax"] = []        #maximum dimension [mm]
    vD_CARE_dic["v_larger0"] = []           #median velocity
    vD_CARE_dic["v_25perc_larger0"] = []    #25% percenticle of the velocity
    vD_CARE_dic["v_75perc_larger0"] = []    #75% percenticle of the velocity
    vD_CARE_dic["v_larger05"] = []           #median velocity
    vD_CARE_dic["v_25perc_larger05"] = []    #25% percenticle of the velocity
    vD_CARE_dic["v_75perc_larger05"] = []    #75% percenticle of the velocity
    vD_CARE_dic["Dmax_wNAN"] = []   #maximum dimension, when including NaN velocity

    with open(hyytiala_data,"rb") as txtfile: #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in
        prop_reader = csv.reader(filter(lambda row: row[0]!='D',txtfile), delimiter=';', quoting=csv.QUOTE_NONNUMERIC, lineterminator=os.linesep) #row[0]!='#': skip the header; quoting avoids '' for formatted string; lineterminator avoids problems with system dependend lineending format https://unix.stackexchange.com/questions/309154/strings-are-missing-after-concatenating-two-or-more-variable-string-in-bash?answertab=active#tab-top
        #for i_row,row_content in enumerate(prop_reader):
        for row_content in prop_reader:
        
            #replace missing values and matlab "nan"'s with python "NaN"s
            for i_entry,entry in enumerate(row_content):
                if entry in ("","nan"):
                    row_content[i_entry] = "NaN"

            vD_CARE_dic["Dmax"] = np.append(vD_CARE_dic["Dmax"],float(row_content[0]))
            
            vD_CARE_dic["Dmax"][vD_CARE_dic["Dmax"]>10] = np.nan #mask out noisy values (>... mm)

            
            vD_CARE_dic["v_larger0"] = np.append(vD_CARE_dic["v_larger0"],float(row_content[1]))
            vD_CARE_dic["v_25perc_larger0"] = np.append(vD_CARE_dic["v_25perc_larger0"],float(row_content[2]))
            vD_CARE_dic["v_75perc_larger0"] = np.append(vD_CARE_dic["v_75perc_larger0"],float(row_content[3]))
            #vD_CARE_dic["Dmax"] = np.append(vD_CARE_dic["Dmax"],float(row_content[3]))
            vD_CARE_dic["v_larger05"] = np.append(vD_CARE_dic["v_larger05"],float(row_content[5]))
            vD_CARE_dic["v_25perc_larger05"] = np.append(vD_CARE_dic["v_25perc_larger05"],float(row_content[6]))
            vD_CARE_dic["v_75perc_larger05"] = np.append(vD_CARE_dic["v_75perc_larger05"],float(row_content[7]))
            #vD_CARE_dic["Dmax_wNAN"] = np.append(vD_CARE_dic["Dmax_wNAN"],float(row_content[3]))
        
    ####
    #END: read in Hyytiala data from Annakaisa
    ####
    
    ####
    #START: read in CARE (??) data from Annakaisa (filtered by w<4m/s)
    ####
    data_dic = dict()
    hyytiala_data_wfilt = "/home/mkarrer/Dokumente/insitu_vterm/Hyytiala_MedianV_perc25_75_withwind.csv"
    data_dic["vD_CAREwindfilt_dic"] = dict() #dictionary which contains arrays of the observed variables from the Hyytiala site
    #b2p1 means the filtering for unrimed applied b_mass>2.1 
    data_dic["vD_CAREwindfilt_dic"]["Dmax"] = []        #maximum dimension [mm]
    data_dic["vD_CAREwindfilt_dic"]["v_b2p1"] = []           #median velocity
    data_dic["vD_CAREwindfilt_dic"]["v_25perc_b2p1"] = []    #25% percenticle of the velocity
    data_dic["vD_CAREwindfilt_dic"]["v_75perc_b2p1"] = []    #75% percenticle of the velocity
    data_dic["vD_CAREwindfilt_dic"]["numpart_b2p1"] = []     #number of particles with b_ms>2.1


    data_dic["vD_CAREwindfilt_dic"]["v_b2p2"] = []           #median velocity
    data_dic["vD_CAREwindfilt_dic"]["v_25perc_b2p2"] = []    #25% percenticle of the velocity
    data_dic["vD_CAREwindfilt_dic"]["v_75perc_b2p2"] = []    #75% percenticle of the velocity
    data_dic["vD_CAREwindfilt_dic"]["numpart_b2p2"] = []     #number of particles with b_ms>2.2

    with open(hyytiala_data_wfilt,"rb") as txtfile: #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in
        prop_reader = csv.reader(filter(lambda row: row[0]!='D',txtfile), delimiter=',', quoting=csv.QUOTE_NONNUMERIC, lineterminator=os.linesep) #row[0]!='#': skip the header; quoting avoids '' for formatted string; lineterminator avoids problems with system dependend lineending format https://unix.stackexchange.com/questions/309154/strings-are-missing-after-concatenating-two-or-more-variable-string-in-bash?answertab=active#tab-top
        #for i_row,row_content in enumerate(prop_reader):
        for row_content in prop_reader:

            #replace missing values and matlab "nan"'s with python "NaN"s
            for i_entry,entry in enumerate(row_content):
                if entry in ("","nan"):
                    row_content[i_entry] = "NaN"

            data_dic["vD_CAREwindfilt_dic"]["Dmax"] = np.append(data_dic["vD_CAREwindfilt_dic"]["Dmax"],float(row_content[0]))
            if row_content[4]<1000:
                data_dic["vD_CAREwindfilt_dic"]["Dmax"][-1] = np.nan #mask out noisy values (>... mm) !ATTENTION: this is now done by modifing the file
            data_dic["vD_CAREwindfilt_dic"]["v_b2p1"] = np.append(data_dic["vD_CAREwindfilt_dic"]["v_b2p1"],float(row_content[1]))
            data_dic["vD_CAREwindfilt_dic"]["v_25perc_b2p1"] = np.append(data_dic["vD_CAREwindfilt_dic"]["v_25perc_b2p1"],float(row_content[2]))
            data_dic["vD_CAREwindfilt_dic"]["v_75perc_b2p1"] = np.append(data_dic["vD_CAREwindfilt_dic"]["v_75perc_b2p1"],float(row_content[3]))
            #data_dic["vD_CAREwindfilt_dic"]["Dmax"] = np.append(data_dic["vD_CAREwindfilt_dic"]["Dmax"],float(row_content[3]))
            data_dic["vD_CAREwindfilt_dic"]["v_b2p2"] = np.append(data_dic["vD_CAREwindfilt_dic"]["v_b2p2"],float(row_content[5]))
            data_dic["vD_CAREwindfilt_dic"]["v_25perc_b2p2"] = np.append(data_dic["vD_CAREwindfilt_dic"]["v_25perc_b2p2"],float(row_content[6]))
            data_dic["vD_CAREwindfilt_dic"]["v_75perc_b2p2"] = np.append(data_dic["vD_CAREwindfilt_dic"]["v_75perc_b2p2"],float(row_content[7]))
            #data_dic["vD_CAREwindfilt_dic"]["Dmax_wNAN"] = np.append(data_dic["vD_CAREwindfilt_dic"]["Dmax_wNAN"],float(row_content[3]))

    ####
    #END: read in Hyytiala data from Annakaisa (wind filtered u>4m/s)
    ####


    #add displayed lines to string to distinguish them in the saved files
    add_displayed_lines2string = '' 
    for key in show_lines.keys():
        if show_lines[key]:
            add_displayed_lines2string+= '_' + key
    if forw_mod:
        particle_types = ["plate","dendrite","column","needle","mixcoldend1","mixcolumndend"]
        #particle_types = ["plate","plate_sideproj","plate_vnoise0p2vcut","plate_vnoise0p4vcut"] #,"plate_vnoisesideprojvcut"]
        #particle_types = ["dendrite","dendrite_vcut","dendrite_vnoisevcut"] #,"dendrite_vnoisesideprojvcut"]
        #particle_types = ["mixdendneedle","mixdendneedle_vcut","mixdendneedle_vnoisevcut"] #,"mixdendneedle_vnoisesideprojvcut"]
        #particle_types = ["needle","needle_vcut","needle_vnoisevcut"] #,"needle_vnoisesideprojvcut"]
        #particle_types = ["column","column_vcut","column_vnoisevcut"] #,"column_vnoisesideprojvcut"]
        #particle_types = ["mixcolumndend","mixcolumndend_sideproj","mixcolumndend_vnoise0p2vcutsideprojvcut","mixcolumndend_vnoise0p4vcut"] #,"mixcolumndend_vcut","mixcolumndend_vnoisevcut"]

        #particle_types = ["column","column_sideproj","column_sideprojvcut"]
        #particle_types = ["column","column_areasideproj"] #,"column_sideproj","column_sideprojvcut"]

        #particle_types = ["mixcolumndend","mixcolumndend_sideproj","mixcolumndend_sideprojvcut"]

    else:
        particle_types = ["plate","dendrite","column","needle","mixcoldend1","mixcolumndend"]
        
    for i_particle_type,particle_type in enumerate(particle_types):
        area_side_proj=False
        
        print "run: ", particle_type

        if "_vcut" in particle_type:
            use_side_proj_Dmax=False
            particle_type = rchop(particle_type, "_vcut") #particle_type.strip("_sideproj")
            v_small=v_small_notuse
            vnoisestd=0
        elif "_vnoise0p2vcutsideprojvcut" in particle_type:
            use_side_proj_Dmax=True
            particle_type = rchop(particle_type, "_vnoise0p2vcutsideprojvcut")
            v_small=v_small_notuse
            vnoisestd=vnoisestdval
        elif "_vnoise0p2vcut" in particle_type:
            use_side_proj_Dmax=False
            particle_type = rchop(particle_type, "_vnoise0p2vcut")
            v_small=v_small_notuse
            vnoisestd=vnoisestdval
        elif "_vnoise0p4vcut" in particle_type:
            use_side_proj_Dmax=False
            particle_type = rchop(particle_type, "_vnoise0p4vcut")
            v_small=v_small_notuse
            vnoisestd=vnoisestdval2
        elif "_vnoisesideprojvcut" in particle_type:
            use_side_proj_Dmax=True
            particle_type = rchop(particle_type, "_vnoisesideprojvcut")
            v_small=v_small_notuse
            vnoisestd=vnoisestdval
        elif "_sideprojvcut" in particle_type:
            use_side_proj_Dmax=True
            particle_type = rchop(particle_type, "_sideprojvcut") #particle_type.strip("_sideproj")
            v_small=v_small_notuse
            vnoisestd=0
        elif "_sideproj" in particle_type:
            use_side_proj_Dmax=True
            particle_type = rchop(particle_type, "_sideproj") #particle_type.strip("_sideproj")
            v_small=0.0
            vnoisestd=0
        elif "_areasideproj" in particle_type: #use a side projected area as in Szyrmer&Zawadzki10 or vonLerber17
            use_side_proj_Dmax=False
            area_side_proj=True
            particle_type = rchop(particle_type, "_areasideproj") #particle_type.strip("_sideproj")
            v_small=0.0
            vnoisestd=0
        else:
            use_side_proj_Dmax=False
            v_small=0.0
            vnoisestd=0

        if particle_type=="Seifert14":
            '''
                            particle_dic["particle_type"] = np.append(particle_dic["particle_type"],particle_type)
                            particle_dic["N_monomer"] = np.append(particle_dic["N_monomer"],row_content[0]) #python indices start with 0
                            particle_dic["mass"] = np.append(particle_dic["mass"],row_content[1])

                            particle_dic["diam"] = np.append(particle_dic["diam"],row_content[3]) #np.sqrt(row_content[4]/np.pi)) #ATTENTION: test take area equivalent diameter

                            particle_dic["area"] = np.append(particle_dic["area"],row_content[4]) #,row_content[2]) -[2] for tumbling [4] for disabled tumbling

                            particle_dic["as_ratio"] = np.append(particle_dic["as_ratio"],row_content[5]) 
                            particle_dic["area_partal10std"] = np.append(particle_dic["area_partal10std"],row_content[6]) 
                            particle_dic["area_partal20std"] = np.append(particle_dic["area_partal20std"],row_content[7])
                            particle_dic["area_partal40std"] = np.append(particle_dic["area_partal40std"],row_content[2])
                            particle_dic["area_partal60std"] = np.append(particle_dic["area_partal60std"],row_content[8])
            '''                
            particle_dic = dict()
            diam_selfsim = np.logspace(-4,np.log10(4e-2),1000)
            particle_dic["diam"] = diam_selfsim
            particle_dic["particle_type"] = particle_type*len(particle_dic["diam"])
            particle_dic["mass"] = 0.038*particle_dic["diam"]**2
            particle_dic["area"] = 0.45*particle_dic["diam"]**2 #based in Field et. al (2008)
            #particle_dic["area"] = 0.2285*particle_dic["diam"]**1.88 #Mitchell (1996)
            
            N_mono_list = np.array([1,2])
            #from IPython.core.debugger import Tracer ; Tracer()()
            
        else:
            #read the properties of the individual particles from the files
            particle_dic,N_mono_list = __tools_for_processing_Jagg.read_particle_prop_files(prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation/",
                                                                                    #prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation_local/tumbling_asratio/", #"/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation/",
                                                                                    D_small_notuse=1e-4, #ATTENTION: particle smaller than D_small_notuse are not considered (e.g. because of resolution issues)
                                                                                    N_small_notuse=1, #ATTENTION:  monomer numbers smaller than N_small_notuse are not considered
                                                                                    grid_res_array = [5e-6,10e-6], #[5e-6,10e-6], #[1e-6,5e-6,10e-6], #array of grid resolutions (defines the files which are read in)
                                                                                    particle_type = particle_type, #define the particle_type
                                                                                    test_with_small_sample = False
                                                                                    )
        #show the dictionary (e.g. to check that it's not empty because the path was wrong)
        print particle_dic

        #get m-D from the assumptions in the aggregate model (this is used at various places)
        if particle_type!="Seifert14":
            a,b,c,d = __tools_for_processing_Jagg.calc_mD_AD_coeffs(particle_type,skip_mono=skip_mono)

        if tumbling:
            particle_dic["area"] = particle_dic["area_partal40std"] #use the tumbled area instead
        if use_side_proj_Dmax:
            particle_dic["diam_obs"] = particle_dic["Dmax_xz"]
        else:
            particle_dic["diam_obs"] = particle_dic["diam"]

            
        #calculate the terminal velocitys of the simulated particles individually for each available velocity model
        for velocity_model in ["HW10","KC05","bohm","mitch_heym"]:
            #from IPython.core.debugger import Tracer ; Tracer()()
            if area_side_proj: #use a side projected area as in Szyrmer&Zawadzki10 or vonLerber17
                particle_dic["vterm_" + velocity_model] = __fallspeed_relations.calc_vterm(velocity_model,particle_dic["mass"],particle_dic["diam"],particle_dic["area_aligned_side"])+np.random.randn(particle_dic["mass"].shape[0])*vnoisestd

            else:
                particle_dic["vterm_" + velocity_model] = __fallspeed_relations.calc_vterm(velocity_model,particle_dic["mass"],particle_dic["diam"],particle_dic["area"])+np.random.randn(particle_dic["mass"].shape[0])*vnoisestd
            #particle_dic["vterm_" + velocity_model] = __fallspeed_relations.calc_vterm(velocity_model,particle_dic["mass"],particle_dic["diam"],particle_dic["area_aligned_side"])+np.random.randn(particle_dic["mass"].shape[0])*vnoisestd
        
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
            '''
            #calculate also percentiles (binned)
            diam_percentile = np.logspace(low_diam_log_detailed,high_diam_log_detailed,50)
            for i_diam,diam_now in enumerate(diam_percentile[:-1]):
                #from IPython.core.debugger import Tracer ; Tracer()()
                
                prop_in_bin = particle_dic[prop][((particle_dic["diam"]>diam_percentile[i_diam]) & (particle_dic["diam"]<diam_percentile[i_diam+1]))]
                for percentile in [25,50,75]:
                    if prop_in_bin.shape[0]>1:
                        fit_dic[prop + "_" + str(percentile) + "perc_Nmono_allagg"] =  np.percentile(prop_in_bin,percentile)
                    else:
                        fit_dic[prop + "_" + str(percentile) + "perc_Nmono_allagg"] = np.nan
            '''
            
            if take_mono_prop_theor and particle_type!="Seifert14" and not skip_mono: #overwrite fit for monomers with theoretical value
            
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
        for i_prop,prop in enumerate(["vterm_bohm","vterm_HW10","vterm_KC05","vterm_mitch_heym"]): #,"vterm_mitch_heym"]): #,"HWJussi","KCJussi"]):
            print "fitting and plotting: ",prop
            #first calculate the exact v-D
            for i_mono,(str_monomer_agg,corr_cat) in enumerate(zip(["1","allagg"],["cloud ice","snow"])):
                if (no_monomers_in_plot or particle_type=="Seifert14") and i_mono==0:
                    continue
                fit_dic[prop +  "_fitted_via_mD_AD_Nmono" + str_monomer_agg] = __fallspeed_relations.calc_vterm(prop[6:],fit_dic["mass_Nmono_" + str_monomer_agg],fit_dic["diam"],fit_dic["area_Nmono_" + str_monomer_agg])
                
            '''
            #calculate velocity also for percentiles of m and A
            for percentile in [25,50,75]:
                fit_dic[prop + "_" + str(percentile) + "_via_mD_AD_Nmono" + str_monomer_agg] = __fallspeed_relations.calc_vterm(prop[6:],fit_dic["mass_" + str(percentile) + "perc_Nmono_allagg"],diam_percentile,fit_dic["area_" + str(percentile) + "perc_Nmono_allagg"]) 
            '''
            #calculate the percentiles of the velocities directly
            particle_dic["vterm_" + velocity_model]
            #calculate also percentiles (binned)
            #diam_percentile = np.logspace(low_diam_log_detailed,high_diam_log_detailed,20)
            diam_percentile = np.arange(0.02e-3,29.82e-3,0.28e-3)
            diam_percentile_center = (diam_percentile[:-1]+diam_percentile[1:])/2.
            #initialize percentiles
            for percentile in [25,50,75]:
                fit_dic[prop + "_" + str(percentile) + "perc_Nmono_allagg"] = np.zeros_like(diam_percentile[:-1])
            for i_diam,diam_now in enumerate(diam_percentile[:-1]):

                prop_in_bin = particle_dic[prop][(((particle_dic["diam_obs"]>diam_percentile[i_diam]) & (particle_dic["diam_obs"]<diam_percentile[i_diam+1])) & (particle_dic["vterm_bohm"]>v_small))]
                for percentile in [25,50,75]:
                    if prop_in_bin.shape[0]>1:
                        fit_dic[prop + "_" + str(percentile) + "perc_Nmono_allagg"][i_diam] =  np.percentile(prop_in_bin,percentile)
                    else:
                        fit_dic[prop + "_" + str(percentile) + "perc_Nmono_allagg"][i_diam] = np.nan
                    #print diam_now,fit_dic[prop + "_" + str(percentile) + "perc_Nmono_allagg"]

            #fit an Atlas type to the "mean v-D" relations (derived from the A-D and m-D fits)
            for i_mono,(str_monomer_agg,corr_cat) in enumerate(zip(["1","allagg"],["cloud ice","snow"])): 
                if (no_monomers_in_plot or particle_type=="Seifert14") and i_mono==0:
                    continue
                ##Atlas-type
                [fitresult,covar] = __tools_for_processing_Jagg.fit_data(fit_dic["diam"],fit_dic[prop +  "_fitted_via_mD_AD_Nmono" + str_monomer_agg],func="Atlas") #, weight=fit_dic["dens" + "_Nmono_" + str_monomer_agg])
                #from IPython.core.debugger import Tracer ; Tracer()()


                [fitresult_Deq,covar_Deq] = __tools_for_processing_Jagg.fit_data(((6.*fit_dic["mass_coeff_Nmono_" + str_monomer_agg][0]*fit_dic["diam"]**fit_dic["mass_coeff_Nmono_" + str_monomer_agg][1]/(np.pi*1000.))**(1./3.)),fit_dic[prop +  "_fitted_via_mD_AD_Nmono" + str_monomer_agg],func="Atlas") #, weight=fit_dic["dens" + "_Nmono_" + str_monomer_agg])
                
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
        markersize=5.0
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
        
        
        #         plate,dendrite,column,needle,mix1,mix2
        linestyles=["-","--",      ":",  "-.",    "-","-"]
        markers=   ["", ""  ,      "" ,   ""  ,   "o","x"]
        ####
        #plot vterm (and its fits)
        ####
        for i_prop,prop in enumerate(["vterm_bohm","vterm_HW10","vterm_KC05","vterm_mitch_heym"]): #,"vterm_mitch_heym"]): #,"HWJussi","KCJussi"]):
            if hydro_model!=all and hydro_model!=prop:
                continue
            
            print "fitting and plotting: ",prop
            
            #define current axis
            i_ax+=1         
            #change the scale (according to xscale_vec and yscale_vec)
            axes[i_ax].set_xscale("log") #axes[i_ax].set_xscale(xscale_vec[i_prop]) 
            axes[i_ax].set_yscale("linear") #axes[i_ax].set_yscale(yscale_vec[i_prop])
            
            #add a line for the monomers
            #fitlines_monomers = axes[i_ax].semilogx(fit_dic["diam"],fit_dic[prop +  "_fitted_via_mD_AD_Nmono1"],c="b",linestyle=["-","--","-.",":"][i_particle_type],linewidth=linewidth,marker='o',markevery=20)
            if not (("_sideproj" in ''.join(particle_types)) or no_monomers_in_plot ):
                if particle_type!="Seifert14":
                    fitlines_monomers = axes[i_ax].semilogx(fit_dic["diam"],fit_dic[prop +  "_fitted_via_mD_AD_Nmono1"],c="b",linestyle=linestyles[i_particle_type],marker=markers[i_particle_type],markevery=50,linewidth=linewidth)

            
            #add a line for the aggregates
            if (not forw_mod) and (not no_agg_in_plot):
                #from IPython.core.debugger import Tracer ; Tracer()()
                print i_particle_type
                fitlines_allagg = axes[i_ax].semilogx(fit_dic["diam"],fit_dic[prop +  "_fitted_via_mD_AD_Nmonoallagg"],c=["g","g","g","g","black","black","g","r"][i_particle_type],linestyle=linestyles[i_particle_type],marker=markers[i_particle_type],markevery=50,markersize=markersize,linewidth=linewidth,fillstyle='none')
                #add 50% percentile for all aggregates
                #if particle_type!="Seifert14":
                #    axes[i_ax].fill_between(diam_percentile[:-1],fit_dic[prop + "_" + str(25) + "perc_Nmono_allagg"],fit_dic[prop + "_" + str(75) + "perc_Nmono_allagg"],color="g",linestyle=["-","--",":",":","--","--","-","-","-."][i_particle_type],alpha=0.3)
                ##add fitted lines
                #monomers
                if show_lines["Atlas_fit"]:
                    #Atlas
                    Atlasfit_monomers, = axes[i_ax].semilogx(fit_dic["diam_full"],fit_dic[prop + "_Nmono_1"],marker='None',c="b",linestyle="--",label="Atlas-type cl. ice "  + particle_type,linewidth=2*linewidth)
                if show_lines["powerlaw_fits"]:
                    #power-law
                    powerlaw_monomers, = axes[i_ax].semilogx(fit_dic["diam_full"], fit_dic[prop + "_Nmono_1_powerlaw"],marker='None',c="b",linestyle="-.",label="powerlaw cl. ice "  + particle_type,linewidth=2*linewidth)
                
                #all aggregates
                if show_lines["Atlas_fit"]:
                    #Atlas
                    Atlasfit_Nmonoallagg, = axes[i_ax].semilogx(fit_dic["diam_full"],fit_dic[prop + "_Nmono_allagg"],marker='None',c=["g","g","g","g","g","black","black","r"][i_particle_type],linestyle="--",label="Atlas-type " + corr_cat  +" " + particle_type,linewidth=2*linewidth)
                if show_lines["powerlaw_fits"]:
                    #power-law
                    powerlaw_handle, = axes[i_ax].semilogx(fit_dic["diam_full"], fit_dic[prop + "_Nmono_allagg" + "_powerlaw"],marker='None',c=["g","g","g","g","g","black","black","r"][i_particle_type],linestyle="-.",label="powerlaw " + corr_cat  + " " + particle_type,linewidth=2*linewidth)
                
                #plot v-D quantiles (here just 50) from Hyytiala data
                #axes[i_ax].plot(vD_hyyt_dic["Dmax"]*1e-3,vD_hyyt_dic["v"],linestyle="-",color="darkorange")
            elif forw_mod  and (not no_agg_in_plot): 
                
                color_forw_mod      =["g" ,"g","g","g" ,"black","black"] #["plate","dendrite","column","needle","rosette","mixcoldend1","mixcolumndend"]
                linestyle_forw_mod= linestyles#["-","--",":","-.",":","-",        "-"]
                marker_forw_mod     = markers#["", ""  ,"" ,""  ,"o","o",      "x"]

                
                #color_forw_mod=    ["g","g" ,"g","g" ,"g" ,"black","black"] #["plate","dendrite","column","needle","rosette","mixcoldend1","mixcolumndend"]
                #linestyle_forw_mod=["-","--",":","-.","--","-."     ,"-"]
                #marker_forw_mod=   ["" ,""  ,"" ,""  ,"d" ,"o"      ,"d"]
                axes[i_ax].semilogx(diam_percentile_center,fit_dic[prop + "_" + str(50) + "perc_Nmono_allagg"],color=color_forw_mod[i_particle_type],linewidth=linewidth,linestyle=linestyle_forw_mod[i_particle_type],marker=marker_forw_mod[i_particle_type],markersize=markersize,fillstyle='none')
                #add a line for the Hyytiala data in case we analyze the sideprojection data
                #add percentile range
                if called_by_main:
                    axes[i_ax].fill_between(diam_percentile[:-1],fit_dic[prop + "_" + str(25) + "perc_Nmono_allagg"],fit_dic[prop + "_" + str(75) + "perc_Nmono_allagg"],color=color_forw_mod[i_particle_type],linestyle=linestyle_forw_mod[i_particle_type],alpha=0.1)
                #plot v-D quantiles from Hyytiala data
                #axes[i_ax].plot(vD_hyyt_dic["Dmax"]*1e-3,vD_hyyt_dic["v"],linestyle="-",color="darkorange")
                #axes[i_ax].fill_between(vD_hyyt_dic["Dmax"]*1e-3,vD_hyyt_dic["v_25perc"],vD_hyyt_dic["v_75perc"],color="darkorange",alpha=0.1)
                #plot v-D quantiles from CARE data
                #axes[i_ax].plot(vD_CARE_dic["Dmax"]*1e-3,vD_CARE_dic["v_larger0"],linestyle="-",color="brown") #,label="PIP CARE(v>0.0m/s)")
                #axes[i_ax].fill_between(vD_CARE_dic["Dmax"]*1e-3,vD_CARE_dic["v_25perc_larger0"],vD_CARE_dic["v_75perc_larger0"],color="brown",alpha=0.05)
                #PIP (u>4m/s) b_mass>2.1
                axes[i_ax].plot(data_dic["vD_CAREwindfilt_dic"]["Dmax"]*1e-3,data_dic["vD_CAREwindfilt_dic"]["v_b2p2"],linestyle="-",color="brown",label="PIP CARE")
                axes[i_ax].fill_between(data_dic["vD_CAREwindfilt_dic"]["Dmax"]*1e-3,data_dic["vD_CAREwindfilt_dic"]["v_25perc_b2p2"],data_dic["vD_CAREwindfilt_dic"]["v_75perc_b2p2"],color="brown",alpha=0.02)
                #add LH74 data and M96 data
                #####
                #Locatelli & Hobbs 
                #####
                #data_dic= dict() #same structure as in compare_geom
                low_diam_log= -4; high_diam_log=np.log10(4e-2)
                diam = np.logspace(low_diam_log,high_diam_log,1000) #set up array of diameters
                data_dic["diam"]=diam
                data_dic["LH74"] = dict()
                data_dic["LH74"]["particle_type"] = ["dendrite","mixed","sideplanes"]
                data_dic["LH74"]["particle_type_label"] = ["dendrite","mixed","side planes"]
                data_dic["LH74"]["mass_coeff_dendrite"] = [0.073*10.**(3*1.4-6.),1.4,2e-3,10e-3]
                data_dic["LH74"]["vterm_coeff_dendrite"] = [0.8*10.**(3*0.16),0.16,2e-3,10e-3]
            
                data_dic["LH74"]["mass_coeff_mixed"] = [0.037*10.**(3*1.9-6.),1.9,1e-3,3e-3]
                data_dic["LH74"]["vterm_coeff_mixed"] = [0.69*10.**(3*0.41),0.41,0.2e-3,3e-3]
            
                data_dic["LH74"]["mass_coeff_sideplanes"] = [0.04*10.**(3*1.4-6.),1.4,0.5e-3,4.0e-3]
                data_dic["LH74"]["vterm_coeff_sideplanes"] = [0.82*10.**(3*0.12),0.12,0.5e-3,4.0e-3]
            
                #calculate mass and vterm array of geometric properties
                for key in data_dic["LH74"]["particle_type"]:
                    for prop in ["mass","vterm"]:
                        data_dic["LH74"][prop + '_' + key]= data_dic["LH74"][prop + "_coeff_" + key][0]*diam**data_dic["LH74"][prop + "_coeff_" + key][1]
                for i_particle_type,particle_type in enumerate(data_dic["LH74"]["particle_type"]):
                    diam_flagged = np.where(((data_dic["diam"]>(data_dic["LH74"]["vterm_coeff_" + particle_type][2])) & (data_dic["diam"]<(data_dic["LH74"]["vterm_coeff_" + particle_type][3]))),data_dic["diam"],np.nan)
                    axes[i_ax].plot(diam_flagged,data_dic["LH74"]["vterm_" + particle_type],linestyle=["--","-.",":","-",":"][i_particle_type],color="orange",label="LH74 " + particle_type)
                #'''
                #######
                #Mitchell 1996 #parameterization for m-D and A-D
                #######
                data_dic["M96"] = dict()
                data_dic["M96"]["particle_type"] = ["mixS3","sideplane","polycrys"]
                data_dic["M96"]["mass_coeff_mixS3"] = [0.0028*10.**(2.*2.1-3.),2.1,800e-6,4500e-6]
                data_dic["M96"]["area_coeff_mixS3"] = [0.2285*10.**(2.*1.88-4.),1.88] #area approximated as assemblage of planar polycrstals

                data_dic["M96"]["mass_coeff_sideplane"] = [0.0033*10.**(2.*2.2-3.),2.2,600e-6,4100e-6] 
                data_dic["M96"]["area_coeff_sideplane"] = [0.2285*10.**(2.*1.88-4.),1.88] #area approximated as assemblage of planar polycrstals

                data_dic["M96"]["mass_coeff_polycrys"] = [0.00739*10.**(2.*2.45-3.),2.45,20e-6,450e-6]
                data_dic["M96"]["area_coeff_polycrys"] = [0.2285*10.**(2.*1.88-4.),1.88]

                #calculate array of geometric properties
                for key in data_dic["M96"]["particle_type"]:
                    for prop in ["mass","area"]:
                        data_dic["M96"][prop + '_' + key]= data_dic["M96"][prop + "_coeff_" + key][0]*diam**data_dic["M96"][prop + "_coeff_" + key][1]
                                
                #calculate v-D
                for key in data_dic["M96"]["particle_type"]:
                    data_dic["M96"]["vterm_" + key] = __fallspeed_relations.calc_vterm("bohm",data_dic["M96"]["mass_" + key],data_dic["diam"],data_dic["M96"]["area_" + key])
                #data_dic["M96"]
                for i_particle_type,particle_type in enumerate(data_dic["M96"]["particle_type"]):
                    diam_flagged = np.where(((data_dic["diam"]>(data_dic["M96"]["mass_coeff_" + particle_type][2])) & (data_dic["diam"]<(data_dic["M96"]["mass_coeff_" + particle_type][3]))),data_dic["diam"],np.nan)
                    axes[i_ax].plot(diam_flagged,data_dic["M96"]["vterm_" + particle_type],linestyle=["-.",":","-",":",":"][i_particle_type],color="r",label="__None") #,label="M96" + particle_type + "_bohm")
                #####
                #Kajikawa 1973 individual crystals 
                #####
                data_dic["K72"] = dict()

                data_dic["K72"]["particle_type"] = ["plate","dendrite"] #,"thickplate"]
                data_dic["K72"]["particle_type_label"] = ["plate","dendrite"] #,"thick plate"]

                data_dic["K72"]["diam"]          = np.arange(0.1,1.9,0.05) #np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]) #[mm]
                data_dic["K72"]["diam"] = data_dic["K72"]["diam"]/1000 #[m]
                data_dic["K72"]["vterm_thickplate"] = np.array([10,15, 20,25, 31, 37 , 44,50, 57,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]) #[cm/s]
                data_dic["K72"]["vterm_thickplate"] = np.array([10,15, 20,25, 31, 37 , 44,50, 57,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]) #[cm/s]
                #data_dic["K72"]["vterm_plate"]      = np.array([np.nan,8,  10,13, 15, 17.5, 20, 22.5, 25,28,30.5,33,35,37.5,40,43,46,48.5,51,53.5,56,58,60,62.5,65,67.5,70,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]) #[cm/s]
                data_dic["K72"]["diam_platelin"]      = np.array([0.15,1.4])/1000.
                data_dic["K72"]["vterm_plate"]      = np.array([8.,70.])
                data_dic["K72"]["vterm_dendrite"]      = np.array([np.nan,np.nan,  np.nan,np.nan, np.nan,np.nan,np.nan, 11.5 ,12.5,13,14,14.9,15.5,16.2,17,17.5,18,18.5,19,19.6,20,20.6,21,21.3,21.5,21.8,22,22.3,22.5,22.8,23,23.3,23.5,23.8,24,24.2]) #[cm/s]

    
                for prop in data_dic["K72"]["particle_type"]:
                    data_dic["K72"]["vterm_"+ prop] = data_dic["K72"]["vterm_"+ prop]/100. #[m/s]
                #data_dic["K72"]
                #for i_particle_type,particle_type in enumerate(["dendrite"]):
                #    axes[i_ax].semilogx(data_dic["K72"]["diam"],data_dic["K72"]["vterm_" + particle_type],linestyle=["--",":",":",":"][i_particle_type],color="b",label="K72 " + data_dic["K72"]["particle_type_label"][i_particle_type])
                axes[i_ax].semilogx(data_dic["K72"]["diam"],data_dic["K72"]["vterm_dendrite"],linestyle="--",color="b",label="K72 dendrite")
                axes[i_ax].semilogx(data_dic["K72"]["diam_platelin"],data_dic["K72"]["vterm_plate"],linestyle="-",color="b",label="K72 plate")

            #make labels
            axes[i_ax].set_xlabel("$D_{max}$ [m]")
            axes[i_ax].set_ylabel("$v_{term}$ [m/s]" ) #TODO: plot also the untis of these properties

            #change the axis
            axes[i_ax].set_xlim([10**low_diam_log_detailed,10**high_diam_log_detailed_plotrange]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
            axes[i_ax].set_ylim([0,2.5]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
            axes[i_ax].grid(b=True,which="both",axis="both")

    for i_ax,ax_now in enumerate(axes):
        for i_particle_type,particle_type in enumerate(particle_types):
            if particle_type=="mixcolumndend":
                particle_type_label="Mix2"
            elif particle_type=="mixcoldend1":
                particle_type_label="Mix1"
            else:
                particle_type_label=particle_type

            if forw_mod:
                ax_now.plot(np.nan,np.nan,color=color_forw_mod[i_particle_type],linestyle=linestyle_forw_mod[i_particle_type],marker=marker_forw_mod[i_particle_type],linewidth=linewidth,label=particle_type_label,markevery=50,markersize=markersize,fillstyle='none')
            else:
                ax_now.plot(np.nan,np.nan,color=["g","g","g","g","black","black","g","r"][i_particle_type],linestyle=linestyles[i_particle_type],marker=markers[i_particle_type],linewidth=linewidth,label=particle_type_label,markevery=50, markersize=markersize,fillstyle='none')
        if not (("_sideproj" in ''.join(particle_types)) or no_monomers_in_plot):
            ax_now.plot(np.nan,np.nan,color="b",label="Nmono=1")
                
        #if not forw_mod:
        #    ax_now.plot(np.nan,np.nan,color="g",label="Nmono>1")
        if forw_mod:
            #ax_now.plot(np.nan,np.nan,linestyle="-",color="darkorange",label="Hyytiala(v>0.5m/s)")
            pass #ax_now.plot(np.nan,np.nan,linestyle="-",color="red",label="PIP CARE(v>0.0m/s)")
        #show legend
        if called_by_main: 
            axes[i_ax].legend(ncol=1,loc='upper left') #, bbox_to_anchor=(1.0, 0.5)) #ncol=1,loc='center left', bbox_to_anchor=(1.0, 0.5))


    return axes


if __name__ == '__main__':

    ##general settings for the plotting
    number_of_plots = 3

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
    axes = read_and_plot(axes,hydro_model="all",forw_mod=forw_mod,tumbling=tumbling,called_by_main=True)

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
