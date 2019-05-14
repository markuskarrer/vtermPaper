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
from matplotlib.ticker import FixedFormatter
from netCDF4 import Dataset
from sys import argv
import re

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
#import other self-defined functions
import __postprocess_McSnow
import __postprocess_SB
import __fallspeed_relations
import __tools_for_processing_Jagg
import __plotting_functions
#from IPython.core.debugger import Tracer ; Tracer()()

'''
this code reads in properties from Jussis aggregate model and McSnow output and creates a 2D histogram (or density estimate) of both + the ratio (McSnow/Jagg) which can be used as a weighting function for the fits
'''


def calc_histogram(diam_particle_array,N_mono_particle_array,low_diam_log=-4, high_diam_log=-1,nbins=40,weights="None"):
    '''
    calculate the histogram
    INPUT:  diam_particle_array: array containing the diameter of all particles
            N_mono_particle_array:array containing the monomer number of all particles
            low_diam_log=-4, high_diam_log=-1,nbins=40 define the array of diameter
            weights: possibility to apply weights (e.g. the multiplicity from McSnow)

    OUTPUT: diam_edges,Nmono_edges: edges of the histogram
            H values of the histogram
    '''
    
    #set up array of diameters (for bin edges and KDE-grid) and the fitting
    diam_edges = np.logspace(low_diam_log,high_diam_log,nbins) #set up array of diameters
    Nmono_edges = np.array([2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,1e3,1e4,1e5])
    
    #calculate the histogram
    if isinstance(weights,str):
        H, xedges, yedges = np.histogram2d(diam_particle_array, N_mono_particle_array, bins=(diam_edges, Nmono_edges))
    else:
        H, xedges, yedges = np.histogram2d(diam_particle_array, N_mono_particle_array, bins=(diam_edges, Nmono_edges), weights=weights)
        
    return diam_edges,Nmono_edges,H

def calc_median_prop_in_hist(diam_particle_array,N_mono_particle_array,area_particle_array,mass_particle_array,a,b,c,d,low_diam_log=-4, high_diam_log=-1,nbins=40,Nmono_edges=np.array([2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,1e3,1e4,1e5]),weights="None"):
    '''
    calculate the median of the area and mass in each bin
    INPUT:  diam_particle_array: array containing the diameter of all particles
            N_mono_particle_array:array containing the monomer number of all particles
            area_particle_array: array containing the area of all particles
            mass_particle_array: array containing the mass of all particles
            a,b,c,d: coefficients of the monomer (m=aD**b; A=c*D**b)
            low_diam_log=-4, high_diam_log=-1,nbins=40 define the array of diameter
            weights: possibility to apply weights (e.g. the multiplicity from McSnow)
            Nmono_edges: edges of monomer number
    '''
    
    #set up array of diameters (for bin edges and KDE-grid) and the fitting
    diam_edges = np.logspace(low_diam_log,high_diam_log,nbins) #set up array of diameters
    #Nmono_edges = np.array([2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,1e3,1e4,1e5])

    H_median_mass = np.zeros([Nmono_edges.shape[0],diam_edges.shape[0]]) #median mass in each bin
    H_median_area = np.zeros([Nmono_edges.shape[0],diam_edges.shape[0]]) #median area in each bin
    H_median_mass_ratio = np.zeros([Nmono_edges.shape[0],diam_edges.shape[0]]) #median mass in each bin (divided by monomer mass)
    H_median_area_ratio = np.zeros([Nmono_edges.shape[0],diam_edges.shape[0]]) #median area in each bin (divided by monomer area)
    for i_Nmono,Nmono_lowbound in enumerate(Nmono_edges[:-1]):
        for i_diam,diam_lowbound in enumerate(diam_edges[:-1]):
            diam_center = 10**((np.log10(diam_edges[i_diam])+np.log10(diam_edges[i_diam+1]))/2.)
            #diam_log_center = 10**((np.log10(diam_edges[:-1])+np.log10(diam_edges[1:]))/2)

            print "diam",diam_center
            in_bin = np.logical_and(np.logical_and(np.logical_and(diam_edges[i_diam]<=diam_particle_array,diam_particle_array<diam_edges[i_diam+1]),Nmono_edges[i_Nmono]<=N_mono_particle_array),N_mono_particle_array<Nmono_edges[i_Nmono+1])
            #print diam_edges[i_diam],diam_edges[i_diam+1],Nmono_edges[i_Nmono],Nmono_edges[i_Nmono+1],sum(in_bin),np.nanmean(mass_particle_array[in_bin]),(a*diam_center**b)#; raw_input("wait")

            
            H_median_mass[i_Nmono,i_diam] = np.nanmedian(mass_particle_array[in_bin])
            H_median_area[i_Nmono,i_diam] = np.nanmedian(area_particle_array[in_bin])
            H_median_mass_ratio[i_Nmono,i_diam] = np.nanmedian(mass_particle_array[in_bin])/(a*diam_center**b)
            print np.nanmedian(mass_particle_array[in_bin]),(a*diam_center**b),"rat",H_median_mass_ratio[i_Nmono,i_diam]
            H_median_area_ratio[i_Nmono,i_diam] = np.nanmedian(area_particle_array[in_bin])/(c*diam_center**d)
            #if sum(in_bin)>10:
            #    from IPython.core.debugger import Tracer ; Tracer()()
    return H_median_mass,H_median_area,H_median_mass_ratio,H_median_area_ratio

def N_D_Dmono_from_MC_and_Jagg(particle_type="plate"):
    
    #this function is necessary so that this script can be run by another script
    #print len(argv); raw_input() #TODO: this trick is conflicting with being called as an imported function
    #if len(argv)>1: #in case this script is called with arguments
    #    noplotting = argv[1] #read a bool wheter there should be plotting
    #else:
    noplotting = False



    def plot_histogram(diam_edges,Nmono_edges,H,axes):
        '''
        plot the histogram
        INPUT:  diam_edges,Nmono_edges: edges of the histogram
                H values of the histogram
                axes: the axes of the figure created in this script
                
        '''
        #self-defined discrete colormap
        # define the colormap and the discrete boundaries (defined by the considered monomer numbers) and set up an array of the same colors (for line-plotting)
        cmap = plt.cm.viridis
        cmaplist = [cmap(i) for i in range(cmap.N)]
        #set zero counts to white
        #cmaplist[0] = (0.0,0.0,0.0,0.0)
        #cmap = colors.LinearSegmentedColormap.from_list('mcm',cmaplist, cmap.N)
        cmap.set_under('black',0.2)

        #plot the histogram of Jagg
        im = axes[i_ax].pcolor(range(0,diam_edges.shape[0]), range(0,Nmono_edges.shape[0]), np.transpose(H)/np.sum(H[:]),cmap=cmap,norm=colors.LogNorm())
        
        #add colorbar
        cbar = fig.colorbar(im,ax=axes[i_ax])
        cbar.set_label("normalized counts / 1")
        axes[i_ax].tick_params(axis='x', which='minor', bottom=False)

        ##
        #set xticks according to boundaries
        xtick_indices = range(0,diam_edges.shape[0],13)
        axes[i_ax].set_xticks(xtick_indices)
        new_labels = [ "%.2e" % l for l in diam_edges[xtick_indices]]
        axes[i_ax].set_xticklabels(new_labels)
        #set yticks according to boundaries
        ytick_indices = range(0,Nmono_edges.shape[0])
        axes[i_ax].set_yticks(ytick_indices)
        new_labels = [ "%d" % tick if i_tick%2 else "" for i_tick,tick in enumerate(Nmono_edges[ytick_indices])]
        axes[i_ax].set_yticklabels(new_labels)
        #set ticklabels
        axes[i_ax].set_xlabel("Diameter D / m")
        axes[i_ax].set_ylabel("Monomer number / 1")
        #switch minor ticks on
        #axes[i_ax].tick_params(axis='x', which='minor', bottom=False)
        #axes[i_ax].minorticks_on() #.tick_params(axis='y', which='minor', bottom=False)

        return axes

    def plot_histogram_ratio(diam_edges,Nmono_edges,H1,H2,axes):
        '''
        plot the histogram
        INPUT:  diam_edges,Nmono_edges: edges of the histogram
                H values of the histogram
                axes: the axes of the figure created in this script
                
        '''
        #calculate ratio
        #H = np.divide(H_MC, H_Jagg, out=np.inf*np.ones_like(H_MC), where=H_Jagg!=0)
        
        #weight according to MCsnows distribution of monomer numbers at each diameter bin (but weigh each diameter bin equally)
        H_weight = np.zeros_like(H_MC)
        H_plot = np.zeros_like(H_MC) #here we are allowed to have nans and need them to indicate in gray where we have unused Jagg-data

        for i_diam in range(H_MC.shape[0]): #first dimension is the diameter
            if i_diam==30:   
                pass #from IPython.core.debugger import Tracer ; Tracer()()
            ###calculate the weighting factor
            #we normalize the weighting factor at each diameter to one : therefore we need the overlap of valid data from Jagg and MC
            sum_valid_Jagg = np.nansum(H_Jagg[i_diam,:][H_MC[i_diam,:]>0])
            sum_valid_MC = np.nansum(H_MC[i_diam,:][H_Jagg[i_diam,:]>0])
            sum_product = np.nansum(H_MC[i_diam,:][H_Jagg[i_diam,:]>0]/H_Jagg[i_diam,:][H_Jagg[i_diam,:]>0])
            if np.nansum(H_Jagg[i_diam,:])>0 and np.nansum(H_MC[i_diam,:][H_Jagg[i_diam,:]>0])>0: #check if there is any bin with an overlap of Jagg and MC data
                H_weight[i_diam,:] = np.divide(H_MC[i_diam,:],H_Jagg[i_diam,:], out=np.ones_like(H_MC[i_diam,:])*0, where=H_Jagg[i_diam,:]!=0)/sum_product #sum_valid_Jagg/sum_valid_MC #np.nansum(H_MC[i_diam,:][H_Jagg[i_diam,:]>0])
                H_plot[i_diam,:] = np.divide(H_MC[i_diam,:],H_Jagg[i_diam,:], out=np.ones_like(H_MC[i_diam,:])*np.nan, where=np.logical_and(H_Jagg[i_diam,:]!=0,H_MC[i_diam,:]!=0))/sum_product #*sum_valid_Jagg/sum_valid_MC #np.nansum(H_MC[i_diam,:][H_Jagg[i_diam,:]>0])
                
            else: 
                H_weight[i_diam,:] = 0
                H_plot[i_diam,:] = np.inf
        
        #from IPython.core.debugger import Tracer ; Tracer()()
        #self-defined discrete colormap
        # define the colormap and the discrete boundaries (defined by the considered monomer numbers) and set up an array of the same colors (for line-plotting)
        cmap = plt.cm.brg #coolwarm
        
        #cmaplist = [cmap(i) for i in range(cmap.N)]
        #set zero counts to white
        #cmaplist[0] = (0.8,0.8,0.8,1.0)
        #cmap = colors.LinearSegmentedColormap.from_list('mcm',cmaplist, cmap.N)
        #cmap.set_bad('black',0.2) #this is masking the zero-values in gray
        #plot the histogram of Jagg
        im = axes[i_ax].pcolor(range(0,diam_edges.shape[0]), range(0,Nmono_edges.shape[0]), np.transpose(H_plot),cmap=cmap) #,norm=colors.LogNorm())

        #plot the nans (representing missing desired values in the aggregate model)
        H_masked = np.nan*np.ones_like(H_plot) #np.ma.masked_where(H_plot,H_Jagg>0) #np.logical_and(H_Jagg>0,H_MC>0))
        H_masked[np.logical_and(H_MC>0,H_Jagg==0)] = 1.
        axes[i_ax].pcolor(range(0,diam_edges.shape[0]), range(0,Nmono_edges.shape[0]), np.transpose(H_masked),hatch='\\\\\\\\\\\\\\\\\\\\',alpha=0)

        #plot the unused data (representing missing desired values in the aggregate model)
        H_masked = np.nan*np.ones_like(H_plot) #np.ma.masked_where(H_plot,H_Jagg>0) #np.logical_and(H_Jagg>0,H_MC>0))
        H_masked[np.logical_and(H_MC==0,H_Jagg>0)] = 1.
        axes[i_ax].pcolor(range(0,diam_edges.shape[0]), range(0,Nmono_edges.shape[0]), np.transpose(H_masked),cmap=plt.cm.Greys_r,alpha=0.2)
        
        #add colorbar
        cbar = fig.colorbar(im,ax=axes[i_ax])
        cbar.set_label("ratio RP/Jagg / 1")
        ##
        #set xticks according to boundaries
        xtick_indices = range(0,diam_edges.shape[0],13)
        axes[i_ax].set_xticks(xtick_indices)
        new_labels = [ "%.2e" % l for l in diam_edges[xtick_indices]]
        axes[i_ax].set_xticklabels(new_labels)
        #set yticks according to boundaries
        ytick_indices = range(0,Nmono_edges.shape[0])
        axes[i_ax].set_yticks(ytick_indices)
        new_labels = [ "%d" % tick if i_tick%2 else "" for i_tick,tick in enumerate(Nmono_edges[ytick_indices])]
        axes[i_ax].set_yticklabels(new_labels)
        #set ticklabels
        axes[i_ax].set_xlabel("Diameter D / m")
        axes[i_ax].set_ylabel("Monomer number / 1")
        #from IPython.core.debugger import Tracer ; Tracer()()
        return axes,H_weight

    ###############
    #aggregate model
    #################

    #define where the txt files with properties of the aggregates are stored
    prop_file_folder = "/data/optimice/Jussis_aggregates/tumbling_asratio/"

    grid_res = 10e-6
    if grid_res==40e-6:
        print "processing particle with resolution: ",grid_res
        sensrun_folder = '/'
    else:
        print "processing particle with resolution: ",grid_res
        sensrun_folder = 'res_'+ str(int(grid_res*1e6)) + 'mum/'# if '/' this will be in /data/optimice/Jussis_aggregates

        #sensrun_folder = '/' #res_10mum/'# if '/' this will be in /data/optimice/Jussis_aggregates


    ###
    #set up the figure
    ###
    number_of_plots = 2 #28

    #optimize the appearance of the plot (figure size, fonts)
    [fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots)


    #read the properties of the aggregates into the particle_dic dictionary
    #particle_type = "column" #,"dendrite","plate","column","rosette"] # ,"bullet"rosette
    print "########################"
    print "#####" + particle_type + "########"
    print "########################"
    particle_dic = __tools_for_processing_Jagg.init_particle_dict() #initialize dictionary which contains information about the type of the pristine crystals (dentrites, needles, plates, ...) the underlying distribution of the   
    for i_file,filename in enumerate(glob.glob(prop_file_folder + sensrun_folder +  particle_type + '*properties.txt')): #loop over all property files (with different mean sizes)
        #read size parameter from filename in order to disregard some size parameter
        m=re.search(prop_file_folder + sensrun_folder + particle_type + '_(.+?)_properties.txt', filename) #search for patter
        size_param_now = float(m.group(1))
        if size_param_now>1000: #ignore all sizeparameter above ... or below ...
            continue
        print "reading: " + filename
        #get m-D from the assumptions in the aggregate model (this is used at various places)
        a,b,c,d = __tools_for_processing_Jagg.calc_mD_AD_coeffs(particle_type)

        with open(filename,"rb") as txtfile: #TODO: this is just one file! #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in
            prop_reader = csv.reader(filter(lambda row: row[0]!='#',txtfile), delimiter=' ', quoting=csv.QUOTE_NONNUMERIC, lineterminator=os.linesep) #row[0]!='#': skip the header; quoting avoids '' for formatted string; lineterminator avoids problems with system dependend lineending format https://unix.stackexchange.com/questions/309154/strings-are-missing-after-concatenating-two-or-more-variable-string-in-bash?answertab=active#tab-top

            #define the monomer numbers to read in
            take_all_Nmono_bool = True
            if take_all_Nmono_bool:
                N_mono_list = np.append(np.append(np.array(range(1,10)),np.array(range(10,100,1))),np.array(range(100,1001,100)))
            else: #select certain monomer numbers
                N_mono_list = np.array([1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100])


            for i_row,row_content in enumerate(prop_reader): #TODO: the row is not any more equal to the monomer number #read the row with N_mono_now (the currently considered monomer number)
                if row_content[0] in N_mono_list: #i_row==0 or i_row==4 or i_row==9:
                    particle_dic["particle_type"] = np.append(particle_dic["particle_type"],particle_type)
                    particle_dic["N_monomer"] = np.append(particle_dic["N_monomer"],row_content[0]) #python indices start with 0
                    particle_dic["mass"] = np.append(particle_dic["mass"],row_content[1])
                    particle_dic["area"] = np.append(particle_dic["area"],row_content[2])
                    particle_dic["diam"] = np.append(particle_dic["diam"],row_content[3])
                    particle_dic["sizeparam_index"] = np.append(particle_dic["sizeparam_index"],i_file) #
    print particle_dic #overview of the dictionary; might be helpful



    ###
    #subplot 1: histogram of Jagg
    ###
    i_ax=0

    #define the overall diameter array
    low_diam_log=-4; high_diam_log=-1;nbins=40

    #calculate the histogram
    diam_edges,Nmono_edges,H_Jagg = calc_histogram(particle_dic["diam"],particle_dic["N_monomer"],low_diam_log=low_diam_log, high_diam_log=high_diam_log,nbins=nbins)
    if not noplotting:
        #plot the histogram
        axes = plot_histogram(diam_edges,Nmono_edges,H_Jagg,axes)
        #add title
        axes[i_ax].set_title("Aggregate model")
        #add number of analyzed aggregates
        axes[i_ax].text(0.1, 1.05, 'N_agg= ' +  str(particle_dic["diam"].shape[0]), horizontalalignment='center', verticalalignment='center', transform=axes[i_ax].transAxes)

    ########
    #MCSNOW
    ########
    directory = "/home/mkarrer/Dokumente/McSnow/MCSNOW/experiments/"
    #experiment = "1d__param_xi1000_nz250_lwc0_ncl0_ssat5_dtc5_nrp500_rm10_rt0_vt3_at2_stick1_dt1_meltt0_multt0_h0-0_ba500"
    experiment = ["1d__param_xi1000_nz250_lwc0_ncl0_ssat5_dtc5_nrp50_rm10_rt0_vt3_at2_stick1_dt1_meltt0_multt0_h0-0_ba500", #low qi #ssat=0.5%
                  "1d__param_xi1000_nz250_lwc0_ncl0_ssat5_dtc5_nrp500_rm10_rt0_vt3_at2_stick1_dt1_meltt0_multt0_h0-0_ba500",#medium qi
                  #"1d__param_xi10000_nz250_lwc0_ncl0_ssat5_dtc5_nrp5000_rm10_rt0_vt3_at2_stick1_dt1_meltt0_multt0_h0-0_ba500", #high qi
                  "1d__param_xi1000_nz250_lwc0_ncl0_ssat10_dtc5_nrp50_rm10_rt0_vt3_at2_stick1_dt1_meltt0_multt0_h0-0_ba500", #low qi #ssat=1.0%
                  "1d__param_xi1000_nz250_lwc0_ncl0_ssat10_dtc5_nrp500_rm10_rt0_vt3_at2_stick1_dt1_meltt0_multt0_h0-0_ba500",#medium qi
                  "1d__param_xi10000_nz250_lwc0_ncl0_ssat10_dtc5_nrp5000_rm10_rt0_vt3_at2_stick1_dt1_meltt0_multt0_h0-0_ba500", #high qi
                  "1d__param_xi1000_nz250_lwc0_ncl0_ssat20_dtc5_nrp50_rm10_rt0_vt3_at2_stick1_dt1_meltt0_multt0_h0-0_ba500", #, #low qi #ssat=2.0%
                  "1d__param_xi1000_nz250_lwc0_ncl0_ssat20_dtc5_nrp500_rm10_rt0_vt3_at2_stick1_dt1_meltt0_multt0_h0-0_ba500",#medium qi
                  "1d__param_xi10000_nz250_lwc0_ncl0_ssat20_dtc5_nrp5000_rm10_rt0_vt3_at2_stick1_dt1_meltt0_multt0_h0-0_ba500", #high qi
                  ]
    #create dictionary for all superparticle properties
    SP = dict()

    for exp in experiment:
        print "reading:",exp
        #load file
        SP_file = Dataset(directory + exp + '/mass2fr_' + str(600).zfill(4) + 'min_avtstep_' + str(5) + '.ncdf',mode='r')

        #if necessary change name of variables
        varlist = SP_file.variables
        for var in varlist: #copy variables to SP dictionary
            if var == "xi": #divide additionally by the total of "real particles" in order to not overweight runs with high number concentration
                if "xi" in SP.keys():
                    SP["xi"] = np.append(SP["xi"],SP_file.variables["xi"]/np.sum(SP_file.variables["xi"]))
                else: #this should be the case for the first exp in experiment
                    SP["xi"] = SP_file.variables["xi"]/np.sum(SP_file.variables["xi"])
                #print np.sum(SP_file.variables["mm"]),SP["mm"][0]; raw_input("wait")
            elif var in SP.keys():
                SP[var] = np.append(SP[var],SP_file.variables[var])
            else: #this should be the case for the first exp in experiment
                SP[var] = SP_file.variables[var]
   
    #get maximum height for plotting
    model_top = np.nanmax(SP['height'])

    
    ###
    #subplot 2: histogram of McSnow
    ###
    i_ax+=1
    
    #calculate the histogram
    #diam_edges,Nmono_edges,H_MC = calc_histogram(SP["diam"][SP["height"]<1000],SP["mm"][SP["height"]<1000],low_diam_log=low_diam_log, high_diam_log=high_diam_log,nbins=nbins, weights=SP["xi"][SP["height"]<1000]) #ATTENTION "SP["mm"] is the monomer number per superparticle
    diam_edges,Nmono_edges,H_MC = calc_histogram(SP["diam"],SP["mm"],low_diam_log=low_diam_log, high_diam_log=high_diam_log,nbins=nbins, weights=SP["xi"]) #ATTENTION "SP["mm"] is the monomer number per superparticle
    if not noplotting:
        #plot the histogram
        #from IPython.core.debugger import Tracer ; Tracer()()
        axes = plot_histogram(diam_edges,Nmono_edges,H_MC,axes)
        #add title
        axes[i_ax].set_title("McSnow")
        #add number of SP on top
        axes[i_ax].text(0.1, 1.05, 'N_SP= ' +  str(SP["diam"].shape[0]), horizontalalignment='center', verticalalignment='center', transform=axes[i_ax].transAxes)
        '''
        ###
        #subplot 3: ratio of Jagg/McSnow counts (this can be used as a weighting function for the fit)
        ###
        i_ax+=1

        #plot the histogram
        axes,H = plot_histogram_ratio(diam_edges,Nmono_edges,H_MC,H_Jagg,axes) #diam_edges,Nmono_edges,H1,H2,axes
        axes[i_ax].set_title("weighting factor for fits: McSnow/Jagg")
        
        ###
        #subplot 3: median mass divided by momomer mass in bins
        ###
        i_ax+=1
        axes[i_ax].set_xscale("log")
        axes[i_ax].set_yscale("log")

        #from IPython.core.debugger import Tracer ; Tracer()()
        Nmono_edges = np.array([2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,1e3,1e4,1e5])
        H_median_mass,H_median_area,H_median_mass_ratio,H_median_area_ratio = calc_median_prop_in_hist(particle_dic["diam"],particle_dic["N_monomer"],particle_dic["mass"],particle_dic["area"],a,b,c,d,low_diam_log=-4, high_diam_log=-1,nbins=40,Nmono_edges=Nmono_edges)
        #from IPython.core.debugger import Tracer ; Tracer()()
        pcol_mass = axes[i_ax].pcolor(diam_edges,Nmono_edges,H_median_mass_ratio,norm=colors.LogNorm())
        #add colorbar
        cbar = fig.colorbar(pcol_mass,ax=axes[i_ax])
        cbar.set_label("median mass ratio / 1")
        axes[i_ax].tick_params(axis='x', which='minor', bottom=False)
        #set ticklabels
        axes[i_ax].set_xlabel("Diameter D / m")
        axes[i_ax].set_ylabel("Monomer number / 1")
        ###
        #subplot 4: median area divided by momomer area in bins
        ###
        i_ax+=1
        axes[i_ax].set_xscale("log")
        axes[i_ax].set_yscale("log")
        
        pcol_area = axes[i_ax].pcolor(diam_edges,Nmono_edges,H_median_area_ratio,norm=colors.LogNorm())
        
        #add colorbar
        cbar = fig.colorbar(pcol_area,ax=axes[i_ax])
        cbar.set_label("median proj. area / m2")
        axes[i_ax].tick_params(axis='x', which='minor', bottom=False)
        #set ticklabels
        axes[i_ax].set_xlabel("Diameter D / m")
        axes[i_ax].set_ylabel("Monomer number / 1")
        '''
        #save the plot (and open it)
        plt.tight_layout()
        dir_save = '/home/mkarrer/Dokumente/plots/Jagg/'
        if not os.path.exists(dir_save): #create direktory if it does not exists
            os.makedirs(dir_save)
        out_filestring = "McSnow_Jagg_comparison_" + particle_type + "_gridres" + str(grid_res)
        plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
        plt.savefig(dir_save + out_filestring + '.png', dpi=100)
        
        print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
        subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
        plt.clf()
        plt.close()
        
        H=np.zeros_like(H_MC)
    return diam_edges,Nmono_edges,H_MC,H_Jagg,H
    
if __name__ == "__main__":
   	# stuff only to run when not called via 'import' here
	#execute this function
	N_D_Dmono_from_MC_and_Jagg()
