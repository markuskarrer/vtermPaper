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
from matplotlib import rc

'''
this code reads in properties from Jussis aggregate model, calculate + plots the following fallspeed
and fits different functions to each property (m-D,A-D,v-D for all fallspeed models)
'''



def get_model_params():
    ####
    #START: get model parameters
    ####
    fit_dic = dict() #initialize fit dictionary


    #McSnow
    fit_dic["mDADvD_dict_MC"] = __setup_mDAD_frommodels.get_model_mDADs(model="MC")
    #SB
    fit_dic["mDADvD_dict_SBcloudice"] = __setup_mDAD_frommodels.get_model_mDADs(model="SBcloudice")
    fit_dic["mDADvD_dict_SBsnow"] = __setup_mDAD_frommodels.get_model_mDADs(model="SBsnow")
    fit_dic["mDADvD_dict_SBcloudiceAtlas"] = __setup_mDAD_frommodels.get_model_mDADs(model="SBcloudice_Atlas")
    fit_dic["mDADvD_dict_SBsnowAtlas"] = __setup_mDAD_frommodels.get_model_mDADs(model="SBsnow_Atlas")
    #P3
    fit_dic["mDADvD_dict_P3"] = __setup_mDAD_frommodels.get_model_mDADs(model="P3")
    #Morrison scheme (from WRF 4.0)
    fit_dic["morr2mom_cloudice"]  = __setup_mDAD_frommodels.get_model_mDADs(model="morr2mom_cloudice")
    fit_dic["morr2mom_snow"]  = __setup_mDAD_frommodels.get_model_mDADs(model="morr2mom_snow")
    #Goddard scheme (from WRF 4.0)
    fit_dic["GSFCsnow"]  = __setup_mDAD_frommodels.get_model_mDADs(model="GSFCsnow")
    fit_dic["thompson_cloudice"]   = __setup_mDAD_frommodels.get_model_mDADs(model="thompson_cloudice")
    fit_dic["thompson_snow"]   = __setup_mDAD_frommodels.get_model_mDADs(model="thompson_snow")

    #calculate the arrays with masses and areas corresponding to the common fit_dic["diam"] and based on the (piecewise) power-law fits
    #set up array of diameters for bin edges, the KDE-grid and the fitting
    low_diam_log=-4; high_diam_log=-1.3979400086720375 #np.log10(3e-2)=-1.5228787452803376 #np.log10(4e-2)=-1.3979400086720375
    diam = np.logspace(low_diam_log,high_diam_log,1000) #set up array of diameters (for bin edges and KDE-grid)
    fit_dic["diam"] = diam
    #get the m,A,v-array corresponding to diam and save them in separate dictionaries
    for dict_now in (fit_dic["mDADvD_dict_MC"],fit_dic["mDADvD_dict_P3"],fit_dic["mDADvD_dict_SBcloudice"],fit_dic["mDADvD_dict_SBsnow"],fit_dic["mDADvD_dict_SBcloudiceAtlas"],fit_dic["mDADvD_dict_SBsnowAtlas"],fit_dic["morr2mom_cloudice"],fit_dic["morr2mom_snow"],fit_dic["GSFCsnow"],fit_dic["thompson_cloudice"],fit_dic["thompson_snow"]):
        dict_now = __setup_mDAD_frommodels.calc_area_mass_vterm_arrays(fit_dic["diam"],dict_now)
        

    ####
    #END: get model parameters
    #for comparing with models 
    ####
    return fit_dic


def plot_model_vterms(ax,fit_dic,scale='xlogylin',show_lines=None,linewidth=1.2):
    #plot the fall speed in the current plot
    
    if show_lines==None:
        show_lines = dict()
        #define which lines to show (ATTENTION: default, can be overwritten by optional arguments below):
        show_lines["SB_mD"] = True #show the m-D relation of the SB scheme
        #from previous model settings
        show_lines["SB_powerlaw"] = True #show also the old fit of Axels power-law
        show_lines["SB_Atlas"] = False #show also the old fit of Axels Atlas-type
        #from other models
        show_lines["P3"] = True
        show_lines["MC"] = False
        show_lines["Morr2mom"] = True
        show_lines["GSFC"] = False #single Moment Godard scheme
        show_lines["THOM"] = False #two-moment Thompson scheme
    
    
    
    if show_lines["SB_powerlaw"]:
        SB_ci_handle,   = ax.plot(fit_dic["diam"],fit_dic["mDADvD_dict_SBcloudice"]["v(D_array)"],color='blue', label="SB cloud ice",linestyle='-',linewidth=linewidth)
        SB_snow_handle, = ax.plot(fit_dic["diam"],fit_dic["mDADvD_dict_SBsnow"]["v(D_array)"]    ,color='green',label="SB snow",     linestyle='-',linewidth=linewidth)
    if show_lines["MC"]:
        MC_handle, = ax.plot(fit_dic["diam"],mDADvD_dict_MC["v_" + velocity_model + "(D_array)"],color='red',label="McSnow(" + velocity_model + ")",linestyle='-',linewidth=linewidth)
    if show_lines["P3"]:
        
        D_min_P3=20e-6     #see do jj = 1,1000 ...d1 = real(jj)*20.*1.e-6 - 10.e-6 in create_p3_lookupTable_1.f90
        D_max_P3=20.*1000.*1e-6       #see do jj = 1,1000 ...d1 = real(jj)*20.*1.e-6 - 10.e-6 in create_p3_lookupTable_1.f90
        fit_dic["mDADvD_dict_P3"]["v_mitch_heym(D_array)"] = np.where(fit_dic["diam"]<D_min_P3,np.nan,fit_dic["mDADvD_dict_P3"]["v_mitch_heym(D_array)"])
        fit_dic["mDADvD_dict_P3"]["v_mitch_heym(D_array)"] = np.where(fit_dic["diam"]>D_max_P3,np.nan,fit_dic["mDADvD_dict_P3"]["v_mitch_heym(D_array)"])
        
        P3_handle, = ax.plot(fit_dic["diam"],fit_dic["mDADvD_dict_P3"]["v_mitch_heym(D_array)"],color='orange',label="P3 unrimed",linestyle='-',linewidth=linewidth) #mitch heym is hardcoded because there are no other options yet
    if show_lines["SB_Atlas"]:
        Axel_iceAtlas_handle, = ax.plot(fit_dic["mDADvD_dict_SBcloudiceAtlas"]["D_max_from_moltenD"],fit_dic["mDADvD_dict_SBcloudiceAtlas"]["v(D_array)"],color='blue',label="SB cloud ice Atlas",linestyle=':')
        Axel_snowAtlas_handle, = ax.plot(fit_dic["mDADvD_dict_SBsnowAtlas"]["D_max_from_moltenD"],fit_dic["mDADvD_dict_SBsnowAtlas"]["v(D_array)"],color='green',label="SB snow Atlas",linestyle=':')
    if show_lines["Morr2mom"]:
        Morr_ci_handle,   = ax.plot(fit_dic["diam"],fit_dic["morr2mom_cloudice"]["v(D_array)"],color='blue', label="Morr cloud ice",linestyle='--',linewidth=linewidth)
        Morr_snow_handle, = ax.plot(fit_dic["diam"],fit_dic["morr2mom_snow"]["v(D_array)"]    ,color='green',label="Morr snow",     linestyle='--',linewidth=linewidth)
    if show_lines["GSFC"]:
        GSFC_snow_handle, = ax.plot(fit_dic["diam"],fit_dic["GSFCsnow"]["v(D_array)"]    ,color='green',label="GSFC snow",     linestyle='-.')
        #TODO: what about GSFC snow
    if show_lines["THOM"]:
        D_max_cloud_ice_thompson=5*200e-6 #see xDx(nbi+1) = 5.0d0*D0s in module_mp_thompson.F
        D_min_snow_thompson=200e-6  #see xDx(1) = D0s*1.0d0 in module_mp_thompson.F
        D_max_snow_thompson=0.02 #see xDx(nbi+1) = 5.0d0*D0s in module_mp_thompson.F

        GSFC_clice_handle, = ax.plot(fit_dic["diam"][fit_dic["diam"]<D_max_cloud_ice_thompson],fit_dic["thompson_cloudice"]["v(D_array)"][fit_dic["diam"]<D_max_cloud_ice_thompson]    ,color='blue',label="THOM cloud ice",     linestyle='-.')
        thompson_snow["v(D_array)"] = np.where(fit_dic["diam"]<D_min_snow_thompson,np.nan,fit_dic["thompson_snow"]["v(D_array)"])
        thompson_snow["v(D_array)"] = np.where(fit_dic["diam"]>D_max_snow_thompson,np.nan,fit_dic["thompson_snow"]["v(D_array)"])

        GSFC_snow_handle, = ax.plot(fit_dic["diam"],fit_dic["thompson_snow"]["v(D_array)"],color='green',label="THOM snow",linestyle='-.')

    #show legend
    ax.legend() 
    #expl1, = ax.plot(np.nan,np.nan,label="model assumptions:",linestyle=""); expl2, = ax.plot(np.nan,np.nan,label="fit to simulated particles:",linestyle="");axes.legend(handles=[expl2,Atlasfit_Nmono1,Atlasfit_Nmonoallagg,expl1,SB_ci_handle,SB_snow_handle,MC_handle,P3_handle]) #uncomment this line for a special order

    #make labels
    ax.set_xlabel("diameter D / m")
    ax.set_ylabel("terminal velocity $v_{term}$ / $m s^{-1}$" ) #TODO: plot also the untis of these properties
    
    #change the axis
    if "xlog" in scale:
        ax.set_xscale('log')
    #else:
    ax.set_xlim([min(fit_dic["diam"]),max(fit_dic["diam"])])

        
    if "ylog" in scale:
        ax.set_yscale('log')
    else:
        ax.set_ylim([0,2.0])

    #add grid
    #ax.grid(which="both")
    
    return ax

if __name__ == '__main__':
    
    fit_dic = get_model_params()
    
    show_lines = dict()
    #define which lines to show (ATTENTION: default, can be overwritten by optional arguments below):
    show_lines["SB_mD"] = True #show the m-D relation of the SB scheme
    #from previous model settings
    show_lines["SB_powerlaw"] = True #show also the old fit of Axels power-law
    show_lines["SB_Atlas"] = False #show also the old fit of Axels Atlas-type
    #from other models
    show_lines["P3"] = True
    show_lines["MC"] = False
    show_lines["Morr2mom"] = True
    show_lines["GSFC"] = False #single Moment Godard scheme
    show_lines["THOM"] = False #two-moment Thompson scheme



    #add displayed lines to string to distinguish them in the saved files
    add_displayed_lines2string = '' 
    for key in show_lines.keys():
        if show_lines[key]:
            add_displayed_lines2string+= '_' + key
    ###
    #initialize the plot
    ###
    number_of_plots = 3 #28

    #optimize the appearance of the plot (figure size, fonts)
    [fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots,aspect_ratio=0.25)
    #optimize the appearance of the plot (figure size, fonts)
    aspect_ratio=1./4.
    legend_fontsize='medium'
    #increase font sizes
    params = {'legend.fontsize': legend_fontsize,
        'axes.labelsize': 'x-large', #size: Either a relative value of 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large' or an absolute font size, e.g., 12
        'axes.titlesize':'x-large',
        'xtick.labelsize':'x-large',
        'ytick.labelsize':'x-large'}
    pylab.rcParams.update(params)
    #define size of colorbar labels
    colbarlabelsize=11

    for i_ax,scale in enumerate(["xlinylin","xlogylin","xlogylog"]):

        axes[i_ax] = plot_model_vterms(axes[i_ax],fit_dic,scale='xlogylin',show_lines=show_lines)
        
        #add grid
        axes[i_ax].grid(which="both")



    dir_save = '/home/mkarrer/Dokumente/plots/'
    out_filestring = "modelcomparison_vD_" +  add_displayed_lines2string

        
    plt.tight_layout()

    fig.savefig(dir_save + out_filestring + 'spec_ax_figure.png')
    fig.savefig(dir_save + out_filestring + 'spec_ax_figure.pdf')
    print 'The pdf is at: ' + dir_save + out_filestring + 'spec_ax_figure.pdf'
    subprocess.Popen(['evince',dir_save + out_filestring + 'spec_ax_figure.pdf'])

    # Save just the portion _inside_ the second axis's boundaries
    save_axis=1
    for i_ax in range(0,len(axes)):#clear all other axes
        if i_ax==save_axis:
            continue
        axes[i_ax].set_xlabel("")
        axes[i_ax].axes.get_xaxis().set_ticklabels([])
    extent = axes[save_axis].get_window_extent().transformed(fig.dpi_scale_trans.inverted())

    fig.savefig(dir_save + 'tmp.pdf',bbox_inches=extent.expanded(1.5, 1.4),dpi=400)

    subprocess.call('cp ' + dir_save + 'tmp.pdf' + ' ' + dir_save + '4paper/' + out_filestring + '.pdf',shell=True)
    plt.clf()
    plt.close()
