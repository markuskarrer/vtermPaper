# coding: utf-8
#import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.pylab as pylab
from matplotlib.ticker import FormatStrFormatter
import os
import subprocess
import csv #to read the txt files
import sys
#import other self-defined functions
import __plotting_functions
import __fallspeed_relations
import __isa
import __postprocess_McSnow
#from IPython.core.debugger import Tracer ; Tracer()()
from matplotlib import rc

'''
this script plots v-D for numerous combinations of m-D and A-D from M96
'''


##general settings for the plotting
number_of_plots = 7

#optimize the appearance of the plot (figure size, fonts)
aspect_ratio=1./3.
legend_fontsize='medium'
#increase font sizes
params = {'legend.fontsize': legend_fontsize,
    'axes.labelsize': 'x-large', #size: Either a relative value of 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large' or an absolute font size, e.g., 12
    'axes.titlesize':'x-large',
    'xtick.labelsize':'x-large',
    'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
#define figure
figsize_height = 1./aspect_ratio*number_of_plots
fig, axes = plt.subplots(nrows=number_of_plots, ncols=1, figsize=(12.5,figsize_height))

#set up array of diameters for displaying the fits
low_diam_log= -4; high_diam_log=np.log10(4e-2)

#set up the diameter array
diam = np.logspace(low_diam_log,high_diam_log,500) #set up array of diameters
diam_log = np.log10(diam)
diam_center = (diam[:-1]+diam[1:])/2 #center of the diameter bins
diam_log_center = 10**((np.log10(diam[:-1])+np.log10(diam[1:]))/2) #logarithmic center of the diameter bins

#also show a sphere?
sphere_dic = dict()
rho_i = 917.6 #define ice density
sphere_dic["particle_type"] = "sphere"
sphere_dic["mass_allagg_coeff"] = [rho_i*np.pi/6.,3]
sphere_dic["area_allagg_coeff"] = [np.pi/4.,2]
for prop in ["mass","area"]:
    sphere_dic[prop] = sphere_dic[prop + "_allagg_coeff"][0]*diam**sphere_dic[prop + "_allagg_coeff"][1]
#calculate v-D
vterm_model = "mitch_heym"
sphere_dic["vterm"] = __fallspeed_relations.calc_vterm(vterm_model,sphere_dic["mass"],diam,sphere_dic["area"],turb_corr="all") #,turb_corr="all" is only for mitch_heym

#use area of Mitchell 1996
M96 = dict()
M96["area_coeff_mixS3"] = [0.2285*10.**(2.*1.88-4.),1.88] #area approximated as assemblage of planar polycrstals
#M96["area_coeff_mixS3"] = [0.2285*10.**(2.*2.2-4.),2.2] #area approximated as assemblage of planar polycrstals #ATTENTION: testing other exponents

##create arrays of m-D A-D relationships
variation_dic = dict()


variation_dic["mass_exp"] = np.arange(1.7,2.31,0.1)
variation_dic["mass_pref"]= np.logspace(-4,-0,20)
variation_dic["area"] = M96["area_coeff_mixS3"][0]*diam**M96["area_coeff_mixS3"][1] #take area from Mitchell 96

allow_mass_fac = 2 #factor by which mass at 100mum can be larger than sphere
allow_vterm = 0.01 #allow velocity at largest diameter

#calculate the arrays
for i_mass_exp,mass_exp in enumerate(variation_dic["mass_exp"]):
    for i_mass_pref,mass_pref in enumerate(variation_dic["mass_pref"]):
        #from IPython.core.debugger import Tracer ; Tracer()()
        #mass_smaller_sphere = (variation_dic["mass_pref"][i_mass_pref]*diam[0]**variation_dic["mass_exp"][i_mass_exp])<(allow_mass_fac*sphere_dic["mass"][0])
        if True: #mass_smaller_sphere:
            variation_dic["mass_" + str(mass_pref) + "_" + str(mass_exp)] =  variation_dic["mass_pref"][i_mass_pref]*diam**variation_dic["mass_exp"][i_mass_exp]
            variation_dic["vterm_" + str(mass_pref) + "_" + str(mass_exp)] = __fallspeed_relations.calc_vterm(vterm_model,variation_dic["mass_" + str(mass_pref) + "_" + str(mass_exp)],diam,variation_dic["area"],turb_corr="all") #,turb_corr="all" is only for mitch_heym

#from IPython.core.debugger import Tracer ; Tracer()()
###
#plot modelled combination of 
###

i_ax=0
for i_mass_exp,mass_exp in enumerate(variation_dic["mass_exp"]):

    #for i_prop,(prop,prop_unit,prop_short) in enumerate(zip(["mass","area","vterm"],["kg","m2","ms-1"],["m","A","vterm"])):
    for i_prop,(prop,prop_unit,prop_short) in enumerate(zip(["vterm"],["ms-1"],["vterm"])):

        if prop=="area":
                axes[i_ax].loglog(diam,variation_dic["area"],linestyle="-",marker="",markevery=20,color="g",label="M96 ")
        if prop in ["mass"]: #plot m/A-D

            #sphere
            axes[i_ax].plot(diam,sphere_dic[prop],linestyle="-",color="gray",label="sphere")

            #variations
            #for i_mass_exp,mass_exp in enumerate(variation_dic["mass_exp"]):
            for i_mass_pref,mass_pref in enumerate(variation_dic["mass_pref"]):
                #mass_smaller_sphere = ((variation_dic["mass_pref"][i_mass_pref]*diam[0]**variation_dic["mass_exp"][i_mass_exp])<(allow_mass_fac*sphere_dic["mass"][0])) && (variation_dic["vterm_" + str(mass_pref) + "_" + str(mass_exp)][-1]>allow_vterm)
                mass_smaller_sphere = (variation_dic["mass_pref"][i_mass_pref]*diam[0]**variation_dic["mass_exp"][i_mass_exp])<(allow_mass_fac*sphere_dic["mass"][0])
                vterm_large_enough  = (variation_dic["vterm_" + str(mass_pref) + "_" + str(mass_exp)][-1]>allow_vterm)

                if mass_smaller_sphere and vterm_large_enough:
                    axes[i_ax].loglog(diam,variation_dic[prop + "_" + str(mass_pref) + "_" + str(mass_exp)],linestyle="-",marker=["","","","","x","x","x"][i_mass_exp],markevery=20,color=["b","r","g","y","orange","cyan","purple"][i_mass_pref],label="{:6.4f}".format(mass_pref) + "D^" +str(mass_exp))


                
        elif prop=="vterm":

            #for scale in ["linear","log"]:
            for scale in ["log"]:

                #plot v-D
                
                #sphere
                axes[i_ax].plot(diam,sphere_dic["vterm"],linestyle="-",color="gray",label="sphere")
                
                #variations
                #for i_mass_exp,mass_exp in enumerate(variation_dic["mass_exp"]):
                for i_mass_pref,mass_pref in enumerate(variation_dic["mass_pref"]):
                    mass_smaller_sphere = (variation_dic["mass_pref"][i_mass_pref]*diam[0]**variation_dic["mass_exp"][i_mass_exp])<(allow_mass_fac*sphere_dic["mass"][0])
                    mass_smaller_sphere = (variation_dic["mass_pref"][i_mass_pref]*diam[0]**variation_dic["mass_exp"][i_mass_exp])<(allow_mass_fac*sphere_dic["mass"][0])
                    vterm_large_enough  = (variation_dic["vterm_" + str(mass_pref) + "_" + str(mass_exp)][-1]>allow_vterm)

                    if mass_smaller_sphere and vterm_large_enough:
                        axes[i_ax].plot(diam,variation_dic[prop + "_" + str(mass_pref) + "_" + str(mass_exp)],linestyle="-",markevery=20,label="{:6.4f}".format(mass_pref) + "D^" +str(mass_exp))
                
    
                #i_ax+=1
            #i_ax+=-2
            
        #do the labelling once in the end
        if prop=="vterm":
            #for scale in ["linear","log"]:
            for scale in ["log"]:
                #labelling
                axes[i_ax].set_xlabel("diameter D / m")
                axes[i_ax].set_ylabel(r"terminal vel. $v_{term}$ / m $s^{-1}$")
                axes[i_ax].set_xscale(scale)
                axes[i_ax].legend(ncol=2,loc='center left', bbox_to_anchor=(1.0, 0.5))
                axes[i_ax].grid(which="both")
                if scale=="linear":
                    axes[i_ax].set_xlim([0,5e-3])
                else:
                    axes[i_ax].set_xlim([1e-4,4e-2])

                axes[i_ax].set_ylim([0.0,2.0])

                i_ax+=1
        else:

            axes[i_ax].grid(which="both")
            axes[i_ax].set_xlim([1e-4,4e-2])
            #labelling
            axes[i_ax].set_xlabel("diameter D / m")
            axes[i_ax].set_ylabel((" / ").join((prop+ ' ' + prop_short,prop_unit)))
            #axes[i_ax].set_xscale(scale)
            axes[i_ax].legend(ncol=2,loc='center left', bbox_to_anchor=(1.0, 0.5))

            i_ax+=1


###########
###save the plot (and open it)
###########
dir_save = '/home/mkarrer/Dokumente/plots/geom_comp/'
out_filestring = "geom_combinations_" + vterm_model


plt.tight_layout()

plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
plt.savefig(dir_save + out_filestring + '.png', dpi=100)
print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])