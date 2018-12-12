# coding: utf-8
#import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.pylab as pylab
import csv #to read the txt files
import os
import subprocess
#import other self-defined functions
import __postprocess_McSnow
import __postprocess_SB
import __fallspeed_relations
#from IPython.core.debugger import Tracer ; Tracer()()

'''
this code reads in properties from Jussis aggregate model and calculate+plots the following fallspeed
'''

#define where the txt files with properties of the aggregates are stored
prop_file_folder = "/home/mkarrer/Dokumente/Jussis_aggregation_model/Leonie_output/"




#read the properties of the aggregates into the particle_dic dictionary
particle_types = ["dendrite","needles"]
for particle_type in particle_types: 
    
    particle_dic = dict() #dictionary which contains information about the type of the pristine crystals (dentrites, needles, plates, ...) the underlying distribution of the pristine crystals and the properties of the aggregates (m,A,D,...)
    particle_dic["particle_type"] = [] #type of the pristine crystals
    particle_dic["N_monomer"] = [] #number of monomers
    particle_dic["mass"] = [] #particle mass
    particle_dic["area"] = [] #particle area
    particle_dic["diam"] = [] #particle diameter
    
    for mean_size in ["50","100","650","1000"]: #loop over the mean size of the distribution of the pristine crystals
        with open(prop_file_folder + particle_type +"_"+ mean_size +"_properties.txt","rb") as txtfile: #TODO: this is just one file! #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in
            prop_reader = csv.reader(txtfile, delimiter=' ', quoting=csv.QUOTE_NONNUMERIC, lineterminator=os.linesep) #quoting avoids '' for formatted string; lineterminator avoids problems with system dependend lineending format https://unix.stackexchange.com/questions/309154/strings-are-missing-after-concatenating-two-or-more-variable-string-in-bash?answertab=active#tab-top
            for i_row,row_content in enumerate(prop_reader):
                particle_dic["particle_type"] = np.append(particle_dic["particle_type"],particle_type)
                particle_dic["N_monomer"] = np.append(particle_dic["N_monomer"],i_row)
                particle_dic["mass"] = np.append(particle_dic["mass"],row_content[0])
                particle_dic["area"] = np.append(particle_dic["area"],row_content[1])
                particle_dic["diam"] = np.append(particle_dic["diam"],row_content[2])

    #calculate the fall speeds for each available velocity model
    for velocity_model in ["HW10"]:
        particle_dic["vterm_" + velocity_model] = __fallspeed_relations.calc_vterm(velocity_model,particle_dic["mass"],particle_dic["diam"],particle_dic["area"])


    print particle_dic,particle_dic["particle_type"].shape

    ###
    #plot the m-D, A-D and v-D relations
    ###

    #increase font sizes
    params = {'legend.fontsize': 'large',
        'figure.figsize': (15, 5),
        'axes.labelsize': 'x-large', #size: Either a relative value of 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large' or an absolute font size, e.g., 12
        'axes.titlesize':'x-large',
        'xtick.labelsize':'x-large',
        'ytick.labelsize':'x-large'}
    pylab.rcParams.update(params)
    #define figure
    number_of_plots = 20
    figsize_height = 6.0/2.0*number_of_plots
    fig, axes = plt.subplots(nrows=number_of_plots, ncols=1, figsize=(8.0,figsize_height))

    ###
    #subplot 1: m vs D
    ###
    #axes2 = axes[0].twinx()
    #plot in loglog


    for i_prop,prop in enumerate(["mass","area","vterm_HW10"]):
        #minval = np.array([])[i_prop]
        #loop for different scales
        for i_scalex,scalex in enumerate(["linear","log"]):
            for i_scaley,scaley in enumerate(["linear","log"]):
                i_ax = i_prop*4+i_scalex*2+i_scaley
                #change the scale (according to the loop)
                if scalex=="log":
                    axes[i_ax].set_xscale(scalex)
                if scaley=="log":
                    axes[i_ax].set_yscale(scaley)

                #make the scatter-splot
                #for particle_type,symbol in zip(particle_types,["s","o"]):
                #select diameter array ,property array and monomer number with just these particle types
                diam_for_this_particle_type = particle_dic["diam"] #[particle_dic["particle_type"]==particle_type]
                prop_for_this_particle_type = particle_dic[prop] #[particle_dic["particle_type"]==particle_type]
                Nmono_for_this_particle_type = particle_dic["N_monomer"] #[particle_dic["particle_type"]==particle_type]
                im = axes[i_ax].scatter(diam_for_this_particle_type,prop_for_this_particle_type,s=1,c=Nmono_for_this_particle_type,rasterized=True,norm=colors.LogNorm())

                #make fast labels
                axes[i_ax].set_xlabel("diameter / m")
                axes[i_ax].set_ylabel(prop)

                if scaley=="linear":
                    axes[i_ax].set_ylim([0,np.array([1e-5,1e-4,1.5])[i_prop]])
                #add colorbar
                cbar = fig.colorbar(im,ax=axes[i_ax])
                cbar.set_label("monomer number")
                # recompute the ax.dataLim
                #axes[i_ax].relim()
                # update ax.viewLim using the new dataLim
                #axes[i_ax].autoscale_view()



    plt.tight_layout()
    dir_save = '/home/mkarrer/Dokumente/plots/'
    if not os.path.exists(dir_save): #create direktory if it does not exists
        os.makedirs(dir_save)
    out_filestring = "Jagg_mD_AD_vD_" + particle_type
    plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
    plt.savefig(dir_save + out_filestring + '.png', dpi=400)
    print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
    subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
    plt.clf()
    plt.cla()
    plt.close()
