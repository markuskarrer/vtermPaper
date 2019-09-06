# coding: utf-8
#import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.pylab as pylab
import os
import subprocess
from scipy.stats import gaussian_kde
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
tumbling=False 
take_mono_prop_theor = True #so far only True is implemented

#show which types of fit
show_lines=dict()
show_lines["powerlaw_fits"] = False
show_lines["highest_diam_for_fit"] = -999 #-999 for all -2 for increasing v range only
show_lines["Atlas_fit"] = True #show the Atlas-type fit for cloud ice & snow
#from previous model settings
show_lines["SB_powerlaw"] = False #show also the old fit of Axels power-law

##general settings for the plotting
number_of_plots = 2

#optimize the appearance of the plot (figure size, fonts)
[fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots,legend_fontsize='medium')


#add displayed lines to string to distinguish them in the saved files
add_displayed_lines2string = '' 
for key in show_lines.keys():
    if show_lines[key]:
        add_displayed_lines2string+= '_' + key

#custom colorbar
cmap = plt.cm.hot #brg
#modify lowest color in colormap
cmaplist = [cmap(i) for i in range(cmap.N)]
cmaplist = cmaplist[30:256-30] #crop colormap to avoid too dark/light scatterer #out of 256 colors
cmap = colors.LinearSegmentedColormap.from_list('mcm',cmaplist, cmap.N)

particle_types = ["column"] #["plate","dendrite","mixdendneedle","needle","column"] #"mixdendneedle""needle","column","plate","dendrite"] #["needle","column","plate","dendrite","bullet","rosette"] # ,"bullet"rosette
for i_particle_type,particle_type in enumerate(particle_types):


    #read the properties of the individual particles from the files
    particle_dic,N_mono_list = __tools_for_processing_Jagg.read_particle_prop_files(prop_file_folder = "/data/optimice/aggregate_model/Jussis_aggregates_bugfixedrotation/",
                                                                            D_small_notuse=1e-4, #ATTENTION: particle smaller than D_small_notuse are not considered (e.g. because of resolution issues)
                                                                            N_small_notuse=1, #ATTENTION:  monomer numbers smaller than N_small_notuse are not considered
                                                                            grid_res_array = [5e-6,10e-6], #[5e-6,10e-6], #[1e-6,5e-6,10e-6], #array of grid resolutions (defines the files which are read in)
                                                                            particle_type = particle_type, #define the particle_type
                                                                            test_with_small_sample = False
                                                                            )
    #show the dictionary (e.g. to check that it's not empty because the path was wrong)
    print particle_dic
    i_ax=0
    for xscale in ["log","linear"]:
        for Dmax_side in ["Dmax_xz"]: #,"Dmax_yz","Dmax_xz_aligned","Dmax_yz_aligned"]:
            #calculate the density with kde
            if xscale=="linear":
                xy=np.vstack([particle_dic["diam"],particle_dic[Dmax_side]/particle_dic["diam"]])
            elif xscale=="log":
                xy=np.vstack([np.log10(particle_dic["diam"]),particle_dic[Dmax_side]/particle_dic["diam"]])
            kde_density = gaussian_kde(xy,bw_method="scott")(xy) #"scott" is the default also

            im = axes[i_ax].scatter(particle_dic["diam"],particle_dic[Dmax_side]/particle_dic["diam"],s=0.1,c=kde_density/max(kde_density),cmap=cmap,rasterized=True)
            axes[i_ax].set_xscale(xscale)
            
            #make labels
            axes[i_ax].set_xlabel("diameter D / m")
            axes[i_ax].set_ylabel(r"$Dmax_{sideproj}$/$Dmax$" + " / m s-1" ) #TODO: plot also the untis of these properties

            #change the axis
            axes[i_ax].set_xlim([1e-4,4e-2]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
            axes[i_ax].set_ylim([0.2,1.0]) 
            axes[i_ax].grid(b=True,which="both",axis="both")
            #add colorbar
            cbar = fig.colorbar(im,ax=axes[i_ax])
            cbar.set_label("density")
            i_ax+=1

#save the plot (and open it)
plt.tight_layout()
dir_save = '/home/mkarrer/Dokumente/plots/Jagg/'
if not os.path.exists(dir_save): #create direktory if it does not exists
    os.makedirs(dir_save)
all_particle_types = "".join(particle_types)
out_filestring = "Dmax_def"

plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
plt.savefig(dir_save + out_filestring + '.png', dpi=200)
print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])