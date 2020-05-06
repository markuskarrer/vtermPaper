import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
from matplotlib import rc
from collections import OrderedDict
from IPython.core.debugger import Tracer ; debug=Tracer()


import processing_Jussis_aggregate_model_plots4paper_aggvDfits
import processing_Jussis_aggregate_model_plots4paper_monovDfits
import compare_geom
#create plot of vterm for different monomer types #inspired by: https://jonchar.net/notebooks/matplotlib-styling/

# Set the global font to be DejaVu Sans, size 10 (or any other sans-serif font of your choice!)
rc('font',**{'family':'sans-serif','sans-serif':['DejaVu Sans'],'size':8})

# Set the font used for MathJax - more on this later
rc('mathtext',**{'default':'regular'})


def stylize_axes(ax, title, xlabel, ylabel, xticks, yticks, xticklabels, yticklabels):
    """Customize axes spines, title, labels, ticks, and ticklabels."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.xaxis.set_tick_params(top='off', direction='out', width=1)
    ax.yaxis.set_tick_params(right='off', direction='out', width=1)
    
    ax.set_title(title)
    
    #ax.set_xlabel(xlabel)
    #ax.set_ylabel(ylabel)
    
    
    #ax.set_xticklabels(xticklabels)
    #ax.set_yticklabels(yticklabels)

    ax.grid(b=True,which="both")
def stylize_axes2(ax, title, xlabel, ylabel, xticks=None, yticks=None, xticklabels=None, yticklabels=None,changeyticks=False,xlims=None,ylims=None):
    """Customize axes spines, title, labels, ticks, and ticklabels."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.xaxis.set_tick_params(top=False, direction='out', width=1)
    ax.yaxis.set_tick_params(right=False, direction='out', width=1)
    
    ax.set_title(title)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    if xticks!=None:    
        ax.set_xticks(xticks)
    if yticks!=None:    
        ax.set_yticks(yticks)
    if xlims!=None:
        ax.set_xlim(xlims)
    if ylims!=None:
        ax.set_ylim(ylims)
    if xticklabels!=None:
        ax.set_xticklabels(xticklabels)
    if yticklabels!=None:
        ax.set_yticklabels(yticklabels)
    
    ax.grid(b=True,which="both")



#define the figure grid
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8,4)) #for bigger labels decrease figsize

#axes[0] = processing_Jussis_aggregate_model_plots4paper_aggvDfits.read_and_plot(axes[0],hydro_model
#    ="vterm_bohm",forw_mod=False)[0]
#read data
data_dic = compare_geom.read_data()
data_dic["fit_dic"]["particle_type"]=["plate","dendrite","column","needle","mixcoldend1","mixcolumndend"]

axes[0] = compare_geom.comp_prop(data_dic,axes[0],prop="spheredensity_Nmono1",get_reldiff=False,show_lit=False)
axes[1] = compare_geom.comp_prop(data_dic,axes[1],prop="spheredensity",get_reldiff=False,show_lit=False)


titles = [' ', ' ', ' ', ' ']
'''
xlims = ((-1, 5),(-1, 5),(-1, 5),(-1, 5))
ylims = ((-1, 5),(-1, 5),(-1, 5),(-1, 5))
bar_ylims = ((-1, 5),(-1, 5),(-1, 5),(-1, 5))
'''
xlabels = [r'$D_{max}$ [m]', r'$D_{max}$ [m]',r'$D_{max}$ [m]', '']
ylabels = [r'$\rho_{sphere}$ [kg/m3]',r'$\rho_{sphere}$ [kg/m3]',r'$\rho_{oblate spheroid}$ [kg/m3]', '']
letters=['a) $N_{mono}=1$', 'b) $N_{mono}=1$', 'c) $N_{mono}>1$','d) $N_{mono}>1$','e)','f)']
xticks = range(1,6)
xticklabels = range(1,6)


legs_handles = []
legs_labels = []
number_of_used_plots= 0
for i, ax in enumerate(axes.flat):
    if i==0:
        stylize_axes2(ax, titles[i], xlabels[i], ylabels[i],  xlims=[1e-4,4e-2],ylims=[1e-0,1e3])
    elif i==1:
        #yticks = np.linspace(0.0, 2.0, 5)
        #yticklabels = yticks
        stylize_axes2(ax, titles[i], xlabels[i], ylabels[i], xlims=[1e-4,4e-2],ylims=[1,1e3]  )
        #yticks = np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], 5)
        #yticklabels = yticks
    elif i==2:
        stylize_axes2(ax, titles[i], xlabels[i], ylabels[i], xlims=[1e-4,4e-2],ylims=[1,1e3]  )
    elif i==3:
        yticks = np.linspace(0.0, 2.0, 5)
        yticklabels = yticks
        stylize_axes(ax, titles[i], xlabels[i], ylabels[i], xticks, yticks, xticklabels, yticklabels)
 # Customize y ticks on a per-axes basis
        yticks = np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], 5)
        yticklabels = yticks
    #stylize_axes(ax, titles[i], xlabels[i], ylabels[i], xticks, yticks, xticklabels, yticklabels)
    ax.text(0.02, 0.98,letters[i],fontsize=14, fontweight='bold',
     horizontalalignment='left',
     verticalalignment='top',
     transform = ax.transAxes)
    #add the legend handles and labels from the current axes
    legs_handles.append(ax.axes.get_legend_handles_labels()[0])
    legs_labels.append(ax.axes.get_legend_handles_labels()[1])
    if i==(len(axes.flat)-1): #joined legend
        #concatenate legend handles and labels
        handles_flat = [item for sublist in legs_handles for item in sublist]
        labels_flat = [item for sublist in legs_labels for item in sublist]        
        
        #create and ordered dictionary such that duplicates are remove duplicates
        by_label = OrderedDict(zip(labels_flat, handles_flat))
        # add a legend for the whole figure
        fig.legend(by_label.values(),by_label.keys(),bbox_to_anchor=(0.5,0), loc="lower center", 
                                bbox_transform=fig.transFigure, ncol=3)
     
###########
###save the plot (and open it)
###########
plt.tight_layout(rect=[0,0.05,1,1]) #rect : tuple (left, bottom, right, top), optional
dir_save = '/home/mkarrer/Dokumente/plots/4paper/'
out_filestring = "densities"

plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
plt.savefig(dir_save + out_filestring + '.png', dpi=100)
print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
