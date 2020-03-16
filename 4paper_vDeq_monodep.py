import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
from matplotlib import rc
#from IPython.core.debugger import Tracer ; Tracer()()

import processing_Jussis_aggregate_model_plots4paper_monodep_binary
#merge the plots of the m-D/A-D plots  #inspired by: https://jonchar.net/notebooks/matplotlib-styling/

# Set the global font to be DejaVu Sans, size 10 (or any other sans-serif font of your choice!)
fontsize=8
rc('font',**{'family':'sans-serif','sans-serif':['DejaVu Sans'],'size':fontsize})

# Set the font used for MathJax - more on this later
rc('mathtext',**{'default':'regular'})

import matplotlib.pylab as pylab
params = {'legend.fontsize': 10,
          'legend.handlelength': 2,
          'figure.figsize': (15, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

def stylize_axes(ax, title, xlabel, ylabel, xticks, yticks, xticklabels, yticklabels):
    """Customize axes spines, title, labels, ticks, and ticklabels."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.xaxis.set_tick_params(top='off', direction='out', width=1)
    ax.yaxis.set_tick_params(right='off', direction='out', width=1)
    
    #ax.set_title(title)
    
    
    #ax.set_xlabel(xlabel)
    #ax.set_ylabel(ylabel)
    
    #ax.set_xticks(xticks)
    #ax.set_yticks(yticks)
    
    #ax.set_xticklabels(xticklabels)
    #ax.set_yticklabels(yticklabels)

plate_and_needle=True

#select monotype
for monotype in ["plate"]: # ,"dendrite","column","needle","rosette"]:
    if plate_and_needle and monotype!="plate":
        continue
    if plate_and_needle:
        #define the figure grid
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16,6))
    else:
        #define the figure grid
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
    vterm_model="bohm"
    #vterm_model="KC05"
    #do all the processing and plotting externally
    #fig,axes,im,tick_locs,N_mono_list = processing_Jussis_aggregate_model_plots4paper_monodep_binary.read_and_plot(fig,[ax],[monotype],plot_vars=["vterm_"+vterm_model],use_Deq=True)

    if plate_and_needle:
        fig,axes[0],im,tick_locs,N_mono_list = processing_Jussis_aggregate_model_plots4paper_monodep_binary.read_and_plot(fig,[axes[0]],["plate"],plot_vars=["vterm_"+vterm_model],use_Deq=True)
        fig,axes[1],im,tick_locs,N_mono_list = processing_Jussis_aggregate_model_plots4paper_monodep_binary.read_and_plot(fig,[axes[1]],["needle"],plot_vars=["vterm_"+vterm_model],use_Deq=True)
    else:
        fig,axes,im,tick_locs,N_mono_list = processing_Jussis_aggregate_model_plots4paper_monodep_binary.read_and_plot(fig,[axes],[monotype],plot_vars=["vterm_"+vterm_model])
    #some labelling
    titles = [' ', ' ', ' ', ' ']
    xlabels = ['x1', 'x2', 'x3', 'x4']
    ylabels = ['y1', 'y2', '', '']
    if plate_and_needle:
        letters=['a) plate', 'b) needle', 'c)','d)','e)','f)']
    else:
        letters=['', 'b)', 'c)','d)','e)','f)']
    xticks = range(1,6) #dummy arguments
    xticklabels = range(1,6) #dummy arguments
    for i, ax in enumerate(axes): #ax.flat):
            
        if plate_and_needle:
            ax=ax[0] #somehow this is still an array of size 1 
        ax.set_xlim([5e-5,5e-3])
        #if i==0:
        #    if plate_and_needle:
        #        ax.legend(loc="upper left", bbox_to_anchor=(0.0,0.9))
        #    else:
        #        ax.legend(loc="upper left")
        #else:
        ax.get_legend().remove()
        if xlabels[i]=='':
            fig.delaxes(ax)
            continue
        ax.set_ylabel("$v_{term}$ [m/s]")
        # Customize y ticks on a per-axes basis
        yticks = np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], 5)
        yticklabels = yticks
        stylize_axes(ax, titles[i], xlabels[i], ylabels[i], xticks, yticks, xticklabels, yticklabels)
        ax.text(0.02, 0.98,letters[i],fontsize=14, fontweight='bold',
            horizontalalignment='left',
            verticalalignment='top',
            transform = ax.transAxes,
            zorder=10)

    ###########
    ###save the plot (and open it)
    ###########
    #plt.tight_layout()
    dir_save = '/home/mkarrer/Dokumente/plots/4paper/'
    if plate_and_needle:
        monotype="plateandneedle"
    if vterm_model!="bohm":
        out_filestring = "vDeq" + monotype
    else:
        out_filestring = "vDeq" + monotype +'_'+ vterm_model

    plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
    plt.savefig(dir_save + out_filestring + '.png', dpi=100)
    print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
    subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
