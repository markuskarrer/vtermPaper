import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
from matplotlib import rc
import pylab
params = {'legend.fontsize': 6,
          'legend.handlelength': 2}
pylab.rcParams.update(params)
#from IPython.core.debugger import Tracer ; Tracer()()

import processing_Jussis_aggregate_model_plots4paper_monodep_binary_Deq
#merge the plots of the m-D/A-D plots  #inspired by: https://jonchar.net/notebooks/matplotlib-styling/

# Set the global font to be DejaVu Sans, size 10 (or any other sans-serif font of your choice!)
rc('font',**{'family':'sans-serif','sans-serif':['DejaVu Sans'],'size':8})

# Set the font used for MathJax - more on this later
rc('mathtext',**{'default':'regular'})
import matplotlib.pylab as pylab
params = {'legend.fontsize': 6,
          'legend.handlelength': 2,
          'figure.figsize': (15, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
def stylize_axes(ax, title, xlabel, ylabel, xticks, yticks, xticklabels, yticklabels):
    """Customize axes spines, title, labels, ticks, and ticklabels."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.xaxis.set_tick_params(top='off', direction='out', width=1)
    ax.yaxis.set_tick_params(right='off', direction='out', width=1)
    
    ax.set_title(title)
    
    #ax.set_xlabel(xlabel)
    #ax.set_ylabel(ylabel)
    
    #ax.set_xticks(xticks)
    #ax.set_yticks(yticks)
    
    #ax.set_xticklabels(xticklabels)
    #ax.set_yticklabels(yticklabels)


#select monotype
for monotype in ["plate"]: #,"dendrite","column","needle","rosette","mixcoldend1","mixcolumndend"]:

    #define the figure grid
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(8,6))

    #do all the processing and plotting externally
    fig,axes,im,tick_locs,N_mono_list = processing_Jussis_aggregate_model_plots4paper_monodep_binary_Deq.read_and_plot(fig,ax,[monotype],plot_vars=["mass","area"])


    titles = [' ', ' ', ' ', ' ']

    xlims = ((-1, 5),(-1, 5),(-1, 5),(-1, 5))
    ylims = ((-1, 5),(-1, 5),(-1, 5),(-1, 5))
    bar_ylims = ((-1, 5),(-1, 5),(-1, 5),(-1, 5))

    xlabels = ['x1', 'x2', 'x3', 'x4']
    ylabels = ['y1', 'y2', '', '']
    letters=['a', 'b', 'c','d','e','f']
    xticks = range(1,6) #dummy arguments
    xticklabels = range(1,6) #dummy arguments

    #from IPython.core.debugger import Tracer ; Tracer()()
    #add one colorbar at right
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.91, 0.2, 0.015, 0.6])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("$N_{mono}$")
    #shift ticks to the middle of the color
    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels(N_mono_list[::3])

    for i, axis in enumerate(axes): #ax.flat):
        if xlabels[i]=='':
            fig.delaxes(axis)
            continue
        # Customize y ticks on a per-axes basis
        yticks = np.linspace(axis.get_ylim()[0], axis.get_ylim()[1], 5)
        yticklabels = yticks
        stylize_axes(axis, titles[i], xlabels[i], ylabels[i], xticks, yticks, xticklabels, yticklabels)
        axis.text(0.02, 0.98,letters[i] + ')',fontsize=14, fontweight='bold',
            horizontalalignment='left',
            verticalalignment='top',
            transform = axis.transAxes)

    ###########
    ###save the plot (and open it)
    ###########
    #plt.tight_layout()
    dir_save = '/home/mkarrer/Dokumente/plots/4paper/'
    out_filestring = "mDAD" + monotype

    plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
    plt.savefig(dir_save + out_filestring + '.png', dpi=100)
    print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
    subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
