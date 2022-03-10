import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
from matplotlib import rc
#from IPython.core.debugger import Tracer ; Tracer()()

import processing_Jussis_aggregate_model_plots4paper_vDfits
#merge the plots of ... for ... #inspired by: https://jonchar.net/notebooks/matplotlib-styling/

# Set the global font to be DejaVu Sans, size 10 (or any other sans-serif font of your choice!)
rc('font',**{'family':'sans-serif','sans-serif':['DejaVu Sans'],'size':8})

# Set the font used for MathJax - more on this later
rc('mathtext',**{'default':'regular'})
import matplotlib.pylab as pylab
params = {'legend.fontsize': 6,
          'legend.handlelength': 2,
          'figure.figsize': (15, 5),
         'axes.labelsize': 30,
         'axes.titlesize':'xx-large',
         'xtick.labelsize':20,
         'ytick.labelsize':60} #THIS DOESNT SEEM TO HAVE AN EFFEECT



def stylize_axes(ax, title, xlabel, ylabel, xticks, yticks, xticklabels, yticklabels):
    """Customize axes spines, title, labels, ticks, and ticklabels."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.xaxis.set_tick_params(top='off', direction='out', width=1)
    ax.yaxis.set_tick_params(right='off', direction='out', width=1)
    
    #ax.set_title(title)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    #ax.set_xticks(xticks)
    #ax.set_yticks(yticks)
    
    #ax.set_xticklabels(xticklabels)
    #ax.set_yticklabels(yticklabels)
    #'''
#define the figure grid
#select monomer type
#monotypes=["plate","dendrite","column","needle","rosette","mixcolumndend"]

for hydro_model in ["vterm_bohm","vterm_HW10","vterm_KC05"]:
    monotypes=["plate","dendrite","column","needle","mixcolumndend"]
    rows = len(monotypes)
    fig, ax = plt.subplots(nrows=rows, ncols=1, figsize=(4,4*rows)) #for bigger labels decrease figsize
    for i_ax,monotype in enumerate(monotypes): # ,"dendrite","column","needle","rosette","mixcoldend1","mixcolumndend"]:


        ###call the functions which plot the subplots

        #do all the processing and plotting externally
        
        #ax[0] = processing_Jussis_aggregate_model_plots4paper_vDfits.read_and_plot(fig,[ax[0]],[monotype],powerlawAtlas="powerlaw",hydro_model
        #    ="vterm_bohm")[0]
        
        ax[i_ax] = processing_Jussis_aggregate_model_plots4paper_vDfits.read_and_plot(fig,[ax[i_ax]],[monotype],powerlawAtlas="Atlas",hydro_model
            =hydro_model,forKamil=True)[0]
        
        #ax[2] = processing_Jussis_aggregate_model_plots4paper_vDfits.read_and_plot(fig,[ax[2]],[monotype],powerlawAtlas="schemes",hydro_model
        #    ="vterm_bohm")[0]
        #from IPython.core.debugger import Tracer ; Tracer()()


        titles = ['', '', '', '']
        '''
        xlims = ((-1, 5),(-1, 5),(-1, 5),(-1, 5))
        ylims = ((-1, 5),(-1, 5),(-1, 5),(-1, 5))
        bar_ylims = ((-1, 5),(-1, 5),(-1, 5),(-1, 5))
        '''
        #xlabels = ['$D_{max} [m]$', '$D_{max} [m]$', '$D_{max} [m]$', '']
        #ylabels = ['$v_{term}$ [m/s]', '$v_{term}$ [m/s]', '$v_{term}$ [m/s]', '']
        
        #letters=['a) '+ monotypes[0], 'b) ' + monotypes[1], 'c) ' ,'d)','e)','f)']
        xticks = range(1,6)
        xticklabels = range(1,6)


        number_of_used_plots= 0
        for i, axes in enumerate(ax):
            #if xlabels[i]=='':
            #    fig.delaxes(axes)
            #    continue
            # Customize y ticks on a per-axes basis
            yticks = np.linspace(axes.get_ylim()[0], axes.get_ylim()[1], 5)
            yticklabels = yticks
            stylize_axes(axes, '', '$D_{max} [m]$', '$v_{term}$ [m/s]', xticks, yticks, xticklabels, yticklabels)
            print i,monotypes[i]
            axes.text(0.02, 0.98,monotypes[i],fontsize=10, fontweight='bold',
             horizontalalignment='left',
             verticalalignment='top',
             transform = axes.transAxes)

    ###########
    ###save the plot (and open it)
    ###########
    plt.tight_layout()
    dir_save = '/home/mkarrer/Dokumente/plots/4Kamilspaper/'
    out_filestring = "vDfits_" + hydro_model 

    plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
    plt.savefig(dir_save + out_filestring + '.png', dpi=100)
    print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
    #subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
