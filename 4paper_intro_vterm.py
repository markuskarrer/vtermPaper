import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
from matplotlib import rc

#own functions
import compare_geom
import compare_fallspeed_parameterizations

#merge the plots of vterm (obs vs. model) for the introduction #inspired by: https://jonchar.net/notebooks/matplotlib-styling/

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

    ax.xaxis.set_tick_params(top=False, direction='out', width=1)
    ax.yaxis.set_tick_params(right=False, direction='out', width=1)
    
    ax.set_title(title)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    #ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    
    #ax.set_xticklabels(xticklabels)
    ax.set_yticklabels(yticklabels)
    
    ax.grid(which="both")



#define the figure grid
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8,3))

#read data
data_dic = compare_geom.read_data()
fit_dic  = compare_fallspeed_parameterizations.get_model_params()

###call the functions which plot the subplots
ax[0] = compare_geom.comp_vterm_intro(data_dic,ax[0])
ax[1] = compare_fallspeed_parameterizations.plot_model_vterms(ax[1],fit_dic)

'''
#show not all lines
show_lines = dict()
show_lines["SB_mD"] = True #show the m-D relation of the SB scheme
#from previous model settings
show_lines["SB_powerlaw"] = True #show also the old fit of Axels power-law
show_lines["SB_Atlas"] = False #show also the old fit of Axels Atlas-type
#from other models
show_lines["P3"] = False
show_lines["MC"] = False
show_lines["Morr2mom"] =False 
show_lines["GSFC"] = False #single Moment Godard scheme
show_lines["THOM"] = False #two-moment Thompson scheme

ax[1] = compare_fallspeed_parameterizations.plot_model_vterms(ax[1],fit_dic,show_lines=show_lines) #show only SB
'''

titles = ['', '', '', '']
'''
xlims = ((-1, 5),(-1, 5),(-1, 5),(-1, 5))
ylims = ((-1, 5),(-1, 5),(-1, 5),(-1, 5))
bar_ylims = ((-1, 5),(-1, 5),(-1, 5),(-1, 5))
'''
xlabels = [r'$D_{max}$ [m]', r'$D_{max}$ [m]', '', '']
ylabels = [r'$v_{term}$ [m/s]',r'$v_{term}$ [m/s]', '', '']
letters=['a) observations', 'b) microphysics schemes', 'c','d','e','f']
xticks = np.logspace(-4,-2,20) #ATTENTION: currently not used
xticklabels = xticks


number_of_used_plots= 0
for i, axes in enumerate(ax.flat):
    if xlabels[i]=='':
        fig.delaxes(axes)
        continue
    # Customize y ticks on a per-axes basis
    yticks = np.linspace(0.0, 2.0, 5)
    yticklabels = yticks
    stylize_axes(axes, titles[i], xlabels[i], ylabels[i], xticks, yticks, xticklabels, yticklabels)
    #add a label for the subplot
    axes.text(0.02, 1.07,letters[i],fontsize=10, fontweight='bold',
     horizontalalignment='left',
     verticalalignment='top',
     transform = axes.transAxes)
    axes.legend(ncol=2,fontsize=7)

###########
###save the plot (and open it)
###########
plt.tight_layout()
dir_save = '/home/mkarrer/Dokumente/plots/4paper/'
out_filestring = "intro_vterm"

plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
plt.savefig(dir_save + out_filestring + '.png', dpi=100)
print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
