import matplotlib.pyplot as plt
import numpy as np
import subprocess
from matplotlib import rc
from collections import OrderedDict
#own functions
import compare_geom
import processing_Jussis_aggregate_model_plots4paper_aggvDfits

from IPython.core.debugger import Tracer ; debug = Tracer()
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
def stylize_axes(ax, title, xlabel, ylabel, xticks=None, yticks=None, xticklabels=None, yticklabels=None,changeyticks=False,xlims=None,ylims=None):
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
    
    ax.grid(which="both")



#define the figure grid
fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(8,8))

#read data
data_dic = compare_geom.read_data()
###call the functions which plot the subplots
ax[0][0] = compare_geom.comp_prop(data_dic,ax[0][0],prop="mass",get_reldiff=False)
ax[0][1] = compare_geom.comp_prop(data_dic,ax[0][1],prop="area",get_reldiff=False)

#ax[1][0] = .compare_geomcomp_vterm_discussion(data_dic,ax[1][0])

#ax[1][0] = compare_geom.comp_vterm_discussion(data_dic,ax[1][0])
ax[1][0] = processing_Jussis_aggregate_model_plots4paper_aggvDfits.read_and_plot([ax[1][0]],hydro_model="vterm_bohm",forw_mod=True,called_by_main=False)[0]
ax[1][1] = processing_Jussis_aggregate_model_plots4paper_aggvDfits.read_and_plot([ax[1][1]],hydro_model="vterm_bohm",forw_mod=True,called_by_main=False)[0]
ax[1][1].set_xscale("linear")
ax[1][0].set_xscale("linear")

titles = ['', '', '', '']
xlabels = [r'$D_{max}$ [m]', r'$D_{max}$ [m]', r'$D_{max,side}$ [mm]', r'$D_{max,side}$ [mm]']
ylabels = [r'm [kg]',r'A [$m^2$]', r'$v_{term}$ [m/s]', r'$v_{term}$ [m/s]']
letters=['a) ', 'b) ', 'c)','d)','e','f']
xticks = np.logspace(-4,-2,20) #ATTENTION: currently not used
xticklabels = xticks


number_of_used_plots= 0
legs_handles = []
legs_labels = []
for i, axes in enumerate(ax.flat):
    if xlabels[i]=='':
        fig.delaxes(axes)
        continue
    # Customize y ticks on a per-axes basis
    yticks = np.linspace(0.0, 2.0, 5)
    yticklabels = yticks
    if i==0:
        stylize_axes(axes, titles[i], xlabels[i], ylabels[i],  xlims=[1e-4,1e-2],ylims=[2e-11,4e-6])
        #axes.grid(which="both") #this is turned off when called twice (here thirs call)
    if i==1:
        stylize_axes(axes, titles[i], xlabels[i], ylabels[i],  xlims=[1e-4,1e-2],ylims=[2e-9,3e-5])
        #axes.grid(which="both") #this is turned off when called twice (here thirs call)
    elif i==2:
        stylize_axes(axes, titles[i], xlabels[i], ylabels[i], yticks=yticks,  xlims=[1e-4,2e-3],ylims=[0.0,1.0])#,xticks=[1e-4,1e-3]) #,yticks=[0.0,0.5]) 

        axes.set_xticklabels(axes.get_xticks()*1000.)#convert to mm
        axes.grid(which="both")
    elif i==3:
        stylize_axes(axes, titles[i], xlabels[i], ylabels[i],yticks=yticks, xlims=[2e-3,1e-2],ylims=[0.3,1.6]) #,yticks=[0.0,0.5]) 
        axes.set_xticklabels(axes.get_xticks()*1000.)#convert to mm
        axes.grid(which="both")
    #add a label for the subplot
    axes.text(0.02, 0.99,letters[i],fontsize=10, fontweight='bold',
     horizontalalignment='left',
     verticalalignment='top',
     transform = axes.transAxes)
    #add the legend handles and labels from the current axes
    legs_handles.append(axes.axes.get_legend_handles_labels()[0])
    legs_labels.append(axes.axes.get_legend_handles_labels()[1])
    if i==3: #special legend for zoom in vterm
        #concatenate legend handles and labels
        handles_flat = [item for sublist in legs_handles for item in sublist]
        labels_flat = [item for sublist in legs_labels for item in sublist]        
        
        #create and ordered dictionary such that duplicates are remove duplicates
        by_label = OrderedDict(zip(labels_flat, handles_flat))
        # add a legend for the whole figure
        fig.legend(by_label.values(),by_label.keys(),bbox_to_anchor=(0.5,0), loc="lower center", 
                                bbox_transform=fig.transFigure, ncol=5)
        #debug()
        #axes.minorticks_off()
###########
###save the plot (and open it)
###########
plt.tight_layout(rect=[0,0.1,1,1]) #rect : tuple (left, bottom, right, top), optional

dir_save = '/home/mkarrer/Dokumente/plots/4paper/'
out_filestring = "discussion_comp"

plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
plt.savefig(dir_save + out_filestring + '.png', dpi=100)
print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
