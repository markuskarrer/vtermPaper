import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
from matplotlib import rc
#from IPython.core.debugger import Tracer ; Tracer()()

#merge the plots of ... for ... #inspired by: https://jonchar.net/notebooks/matplotlib-styling/

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
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    
    ax.set_xticklabels(xticklabels)
    ax.set_yticklabels(yticklabels)



#define the figure grid
fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(8,6))

'''
y = data.mean()
y_all = data.values
x = np.arange(len(means))
error = data.std()


xlims = (-1, 5)
ylims = (-5, 15)
bar_ylims = (0, 15)
'''

###call the functions which plot the subplots

ax[0][0].plot(range(1,6),range(1,6))
ax[0][1].plot(range(1,6),range(1,6))
#ax[1][1].plot(range(1,6),range(1,6))

'''
custom_lineplot(ax[0][0], x, y, error, xlims, ylims)
custom_scatterplot(ax[0][1], x, y, error, xlims, ylims)
custom_barchart(ax[1][0], x, y, error, xlims, bar_ylims, error_kw)
custom_boxplot(ax[1][1], x, y_all, error, xlims, ylims)
'''

titles = ['', '', '', '']
'''
xlims = ((-1, 5),(-1, 5),(-1, 5),(-1, 5))
ylims = ((-1, 5),(-1, 5),(-1, 5),(-1, 5))
bar_ylims = ((-1, 5),(-1, 5),(-1, 5),(-1, 5))
'''
xlabels = ['x1', 'x2', '', '']
ylabels = ['y1', 'y2', '', '']
letters=['a)', 'b)', 'c)','d)','e)','f)']
xticks = range(1,6)
xticklabels = range(1,6)


number_of_used_plots= 0
for i, axes in enumerate(ax.flat):
    if xlabels[i]=='':
        fig.delaxes(axes)
        continue
    # Customize y ticks on a per-axes basis
    yticks = np.linspace(axes.get_ylim()[0], axes.get_ylim()[1], 5)
    yticklabels = yticks
    stylize_axes(axes, titles[i], xlabels[i], ylabels[i], xticks, yticks, xticklabels, yticklabels)
    axes.text(0.02, 0.98,letters[i],fontsize=14, fontweight='bold',
     horizontalalignment='left',
     verticalalignment='top',
     transform = axes.transAxes)

###########
###save the plot (and open it)
###########
plt.tight_layout()
dir_save = '/home/mkarrer/Dokumente/plots/4paper/'
out_filestring = "template"

plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
plt.savefig(dir_save + out_filestring + '.png', dpi=100)
print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])