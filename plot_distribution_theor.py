'''
this script plots some distributions
'''
#from IPython.core.debugger import Tracer ; Tracer()()

#import modules
from matplotlib import pyplot as plt
import numpy as np
import subprocess
import sys
import os

#set up figure
number_of_plots = 1

import matplotlib.pylab as pylab

#increase font sizes
params = {'legend.fontsize': 'small',
    'figure.figsize': (15, 5),
    'axes.labelsize': 'x-large', #size: Either a relative value of 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large' or an absolute font size, e.g., 12
    'axes.titlesize':'x-large',
    'xtick.labelsize':'x-large',
    'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
#define figure
figsize_height = 6.0/2.0*number_of_plots
fig, axes = plt.subplots(nrows=number_of_plots, ncols=1, figsize=(8.0,figsize_height))

number_of_lines=15
colormap = plt.cm.nipy_spectral
colors = [colormap(i) for i in np.linspace(0, 1,number_of_lines)]
axes.set_prop_cycle('color', colors)

#define the distribution and some properties
distribution = 'exponential'
mean_size_array=np.logspace(np.log10(200e-6),np.log10(10000e-6),500) #array of mean sizes
mean_size_array2=np.logspace(np.log10(50e-6),np.log10(200e-6),200) #array of mean sizes

size_lims = [1e-4,3e-3] #cut-offs of the size distribution in m

#create diameter array
D_array= np.logspace(np.log10(size_lims[0]),np.log10(size_lims[1]),1000)
if distribution=='exponential':
    
    axes.plot(np.nan,np.nan,label="i_size, mean_size",linestyle='None')
    for i_mean,mean_size in enumerate(mean_size_array2):
        if (not (i_mean%50)==0) and (not i_mean==(mean_size_array.shape[0]-1)): #i_mean should count all elements but not all are plotted
            continue
        lam = 1./mean_size
        f_dist = lam*np.exp(-lam*D_array) #calculate distribution
        axes.semilogx(D_array,f_dist,label=str("{:03d}       {:.3e}".format(i_mean+1000,mean_size)))
        
    for i_mean,mean_size in enumerate(mean_size_array):
        if (not (i_mean%50)==0) and (not i_mean==(mean_size_array.shape[0]-1)): #i_mean should count all elements but not all are plotted
            continue
        lam = 1./mean_size
        f_dist = lam*np.exp(-lam*D_array) #calculate distribution
        axes.semilogx(D_array,f_dist,label=str("{:03d}       {:.3e}".format(i_mean,mean_size)))


axes.set_ylabel("number / 1")
axes.set_xlabel("diameter / m")
#create legend        
axes.legend()

plt.tight_layout()
#save figure
dir_save = '../plots/'
out_filestring = "size_dist_" + distribution
plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
plt.savefig(dir_save + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
plt.clf()
plt.close()