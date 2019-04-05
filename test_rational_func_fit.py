# coding: utf-8
#import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.pylab as pylab
import csv #to read the txt files
import os
import sys
import glob #to get filename to list
import subprocess
import itertools#to read only certain lines of the txt files
from matplotlib.colors import LogNorm
#import other self-defined functions
import __tools_for_processing_Jagg
import __plotting_functions
#from IPython.core.debugger import Tracer ; Tracer()()

'''
this code tests the rational function fit
'''

method = "rational_powerlaw" # "polynom","rational""polynom1_powerlaw" "rational_powerlaw"

#generate some x- and y-data
x1_orig = np.linspace(0,4,70)
x2_orig = np.linspace(0,3,70)
#setup x1 and x2 so that they contain all combinations of x1 vs x2
x1 = np.meshgrid(x1_orig,x2_orig)[0].flatten()
x2 = np.meshgrid(x1_orig,x2_orig)[1].flatten()


#spanned grid from D-N
x1_x2_grid =  np.meshgrid(x1_orig,x2_orig)
y = (1.-(x1*x2)**0.3)/(1.+0.01*x1)
y_noisy = y + np.random.normal(0,0.1,y.shape[0])

#remove parts with small

#keep_only = np.logical_and(0.5<x2/x1,x2/x1<1.5)
keep_only = (x1>=0)
x1 = x1[keep_only]
x2 = x2[keep_only]
y_noisy = y_noisy[keep_only]



#perform fit
if method=="polynom":
    fitcoeffs = __tools_for_processing_Jagg.fit_2D_rational_function(         [x1,x2],y_noisy,func='polynom1_1_pade_fixp00')
elif method=="rational":
    fitcoeffs = __tools_for_processing_Jagg.fit_2D_rational_function(         [x1,x2],y_noisy,func='rational1_1_pade_fixp00_and_fixq00')
elif method=="polynom1_powerlaw":
    fitcoeffs = __tools_for_processing_Jagg.fit_2D_rational_powerlaw_function([x1,x2],y_noisy,func='polynom1_powerlaw_fixp00_and_fixq00')
elif method=="rational_powerlaw":
    fitcoeffs = __tools_for_processing_Jagg.fit_2D_rational_powerlaw_function([x1,x2],y_noisy,func='rational1_1_powerlaw_fixp00_and_fixq00')
#from IPython.core.debugger import Tracer ; Tracer()()

#dictionary for fits as in processing_Jussis_...
fit_dic = dict()
fit_dic["coeffs"] = fitcoeffs

if method=="polynom":
    prop_fitted = __tools_for_processing_Jagg.polynom1_1_pade_fixp00(x1_x2_grid,fit_dic["coeffs"][0],fit_dic["coeffs"][1],fit_dic["coeffs"][2])
elif method=="rational":
    prop_fitted = __tools_for_processing_Jagg.rational1_1_pade_fixp00_and_fixq00(x1_x2_grid,fit_dic["coeffs"][0],fit_dic["coeffs"][1],fit_dic["coeffs"][2],fit_dic["coeffs"][3],fit_dic["coeffs"][4],fit_dic["coeffs"][5])
elif method=="polynom1_powerlaw":
    prop_fitted = __tools_for_processing_Jagg.polynom1_powerlaw_fixp00_and_fixq00(x1_x2_grid,fit_dic["coeffs"][0],fit_dic["coeffs"][1],fit_dic["coeffs"][2],fit_dic["coeffs"][3])    
elif method=="rational_powerlaw":
    prop_fitted = __tools_for_processing_Jagg.rational1_1_powerlaw_fixp00_and_fixq00(x1_x2_grid,fit_dic["coeffs"][0],fit_dic["coeffs"][1],fit_dic["coeffs"][2],fit_dic["coeffs"][3],fit_dic["coeffs"][4],fit_dic["coeffs"][5],fit_dic["coeffs"][6],fit_dic["coeffs"][7],fit_dic["coeffs"][8],fit_dic["coeffs"][9],fit_dic["coeffs"][10],fit_dic["coeffs"][11]) #,fit_dic["coeffs"][12])


#self-defined discrete colormap
# define the colormap and the discrete boundaries (defined by the considered monomer numbers) and set up an array of the same colors (for line-plotting)
cmap = plt.cm.viridis #brg
#modify lowest color in colormap
cmaplist = [cmap(i) for i in range(cmap.N)]
#cmaplist = cmaplist[80:] #crop colormap to avoid to light scatterer
cmaplist[0] = (1.0,0.0,0.0,1.0)
cmap = colors.LinearSegmentedColormap.from_list('mcm',cmaplist, cmap.N)
bounds = x2_orig
norm = colors.BoundaryNorm(bounds, cmap.N)
usecolors = cmap(np.linspace(0,1,x2_orig.shape[0])) #get colors from colormap to make consistent line colors for the KDE-lines #pylab.cm.Greens(np.linspace(0,1,N_mono_list.shape[0]))

#plot the raw_data+fit:
ax1=plt.subplot(1, 1, 1)
raw_plot = ax1.scatter(x1,y_noisy,c=x2,s=1,norm=norm,cmap=cmap)
plt.colorbar(raw_plot)

#from IPython.core.debugger import Tracer ; Tracer()()
for i_x2 in range(0,x2_orig.shape[0]):
    if not i_x2%10==0:
        continue
    #from IPython.core.debugger import Tracer ; Tracer()()
    fitted_plot = ax1.plot(x1_x2_grid[0][i_x2],prop_fitted[i_x2],linestyle='--',c=usecolors[i_x2]) #,color=raw_plot[0].get_color())
    #fitted_plot = ax1.plot(x1,y_noisy,linestyle='--',c=usecolors[i_x2]) #,color=raw_plot[0].get_color())
plt.show()



