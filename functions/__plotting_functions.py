
# coding: utf-8
#import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# define some plotting functions

def plot1Dhist(x,bins,normed=1,facecolor='blue',alpha=1, 
               xlabel='',ylabel='counts',title='',logflag=0,axlims=[0]):
    '''
    input: x->main variable of the histogram, bins, normed-> normalize flag, number_of_bins, facecolor, 
    alpha-> transparency value, xlabel,ylabel,title, logflag-> 0: no 1: logx 2:logy 3:loglog
    '''
    if logflag==1 or logflag==3:
        plt.xscale('log')
    if logflag==2 or logflag==3:
            logy=True
    else: logy=False
    
    # the histogram of the data
    n, bins, patches = plt.hist(x, bins=bins, normed=normed, facecolor=facecolor, alpha=alpha, log=logy)
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.axis([10**-5, 1, 1, 10**6])
    if len(axlims)==4:
        plt.axis(axlims)
    plt.grid(True)

def plotbar(x,counts,facecolor='blue',alpha=1, 
               xlabel='',ylabel='counts',title='',logflag=0,axlims=[0]):
    '''
    input: x-> x-axes of bar-plot of the histogram, counts-> corresponding counts to x-axis, facecolor, 
    alpha-> transparency value, xlabel,ylabel,title, logflag-> 0: no 1: logx 2:logy 3:loglog
    '''
    if logflag==1 or logflag==3:
        plt.xscale('log')
    if logflag==2 or logflag==3:
            logy=True
    else: logy=False
    
    # the histogram of the data
    plt.bar(x[:-1], counts[:-1], np.diff(x) , facecolor=facecolor, alpha=alpha, log=logy)
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if len(axlims)==4:
        plt.axis(axlims)
    plt.grid(True)
    plt.rcParams.update({'font.size': 20})

def pcol_height_diam(ax,diam,heights,var,mincol=-999,maxcol=-999,logcol=0,cmap='viridis_r',collabel='',ticks=np.array([-999])):
    '''
    create a 2D-histogram with axis diameter & height : min&max, ticks and label of colorbar can be specified
    INPUT:  ax: axis handle (subplots or similar should be organized by calling script)
            diam: array of diameters
            height: array of heights
            var: variable which should be plotted now
            
            mincol,maxcol: minimum and maximum colors of the colorbar  (can be done automatically)
            logcol: 0-> linear colorbar 1-> logarithmic colorbar
            cmap: define colormap  (can be done automatically)
            collabel: label for colorbar
            ticks: specifiy ticks of colorbar (can be done automatically)
    '''

    if mincol==-999: #set range of colorbar automatically
        pcol = ax.pcolormesh(diam,heights,var,rasterized=True,cmap=cmap) #norm=norm)
        if logcol==0:
            pcol = ax.pcolormesh(diam,heights,var,rasterized=True,cmap=cmap) #,norm=colors.LogNorm(vmin=mincol,vmax=maxcol))
        elif logcol==1:
            if mincol<=0: print "Error: logarithmic colorbar and non-positive and mincol<= 0 is not possible"
            pcol = ax.pcolormesh(diam,heights,var,rasterized=True,cmap=cmap,norm=colors.LogNorm())
    else: #set range of colorbar via mincol, maxcol
        if logcol==0:
            pcol = ax.pcolormesh(diam,heights,var,rasterized=True,cmap=cmap,vmin=mincol,vmax=maxcol) #,norm=colors.LogNorm(vmin=mincol,vmax=maxcol))
        elif logcol==1:
            if mincol<=0: print "Error: logarithmic colorbar and non-positive and mincol<= 0 is not possible"
            pcol = ax.pcolormesh(diam,heights,var,rasterized=True,cmap=cmap,norm=colors.LogNorm(vmin=mincol,vmax=maxcol)) #,norm=colors.LogNorm(vmin=mincol,vmax=maxcol))
        
    #modify axis
    ax.set_xscale("log")
    #TODO: give this as an input
    ax.set_xlim([10**-5,10**-1])
    #plot colorbar
    
    
    if ticks[0]==-999: #automatic ticks
        col = plt.colorbar(pcol)
    else:#predefined ticks
        col = plt.colorbar(pcol,ticks=ticks)

    col.set_label(collabel, rotation=90)
    #set labels
    ax.set_xlabel('diameter / m ')
    ax.set_ylabel('height / m ')
    return plt
def plot_pamtra_out(ax,pamData):
    '''
    plot pamtra output
    INPUT: pamData: dictionary with PAMTRA variables
    '''

    #create figure for plotting of pamtra output

    #plot X-,Ka-,W-Band in three different colors
    ax.plot(pamData["Ze"][:,0],pamData["height"], 'b-')
    ax.plot(pamData["Ze"][:,1],pamData["height"], 'r-')
    ax.plot(pamData["Ze"][:,2],pamData["height"], 'g-')

    #set range
    ax.set_xlim([-20,30]) #range of Ze
    #TODO: set height flexible
    ax.set_ylim([0,5000])
    ax.set_xlabel("reflectivity / dBz")
    ax.set_ylabel("height / m")
    
    #add legend
    ax.legend(["9.6GHz","35.5GHz","95GHz"])

    return plt