
# coding: utf-8
#import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#import other self-defined functions
import __postprocess_McSnow

#from IPython.core.debugger import Tracer ; Tracer()()


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
def plot_pamtra_Ze(ax,pamData):
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
    ax.set_xlim([-20,50]) #range of Ze
    #TODO: set height flexible
    ax.set_ylim([0,pamData["height"][-1]])
    ax.set_xlabel("reflectivity / dBz")
    ax.set_ylabel("height / m")
    
    #add legend
    ax.legend(["9.6GHz","35.5GHz","95GHz"])

    return plt

def plot_pamtra_spectrogram(ax,pamData,freq=35.5,cmap='viridis_r'):
    '''
    plot sprectogram from pamtra output
    INPUT: pamData: dictionary with PAMTRA variables
    freq: frequency, from which the spectrum should be used
    cmap: name of colormap
    '''

    #choose spectrum of right frequency
    freqindex = np.where(abs(pamData["frequency"]-freq)<0.01)[0] #we need some tolerance because the frequency is not always saved correctly
    if freqindex.size==0: #give error if frequency is not in file
        print "frequency: ", freq ,"is not in the PAMTRA output file"
        print "choose between: ",pamData["frequency"]," or change frequencies in run_pamtra during runtime"
        sys.exit(0)

    #get axis and spectrogram data from pamData
    min_shown = -30 #minimum value shown in the diagram
    Radar_Spectrum = np.squeeze(pamData["Radar_Spectrum"][:,freqindex,:]) #dimensions [height, nfft]
    Radar_Spectrum_masked = np.ma.array(Radar_Spectrum,mask=Radar_Spectrum<min_shown) #do not show small values
    Radar_Velocity = np.squeeze(pamData["Radar_Velocity"][freqindex,:])  #dimension [nfft]
    height = np.squeeze(pamData["height"])
    
    #plot
    pcol = ax.pcolormesh(-Radar_Velocity,height,Radar_Spectrum_masked,rasterized=True,cmap=cmap,vmin=min_shown) #norm=norm)

    #set range
    ax.set_xlim([-3,1]) #range of Ze

    #plot labels and create colorbar
    ax.set_xlabel("Doppler velocity / m s-1")
    ax.set_ylabel("height / m")
    col = plt.colorbar(pcol)
    col.set_label("spectral power / dB", rotation=90)
    #plot label of used frequency
    plt.text(0.4, 4500, str(freq) + 'GHz', fontsize=12)
    
    return plt

def plot_McSnows_vt_in_spectrogram(ax,pamData,SP,experiment,cmap='viridis_r'):
    '''
    plot number of RP in spectrogram bins
    INPUT: pamData: dictionary with PAMTRA variables
    SP dictionary with McSnows variables for each SP
    experiment: string describing the experiment (needed for calculating Vbox)
    cmap: name of colormap
    '''

    #get axis and spectrogram data from pamData
    Radar_Velocity = np.squeeze(pamData["Radar_Velocity"][0,:])  #dimension [nfft]
    height = np.squeeze(pamData["height"])
    binned_val=dict()
    binned_val["counts_spectro"] = np.zeros([len(height)-1,len(Radar_Velocity)]) #multiplicity weighted fall speed
    #get the number of particles (SP*sp_multiplicity) at each bin (binned by height and sp_diameter)
    for i in range(0,len(height)-1):
        for j in range(0,len(Radar_Velocity)-1):
                    condition_in_bin = np.logical_and(
                                    np.logical_and(Radar_Velocity[j]<=SP["vt"],SP["vt"]<Radar_Velocity[j+1]),
                                    np.logical_and(height[i]<=SP["height"],SP["height"]<height[i+1]),
                                    )
                    binned_val["counts_spectro"][i,j] = np.sum(np.where(condition_in_bin,SP["xi"],0))
    
    #calculate volume of box to derive number density later
    zres = height[1]-height[0]
    Vbox = __postprocess_McSnow.calculate_Vbox(experiment,zres)
    numdens_spectro=binned_val["counts_spectro"]/Vbox
    #plot
    masked_counts_spectr = np.ma.array(numdens_spectro,mask=numdens_spectro<=0)
    pcol = ax.pcolormesh(-Radar_Velocity,height,masked_counts_spectr,rasterized=True,cmap=cmap,norm=colors.LogNorm(vmin=10,vmax=10**5)) #norm=norm)
    #from IPython.core.debugger import Tracer ; Tracer()()

    #set range
    ax.set_xlim([-3,1]) #range of Ze

    #plot labels and create colorbar
    ax.set_xlabel("fall speed / m s-1")
    ax.set_ylabel("height / m")
    col = plt.colorbar(pcol)
    col.set_label("number density / m-3", rotation=90)
    return plt