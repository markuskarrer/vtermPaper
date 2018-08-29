
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
    ax.set_ylim([0,heights[-1]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
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
def plot_pamtra_Ze(ax,pamData,linestyle='-',marker=' '):
    '''
    plot pamtra output
    INPUT: pamData: dictionary with PAMTRA variables
    '''

    #create figure for plotting of pamtra output

    #plot X-,Ka-,W-Band in three different colors
    #from IPython.core.debugger import Tracer ; Tracer()()
    for i in range(0,len(pamData["frequency"])):
        if linestyle=='-': #just create labels for the first series of frequencies (which should have '-' as a linestyle)
            ax.plot(pamData["Ze"][:,i],pamData["height"],color=np.array(['b','r','g'])[i],linestyle=linestyle,marker=marker,markerfacecolor='None',markevery=5,label='{:5.1f}GHz'.format(pamData["frequency"][i]))
        else:
            ax.plot(pamData["Ze"][:,i],pamData["height"],color=np.array(['b','r','g'])[i],linestyle=linestyle,marker=marker,markerfacecolor='None',markevery=5)

    #set range
    ax.set_xlim([-20,55]) #range of Ze
    #TODO: set height flexible
    ax.set_ylim([0,pamData["height"][-1]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    ax.set_xlabel("reflectivity / dBz")
    ax.set_ylabel("height / m")

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
    ax.set_ylim([0,pamData["height"][-1]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
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
    
    #change height limits and ticks
    ax.set_ylim([0,height[-1]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))

    #plot labels and create colorbar
    ax.set_xlabel("fall speed / m s-1")
    ax.set_ylabel("height / m")
    col = plt.colorbar(pcol)
    col.set_label("number density / m-3", rotation=90)
    return plt

def plot_twomom_moments(ax,ax2,twomom,i_timestep):
    '''
    plot moments (number density and mass density) over height
    INPUT:  ax: first x-axis
            ax2: second x-axis
            twomom: dictionary with (McSnows) two-moment variables
    '''
    #plot mixing ratios
    #from IPython.core.debugger import Tracer ; Tracer()()
    plotqc = ax.semilogx(twomom['qc'][i_timestep,:],twomom['heights'])
    plotqr = ax.semilogx(twomom['qr'][i_timestep,:],twomom['heights'])

    plotqi = ax.semilogx(twomom['qi'][i_timestep,:],twomom['heights'])
    plotqs = ax.semilogx(twomom['qs'][i_timestep,:],twomom['heights'])
    plotqg = ax.semilogx(twomom['qg'][i_timestep,:],twomom['heights'])
    plotqh = ax.semilogx(twomom['qh'][i_timestep,:],twomom['heights'])

    #plot number concentrations
    plotqnc = ax2.semilogx(twomom['qnc'][i_timestep,:],twomom['heights'],'--')
    plotqnr = ax2.semilogx(twomom['qnr'][i_timestep,:],twomom['heights'],'--')

    plotqni = ax2.semilogx(twomom['qni'][i_timestep,:],twomom['heights'],'--')
    plotqns = ax2.semilogx(twomom['qns'][i_timestep,:],twomom['heights'],'--')
    plotqng = ax2.semilogx(twomom['qng'][i_timestep,:],twomom['heights'],'--')
    plotqnh = ax2.semilogx(twomom['qnh'][i_timestep,:],twomom['heights'],'--')
    
    #change height limits
    ax.set_ylim([0,twomom["heights"][0]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    #add labels and legend
    ax.set_xlabel("mixing ratio (solid) / kg m-3")
    ax2.set_xlabel("number concentration (dashed) / m-3")
    ax.set_ylabel("height / m")
    ax.set_xlim([ax.get_xlim()[0],ax.get_xlim()[1]*100])
    ax2.set_xlim([ax2.get_xlim()[0],ax2.get_xlim()[1]*1000])

    ax.legend(["cloud w.","rain","cloud ice","snow","graupel","hail"],loc='center right') #,loc='center left', bbox_to_anchor=(1, 0.5)) #position: the "center left" of the box is at (x=1,y=0.5, in relative coordinates of the whole plot)
    
    return plt
def plot_atmo(ax,ax2,atmo):
    '''
    plot atmospheric variables in one panel
    INPUT:  ax: first x-axis
            ax2: second x-axis
            atmo: dictionary with atmospheric variables
    '''
    #plot T, RHw and RHi in one panel
    plottemp = ax.plot(atmo['T']-273.15,atmo['z'],color='r')
    plotRHw = ax2.plot(atmo["rh"],atmo['z'],color='b')
    plotRHi = ax2.plot(atmo["rh"]*atmo["psatw"]/atmo["psati"],atmo['z'],color='b',linestyle='--')
    #change height limits
    ax.set_ylim([0,atmo['z'][-1]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    #add labels and legend
    ax.set_xlabel("temperature / $^\circ$C")
    ax2.set_xlabel("relative humidity / %")
    ax.set_ylabel("height / m")
    #create space for legend
    ax.set_xlim([ax.get_xlim()[0],ax.get_xlim()[1]+0])
    ax2.set_xlim([ax2.get_xlim()[0],ax2.get_xlim()[1]+15])
    # added these three lines
    handles = plottemp+plotRHw+plotRHi
    ax.legend(handles,["T","RHw","RHi"],loc='upper right')
    #from IPython.core.debugger import Tracer ; Tracer()()
    return plt