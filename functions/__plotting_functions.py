
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
    
def plot1Dhistline(ax,x,counts,xlabel='',ylabel='counts',title='',linelabel='',logflag=0,xlims='None',ylims='None',linestyle='-',color='r'):
    '''
    input: x-> x-axes of bar-plot of the histogram, counts-> corresponding counts to x-axis, facecolor, 
    alpha-> transparency value, xlabel,ylabel,title, logflag-> 0: no 1: logx 2:logy 3:loglog
    '''
    
    # the histogram of the data
    if logflag==0:
        ax.plot(x,counts,label=linelabel,linestyle=linestyle,color=color)
    elif logflag==1:
        ax.semilogx(x,counts,label=linelabel,linestyle=linestyle,color=color)
    elif logflag==2:
        ax.semilogy(x,counts,label=linelabel,linestyle=linestyle,color=color)
    elif logflag==3:
        ax.loglog(x, counts,label=linelabel,linestyle=linestyle,color=color) #, np.diff(x) , facecolor=facecolor, alpha=alpha, log=logy)
    
    if not ylims=='None':
        ax.set_ylim([ylims[0],ylims[1]])
    if not xlims=='None':
        ax.set_xlim([xlims[0],xlims[1]])    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    return ax

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
    for i in range(0,len(pamData["frequency"])): #all available frequencies in the pamtra files are plotted
        if linestyle=='-': #just create labels for the first series of frequencies (which should have '-' as a linestyle)
            ax.plot(pamData["Ze"][:,i],pamData["height"],color=np.array(['b','r','g'])[i],linestyle=linestyle,marker=marker,markerfacecolor='None',markevery=20,label='{:5.1f}GHz'.format(pamData["frequency"][i]))
        else:
            ax.plot(pamData["Ze"][:,i],pamData["height"],color=np.array(['b','r','g'])[i],linestyle=linestyle,marker=marker,markerfacecolor='None',markevery=5)

    #set range
    ax.set_xlim([-40,55]) #range of Ze
    ax.set_ylim([0,pamData["height"][-1]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    ax.set_xlabel("reflectivity / dBz")
    ax.set_ylabel("height / m")

    return plt

def plot_pamtra_highermoments(ax,pamData,linestyle='-',marker=' '):
    '''
    plot pamtra output
    INPUT: pamData: dictionary with PAMTRA variables
    '''

    #create figure for plotting of pamtra output

    #plot X-,Ka-,W-Band in three different colors
    #from IPython.core.debugger import Tracer ; Tracer()()
    for i in range(0,len(pamData["frequency"])): #all available frequencies in the pamtra files are plotted
        if linestyle=='-': #just create labels for the first series of frequencies (which should have '-' as a linestyle)
            ax.plot(pamData["Radar_SpectrumWidth"][:,i],pamData["height"],color=np.array(['b','r','g'])[i],linestyle=linestyle,marker=marker,markerfacecolor='None',markevery=20,label='{:5.1f}GHz'.format(pamData["frequency"][i]))
        else:
            ax.plot(pamData["Radar_SpectrumWidth"][:,i],pamData["height"],color=np.array(['b','r','g'])[i],linestyle=linestyle,marker=marker,markerfacecolor='None',markevery=5)

    #set range
    #ax.set_xlim([-40,55]) #range of Ze
    ax.set_ylim([0,pamData["height"][-1]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    ax.set_xlabel("spectrum width / m s-1")
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

def plot_twomom_moments(ax,ax2,twomom,i_timestep,add_Dmean=True):
    '''
    plot moments (number density and mass density) over height
    INPUT:  ax: first x-axis
            ax2: second x-axis
            twomom: dictionary with (McSnows) two-moment variables
            add_Dmean (boolean): add also Dmean into plot
    '''

    import matplotlib.ticker
    #plot mixing ratios
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    markerlist = ["o","^","s","*","x","d"]
    
    axisq = [];axisDmean = []
    for i,cat in enumerate(['c','i','r','s','g','h']):
        axisq += ax.semilogx(twomom['q' + cat][i_timestep,:],twomom['heights'],color='r',marker=markerlist[i],markevery=20,markerfacecolor='None')
        if add_Dmean:
            axisDmean += ax2.semilogx(twomom['D_mean_' + cat][i_timestep,:],twomom['heights'],color='k',linestyle='',marker=markerlist[i],markevery=20,markerfacecolor='None')

    #plot number concentrations
    axisqn = []
    for i,cat in enumerate(['c','i','r','s','g','h']):
        axisqn += ax2.semilogx(twomom['qn' + cat][i_timestep,:],twomom['heights'],color='k',marker=markerlist[i],markevery=20,markerfacecolor='None')

    #change height limits
    ax.set_ylim([0,twomom["heights"][0]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    #add labels and legend
    ax.set_xlabel("mixing ratio / kg m-3", color="r")
    ax.tick_params(axis='x', colors='r')    
    ax2.set_xlabel("mean diameter (symbols only) / m     number concentration / m-3",color="k")
    ax.set_ylabel("height / m")
    #set a fix range of the x-axis
    xmin=-9;xmax=1
    x2min=-5;x2max=9
    ax.set_xlim([10**xmin,10**xmax]) #[ax.get_xlim()[0],ax.get_xlim()[1]*100])
    ax2.set_xlim([10**x2min,10**x2max]) #[ax2.get_xlim()[0],ax2.get_xlim()[1]*1000])
    #set xticks
    ax.set_xticks(np.logspace(xmin,xmax,xmax-xmin+1))
    ax2.set_xticks(np.logspace(x2min,x2max,x2max-x2min+1))
    #remove every second ticklabel
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    for label in ax2.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    
    #add legend
    ax.legend(axisqn,["cloud w.","cloud ice","rain","snow","graupel","hail"],loc='center right') #,loc='center left', bbox_to_anchor=(1, 0.5)) #position: the "center left" of the box is at (x=1,y=0.5, in relative coordinates of the whole plot)
    
    return plt

def plot_fluxes(ax,ax2,twomom,hei2massdens,i_timestep,mass_num_flag=2):
    '''
    plot number and mass fluxes of SB and McSnow over height
    INPUT:  ax: first x-axis
            ax2: second x-axis
            twomom: dictionary with (McSnows) two-moment variables
            hei2massdens: dictionary with McSnow variables
            i_timestep: timestep used (the dicts contain all timesteps of the simulation output; i_timestep is defined by governing scripts)
            mass_num_flag: boolean which determines if the mass or the number flux should be plotted
    '''
    import matplotlib.ticker
    #prop_cycle = plt.rcParams['axes.prop_cycle']
    #colors = prop_cycle.by_key()['color']
    #colors[-1] = 'red'
    colors=["blue","green","red"]
    markerlist = ["^","s","*","x","d",""]
    SB_specs = ['i','s','_all'] #specify here what categories should be used; all would be ['i','r','s','g','h']
    MC_specs = ["_mm1","_unr",""] #specify here what categories should be used; all would be ["_mm1","_unr","_grp","_liq","_rimed",""]#
    labellistSB = ["cloud ice","snow","SB all"] #check consistency with SB_specs; all would be ["cloud ice","rain","snow","SB graupel","hail"]
    labellistMC = ["pristine","unrimed agg.","MC all"] #check consistency with MC_specs; all would be ["pristine","unrimed agg.","MC graupel","liquid","rimed","MC all"]
    #sum up all Sb categories
    twomom['qn_all'] = np.zeros_like(twomom['qni']);    twomom['q_all'] = np.zeros_like(twomom['qni']);     twomom['fn_all'] = np.zeros_like(twomom['qni']);     twomom['f_all'] = np.zeros_like(twomom['qni'])
    for i,cat in enumerate(['i','s']): #['i','s','g','h','r']):#TODO: include other categories when analyzing rimed cases 
        twomom['qn_all']+= twomom['qn' + cat]
        twomom['q_all']+= twomom['q' + cat]
        twomom['fn_all']+= twomom['fn' + cat]
        twomom['f_all']+= twomom['f' + cat]
    McSNow_plot_only_every = 1 #plot not every height point of the noisy McSnow fluxes (for better visibility)
    if mass_num_flag==0:  #plot number flux
        axisqn = []
        print twomom.keys()
        for i,cat in enumerate(SB_specs):
            if any(twomom['fn' + cat][i_timestep,:])>0:
                axisqn += ax2.semilogx(twomom['fn' + cat][i_timestep,:],twomom['heights'],color=colors[i],label=labellistSB[i])
        for j,cat in enumerate(MC_specs): #choose "" as an entry to get all summed up
            if any(hei2massdens["Fn" + cat][:])>0:
                axisqn += ax2.semilogx(hei2massdens["Fn" + cat][:-1:McSNow_plot_only_every],hei2massdens['z'][:-1:McSNow_plot_only_every],color=colors[j],label=labellistMC[j],linestyle='--')
    if mass_num_flag==1:     #plot mass flux
        axisq = [];axisDmean = []
        if mass_num_flag==1:
            color='r'
        elif mass_num_flag==2:
            color='r'
        for i,cat in enumerate(SB_specs):
            if any(twomom['f' + cat][i_timestep,:])>0:
                axisq += ax.semilogx(twomom['f' + cat][i_timestep,:],twomom['heights'],color=colors[i],label=labellistSB[i])
        for j,cat in enumerate(MC_specs): #choose "" as an entry to get all summed up            
            if any(hei2massdens["Fm" + cat][:])>0:
                axisq += ax.semilogx(hei2massdens["Fm" + cat][:-1:McSNow_plot_only_every],hei2massdens['z'][:-1:McSNow_plot_only_every],color=colors[j],label=labellistMC[j],linestyle='--')    
    elif mass_num_flag==2:
        for i,cat in enumerate(SB_specs):
            if any(twomom['f' + cat][i_timestep,:])>0:
                axisq += ax.semilogx(twomom['f' + cat][i_timestep,:],twomom['heights'],color=color,marker=markerlist[i],markevery=20,markerfacecolor='None',label=labellistSB[i])

        for i,cat in enumerate(MC_specs): #choose "" as an entry to get all summed up            
            if any(hei2massdens["Fm" + cat][:])>0:
                axisq += ax2.semilogx(hei2massdens["Md" + cat][:],hei2massdens['z'],color='b',linestyle='--',marker=markerlist[i],markevery=20,markerfacecolor='None',label=labellistMC[i])
        for i,cat in enumerate(SB_specs):
            if any(twomom['fn' + cat][i_timestep,:])>0:
                axisqn += ax2.semilogx(twomom['fn' + cat][i_timestep,:],twomom['heights'],color='k',marker=markerlist[i],markevery=20,markerfacecolor='None',label=labellistSB[i])
        for i,cat in enumerate(MC_specs): #choose "" as an entry to get all summed up
            if any(hei2massdens["Fn" + cat][:])>0:
                axisqn += ax2.semilogx(hei2massdens["Fn" + cat][:],hei2massdens['z'],color='g',linestyle='--',marker=markerlist[i],markevery=20,markerfacecolor='None',label=labellistMC[i])
    #from IPython.core.debugger import Tracer ; Tracer()()
    #change height limits and set y-label
    ax.set_ylim([0,twomom["heights"][0]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    ax.set_ylabel("height / m")
    
    
    #add labels and legend
    mass_min=-8; mass_max=0
    num_min=2; num_max=6
    if mass_num_flag==2:
        ax.set_xlabel("mass flux density / kg m-2 s-1", color="r")
        ax.tick_params(axis='x', colors='r')    
        ax2.set_xlabel("number flux density / m-2 s-1",color="k")
        #set a fix range of the x-axis
        xmin=mass_min;xmax=mass_max
        x2min=num_min;x2max=num_max
    elif mass_num_flag==0:
        ax.set_xlabel("number flux density / m-2 s-1",color="k")
        xmin=num_min;xmax=num_max
    elif mass_num_flag==1:
        ax.set_xlabel("mass flux density / kg m-2 s-1",color="k")
        xmin=mass_min;xmax=mass_max
    else:
        print "error: mass_num_flag in plot_fluxes must be in [0,1,2]"


    #set a fix range of the x-axis
    ax.set_xlim([10**xmin,10**xmax]) #[ax.get_xlim()[0],ax.get_xlim()[1]*100])
    ax.set_xticks(np.logspace(xmin,xmax,xmax-xmin+1)) #set xticks
    #remove every second ticklabel
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    if mass_num_flag==2: #do the same for the second x-axis if existent
        ax2.set_xlim([10**x2min,10**x2max]) #[ax2.get_xlim()[0],ax2.get_xlim()[1]*1000])
        ax2.set_xticks(np.logspace(x2min,x2max,x2max-x2min+1))
        for label in ax2.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
    
    #add legend
    if mass_num_flag in [0,2]:
        ax.legend(loc='upper right') #,loc='center left', bbox_to_anchor=(1, 0.5)) #position: the "center left" of the box is at (x=1,y=0.5, in relative coordinates of the whole plot)
    else: #if there is no qn in the plot label q
        ax.legend(loc='upper right') #,loc='center left', bbox_to_anchor=(1, 0.5)) #position: the "center left" of the box is at (x=1,y=0.5, in relative coordinates of the whole plot)
    
    return plt

def plot_moments(ax,ax2,twomom,hei2massdens,i_timestep,mass_num_flag=2):
    '''
    plot number and mass concentration of SB and McSnow over height
    INPUT:  ax: first x-axis
            ax2: second x-axis
            twomom: dictionary with (McSnows) two-moment variables
            hei2massdens: dictionary with McSnow variables
            i_timestep: timestep used (the dicts contain all timesteps of the simulation output; i_timestep is defined by governing scripts)
            mass_num_flag: boolean which determines if the mass or the number flux should be plotted
    '''
    import matplotlib.ticker
    #prop_cycle = plt.rcParams['axes.prop_cycle']
    #colors = prop_cycle.by_key()['color']
    colors=["blue","green","red"]
    markerlist = ["^","s","*","x","d",""]
    SB_specs = ['i','s','_all'] #specify here what categories should be used; all would be ['i','r','s','g','h']
    MC_specs = ["_mm1","_unr",""] #specify here what categories should be used; all would be ["_mm1","_unr","_grp","_liq","_rimed",""]#
    labellistSB = ["cloud ice","snow","SB all"] #check consistency with SB_specs; all would be ["cloud ice","rain","snow","SB graupel","hail"]
    labellistMC = ["pristine","unrimed agg.","MC all"] #check consistency with MC_specs; all would be ["pristine","unrimed agg.","MC graupel","liquid","rimed","MC all"]
    #sum up all Sb categories
    twomom['qn_all'] = np.zeros_like(twomom['qni']);    twomom['q_all'] = np.zeros_like(twomom['qni']);     twomom['fn_all'] = np.zeros_like(twomom['qni']);     twomom['f_all'] = np.zeros_like(twomom['qni'])
    for i,cat in enumerate(['i','s']): #['i','s','g','h','r']):#TODO: include other categories when analyzing rimed cases 
        twomom['qn_all']+= twomom['qn' + cat]
        twomom['q_all']+= twomom['q' + cat]
        
    McSNow_plot_only_every = 4 #plot not every height point of the noisy McSnow fluxes (for better visibility)
    if mass_num_flag==0:  #plot number concentrations
        axisqn = []
        for i,cat in enumerate(SB_specs):
            if ('qn' + cat) in twomom.keys() and any(twomom['qn' + cat][i_timestep,:])>0:
                axisqn += ax2.semilogx(twomom['qn' + cat][i_timestep,:],twomom['heights'],color=colors[i],label=labellistSB[i])
                if ('qn' + cat + '_std') in twomom.keys(): #add +- standard deviation in plot if available
                    ax2.semilogx(twomom['qn' + cat][i_timestep,:]-twomom['qn' + cat + '_std'][i_timestep,:],twomom['heights'],color=colors[i],linestyle='--',label='_dontshowinlegend')
                    ax2.semilogx(twomom['qn' + cat][i_timestep,:]+twomom['qn' + cat + '_std'][i_timestep,:],twomom['heights'],color=colors[i],linestyle='--',label='_dontshowinlegend')
        for j,cat in enumerate(MC_specs): #choose "" as an entry to get all summed up
            if ("Nd" + cat) in hei2massdens.keys() and any(hei2massdens["Nd" + cat][:])>0:
                axisqn += ax2.semilogx(hei2massdens["Nd" + cat][::McSNow_plot_only_every],hei2massdens['z'][::McSNow_plot_only_every],color=colors[j],label=labellistMC[j],linestyle='--')
    if mass_num_flag==1:     #plot mixing ratios
        axisq = [];axisDmean = []
        #if mass_num_flag==1:
        #    color='r'
        #elif mass_num_flag==2:
        #    color='r'
        for i,cat in enumerate(SB_specs):
            if ("q" + cat) in twomom.keys() and any(twomom['q' + cat][i_timestep,:])>0:
                axisq += ax.semilogx(twomom['q' + cat][i_timestep,:],twomom['heights'],color=colors[i],label=labellistSB[i])
                if ('q' + cat + '_std') in twomom.keys(): #add +- standard deviation in plot if available
                    ax2.semilogx(twomom['q' + cat][i_timestep,:]-twomom['q' + cat + '_std'][i_timestep,:],twomom['heights'],color=colors[i],linestyle='--',label='_dontshowinlegend')
                    ax2.semilogx(twomom['q' + cat][i_timestep,:]+twomom['q' + cat + '_std'][i_timestep,:],twomom['heights'],color=colors[i],linestyle='--',label='_dontshowinlegend')
                    #print "q" + cat,twomom['q' + cat]
        for j,cat in enumerate(MC_specs): #choose "" as an entry to get all summed up            
            if ("Md" + cat) in hei2massdens.keys() and any(hei2massdens["Md" + cat][:])>0:
                axisq += ax.semilogx(hei2massdens["Md" + cat][::McSNow_plot_only_every],hei2massdens['z'][::McSNow_plot_only_every],color=colors[j],label=labellistMC[j],linestyle='--')
    elif mass_num_flag==2:
        for i,cat in enumerate(SB_specs):
            if any(twomom['f' + cat][i_timestep,:])>0:
                axisq += ax.semilogx(twomom['f' + cat][i_timestep,:],twomom['heights'],color=color,marker=markerlist[i],markevery=20,markerfacecolor='None',label=labellistSB[i])
        for i,cat in enumerate(MC_specs): #choose "" as an entry to get all summed up            
            if any(hei2massdens["Md" + cat][:])>0:
                axisq += ax2.semilogx(hei2massdens["Md" + cat][:],hei2massdens['z'],color='b',linestyle='--',marker=markerlist[i],markevery=20,markerfacecolor='None',label=labellistMC[i])
        for i,cat in enumerate(SB_specs):
            if any(twomom['fn' + cat][i_timestep,:])>0:
                axisqn += ax2.semilogx(twomom['fn' + cat][i_timestep,:],twomom['heights'],color='k',marker=markerlist[i],markevery=20,markerfacecolor='None',label=labellistSB[i])
        for i,cat in enumerate(MC_specs): #choose "" as an entry to get all summed up
            if any(hei2massdens["Nd" + cat][:])>0:
                axisqn += ax2.semilogx(hei2massdens["Nd" + cat][:],hei2massdens['z'],color='g',linestyle='--',marker=markerlist[i],markevery=20,markerfacecolor='None',label=labellistMC[i])
    #change height limits and set y-label
    ax.set_ylim([0,twomom["heights"][0]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    ax.set_ylabel("height / m")
    
    
    #add labels and legend
    mass_min=-8; mass_max=-1
    num_min=0; num_max=6
    if mass_num_flag==2:
        ax.set_xlabel("mass density / kg m-3", color="r")
        ax.tick_params(axis='x', colors='r')    
        ax2.set_xlabel("number density / m-3",color="k")
        #set a fix range of the x-axis
        xmin=mass_min;xmax=mass_max
        x2min=num_min;x2max=num_max
    elif mass_num_flag==0:
        ax.set_xlabel("number density / m-3",color="k")
        xmin=num_min;xmax=num_max
    elif mass_num_flag==1:
        ax.set_xlabel("mass density / kg m-3",color="k")
        xmin=mass_min;xmax=mass_max
    else:
        print "error: mass_num_flag in plot_fluxes must be in [0,1,2]"



    #set a fix range of the x-axis
    ax.set_xlim([10**xmin,10**xmax]) #[ax.get_xlim()[0],ax.get_xlim()[1]*100])
    ax.set_xticks(np.logspace(xmin,xmax,xmax-xmin+1)) #set xticks
    #remove every second ticklabel
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    if mass_num_flag==2: #do the same for the second x-axis if existent
        ax2.set_xlim([10**x2min,10**x2max]) #[ax2.get_xlim()[0],ax2.get_xlim()[1]*1000])
        ax2.set_xticks(np.logspace(x2min,x2max,x2max-x2min+1))
        for label in ax2.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
    
    labels = ax2.get_xticklabels()
    #add legend
    ax.legend(loc='upper right') #,loc='center left', bbox_to_anchor=(1, 0.5)) #position: the "center left" of the box is at (x=1,y=0.5, in relative coordinates of the whole plot)
    
    return plt

def plot_normmix(ax,ax2,twomom,hei2massdens,i_timestep):
    '''
    plot normalized mixing ratio/ mean particle mass (ratio of mass density and number density) of SB and McSnow over height
    INPUT:  ax: first x-axis
            ax2: second x-axis
            twomom: dictionary with (McSnows) two-moment variables
            hei2massdens: dictionary with McSnow variables
            i_timestep: timestep used (the dicts contain all timesteps of the simulation output; i_timestep is defined by governing scripts)
    '''
    import matplotlib.ticker
    #prop_cycle = plt.rcParams['axes.prop_cycle']
    #colors = prop_cycle.by_key()['color']
    #colors[-1] = 'red'
    colors=["blue","green","red"]
    markerlist = ["^","s","*","x","d",""]
    SB_specs = ['i','s','_all'] #specify here what categories should be used; all would be ['i','r','s','g','h']
    MC_specs = ["_mm1","_unr",""] #specify here what categories should be used; all would be ["_mm1","_unr","_grp","_liq","_rimed",""]#
    labellistSB = ["cloud ice","snow","SB all"] #check consistency with SB_specs; all would be ["cloud ice","rain","snow","SB graupel","hail"]
    labellistMC = ["pristine","unrimed agg.","MC all"] #check consistency with MC_specs; all would be ["pristine","unrimed agg.","MC graupel","liquid","rimed","MC all"]
    #sum up all SB categories
    twomom['qn_all'] = np.zeros_like(twomom['qni']); twomom['q_all'] = np.zeros_like(twomom['qni']);
    for i,cat in enumerate(['i','s']): #['i','s','g','h','r']):#TODO: include other categories when analyzing rimed cases 
        twomom['qn_all']+= twomom['qn' + cat]
        twomom['q_all']+= twomom['q' + cat]
    #calculate normalized mixing ratios
    #first for the sum of all categories
    twomom['nq_all'] = twomom['q_all']/twomom['qn_all']
    hei2massdens['nq_all'] = hei2massdens['Md']/hei2massdens['Nd']
    #second for all categories individually
    for i,cat in enumerate(SB_specs):
        twomom['nq' + cat] = twomom['q' + cat]/twomom['qn' + cat]
    for j,cat in enumerate(MC_specs):
        hei2massdens['nq' + cat] = hei2massdens["Md" + cat]/hei2massdens["Nd" + cat]
    #from IPython.core.debugger import Tracer ; Tracer()()

    axisnq = [];axisDmean = [] #initialize handles
    McSNow_plot_only_every = 1 #plot not every height point of the noisy McSnow fluxes (for better visibility)
    for i,cat in enumerate(SB_specs): #['i','s','g','h','r','_all']):#TODO: include other categories when analyzing rimed cases 
        if ('nq' + cat) in twomom.keys() and any(twomom['nq' + cat][i_timestep,:])>0:
            axisnq += ax2.semilogx(twomom['nq' + cat][i_timestep,:],twomom['heights'],color=colors[i],label=labellistSB[i])
    for j,cat in enumerate(MC_specs):  #["_mm1","_unr","_grp","_liq","_rimed",""]): #choose "" as an entry to get all summed up #TODO: include other categories when analyzing rimed cases 
        if ("nq" + cat) in hei2massdens.keys() and any(hei2massdens["nq" + cat][:])>0:
            axisnq += ax2.semilogx(hei2massdens["nq" + cat][::McSNow_plot_only_every],hei2massdens['z'][::McSNow_plot_only_every],color=colors[j],label=labellistMC[j],linestyle='--')

    #change height limits and set y-label
    ax.set_ylim([0,twomom["heights"][0]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    ax.set_ylabel("height / m")
    
    
    #add labels and legend
    nmass_min=-10; nmass_max=-4
    ax.set_xlabel("mean particle mass / kg") #"normalized mixing ratio / kg",color="k")


    #set a fix range of the x-axis
    ax.set_xlim([10**nmass_min,10**nmass_max]) #[ax.get_xlim()[0],ax.get_xlim()[1]*100])
    ax.set_xticks(np.logspace(nmass_min,nmass_max,nmass_max-nmass_min+1)) #set xticks
    '''
    #remove every second ticklabel
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    if mass_num_flag==2: #do the same for the second x-axis if existent
        ax2.set_xlim([10**x2min,10**x2max]) #[ax2.get_xlim()[0],ax2.get_xlim()[1]*1000])
        ax2.set_xticks(np.logspace(x2min,x2max,x2max-x2min+1))
        for label in ax2.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
    '''
    #labels = ax2.get_xticklabels()
    #add legend
    ax.legend(loc='upper right') #,loc='center left', bbox_to_anchor=(1, 0.5)) #position: the "center left" of the box is at (x=1,y=0.5, in relative coordinates of the whole plot)

    return plt

def plot_atmo(ax,ax2,atmo):
    '''
    plot atmospheric variables in one panel
    INPUT:  ax: first x-axis
            ax2: second x-axis
            atmo: dictionary with atmospheric variables
    '''
    #plot T, RHw and RHi in one panel
    handles = []
    if 'T' in atmo.keys():
        handles += ax.plot(atmo['T']-273.15,atmo['z'],color='r',label='T')
        ax.set_xlabel("temperature / $^\circ$C")
        if 'T_std' in atmo.keys():
            handles += ax.plot(atmo['T']-273.15-atmo['T_std'],atmo['z'],color='r',linestyle='--',label='__none')
            handles += ax.plot(atmo['T']-273.15-atmo['T_std'],atmo['z'],color='r',linestyle='--',label='__none')
    if 'rh' in atmo.keys():
        handles += ax2.plot(atmo["rh"],atmo['z'],color='b',label='RHw')
        ax2.plot(np.ones(atmo['z'].shape)*100,atmo['z'],color='grey',linestyle='--',label='_dontshowinlegend')
        ax2.set_xlabel("relative humidity / %")
        if 'rh_std' in atmo.keys():
            handles += ax2.plot(atmo["rh"]+atmo["rh_std"],atmo['z'],color='b',linestyle='--',label='__none')
            handles += ax2.plot(atmo["rh"]-atmo["rh_std"],atmo['z'],color='b',linestyle='--',label='__none')
    #plot rhi (if rh,psatw and psati present calculate it from them) elif rhi is there take it directly
    #if all (key in atmo.keys() for key in ("rh","psatw","psati")):
    #    handles += ax2.plot(atmo["rh"]*atmo["psatw"]/atmo["psati"],atmo['z'],color='b',linestyle='--',label='RHi')
    if 'rhi' in atmo.keys():
        handles += ax2.plot(atmo["rhi"],atmo['z'],color='g',linestyle='-',label='RHi')
        if 'rhi_std' in atmo.keys():
            handles += ax2.plot(atmo["rhi"]+atmo["rhi_std"],atmo['z'],color='g',linestyle='--',label='__none')
            handles += ax2.plot(atmo["rhi"]-atmo["rhi_std"],atmo['z'],color='g',linestyle='--',label='__none')
    #change height limits
    ax.set_ylim([0,np.nanmax(atmo['z'])])#atmo['z'][-1]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    #add labels and legend
    ax.set_ylabel("height / m")
    #create space for legend
    ax.set_xlim([ax.get_xlim()[0],ax.get_xlim()[1]+0])
    low_RH=50 #define lower end of RH by hand
    ax2.set_xlim([low_RH,ax2.get_xlim()[1]+15])
    # added these three lines
    labs = [h.get_label() for h in handles]
    ax.legend(handles,labs,loc='upper right')
    #from IPython.core.debugger import Tracer ; Tracer()()
    return plt
    
def pcolor_timeseries(fig,ax,time,height,var,varname='',time_formatted='None',lin_log=0,unit=''):
    '''
    plot timeseries (written for plotting ICON-meteogram)
    INPUT:  time-> timevector
            height-> heightvector
            var -> variable which should be plotted
    '''
    import matplotlib.colors as colors
    
    #####
    #variable specific settings
    #####
    #default settings #these are overwritten after "specifications for variable groups or individual variables"
    lin_log=0 #lin is default
    cmap = "viridis" #take viridis as default colorbar
    #if lin_log==0:
    varmin=var.min()
    varmax=var.max()
    #elif lin_log==1:
    #    varmin=max(var.min(),1e-100) #ATTENTION:max(..,1e-100) avoids only the crash of the script; define lower limit for variables for which zeros occur below
    #    varmax=var.max()
        
    ###
    #specifications for variable groups or individual variables
    ###
    if varname in ("qc","qr","qi","qs","qg","qh"): #these are logscaled
        lin_log = 1 #set y-axis to logscale
        #set limits
        varmin = 1e-8
        varmax = 1e-2
    if varname in ("qnc","qnr"): #these are logscaled
        lin_log = 1 #set y-axis to logscale
        #set limits
        varmin = 1e0
        varmax = 1e12
    if varname in ("qni","qns","qng","qnh"):#these are logscaled
        lin_log = 1 #set y-axis to logscale
        #set limits
        varmin = 1e0
        varmax = 1e6
    if varname=="rhi":
        #set limits
        varmin = 80
        varmax = 120
        cmap = "bwr"
    #normalize according to lin_log
    if lin_log==0:
        norm = colors.Normalize(vmin=varmin,vmax=varmax)
    elif lin_log==1:
        norm = colors.LogNorm(vmin=varmin,vmax=varmax)
    
    #plot
    ax_now = ax.pcolormesh(time,height,var,rasterized=True,norm=norm,cmap=cmap)#,vmin=mincol,vmax=maxcol,rasterized=True,norm=norm,cmap=new_cmap)
    #create colorbar
    cbar = plt.colorbar(ax_now,ax=ax) #, format='%d')
    cbar.ax.set_ylabel(varname + " / " + unit)
    #ticklabels
    #derive location of ticks dynamically depending on timeranges
    if np.sum((time_formatted.minute==0) & (time_formatted.second==0))>15: #label every 2 hours for long timeranges (f.e. 15 hours); specify hours by the last number
        ticklocations = np.where((time_formatted.minute==0) & (time_formatted.second==0) & ((time_formatted.hour%2)==0))[0]
    elif 2<np.sum((time_formatted.minute==0) & (time_formatted.second==0))<15: #label every hour for delta_t_hour = [2,15]
        ticklocations = np.where((time_formatted.minute==0) & (time_formatted.second==0))[0]
    else:
        print "not implemented x-axis labelling for short timeranges"
        sys.exit(1)
    #set location of ticks
    ax.set_xticks(time[ticklocations])
    #generate tick labels and set them
    timestrings = [ "{:02d}:{:02d}".format(time_formatted[i].hour,time_formatted[i].second) for i in ticklocations]
    ax.set_xticklabels(timestrings) #"{:2.0f}:{:2.0f}".format(time[ticklocations].hour.values,time[ticklocations].minute.values))
    #labelling
    ax.set_xlabel("time / UTC")
    ax.set_ylabel("height / m")

    return ax