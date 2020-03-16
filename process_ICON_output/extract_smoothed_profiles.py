# -*- coding: utf-8 -*-
"""
read profiles from ICON simulation and derive smoothed profiles for semi-idealized runs
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
plt.close('all')
import matplotlib.dates as md
import numpy as np
import matplotlib.colors as colors
import pandas as pd
from scipy.special import gamma
from IPython.core.debugger import Tracer; debug = Tracer()
plt.close('all')
import sys
sys.path.append("../functions/")
import __general_utilities
# Define some paths and filenames

#date='20190122'
date='20181223'

for date in ['20181101','20181110','20181111','20181112','20181124','20181126','20181128','20181129','20181201','20181202','20181205','20181207','20181219','20181221','20181223','20190105','20190110','20190113','20190116','20190122']:
#for date in ['','','','','','','','','']:
    runsFolder = '/data/inscape/icon/experiments/juelich/testbed/testbed_' + date + '/' #/data/optimice/pamtra_runs/'
    runFld = runsFolder + '/'

    meteogramName = 'METEOGRAM_patch001_' + date + '_joyce.nc' #'METEOGRAM_patch001_joyce.nc'
    meteogramFile = runFld + meteogramName

    # Define Plotting Function
    def plot_variable(x,y,v,axes,
                      xlab=None,ylab=None,vlab=None,title=None,
                      vmin=None,vmax=None,xlim=None,ylim=None,
                      loglin="log",cmap='jet'):
        if loglin=="log":
            mesh = axes.pcolormesh(x,y,v,norm=colors.LogNorm(vmin=vmin, vmax=vmax),cmap=cmap)
        elif loglin=="lin":
            mesh = axes.pcolormesh(x,y,v,vmin=vmin,vmax=vmax,cmap=cmap)
        if title is not None:
            axes.text(0.1,0.9,title,transform=axes.transAxes,weight='black',
                      bbox=dict(facecolor='white'))
        if vlab is not None:
            plt.colorbar(mesh,label=vlab,ax=axes)
        else:
            plt.colorbar(mesh,ax=axes)
        if xlab is not None:
            axes.set_xlabel(xlab)
        if ylab is not None:
            axes.set_ylabel(ylab)
        axes.set_xlim(xlim)
        axes.set_ylim(ylim)

    def plot_variable_profile(x,y,v,axes,
                      xlab=None,ylab=None,vlab=None,title=None,
                      vmin=None,vmax=None,xlim=None,ylim=None,
                      loglin="log",i_time=0):
        #debug()
        if loglin=="log":
            mesh = axes.semilogx(v[:,i_time],y[:,i_time])
        elif loglin=="lin":
            mesh = axes.plot(v[:,i_time],y[:,i_time])
        if title is not None:
            axes.text(0.1,0.9,title,transform=axes.transAxes,weight='black',
                      bbox=dict(facecolor='white'))
        if xlab is not None:
            axes.set_xlabel(xlab)
        if ylab is not None:
            axes.set_ylabel(ylab)
        axes.set_xlim(xlim)
        axes.set_ylim(ylim)
    versus = -1 # Top Down
    versus =  1 # Bottom Up

    xfmt = md.DateFormatter('%H')
    ylim=(0,15)
    xDataLim = -1
    num_plot=12
    figsize=(18,12*num_plot/5.)

    # Open the netcdf meteogram file
    meteogramDataset = Dataset(meteogramFile)
    meteogramVars = meteogramDataset.variables

    # Get and Times (X axis)
    times = meteogramVars['time'] # seconds since 2015-11-24 02:00:03 proplectic gregorian
    units = times.units.split('since')[0]
    basetime = pd.to_datetime(times.units.split('since')[-1])
    #debug()
    dtimes = pd.to_timedelta(times[:].data,unit=str(units)[0])

    # Extract Heights
    H = (meteogramVars['height_2'][:]*0.001) # Kilometers
    NH = len(H)
    # Reshape times
    tt = (np.tile( (basetime + dtimes),(NH,1) ))[:,:xDataLim]
    Nt = tt.shape[1]
    # Reshape Heights
    H = np.tile(H,(Nt,1)).T

    #hydrometeors
    QNI = meteogramVars['QNI'][:xDataLim,:].T 
    QNS = meteogramVars['QNS'][:xDataLim,:].T 
    QNG = meteogramVars['QNG'][:xDataLim,:].T 
    QNR = meteogramVars['QNR'][:xDataLim,:].T 
    QI  = meteogramVars['QI'][:xDataLim,:].T 
    QIQNI  = meteogramVars['QI'][:xDataLim,:].T/meteogramVars['QNI'][:xDataLim,:].T 
    QS  = meteogramVars['QS'][:xDataLim,:].T 
    QSQNS  = meteogramVars['QS'][:xDataLim,:].T/meteogramVars['QNS'][:xDataLim,:].T  
    QG  = meteogramVars['QG'][:xDataLim,:].T 
    QR  = meteogramVars['QR'][:xDataLim,:].T 

    #atmosphere
    T  = meteogramVars['T'][:xDataLim,:].T 
    P  = meteogramVars['P'][:xDataLim,:].T 
    REL_HUM  = meteogramVars['REL_HUM'][:xDataLim,:].T 
    #derive rhi
    atmo=dict()
    atmo["rh"] = REL_HUM
    atmo["T"] = T
    rhi = __general_utilities.calc_rhi(atmo)

    #f, axs = plt.subplots(6,2,sharex=True,sharey=True,figsize=figsize)
    f, axs = plt.subplots(num_plot,2,figsize=figsize)
    vmin_qn=1e-2
    vmax_qn=1e6
    vmin_q=1e-10
    vmax_q=1e-2
    vmin_qqn=1e-10
    vmax_qqn=1e-4
    plot_variable(tt,H,QNI,axs[0,0],ylab='H    [km]',title='QNI',ylim=[0,15],vmin=vmin_qn,vmax=vmax_qn)
    plot_variable(tt,H,QNS,axs[0,1],title='QNS',ylim=[0,15],vmin=vmin_qn,vmax=vmax_qn,vlab='# / kg')
    plot_variable(tt,H,QI,axs[1,0],ylab='H    [km]',title='QI',ylim=[0,15],vmin=vmin_q,vmax=vmax_q)
    plot_variable(tt,H,QS,axs[1,1],title='QS',ylim=[0,15],vmin=vmin_q,vmax=vmax_q,vlab='kg / kg')
    plot_variable(tt,H,QIQNI,axs[2,0],title='QI/QNI',ylim=[0,15],vmin=vmin_qqn,vmax=vmax_qqn,vlab='kg')
    plot_variable(tt,H,QSQNS,axs[2,1],title='QS/QNS',ylim=[0,15],vmin=vmin_qqn,vmax=vmax_qqn,vlab='kg')
    plot_variable(tt,H,T-273.15,axs[3,0],title='T',ylim=[0,15],vmin=-40,vmax=20,loglin="lin") 
    plot_variable(tt,H,P,axs[3,1],title='P',ylim=[0,15],loglin="lin") 
    plot_variable(tt,H,REL_HUM,axs[4,0],title='RH',ylim=[0,15],loglin="lin") 
    plot_variable(tt,H,rhi,axs[4,1],title='RHi',ylim=[0,15],loglin="lin",vmin=80,vmax=120,cmap="bwr")
    #if date=="20190122": 
    #    i_time0=4000 #for 20190122
    #    i_time1=6000
    #elif date=="20181223": 
    #    i_time0=3000
    #    i_time1=5000
    #else:
    QSQNS_enoughQS = QSQNS
    QSQNS_enoughQS[QS<1e-5]=0.0
    i_time0=np.unravel_index(np.argmax(QS,axis=None),QS.shape)[1] #get time with maximum qs
    i_time1=np.unravel_index(np.argmax(QSQNS_enoughQS,axis=None),QS.shape)[1] #get time with maximum qs/qns
    for i_plot,i_time in enumerate([i_time0,i_time1]):
        #i_time = i_time1 #select time where we take the profile #i_time1 is '2019-01-22T15:00:00.000000000'
        plot_variable_profile(tt,H,REL_HUM,axs[5,i_plot],title='t(max(qs))\nRH ' + str(tt[0,i_time]),ylim=[0,15],loglin="lin",i_time=i_time,xlim=[90,110]) 
        plot_variable_profile(tt,H,rhi,axs[5,i_plot],title='t(max(qs/qns))\nRH' + str(tt[0,i_time]),ylim=[0,15],loglin="lin",i_time=i_time,xlim=[90,110]) 
    #QI
    for i_plot,i_time in enumerate([i_time0,i_time1]):
        plot_variable_profile(tt,H,QI,axs[6,i_plot],title='QI ' + str(tt[0,i_time]),ylim=[0,15],i_time=i_time,xlim=[vmin_q,vmax_q]) 
    for i_plot,i_time in enumerate([i_time0,i_time1]):
        plot_variable_profile(tt,H,QNI,axs[7,i_plot],title='QNI ' + str(tt[0,i_time]),ylim=[0,15],i_time=i_time,xlim=[vmin_qn,vmax_qn]) 
    #QS
    for i_plot,i_time in enumerate([i_time0,i_time1]):
        plot_variable_profile(tt,H,QS,axs[8,i_plot],title='QS ' + str(tt[0,i_time]),ylim=[0,15],i_time=i_time,xlim=[vmin_q,vmax_q]) 
    for i_plot,i_time in enumerate([i_time0,i_time1]):
        plot_variable_profile(tt,H,QNS,axs[9,i_plot],title='QNS ' + str(tt[0,i_time]),ylim=[0,15],i_time=i_time,xlim=[vmin_qn,vmax_qn]) 
    #mean mass Q/QN
    for i_plot,i_time in enumerate([i_time0,i_time1]):
        plot_variable_profile(tt,H,QI/QNI,axs[10,i_plot],title='QI/QNI ' + str(tt[0,i_time]),ylim=[0,15],i_time=i_time,xlim=[vmin_qqn,vmax_qqn]) 
    for i_plot,i_time in enumerate([i_time0,i_time1]):
        plot_variable_profile(tt,H,QS/QNS,axs[11,i_plot],title='QS/QNS ' + str(tt[0,i_time]),ylim=[0,15],i_time=i_time,xlim=[vmin_qqn,vmax_qqn]) 

    axs[2,0].xaxis.set_major_formatter(xfmt)
    axs[2,1].xaxis.set_major_formatter(xfmt)
    f.tight_layout()
    f.savefig('/home/mkarrer/Dokumente/plots/ICONprofiles/Hydro_content_'+ date+'.png')
    print "plot is at '/home/mkarrer/Dokumente/plots/ICONprofiles/Hydro_content_"+ date +".png"
