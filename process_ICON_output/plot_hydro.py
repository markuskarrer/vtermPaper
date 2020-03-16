# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 10:01:51 2018

@author: dori
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
plt.close('all')
import matplotlib.dates as md
import numpy as np
import pandas as pd
from scipy.special import gamma
from IPython.core.debugger import Tracer; debug = Tracer()
plt.close('all')

# Library of parameters
ice_cosmo5={'ams':1.58783,'bms':2.56,'mu':1.564,'gam':0.8547}
ice_hexag={'ams':150.169,'bms':3.31,'mu':2.31,'gam':1.1033}
snowSBB={'ams':0.038,'bms':2.0,'mu':1.0,'gam':1.0}
snowSBBnosphere={'ams':0.038,'bms':2.0,'mu':1.0,'gam':0.66666}
snow_cosmo5={'ams':0.146,'bms':2.1978,'mu':1.1978,'gam':1.0989}
graupelhail_cosmo5={'ams':500.86,'bms':3.18,'mu':5.37,'gam':1.06}
rainSBB={'ams':np.pi/6.0,'bms':3.0,'mu':2.0,'gam':1.0}

# Define some paths and filenames
date='20190122'
runsFolder = '/data/inscape/icon/experiments/juelich/testbed/testbed_' + date + '/' #/data/optimice/pamtra_runs/'
#runName = '20151124_snowcosmo'
#runName = 'newicon'
#runName = '20151124_noagg'
#runName = '20151124_operational'
#runName = 'radcoupl'
#runName = 'nosnowselfcoll'
#runName = 'hexaplates'
#runName = 'xmin10-14'
#runName = 'stickeff/stickeff4'
#runName = 'stickeff/stickeff5'
#runName = 'stickeff/ref'
#runName = 'atlasvel'
runName = ''
runFld = runsFolder + runName + '/'

ice=ice_hexag#ice_cosmo5
snow=snowSBBnosphere#snowSBB#snow_cosmo5#
graupel=graupelhail_cosmo5
rain=rainSBB

meteogramName = 'METEOGRAM_patch001_' + date + '_joyce.nc' #'METEOGRAM_patch001_joyce.nc'
#meteogramName = '1d_vars_DOM01.nc'
meteogramFile = runFld + meteogramName

runs = ['all_hydro_mom.nc', 'no_snow_mom.nc', 'only_ice_mom.nc', 'only_liquid_mom.nc', 'only_snow_mom.nc', 'only_graupel_hail_mom.nc']
titles = ['all Hydrometeors', 'No Snow', 'Only Ice', 'Only liquid (cloud drops and rain)', 'only Snow', 'only Graupel and Hail']

runTitles=dict(zip(runs,titles))



def mgammaLambda13(N,q,params):
    mu = params['mu']
    gam = params['gam']
    ams = params['ams']
    bms = params['bms']
    
    temp2 = gamma((mu + bms + 1.0)/gam)
    temp3 = gamma((mu + 1.0)/gam)
    Lambda = (ams/q*N*temp2/temp3)**(gam/bms)
    N0 = gam*N/temp3*Lambda**((mu + 1.0)/gam)
    return Lambda, N0

# Define Plotting Function
def plot_variable(x,y,v,axes,
                  xlab=None,ylab=None,vlab=None,title=None,
                  vmin=None,vmax=None,xlim=None,ylim=None):
    mesh = axes.pcolormesh(x,y,v,vmin=vmin,vmax=vmax,cmap='jet')
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

versus = -1 # Top Down
versus =  1 # Bottom Up

xfmt = md.DateFormatter('%H')
ylim=(0,15)
xDataLim = -1
figsize=(18,12)

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

QNI = meteogramVars['QNI'][:xDataLim,:].T 
QNS = meteogramVars['QNS'][:xDataLim,:].T 
QNG = meteogramVars['QNG'][:xDataLim,:].T 
QNR = meteogramVars['QNR'][:xDataLim,:].T 
QI  = meteogramVars['QI'][:xDataLim,:].T 
QS  = meteogramVars['QS'][:xDataLim,:].T 
QG  = meteogramVars['QG'][:xDataLim,:].T 
QR  = meteogramVars['QR'][:xDataLim,:].T 

LamI, N0I = mgammaLambda13(QNI,QI,ice)
LamS, N0S = mgammaLambda13(QNS,QS,snow)
LamG, N0G = mgammaLambda13(QNG,QG,graupel)
LamR, N0R = mgammaLambda13(QNR,QR,rain)


#f, ((ax11,ax12,ax13,ax14),(ax21,ax22,ax23,ax24),(ax31,ax32,ax33,ax34)) = plt.subplots(3,4,sharex=True,sharey=True,figsize=figsize)
#f, ((ax11,ax12),(ax21,ax22),(ax31,ax32)) = plt.subplots(3,4,sharex=True,sharey=True,figsize=figsize)
f, axs = plt.subplots(3,2,sharex=True,sharey=True,figsize=figsize)
plot_variable(tt,H,QNI/1000.,axs[0,0],ylab='H    [km]',title='QNI',ylim=[0,15],vmin=0.0,vmax=1000)
plot_variable(tt,H,QNS/1000.,axs[0,1],title='QNS',ylim=[0,15],vmin=0.0,vmax=80,vlab='thousands / kg')
#plot_variable(tt,H,QNG/1000.,ax13,title='QNG',ylim=[0,15],vmin=0.0,vmax=8)
#plot_variable(tt,H,QNR/1000.,ax14,title='QNR',ylim=[0,15],vmin=0.0,vmax=20,vlab='thousands / m$^3$')
plot_variable(tt,H,1000.*QI,axs[1,0],ylab='H    [km]',title='QI',ylim=[0,15],vmin=0.0,vmax=0.25)
plot_variable(tt,H,1000.*QS,axs[1,1],title='QS',ylim=[0,15],vmin=0.0,vmax=1.0,vlab='g / kg')
#plot_variable(tt,H,1000.*QG,ax23,title='QG',ylim=[0,15],vmin=0.0,vmax=0.5)
#plot_variable(tt,H,1000.*QR,ax24,title='QR',ylim=[0,15],vmin=0.0,vmax=0.05,vlab='g/kg')
plot_variable(tt,H,3670/LamI,axs[2,0],xlab='time',ylab='H    [km]',title='D$_0$I',ylim=[0,15],vmin=0.0,vmax=20)
plot_variable(tt,H,3670/LamS,axs[2,1],xlab='time',title='D$_0$S',ylim=[0,15],vmin=0.0,vmax=20,vlab='mm')
#plot_variable(tt,H,3670/LamG,ax33,title='D$_0$G',ylim=[0,15],vmin=0.0,vmax=20)
#plot_variable(tt,H,3670/LamR,ax34,title='D$_0$R',ylim=[0,15],vmin=0.0,vmax=20,vlab='mm')
axs[2,0].xaxis.set_major_formatter(xfmt)
axs[2,1].xaxis.set_major_formatter(xfmt)
#ax33.xaxis.set_major_formatter(xfmt)
#ax34.xaxis.set_major_formatter(xfmt)
f.tight_layout()
f.savefig('/home/mkarrer/Dokumente/plots/ICONprofiles/Hydro_content_'+ date+'.png')
print "plot is at '/home/mkarrer/Dokumente/plots/ICONprofiles/Hydro_content_"+ date +".png"
#plt.close('all')
