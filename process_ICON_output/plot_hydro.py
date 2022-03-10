# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 10:01:51 2018

@author: dori
"""

from IPython.core.debugger import Tracer; debug = Tracer()
from netCDF4 import Dataset
import matplotlib.pyplot as plt
plt.close('all')
import matplotlib.dates as md
import matplotlib.colors as colors
import numpy as np
import pandas as pd
from scipy.special import gamma
from IPython.core.debugger import Tracer; debug = Tracer()
plt.close('all')
import sys
sys.path.append("../functions/")
import __general_utilities

import argparse
parser =  argparse.ArgumentParser(description='plot hydrometeors')
parser.add_argument('--date', nargs=1, help='gimme datestring in the format YYYYMMDD')
parser.add_argument('-exp','--experiment', nargs=1,
                         help='gimme the name of the sensitivity experiment')

# Library of parameters
ice_cosmo5={'ams':1.58783,'bms':2.56,'mu':1.564,'gam':0.8547}
ice_column_narrow={'ams':0.046,'bms':2.07,'mu':5.21,'gam':0.69}
ice_hexag={'ams':150.169,'bms':3.31,'mu':2.31,'gam':1.1033}
snowSBB={'ams':0.038,'bms':2.0,'mu':1.0,'gam':1.0}
snowSBBnosphere={'ams':0.038,'bms':2.0,'mu':1.0,'gam':0.66666}
snow_mix2={'ams':0.017,'bms':1.95,'mu':0.95,'gam':0.65}
snow_cosmo5={'ams':0.146,'bms':2.1978,'mu':1.1978,'gam':1.0989}
graupelhail_cosmo5={'ams':500.86,'bms':3.18,'mu':5.37,'gam':1.06}
rainSBB={'ams':np.pi/6.0,'bms':3.0,'mu':2.0,'gam':1.0}

args = parser.parse_args()
experiment = args.experiment[0]
date = args.date[0]

# Define some paths and filenames
#date='20181110'
#experiment = 'default'
#experiment = 'colMix2_Akernel'
if experiment=="olddefault":
    oldruns = True
else:
    oldruns = False
if oldruns:
    runFld = '/data/inscape/icon/experiments/juelich/testbed/testbed_' + date + '/' #/data/optimice/pamtra_runs/'
else:
    runFld = '/data/optimice/ICON_output/' + date + '_110km_' + experiment+ '/' 
#runsFolder = '/data/optimice/ICON_output/' + date + '_110km_default' + '/' 
#runFld = runsFolder + experiment + '/'

if experiment in  ['default','olddefault']:
    ice=ice_cosmo5
    snow=snowSBBnosphere #snowSBB#snow_cosmo5#
    graupel=graupelhail_cosmo5
    rain=rainSBB
elif experiment == 'colMix2_Akernel':
    ice=ice_column_narrow
    snow=snow_mix2
    graupel=graupelhail_cosmo5
    rain=rainSBB

#meteogramName = 'METEOGRAM_patch001_' + date + '_joyce.nc' #'METEOGRAM_patch001_joyce.nc'
if oldruns:
    meteogramName = 'METEOGRAM_patch001_' + date + '_joyce.nc'
else:
    meteogramName = 'METEOGRAM_patch001_joyce.nc'
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
                  vmin=None,vmax=None,xlim=None,ylim=None,log=True):
    if log:
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        mesh = axes.pcolormesh(x,y,v,norm=norm,cmap='jet')
    else:
        mesh = axes.pcolormesh(x,y,v,vmin=vmin, vmax=vmax,cmap='jet')
        
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

QNI = np.ma.masked_array(QNI,mask=QNI<1e-4)
QI = np.ma.masked_array(QI,mask=QNI<1e-4)
QNS = np.ma.masked_array(QNS,mask=QNS<1e-4)
QS = np.ma.masked_array(QS,mask=QNS<1e-4)
QNG = np.ma.masked_array(QNG,mask=QNG<1e-4)
QG = np.ma.masked_array(QG,mask=QNG<1e-4)
QNR = np.ma.masked_array(QNR,mask=QNR<1e-4)
QR = np.ma.masked_array(QR,mask=QNR<1e-4)

LamI, N0I = mgammaLambda13(QNI,QI,ice)
#LamI = np.ma.masked_array(LamI,mask=QNI<1e-4)
LamS, N0S = mgammaLambda13(QNS,QS,snow)
LamG, N0G = mgammaLambda13(QNG,QG,graupel)
LamR, N0R = mgammaLambda13(QNR,QR,rain)


#f, ((ax11,ax12,ax13,ax14),(ax21,ax22,ax23,ax24),(ax31,ax32,ax33,ax34)) = plt.subplots(3,4,sharex=True,sharey=True,figsize=figsize)
f, ((ax11,ax12,ax13,ax14),(ax21,ax22,ax23,ax24),(ax31,ax32,ax33,ax34)) = plt.subplots(3,4,sharex=True,sharey=True,figsize=figsize)
#f, ((ax11,ax12),(ax21,ax22),(ax31,ax32)) = plt.subplots(3,4,sharex=True,sharey=True,figsize=figsize)
#f, axs = plt.subplots(3,2,sharex=True,sharey=True,figsize=figsize)
qnmin=1e-4; qnmax=1e6
qmin=1e-10; qmax=1e-2
plot_variable(tt,H,QNI,ax11,ylab='H    [km]',title='QNI',ylim=[0,15],vmin=qnmin,vmax=qnmax)
plot_variable(tt,H,QNS,ax12,title='QNS',ylim=[0,15],vmin=qnmin,vmax=qnmax) #,vlab='1 / kg')
plot_variable(tt,H,QNG,ax13,title='QNG',ylim=[0,15],vmin=qnmin,vmax=qnmax)
plot_variable(tt,H,QNR,ax14,title='QNR',ylim=[0,15],vmin=qnmin,vmax=qnmax,vlab='1 / kg')
plot_variable(tt,H,QI,ax21,ylab='H    [km]',title='QI',ylim=[0,15],vmin=qmin,vmax=qmax)
plot_variable(tt,H,QS,ax22,title='QS',ylim=[0,15],vmin=qmin,vmax=qmax) #,vlab='kg / kg')
plot_variable(tt,H,QG,ax23,title='QG',ylim=[0,15],vmin=qmin,vmax=qmax)
plot_variable(tt,H,QR,ax24,title='QR',ylim=[0,15],vmin=qmin,vmax=qmax,vlab='kg /kg')
print(LamI)
plot_variable(tt,H,3670/LamI,ax31,xlab='time',ylab='H    [km]',title='D$_0$I',ylim=[0,15],vmin=0.0,vmax=5,log=False)
plot_variable(tt,H,3670/LamS,ax32,xlab='time',title='D$_0$S',ylim=[0,15],vmin=0.0,vmax=20,vlab='mm',log=False)
plot_variable(tt,H,3670/LamG,ax33,title='D$_0$G',ylim=[0,15],vmin=0.0,vmax=20,log=False)
plot_variable(tt,H,3670/LamR,ax34,title='D$_0$R',ylim=[0,15],vmin=0.0,vmax=20,vlab='mm',log=False)
ax31.xaxis.set_major_formatter(xfmt)
ax32.xaxis.set_major_formatter(xfmt)
ax33.xaxis.set_major_formatter(xfmt)
ax34.xaxis.set_major_formatter(xfmt)
f.tight_layout()
f.savefig('/home/mkarrer/Dokumente/plots/ICONprofiles/Hydro_content_'+ date+ experiment +'.png')
print("plot is at '/home/mkarrer/Dokumente/plots/ICONprofiles/Hydro_content_"+ date + experiment+ ".png")
plt.close('all')

#plot some atmospheric variables
f, ((ax11,ax12),(ax21,ax22)) = plt.subplots(2,2,figsize=figsize)
REL_HUM = meteogramVars['REL_HUM'][:xDataLim,:].T 
T = meteogramVars['T'][:xDataLim,:].T 
P = meteogramVars['P'][:xDataLim,:].T 
rain = meteogramVars['RAIN_GSP'][:]
snow = meteogramVars['SNOW_GSP'][:]
#derive rhi
atmo=dict()
atmo["rh"] = REL_HUM
atmo["T"] = T
rhi = __general_utilities.calc_rhi(atmo)

plot_variable(tt,H,T-273.15,ax11,xlab='time',ylab='H    [km]',title='T',ylim=[0,15],vmin=-40.0,vmax=20,log=False)
ax12.plot(basetime+dtimes,rain+snow)
ax12.set_ylabel("accumulated precip. kg/m2")
plot_variable(tt,H,REL_HUM,ax21,xlab='time',ylab='H    [km]',title='RH',ylim=[0,15],vmin=0.0,vmax=100,log=False)
plot_variable(tt,H,rhi,ax22,xlab='time',ylab='H    [km]',title='RHi',ylim=[0,15],vmin=80,vmax=120,log=False)
f.savefig('/home/mkarrer/Dokumente/plots/ICONprofiles/Atmos_'+ date+ experiment +'.png')
print("plot is at '/home/mkarrer/Dokumente/plots/ICONprofiles/Atmos_"+ date + experiment+ ".png")
