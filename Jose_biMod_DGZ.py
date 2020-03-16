#!/usr/bin/env python
# coding: utf-8
from IPython.core.debugger import Tracer ; debug=Tracer()
# In[ ]:




# In[1]:


import pyPamtra

import xarray as xr
import numpy as np
import pandas as pd
#from scipy.special import gamma
import scipy.special as spe
#import os


import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('seaborn')
import __fallspeed_relations

# In[2]:


def normedGamma(D, D0, mu, Nw):
    
    Lambda = (3.67 + mu)/D0
    f_mu = (6/3.67**(4)) * ((3.67 + mu)**(mu + 4)/spe.gamma(mu + 4))
#     N_D = Nw * f_mu * (D/D0)**(mu) * np.exp(-Lambda*D)
    
    try:
        N_D = np.ones((len(D0),len(D)))
        for i in range(len(D0)):
            N_D[i] = Nw[i] * f_mu[i] * (D/D0[i])**(mu[i]) * np.exp(-Lambda[i]*D)
            
    except:
            N_D = Nw * f_mu * (D/D0)**(mu) * np.exp(-Lambda*D)
    
    return N_D


# In[3]:


def gammaFunc(D, mu, N0, Lambda):
    
    N_D = N0 * D**(mu) * np.exp(-Lambda * D)
    
    return N_D


# In[4]:


def getN0(Nmax, DNmax, Lambda):
    
    N0 = Nmax / (DNmax**(DNmax*Lambda) * np.exp(-DNmax*Lambda))
    
    return N0


# In[5]:


def getmu(DNmax, Lambda):
    
    mu = DNmax*Lambda
    
    return mu


# In[6]:


#(mitchell 1996)
massAreaParam = {'Hexagonal_plates1':{'alpha':0.00739,'beta':2.45, 'gamma':0.24, 'sigma':1.85, 'D_range':(15,100)},
                 'Hexagonal_plates2':{'alpha':0.00739,'beta':2.45, 'gamma':0.65, 'sigma':2.00, 'D_range':(50,3000)},
                 
                 'Crystal_with_sector_like_branches1':{'alpha':0.00614,'beta':2.42, 'gamma':0.24, 'sigma':1.85, 'D_range':(10,40)},
                 'Crystal_with_sector_like_branches2':{'alpha':0.00142,'beta':2.02, 'gamma':0.55, 'sigma':1.97, 'D_range':(500,3000)}, #50,3000)},
                 
                 'Broad_branched_crystal1':{'alpha':0.00583,'beta':2.42, 'gamma':0.24, 'sigma':1.85, 'D_range':(10,100)},
                 'Broad_branched_crystal2':{'alpha':0.000516,'beta':1.80, 'gamma':0.21, 'sigma':1.76, 'D_range':(50,3000)},
                 
                 'Stellar_crystal_with_broad_arms1':{'alpha':0.00583,'beta':2.42, 'gamma':0.24, 'sigma':1.85, 'D_range':(10,90)},
                 'Stellar_crystal_with_broad_arms2':{'alpha':0.000270,'beta':1.67, 'gamma':0.11, 'sigma':1.63, 'D_range':(50,3000)},
                 
                 'Aggregates_of_side_planes':{'alpha':0.0033,'beta':2.2, 'gamma':0.2285, 'sigma':1.88, 'D_range':(300,4100)},
                }


# In[7]:


binsNumber = 200 
selHgtLev = 3360, 3390


# In[8]:


#---mode_1 definitions---#

particle_M1 = 'Aggregates_of_side_planes'
Drange_M1 = massAreaParam[particle_M1]['D_range']    
alpha_M1 = massAreaParam[particle_M1]['alpha']
beta_M1 = massAreaParam[particle_M1]['beta']
gamma_M1 = massAreaParam[particle_M1]['gamma']
sigma_M1 = massAreaParam[particle_M1]['sigma']

#min and max size of the crystals
#binSize_M1 = 20*1e-6
minSize_M1 = Drange_M1[0]*1e-6
maxSize_M1 = Drange_M1[1]*1e-6

#Ds_M1 = np.arange(minSize_M1, maxSize_M2, binSize_M1)
Ds_M1 = np.linspace(minSize_M1, maxSize_M1, binsNumber)
#centerBin_M1 = Ds_M1[:-1] + 0.5 * Ds_M1[:-1]
centerBin_M1 = Ds_M1[:-1] + np.diff(Ds_M1) * 0.5
binSize_M1 = np.abs(Ds_M1[0]-Ds_M1[1])

mu_M1 = 9
N0_M1 = 0 #1.5e37
Lambda_M1 = 9e3#5e3#9e3

psd_M1 = gammaFunc(centerBin_M1, mu_M1, N0_M1, Lambda_M1)                   



# In[9]:


#---mode_2 definitions---#

#particle_M2 = 'Hexagonal_plates2'
#Drange_M2 = massAreaParam[particle_M2]['D_range']    
#alpha_M2 = massAreaParam[particle_M2]['alpha']
#beta_M2 = massAreaParam[particle_M2]['beta']
#gamma_M2 = massAreaParam[particle_M2]['gamma']
#sigma_M2 = massAreaParam[particle_M2]['sigma']
#
##min and max size of the crystals
#minSize_M2 = Drange_M2[0]*1e-6
#maxSize_M2 = Drange_M2[1]*1e-6
#
#Ds_M2 = np.linspace(minSize_M2, maxSize_M2, binsNumber)
#centerBin_M2 = Ds_M2[:-1] + np.diff(Ds_M2) * 0.5
#binSize_M2 = np.abs(Ds_M2[0]-Ds_M2[1])
#
#Lambda_M2 = 1.5e4
#Nmax_M2 = 6e7
#DNmax_M2 = 2e-4
#
#N0_M2 = getN0(Nmax_M2, DNmax_M2, Lambda_M2)
#mu_M2 = getmu(DNmax_M2, Lambda_M2)
#
#psd_M2 = gammaFunc(centerBin_M2, mu_M2, N0_M2, Lambda_M2)
#
#print 'N0_M2 {0:1.2e}'.format(N0_M2), 'Mu '+ str(mu_M2)
#
#
# In[10]:


# #---mode_2 definitions---#

particle_M2 = 'Crystal_with_sector_like_branches2'
Drange_M2 = massAreaParam[particle_M2]['D_range']    
alpha_M2 = massAreaParam[particle_M2]['alpha']
beta_M2 = massAreaParam[particle_M2]['beta']
gamma_M2 = massAreaParam[particle_M2]['gamma']
sigma_M2 = massAreaParam[particle_M2]['sigma']

#min and max size of the crystals
minSize_M2 = Drange_M2[0]*1e-6
maxSize_M2 = Drange_M2[1]*1e-6

Ds_M2 = np.linspace(minSize_M2, maxSize_M2, binsNumber)
centerBin_M2 = Ds_M2[:-1] + np.diff(Ds_M2) * 0.5
binSize_M2 = np.abs(Ds_M2[0]-Ds_M2[1])

Nmax_M2 = 5e7
DNmax_M2 = 0.5e-4
Lambda_M2 = 6e3

N0_M2 = getN0(Nmax_M2, DNmax_M2, Lambda_M2)
mu_M2 = getmu(DNmax_M2, Lambda_M2)

psd_M2 = gammaFunc(centerBin_M2, mu_M2, N0_M2, Lambda_M2)

print 'N0_M2 {0:1.2e}'.format(N0_M2), 'Mu '+ str(mu_M2)


# In[11]:


#--- plotting PSDs ---#
fs = 14
plt.semilogy(centerBin_M1, psd_M1, lw = 3,
             label='M1 {0:1.3} [# m -3]'.format(np.sum(psd_M1*binSize_M1)))
plt.semilogy(centerBin_M2, psd_M2, lw =3, 
             label='M2 {0:1.3} [# m -3]'.format(np.sum(psd_M2*binSize_M2)))

plt.xlim(0,7e-3)
plt.ylim(1e-6,1e9)
plt.ylabel('# m -4' ,fontsize=fs)
plt.xlabel('diameter [m]', fontsize=fs)
plt.legend(fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)


#plt.show()


# In[12]:


pam = pyPamtra.pyPamtra()
pam = pyPamtra.importer.createUsStandardProfile(pam, hgt_lev=selHgtLev)

scatteringCalc = 'tmatrix'
#scatteringCalc = 'ss-rayleigh-gans'

pam.df.addHydrometeor(('mode1', -99, -1, -99,
                       -99, -99, -99, -99,
                       0, binsNumber-1, 
                       None, -99, -99, -99, -99,
                       -99, -99, scatteringCalc,
                       'heymsfield10_particles', -99)) #'corPowerLaw_4.49_0.25'               

pam.df.addHydrometeor(('mode2', -99, -1, -99,
                       -99, -99, -99, -99,
                       0, binsNumber-1, 
                       None, -99, -99, -99, -99,
                       -99, -99, scatteringCalc,
                       'heymsfield10_particles', -99)) #'corPowerLaw_4.49_0.25'               


# In[13]:


#--- pamtra settings ---#
pam.p["airturb"][:] = 0
pam.p["wind_w"][:] =0
pam.p['temp_lev'][:]=259
pam.p['press_lev'][:]=58517

pam.nmlSet["passive"] = False  # Activate this for Microwave radiometer
pam.nmlSet["radar_mode"] = "spectrum"
pam.nmlSet["radar_pnoise0"]=-100 #set radar noise to an arbitrary low number
#pam.nmlSet["radar_noise_distance_factor"] = 2.0
pam.nmlSet["radar_save_noise_corrected_spectra"]=  False
pam.nmlSet["randomseed"]=  10
pam.nmlSet['radar_airmotion'] = True
pam.nmlSet['radar_airmotion_model'] = 'constant'
pam.nmlSet['radar_airmotion_vmax'] = 0.0
pam.nmlSet['radar_airmotion_vmin'] = 0.0
pam.nmlSet["radar_nfft"] = 512
pam.nmlSet['radar_aliasing_nyquist_interv'] = 1
pam.nmlSet['radar_no_Ave'] = 20
pam.nmlSet["save_psd"] = True
pam.nmlSet["hydro_adaptive_grid"]=False
pam.nmlSet["radar_max_V"] = 3.
pam.nmlSet["radar_min_V"] = -3.
pam.nmlSet["hydro_fullspec"] = True

#polarimetric calculation
pam.nmlSet["radar_polarisation"] = "NN,HH,VV"#,HV,VH"
pam.df.addFullSpectra()


# In[14]:


#elv = 30
elv = 90
canting = 90 - elv

#--- loading mode_1 settings ---#
a_M1 = alpha_M1*(10**(-3+2*beta_M1)) #transformation from cgs to SI units
g_M1 = gamma_M1*(10**(-4+2*sigma_M1)) #transformation from cgs to SI units

as_ratio_M1 = 0.4
mass_M1 = a_M1*centerBin_M1**beta_M1
vol_M1 = (np.pi/6.)*as_ratio_M1*(centerBin_M1)**3.
rho_M1 = mass_M1/vol_M1

area_M1 = g_M1*centerBin_M1**sigma_M1

pam.df.dataFullSpec['canting'][0,0,0,0,:] = canting
pam.df.dataFullSpec['d_bound_ds'][0,0,0,0,:] = Ds_M1
pam.df.dataFullSpec['as_ratio'][0,0,0,0,:] = as_ratio_M1
pam.df.dataFullSpec['d_ds'][0,0,0,0,:] = centerBin_M1
pam.df.dataFullSpec['n_ds'][0,0,0,0,:] = psd_M1*binSize_M1
pam.df.dataFullSpec['area_ds'][0,0,0,0,:] = area_M1
pam.df.dataFullSpec['mass_ds'][0,0,0,0,:] = mass_M1 # needed for SS and TM
#pam.df.dataFullSpec['rho_ds'][0,0,0,0,:] = 200 ##SS
pam.df.dataFullSpec['rho_ds'][0,0,0,0,:] = rho_M1 ##TM


#--- loading mode_2 settings ---#
a_M2 = alpha_M2*(10**(-3+2*beta_M2)) #transformation from cgs to SI units 
g_M2 = gamma_M2*(10**(-4+2*sigma_M2)) #transformation from cgs to SI units

as_ratio_M2 = 0.2
mass_M2 = a_M2*centerBin_M2**beta_M2
vol_M2 = (np.pi/6.)*as_ratio_M2*(centerBin_M2)**3.
rho_M2 = mass_M2/vol_M2

area_M2 = g_M2*centerBin_M2**sigma_M2
print "area ratio:", area_M2/(np.pi/4.*centerBin_M2**2.)
print "eff density:", mass_M2/(np.pi/6.*centerBin_M2**3.)
print "D", Ds_M2
vterm = __fallspeed_relations.calc_vterm("HW10",mass_M2,centerBin_M2,area_M2)
plt.clf()
plt.semilogx(Ds_M2[:-1],vterm)
plt.show()

#raw_input(); import sys; sys.exit()
pam.df.dataFullSpec['canting'][0,0,0,1,:] = canting
pam.df.dataFullSpec['d_bound_ds'][0,0,0,1,:] = Ds_M2
pam.df.dataFullSpec['as_ratio'][0,0,0,1,:] = as_ratio_M2
pam.df.dataFullSpec['d_ds'][0,0,0,1,:] = centerBin_M2
pam.df.dataFullSpec['n_ds'][0,0,0,1,:] = psd_M2*binSize_M2
pam.df.dataFullSpec['area_ds'][0,0,0,1,:] = area_M2
pam.df.dataFullSpec['mass_ds'][0,0,0,1,:] = mass_M2 # needed for SS and TM
#pam.df.dataFullSpec['rho_ds'][0,0,0,1,:] = 200 ##SS
pam.df.dataFullSpec['rho_ds'][0,0,0,1,:] = rho_M2 ##TM

plt.clf()
#debug()
plt.semilogx(pam.df.dataFullSpec['d_ds'][0,0,0,1,:],pam.df.dataFullSpec['n_ds'][0,0,0,1,:])
print "n_ds vs d_ds"
plt.show(); sys.exit()

#debug()
# In[15]:


#---comparison plots 

dicData = {'area':{'data':area_M1, 'units':'m^2'},
           'mass':{'data':mass_M1, 'units':'kg'},
           'volume':{'data':vol_M1, 'units':'m^3'},
           'density':{'data':rho_M1*1e-3, 'units':'g cm^-3'},
          }

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(18,12))
flatAxes = axes.flatten()

for i, variable in enumerate(sorted(dicData.keys())):

    flatAxes[i].plot(centerBin_M1, dicData[variable]['data'])
    flatAxes[i].set_xlabel('diamiter [m]')
    flatAxes[i].set_ylabel(dicData[variable]['units'])
    flatAxes[i].set_title(variable)
    
plt.show()


# In[16]:


pam.set["verbose"] = 0
pam.runParallelPamtra([9.5, 35, 95], pp_deltaX=1, 
                      pp_deltaY=1, pp_deltaF=1,
                      pp_local_workers='auto')


# In[17]:


#--- pamtra output ---#
pam_output = dict()

pam_output["spectra_X_NN"] = pam.r["radar_spectra"][0,0,0,0,0,:]
pam_output["spectra_X_HH"] = pam.r["radar_spectra"][0,0,0,0,1,:]
pam_output["spectra_X_VV"] = pam.r["radar_spectra"][0,0,0,0,2,:]

pam_output["spectra_Ka_NN"] = pam.r["radar_spectra"][0,0,0,1,0,:]
pam_output["spectra_Ka_HH"] = pam.r["radar_spectra"][0,0,0,1,1,:]
pam_output["spectra_Ka_VV"] = pam.r["radar_spectra"][0,0,0,1,2,:]

pam_output["spectra_W_NN"] = pam.r["radar_spectra"][0,0,0,2,0,:]
pam_output["spectra_W_HH"] = pam.r["radar_spectra"][0,0,0,2,1,:]
pam_output["spectra_W_VV"] = pam.r["radar_spectra"][0,0,0,2,2,:]

pam_output["velbins_X"] = pam.r["radar_vel"][0,:]
pam_output["velbins_Ka"] = pam.r["radar_vel"][1,:]
pam_output["velbins_W"] = pam.r["radar_vel"][2,:]

pam_output["Ze_X"] = pam.r["Ze"][0,0,0,0,0,0]
pam_output["Ze_Ka"] = pam.r["Ze"][0,0,0,1,0,0]
pam_output["Ze_W"] = pam.r["Ze"][0,0,0,2,0,0]

pam_output['MDV_Ka'] = pam.r["radar_moments"][0,0,0,1,0,0,0]
pam_output['SW_Ka'] = pam.r["radar_moments"][0,0,0,1,0,0,1]
pam_output["SK_Ka"] = pam.r["radar_moments"][0,0,0,1,0,0,2]

pam_output['radar_hgt'] = pam.r['radar_hgt'][0,0,:]


# In[18]:


print pam_output['Ze_X'], pam_output['Ze_Ka'], pam_output['Ze_W'] 


# In[19]:


print pam_output['MDV_Ka'], pam_output['SW_Ka'], -pam_output['SK_Ka']


# In[20]:


selObsSpec = xr.open_dataset('/home/jdias/Projects/radarCode/tripexPolPro/selSpecDT.nc')

selObsSpec.xSpec.plot(lw=4, label='radar_X', color='C0')
selObsSpec.kaSpec.plot(lw=4, label='radar_Ka', color='C1')
selObsSpec.wSpec.plot(lw=4, label='radar_W', color='C2')

# for i, spec in enumerate(['spectra_X_HH','spectra_Ka_HH','spectra_W_HH']):

#     plt.plot(-pam_output['velbins_Ka'], 
#               pam_output[spec], label=spec + ' Cant: {elv}'.format(elv=canting))

plt.title('temperature: {temp}'.format(temp=np.round(selObsSpec.temp.values,2)))
    
plt.xlim(-2, 0.5)
plt.ylim(-50, -5)
plt.ylabel('SP [dB]')
plt.xlabel('DV [m s-1]')
plt.legend()        
plt.show()


# In[21]:


selObsSpec = xr.open_dataset('/home/jdias/Projects/radarCode/tripexPolPro/selSpecDT.nc')
plot = selObsSpec.kaSpec.plot(lw=4, label='radar_Ka')
#selObsSpec.wSpec.plot()
for i, spec in enumerate(['spectra_X_VV','spectra_Ka_VV','spectra_W_VV']):

    plot = plt.plot(-pam_output['velbins_Ka'], 
                     pam_output[spec], lw=3,
                    label=spec + ' Cant: {elv}'.format(elv=canting))
    
plt.xlim(-2, 0.5)
plt.ylim(-50,-5)
plt.ylabel('SP [dB]', fontsize=fs)
plt.xlabel('DV [m s-1]', fontsize=fs)
plt.legend(fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)     
plt.show()


# In[22]:


#selObsSpec.wSpec.plot()
specZDR = pam_output['spectra_W_HH'] - pam_output['spectra_W_VV']
plot = plt.plot(-pam_output['velbins_W'], 
                specZDR, label='Zdr W Cant: {elv}'.format(elv=canting))
    
plt.xlim(-2, 0.5)
plt.ylim(-0.2, 6)
plt.ylabel('Zdr [dB]')
plt.xlabel('DV [m s-1]')
plt.legend()        
plt.show()


# In[ ]:




