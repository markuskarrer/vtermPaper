# coding: utf-8
#import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.pylab as pylab
import os
import subprocess
#import other self-defined functions
import __postprocess_McSnow
import __postprocess_SB
#from IPython.core.debugger import Tracer ; Tracer()()

'''
plot the mass over diameter for the SB and the McSnow scheme
'''

#set up diameter array
n_steps = 1001
min_D = 1e-5; max_D = 5e-2
D_array = np.logspace(np.log10(min_D),np.log10(max_D),n_steps)

#get parameter from McSnow
mth,unr_alf,unr_bet,rhoi,rhol,Dth = __postprocess_McSnow.return_parameter_mD_AD_rel()
#get parameters from the SB-categories
cloud_water,rain,cloud_ice,snow,graupel,hail = __postprocess_SB.init_particles()
cloud_water,rain,cloud_ice,snow,graupel,hail = __postprocess_SB.convert_Nm_to_ND(cloud_water,rain,cloud_ice,snow,graupel,hail)

#initialize array with masses
m_MC = np.zeros(D_array.shape)
m_SB_ice = np.zeros(D_array.shape)
m_SB_snow = np.zeros(D_array.shape)

#calculate the masses for the different m-D regions of McSnow
for i, D in enumerate(D_array):
    if D<=Dth:
        m_MC[i] = np.pi/6.*rhoi*D**3
    elif D>Dth:
        m_MC[i] = unr_alf*D**unr_bet
#calculate the masses for the SB-categories (we do not need a case distinction as for McSnow)
m_SB_ice = cloud_ice.a_ms*D_array**cloud_ice.b_ms
m_SB_snow = snow.a_ms*D_array**snow.b_ms
#apply size limits
m_SB_ice[m_SB_ice<cloud_ice.xmin] = np.nan
m_SB_ice[m_SB_ice>cloud_ice.xmax] = np.nan
m_SB_snow[m_SB_snow<snow.xmin] = np.nan
m_SB_snow[m_SB_snow>snow.xmax] = np.nan
        
#increase font sizes
params = {'legend.fontsize': 'large',
    'figure.figsize': (15, 5),
    'axes.labelsize': 'x-large', #size: Either a relative value of 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large' or an absolute font size, e.g., 12
    'axes.titlesize':'x-large',
    'xtick.labelsize':'x-large',
    'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

#define figure
number_of_plots = 2
figsize_height = 6.0/2.0*number_of_plots
fig, axes = plt.subplots(nrows=number_of_plots, ncols=1, figsize=(8.0,figsize_height))
###
#subplot 1: m vs D
###
axes2 = axes[0].twinx()
#plot in loglog
axes[0].loglog(D_array,m_MC,linestyle='--',label="McSnow")
axes[0].loglog(D_array,m_SB_ice,color='blue',label="SB cloud ice")
axes[0].loglog(D_array,m_SB_snow,color='green',label="SB snow")
axes2.semilogx(D_array,(m_SB_ice-m_MC)/m_MC,color='red',label="|SB cloud ice - McSnow|/MC_snow")
#define labels
axes[0].set_xlabel("diameter / m")
axes2.set_xlabel("diameter / m")
axes[0].set_ylabel("mass / kg")
axes2.set_ylabel("relative diff. of \n SB cloud ice to McSnow")
#set limits
axes2.set_ylim([-1,3])
#show legend
axes[0].legend()
#axes2.legend()
plt.tight_layout()
###
#calculate D for subplot 2: D vs m
###
#set up mass array
n_steps = 1001
min_m = 1e-13; max_m = 1e-4
m_array = np.logspace(np.log10(min_m),np.log10(max_m),n_steps)


#initialize array with masses
D_MC = np.zeros(m_array.shape)
D_SB_ice = np.zeros(m_array.shape)
D_SB_snow = np.zeros(m_array.shape)

#calculate the masses for the different m-D regions of McSnow
for i, m in enumerate(m_array):
    if m<=mth:
        D_MC[i] = (6.*m/(np.pi*rhoi))**(1./3.)
    elif m>mth:
        D_MC[i] = (m/unr_alf)**(1./unr_bet)
#calculate the masses for the SB-categories (we do not need a case distinction as for McSnow)
D_SB_ice = (m_array/cloud_ice.a_ms)**(1./cloud_ice.b_ms)
D_SB_snow = (m_array/snow.a_ms)**(1./snow.b_ms)
#apply size limits
D_SB_ice[D_SB_ice<cloud_ice.Dmin] = np.nan
D_SB_ice[D_SB_ice>cloud_ice.Dmax] = np.nan
D_SB_snow[D_SB_snow<snow.Dmin] = np.nan
D_SB_snow[D_SB_snow>snow.Dmax] = np.nan
#plot subplot 2
axes22 = axes[1].twinx()
#plot in loglog
axes[1].loglog(m_array,D_MC,linestyle='--',label="McSnow")
axes[1].loglog(m_array,D_SB_ice,color='blue',label="SB cloud ice")
axes[1].loglog(m_array,D_SB_snow,color='green',label="SB snow")
axes22.semilogx(m_array,(D_SB_ice-D_MC)/D_MC,color='red',label="|SB cloud ice - McSnow|/MC_snow")
#define labels
axes[1].set_xlabel("mass / kg")
axes22.set_xlabel("mass / kg")
axes[1].set_ylabel("diameter / m")
axes22.set_ylabel("relative diff. of \n SB cloud ice to McSnow")
#set limits
#axes2.set_ylim([-1,5])
#show legend
axes[1].legend()
#axes2.legend()
plt.tight_layout()
dir_save = '/home/mkarrer/Dokumente/plots/'
if not os.path.exists(dir_save): #create direktory if it does not exists
    os.makedirs(dir_save)
out_filestring = "mD_rel"
plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
plt.savefig(dir_save + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])