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
min_D = 1e-5; max_D = 3e-2
D_array = np.logspace(np.log10(min_D),np.log10(max_D),n_steps)

#get parameter from McSnow in the state of Brdar&Seifert2018
mth,unr_alf,unr_bet,rhoi,rhol,Dth,unr_sig,unr_gam,sph_sig,sph_gam =  __postprocess_McSnow.return_parameter_mD_AD_rel("1d_xi") #the argument is a workaround to get the standard McSnow settings
sph_alf=np.pi/6.*rhoi
sph_bet=3.
print "used constants (control MC)"
print "mth,               unr_alf,    unr_bet,rhoi,  rhol,       Dth,     unr_sig,     unr_gam,  sph_sig,    sph_gam"
print mth,unr_alf,unr_bet,rhoi,rhol,Dth,unr_sig,unr_gam,sph_sig,sph_gam

#add other used parameter here
unr_alf_Jaggdent = 0.01243; unr_bet_Jaggdent = 2.00000 #found commented in McSnows mo_mass_diam.f90
unr_sig_Jaggdent = 0.05625; unr_gam_Jaggdent = 1.81000
print "##################"
print "added parameter"
print "unr_alf_Jaggdent,unr_bet_Jaggdent,unr_sig_Jaggdent,unr_gam_Jaggdent"
print unr_alf_Jaggdent,"         ",unr_bet_Jaggdent,"             ",unr_sig_Jaggdent,"           ",unr_gam_Jaggdent
#calculate Dth and mth for the added coefficients
D_th_Jaggdent=(sph_alf/unr_alf_Jaggdent)**(1./(unr_bet_Jaggdent-3.))
m_th_Jaggdent=sph_alf*D_th_Jaggdent**3
print "D_th_Jaggdent,m_th_Jaggdent"
print D_th_Jaggdent,m_th_Jaggdent

#get parameters from the SB-categories
cloud_water,rain,cloud_ice,snow,graupel,hail = __postprocess_SB.init_particles()
cloud_water,rain,cloud_ice,snow,graupel,hail = __postprocess_SB.convert_Nm_to_ND(cloud_water,rain,cloud_ice,snow,graupel,hail)

#initialize array with masses
m_MC = np.zeros(D_array.shape)

#calculate the masses for the different m-D regions of McSnow
for i, D in enumerate(D_array):
    if D<=Dth:
        m_MC[i] = np.pi/6.*rhoi*D**3
    elif D>Dth:
        m_MC[i] = unr_alf*D**unr_bet
#calculate the masses for the SB-categories (we do not need a case distinction as for McSnow)
m_SB_ice = cloud_ice.a_ms*D_array**cloud_ice.b_ms
m_SB_snow = snow.a_ms*D_array**snow.b_ms
m_Jaggdent = unr_alf_Jaggdent*D_array**unr_bet_Jaggdent

#increase font sizes
params = {'legend.fontsize': 'large',
    'figure.figsize': (15, 5),
    'axes.labelsize': 'x-large', #size: Either a relative value of 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large' or an absolute font size, e.g., 12
    'axes.titlesize':'x-large',
    'xtick.labelsize':'x-large',
    'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

#define figure
number_of_plots = 6
figsize_height = 6.0/2.0*number_of_plots
fig, axes = plt.subplots(nrows=number_of_plots, ncols=1, figsize=(8.0,figsize_height))
###
#subplot 1: m vs D
###
#axes2 = axes[0].twinx()
#plot in loglog
axes[0].loglog(D_array,m_MC,linestyle='--',label="McSnow")
axes[0].loglog(D_array,m_SB_ice,color='blue',label="SB cloud ice")
axes[0].loglog(D_array,m_SB_snow,color='green',label="SB snow")
axes[0].loglog(D_array[D_array>D_th_Jaggdent],m_Jaggdent[D_array>D_th_Jaggdent],color='red',linestyle='-.',label="Jaggdent")
#axes2.semilogx(D_array,(m_SB_ice-m_MC)/m_MC,color='red',label="|SB cloud ice - McSnow|/MC_snow")
#define labels
axes[0].set_xlabel("diameter / m")
#axes2.set_xlabel("diameter / m")
axes[0].set_ylabel("mass / kg")
#axes2.set_ylabel("relative diff. of \n SB cloud ice to McSnow")
#set limits
#axes2.set_ylim([-1,3])
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

#calculate the masses for the different m-D regions of McSnow
for i, m in enumerate(m_array):
    if m<=mth:
        D_MC[i] = (m/(sph_alf))**(1./sph_bet)
    elif m>mth:
        D_MC[i] = (m/unr_alf)**(1./unr_bet)
#calculate the masses for the SB-categories (we do not need a case distinction as for McSnow)
D_SB_ice = (m_array/cloud_ice.a_ms)**(1./cloud_ice.b_ms)
D_SB_snow = (m_array/snow.a_ms)**(1./snow.b_ms)
D_Jaggdent = (m_array/unr_alf_Jaggdent)**(1./unr_bet_Jaggdent);
#from IPython.core.debugger import Tracer ; Tracer()()

#plot subplot 2
#axes22 = axes[1].twinx()
#plot in loglog
axes[1].loglog(m_array,D_MC,linestyle='--',label="McSnow")
axes[1].loglog(m_array,D_SB_ice,color='blue',label="SB cloud ice")
axes[1].loglog(m_array,D_SB_snow,color='green',label="SB snow")
axes[1].loglog(m_array[m_array>m_th_Jaggdent],D_Jaggdent[m_array>m_th_Jaggdent],color='red',linestyle='-.',label="Jaggdent")

#axes22.semilogx(m_array,(D_SB_ice-D_MC)/D_MC,color='red',label="|SB cloud ice - McSnow|/MC_snow")
#define labels
axes[1].set_xlabel("mass / kg")
#axes22.set_xlabel("mass / kg")
axes[1].set_ylabel("diameter / m")
#axes22.set_ylabel("relative diff. of \n SB cloud ice to McSnow")
#set limits
#axes2.set_ylim([-1,5])
#show legend
axes[1].legend()
#axes2.legend()
plt.tight_layout()

###
#subplot 3: A vs D
###
#axes2 = axes[0].twinx()
#plot in loglog
A_MC = np.zeros(m_array.shape)

for i, D in enumerate(D_array):
    if D<=Dth:
        A_MC[i] = sph_sig*D**sph_gam
    elif D>Dth:
        A_MC[i] = unr_sig*D**unr_gam
A_Jaggdent = unr_sig_Jaggdent*D_array**unr_gam_Jaggdent
axes[2].loglog(D_array,A_MC,linestyle='--',label="McSnow")
#axes[2].loglog(D_array,m_SB_ice,color='blue',label="SB cloud ice")
#axes[2].loglog(D_array,m_SB_snow,color='green',label="SB snow")
axes[2].loglog(D_array[D_array>D_th_Jaggdent],A_Jaggdent[D_array>D_th_Jaggdent],color='red',linestyle='-.',label="Jaggdent")
#axes2.semilogx(D_array,(m_SB_ice-m_MC)/m_MC,color='red',label="|SB cloud ice - McSnow|/MC_snow")
#define labels
axes[2].set_xlabel("diameter / m")
#axes2.set_xlabel("diameter / m")
axes[2].set_ylabel("area / m2")
#axes2.set_ylabel("relative diff. of \n SB cloud ice to McSnow")
#set limits
#axes2.set_ylim([-1,3])
#show legend
axes[2].legend()
#axes2.legend()
plt.tight_layout()

###
#subplot 4: A vs m
###
#axes2 = axes[0].twinx()
#plot in loglog
Am_MC = np.zeros(m_array.shape)

for i, m in enumerate(m_array):
    if m<=mth:
        Am_MC[i] = sph_sig*(m/sph_alf)**(sph_gam/sph_bet)
    elif m>mth:
        Am_MC[i] = unr_sig*(m/unr_alf)**(unr_gam/unr_bet)
Am_Jaggdent = unr_sig_Jaggdent*(m_array/unr_alf_Jaggdent)**(unr_gam_Jaggdent/unr_bet_Jaggdent)
axes[3].loglog(m_array,Am_MC,linestyle='--',label="McSnow")
#axes[3].loglog(D_array,m_SB_ice,color='blue',label="SB cloud ice")
#axes[3].loglog(D_array,m_SB_snow,color='green',label="SB snow")
axes[3].loglog(m_array[m_array>m_th_Jaggdent],Am_Jaggdent[m_array>m_th_Jaggdent],color='red',linestyle='-.',label="Jaggdent")
#axes2.semilogx(D_array,(m_SB_ice-m_MC)/m_MC,color='red',label="|SB cloud ice - McSnow|/MC_snow")
#define labels
axes[3].set_xlabel("mass / kg")
#axes2.set_xlabel("diameter / m")
axes[3].set_ylabel("area / m2")
#axes2.set_ylabel("relative diff. of \n SB cloud ice to McSnow")
#set limits
#axes[3].set_ylim([0,2])
#show legend
axes[3].legend()
#axes2.legend()
plt.tight_layout()

###
#subplot 5: v vs D
###
#axes2 = axes[0].twinx()
#plot in loglog
A_MC = np.zeros(m_array.shape)
v_SB_ice = np.zeros(m_array.shape)
v_SB_snow = np.zeros(m_array.shape)

#parameters from SB-scheme
a_vel_icecosmo5 = 2.77e1; b_vel_icecosmo5 = 0.215790
a_vel_ulis_ice = 2.60e1; b_vel_ulis_ice = 0.215790
a_vel_snowSBB = 8.294000; b_vel_snowSBB = 0.125000
a_vel_icecosmo5nonsphere = 1.860; b_vel_icecosmo5nonsphere = 1.872; c_vel_icecosmo5nonsphere = 9.325e2 #fits to Ulis ice (despite the name) #this is as a function of the molten diameter
a_vel_snowSBBnonsphere = 1.271; b_vel_snowSBBnonsphere = 1.252; c_vel_snowSBBnonsphere = 3.698e3


v_SB_ice_cosmo5 = a_vel_icecosmo5*(cloud_ice.a_ms*D_array**cloud_ice.b_ms)**(b_vel_icecosmo5)
v_SB_ulis_ice = a_vel_ulis_ice*(cloud_ice.a_ms*D_array**cloud_ice.b_ms)**(b_vel_ulis_ice)
v_SB_snow = a_vel_snowSBB*(snow.a_ms*D_array**snow.b_ms)**(b_vel_snowSBB)
#define a new array for the molten diameter (which goes into the Atlas-type formulation)
D_array_molten = np.logspace(-12,-1,1001)#(cloud_ice.a_ms*D_array**cloud_ice.b_ms*6./(1000.*np.pi))**(1./.3)
v_SB_icecosmo5nonsphere = a_vel_icecosmo5nonsphere-b_vel_icecosmo5nonsphere*np.exp(-c_vel_icecosmo5nonsphere*D_array_molten) #these are formulated as a melted diameter!!
v_SB_snowSBBnonsphere = a_vel_snowSBBnonsphere-b_vel_snowSBBnonsphere*np.exp(-c_vel_snowSBBnonsphere*D_array_molten)

#get the diameter as a maximum diamension from the molten diameter
m_from_moltenD = 1000.*np.pi/6.*D_array_molten**3
D_max_from_moltenD_ice = cloud_ice.a_geo*m_from_moltenD**cloud_ice.b_geo
D_max_from_moltenD_snow = snow.a_geo*m_from_moltenD**snow.b_geo

#m_for_Darray_ulis_ice #this is m(D) for the v(m)
#from IPython.core.debugger import Tracer ; Tracer()()
#axes[4].loglog(D_array,A_MC,linestyle='--',label="McSnow")
axes[4].semilogx(D_array,v_SB_ice_cosmo5,color='blue',label="SB cloud ice_cosmo5")
axes[4].semilogx(D_array,v_SB_ulis_ice,color='blue',linestyle="--",label="SB cloud ulis_ice")
axes[4].semilogx(D_array,v_SB_snow,color='green',label="SB snow")
axes[4].semilogx(D_max_from_moltenD_ice,v_SB_icecosmo5nonsphere,color='blue',linestyle="-.",label="SB icecosmo5nonsphere")
axes[4].semilogx(D_max_from_moltenD_snow,v_SB_snowSBBnonsphere,color='green',linestyle="-.",label="SB snowSBBnonsphere")
#axes2.semilogx(D_array,(m_SB_ice-m_MC)/m_MC,color='red',label="|SB cloud ice - McSnow|/MC_snow")
#define labels
axes[4].set_xlabel("diameter / m")
#axes2.set_xlabel("diameter / m")
axes[4].set_ylabel("fall speed / kg")
#axes2.set_ylabel("relative diff. of \n SB cloud ice to McSnow")
#set limits
axes[4].set_xlim([1e-5,5e-2])
axes[4].set_ylim([0,2])
#show legend
axes[4].legend()
#axes2.legend()
plt.tight_layout()

###
#subplot 6: v vs m
###
#axes2 = axes[0].twinx()
#plot in loglog
Am_MC = np.zeros(m_array.shape)
#uncommented because we dont need them again
#a_vel_icecosmo5 = 2.77e1; b_vel_icecosmo5 = 0.215790
#a_vel_ulis_ice = 2.60e1; b_vel_ulis_ice = 0.215790
#a_vel_snowSBB = 8.294000; b_vel_snowSBB = 0.125000
#a_vel_icecosmo5nonsphere = 1.860; b_vel_icecosmo5nonsphere = 1.872; c_vel_icecosmo5nonsphere = 9.325e2 #fits to Ulis ice (despite the name)
#a_vel_snowSBBnonsphere = 1.271; b_vel_snowSBBnonsphere = 1.252; c_vel_snowSBBnonsphere = 9.325e2 #3.698e3 
###
#power law
###
#ice
vm_SB_ice_cosmo5 = a_vel_icecosmo5*m_array**b_vel_icecosmo5
vm_SB_ulis_ice = a_vel_ulis_ice*m_array**b_vel_ulis_ice
#snow
vm_SB_snow = a_vel_snowSBB*m_array**b_vel_snowSBB
###
#Atlas-type
###
vm_SB_icecosmo5nonsphere = a_vel_icecosmo5nonsphere-b_vel_icecosmo5nonsphere*np.exp(-c_vel_icecosmo5nonsphere*D_array) #these are formulated as a molten diameter!!
vm_SB_snowSBBnonsphere = a_vel_snowSBBnonsphere-b_vel_snowSBBnonsphere*np.exp(-c_vel_snowSBBnonsphere*D_array)

#get mass from molten diameter
m_from_moltenD = 1000.*np.pi/6.*D_array**3
#axes[5].loglog(D_array,A_MC,linestyle='--',label="McSnow")
axes[5].semilogx(m_array,vm_SB_ice_cosmo5,color='blue',label="SB ice_cosmo5")
axes[5].semilogx(m_array,vm_SB_ulis_ice,color='blue',linestyle='--',label="SB ulis_ice")
axes[5].semilogx(m_array,vm_SB_snow,color='green',label="SB snow")
axes[5].semilogx(m_from_moltenD,vm_SB_icecosmo5nonsphere,color='blue',linestyle="-.",label="SB icecosmo5nonsphere")
axes[5].semilogx(m_from_moltenD,vm_SB_snowSBBnonsphere,color='green',linestyle="-.",label="SB snowSBBnonsphere")
#axes2.semilogx(D_array,(m_SB_ice-m_MC)/m_MC,color='red',label="|SB cloud ice - McSnow|/MC_snow")
#define labels
axes[5].set_xlabel("mass / m")
#axes2.set_xlabel("diameter / m")
axes[5].set_ylabel("fall speed / kg")
#axes2.set_ylabel("relative diff. of \n SB cloud ice to McSnow")
#set limits
axes[5].set_xlim([1e-12,1e-2])
axes[5].set_ylim([0,2])
#show legend
axes[5].legend()
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