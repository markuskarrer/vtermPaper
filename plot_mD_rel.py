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
import __fallspeed_relations
import __setup_mDAD_frommodels
#from IPython.core.debugger import Tracer ; Tracer()()

'''
plot the mass over diameter for the SB and the McSnow scheme
'''

#get shared diameter array
D_array = __setup_mDAD_frommodels.setup_diam_arrays()


#get model parameters
#McSnow
mDADvD_dict_MC = __setup_mDAD_frommodels.get_model_mDADs(model="MC")
mDADvD_dict_MCJagg_old = __setup_mDAD_frommodels.get_model_mDADs(model="MCJagg_old")
#SB
mDADvD_dict_SBcloudice = __setup_mDAD_frommodels.get_model_mDADs(model="SBcloudice")
mDADvD_dict_SBcloudice_uli = __setup_mDAD_frommodels.get_model_mDADs(model="SBcloudice_uli")
mDADvD_dict_SBcloudice_Atlas = __setup_mDAD_frommodels.get_model_mDADs(model="SBcloudice_Atlas")
mDADvD_dict_SBsnow = __setup_mDAD_frommodels.get_model_mDADs(model="SBsnow")
mDADvD_dict_SBsnow_Atlas = __setup_mDAD_frommodels.get_model_mDADs(model="SBsnow_Atlas")
#P3
mDADvD_dict_P3 = __setup_mDAD_frommodels.get_model_mDADs(model="P3")
#from Jussis aggregate model
mDADvD_dict_Jneedles = __setup_mDAD_frommodels.get_model_mDADs(model="Jussis_needles")

#calculate the arrays with masses and areas corresponding to the common D_array and based on the (piecewise) power-law fits

for dict_now in (mDADvD_dict_MC,mDADvD_dict_MCJagg_old,mDADvD_dict_P3,mDADvD_dict_SBcloudice,mDADvD_dict_SBcloudice_uli,
                 mDADvD_dict_SBcloudice_Atlas,mDADvD_dict_SBsnow,mDADvD_dict_SBsnow_Atlas,mDADvD_dict_Jneedles):
    dict_now = __setup_mDAD_frommodels.calc_area_mass_vterm_arrays(D_array,dict_now)

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
#axes[0].loglog(D_array,mDADvD_dict_MCJagg_old["m(D_array)"],color='red',linestyle='-.',label="Jaggdent")
axes[0].loglog(D_array,mDADvD_dict_MC["m(D_array)"],linestyle='--',label="McSnow")
axes[0].loglog(D_array,mDADvD_dict_SBcloudice["m(D_array)"],color='blue',label="SB cloud ice")
axes[0].loglog(D_array,mDADvD_dict_SBsnow["m(D_array)"],color='green',label="SB snow")
axes[0].loglog(D_array,mDADvD_dict_P3["m(D_array)"],color='orange',linestyle='-.',label="P3unrimed")
axes[0].loglog(D_array,mDADvD_dict_Jneedles["m(D_array)"],color='red',linestyle='-.',label="Jussi's needles")
#axes[0].loglog(D_array[D_array>mDADvD_dict_P3["Dth1"]],mDADvD_dict_P3["m(D_array)"][D_array>mDADvD_dict_P3["Dth1"]],color='orange',linestyle='-.',label="P3unrimedagg")

#axes2.semilogx(D_array,(mDADvD_dict_SBcloudice["m(D_array)"]-mDADvD_dict_MC["m(D_array)"])/mDADvD_dict_MC["m(D_array)"],color='red',label="|SB cloud ice - McSnow|/MC_snow")
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
'''
###
#calculate D for subplot 2: D vs m
###
#set up mass array
n_steps = 1001
min_m = 1e-13; max_m = 1e-4
m_array = np.logspace(np.log10(min_m),np.log10(max_m),n_steps)


#initialize array with masses
D_MC = np.zeros(m_array.shape)
D_P3 = np.zeros(m_array.shape)

#calculate the diameters for the different m-D regions of McSnow
for i, m in enumerate(m_array):
    if m<=mDADvD_dict_MC["mth1"]:
        D_MC[i] = (m/(mDADvD_dict_MC["a1"]))**(1./mDADvD_dict_MC["b1"])
    elif m>mDADvD_dict_MC["mth1"]:
        D_MC[i] = (m/mDADvD_dict_MC["a2"])**(1./mDADvD_dict_MC["b2"])
#calculate the diameters for the different m-D regions of P3
for i, m in enumerate(m_array):
    if m<=mDADvD_dict_P3["mth1"]:
        D_P3[i] = (m/(mDADvD_dict_P3["a1"]))**(1./mDADvD_dict_P3["b1"])
    elif m>mDADvD_dict_P3["mth1"]:
        D_P3[i] = (m/(mDADvD_dict_P3["a2"]))**(1./mDADvD_dict_P3["b2"])

#calculate the masses for the SB-categories (we do not need a case distinction as for McSnow)
D_SB_ice = mDADvD_dict_SBcloudice["am1"]*m_array**mDADvD_dict_SBcloudice["bm1"]#(m_array/mDADvD_dict_SBcloudice["a1"])**(1./mDADvD_dict_SBcloudice["b1"])
D_SB_snow = (m_array/mDADvD_dict_SBsnow['a1'])**(1./mDADvD_dict_SBsnow['b1'])
D_Jaggdent = (m_array/mDADvD_dict_MCJagg_old["a2"])**(1./mDADvD_dict_MCJagg_old["b2"]);

#plot subplot 2
#axes22 = axes[1].twinx()
#plot in loglog
axes[1].loglog(m_array,D_MC,linestyle='--',label="McSnow")
axes[1].loglog(m_array,D_SB_ice,color='blue',label="SB cloud ice")
axes[1].loglog(m_array,D_SB_snow,color='green',label="SB snow")
axes[1].loglog(m_array[m_array>mDADvD_dict_MCJagg_old["mth1"]],D_Jaggdent[m_array>mDADvD_dict_MCJagg_old["mth1"]],color='red',linestyle='-.',label="Jaggdent")
axes[1].loglog(m_array,D_P3,color='orange',label="P3 unrimed",linestyle='-.')

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
'''
###
#subplot 3: A vs D
###

#axes[2].loglog(D_array,mDADvD_dict_MCJagg_old["A(D_array)"],color='red',linestyle='-.',label="Jaggdent")
axes[2].loglog(D_array,mDADvD_dict_MC["A(D_array)"],linestyle='--',label="McSnow")
axes[2].loglog(D_array,mDADvD_dict_P3["A(D_array)"],color='orange',linestyle='-.',label="P3unrimed")
axes[2].loglog(D_array,mDADvD_dict_Jneedles["A(D_array)"],color='red',linestyle='-.',label="Jussi's needles")

#define labels
axes[2].set_xlabel("diameter / m")
#axes2.set_xlabel("diameter / m")
axes[2].set_ylabel("area / m2")
#set limits
#axes2.set_ylim([-1,3])
#show legend
axes[2].legend()
#axes2.legend()
plt.tight_layout()
'''
###
#subplot 4: A vs m
###


#axes2 = axes[0].twinx()
#plot in loglog
Am_MC = np.zeros(m_array.shape)
Am_P3 = np.zeros(m_array.shape)

for i, m in enumerate(m_array):
    if m<=mDADvD_dict_MC["mth1"]:
        Am_MC[i] = mDADvD_dict_MC["c1"]*(m/mDADvD_dict_MC["a1"])**(mDADvD_dict_MC["d1"]/mDADvD_dict_MC["b1"])
    elif m>mDADvD_dict_MC["mth1"]:
        Am_MC[i] = mDADvD_dict_MC["c2"]*(m/mDADvD_dict_MC["a2"])**(mDADvD_dict_MC["d2"]/mDADvD_dict_MC["b2"])
for i, m in enumerate(m_array):
    if m<=mDADvD_dict_P3["mth1"]:
        Am_P3[i] = mDADvD_dict_P3["c1"]*(m/mDADvD_dict_P3["a1"])**(mDADvD_dict_P3["d1"]/mDADvD_dict_P3["b1"])
    elif m>mDADvD_dict_P3["mth1"]:
        Am_P3[i] = mDADvD_dict_P3["c2"]*(m/mDADvD_dict_P3["a2"])**(mDADvD_dict_P3["d2"]/mDADvD_dict_P3["b2"])
        
        
AmDADvD_dict_MCJagg_old["m(D_array)"] = mDADvD_dict_MCJagg_old["c2"]*(m_array/mDADvD_dict_MCJagg_old["a2"])**(mDADvD_dict_MCJagg_old["d2"]/mDADvD_dict_MCJagg_old["b2"])

axes[3].loglog(m_array,Am_MC,linestyle='--',label="McSnow")
#axes[3].loglog(D_array,mDADvD_dict_SBcloudice["m(D_array)"],color='blue',label="SB cloud ice")
#axes[3].loglog(D_array,mDADvD_dict_SBsnow["m(D_array)"],color='green',label="SB snow")
axes[3].loglog(m_array[m_array>mDADvD_dict_MCJagg_old["mth1"]],AmDADvD_dict_MCJagg_old["m(D_array)"][m_array>mDADvD_dict_MCJagg_old["mth1"]],color='red',linestyle='-.',label="Jaggdent")
#axes2.semilogx(D_array,(mDADvD_dict_SBcloudice["m(D_array)"]-mDADvD_dict_MC["m(D_array)"])/mDADvD_dict_MC["m(D_array)"],color='red',label="|SB cloud ice - McSnow|/MC_snow")
axes[3].loglog(m_array,Am_P3,linestyle='-.',color='orange',label="P3unrimed")

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

#calculate v vs D for unrimed McSnow (with constants from Brdar&Seifert2018 and KC05)
v_MC_KC05 = __fallspeed_relations.calc_vterm("KC05",mDADvD_dict_MC["m(D_array)"],D_array,mDADvD_dict_MC["A(D_array)"])
#calculate v vs D for unrimed P3
v_P3_mitchheym = __fallspeed_relations.calc_vterm("mitch_heym",mDADvD_dict_P3["m(D_array)"],D_array,mDADvD_dict_P3["A(D_array)"])
'''
###
#subplot 5: v vs D
###

axes[4].semilogx(D_array,mDADvD_dict_SBcloudice["v(D_array)"],color='blue',label="SB cloud ice_cosmo5")
#axes[4].semilogx(mDADvD_dict_SBcloudice_Atlas["D_max_from_moltenD"],mDADvD_dict_SBcloudice_Atlas["v(D_array)"],color='blue',linestyle="-.",label="SB icecosmo5nonsphere")
#axes[4].semilogx(D_array,mDADvD_dict_SBcloudice_uli["v(D_array)"],color='blue',linestyle="--",label="SB ulis_ice")

axes[4].semilogx(D_array,mDADvD_dict_SBsnow["v(D_array)"],color='green',label="SB snow")
#axes[4].semilogx(mDADvD_dict_SBsnow_Atlas["D_max_from_moltenD"],mDADvD_dict_SBsnow_Atlas["v(D_array)"],color='green',linestyle="-.",label="SB snowSBBnonsphere")
#axes[4].semilogx(D_array,v_SB_ulis_ice,color='blue',linestyle="--",label="SB cloud ulis_ice")
axes[4].semilogx(D_array,mDADvD_dict_MC["v_HW10(D_array)"],color='red',label="McSnow(HW10)")
axes[4].semilogx(D_array,mDADvD_dict_P3["v_mitch_heym(D_array)"],color='orange',label="P3unrimed",linestyle='-.')


#axes2.semilogx(D_array,(mDADvD_dict_SBcloudice["m(D_array)"]-mDADvD_dict_MC["m(D_array)"])/mDADvD_dict_MC["m(D_array)"],color='red',label="|SB cloud ice - McSnow|/MC_snow")
#define labels
axes[4].set_xlabel("diameter / m")
#axes2.set_xlabel("diameter / m")
axes[4].set_ylabel("terminal velocity / m s-1")
#axes2.set_ylabel("relative diff. of \n SB cloud ice to McSnow")
#set limits
axes[4].set_xlim([1e-4,1e-2])#[1e-5,5e-2])
axes[4].set_ylim([0,2])
#show legend
axes[4].legend()
#axes2.legend()
plt.tight_layout()
'''
###
#subplot 6: v vs m
###
#axes2 = axes[0].twinx()
#plot in loglog
#mDADvD_dict_MC["m(D_array)"] = np.zeros(m_array.shape)

###
#power law
###
#ice
#vm_SB_cloud_cloud_ice_cosmo5 = mDADvD_dict_SBcloudice["aterm"]*m_array**mDADvD_dict_SBcloudice["bterm"]
#vm_SB_ulis_ice = mDADvD_dict_SBcloudice_uli["aterm"]*m_array**mDADvD_dict_SBcloudice_uli["bterm"]
#snow
mDADvD_dict_SBsnow["v(D_array)"] = mDADvD_dict_SBsnow["aterm"]*m_array**mDADvD_dict_SBsnow["bterm"]
###
#Atlas-type
###
mDADvD_dict_SBcloudice_Atlas["v(D_array)"] = mDADvD_dict_SBcloudice_Atlas["aterm"]-mDADvD_dict_SBcloudice_Atlas["bterm"]*np.exp(-mDADvD_dict_SBcloudice_Atlas["cterm"]*D_array) #these are formulated as a molten diameter!!
mDADvD_dict_SBsnow_Atlas["v(D_array)"] = mDADvD_dict_SBsnow_Atlas["aterm"]-mDADvD_dict_SBsnow_Atlas["bterm"]*np.exp(-mDADvD_dict_SBsnow_Atlas["cterm"]*D_array)

#get mass from molten diameter
m_from_moltenD = 1000.*np.pi/6.*D_array**3
#axes[5].loglog(D_array,mDADvD_dict_MC["A(D_array)"],linestyle='--',label="McSnow")
axes[5].semilogx(m_array,mDADvD_dict_SBcloudice["v(D_array)"],color='blue',label="SB ice_cosmo5")
#axes[5].semilogx(m_array,vm_SB_ulis_ice,color='blue',linestyle='--',label="SB ulis_ice")
axes[5].semilogx(m_array,mDADvD_dict_SBsnow["v(D_array)"],color='green',label="SB snow")
#axes[5].semilogx(m_from_moltenD,vmDADvD_dict_SBcloudice["m(D_array)"]cosmo5nonsphere,color='blue',linestyle="-.",label="SB icecosmo5nonsphere")
#axes[5].semilogx(m_from_moltenD,vmDADvD_dict_SBsnow["m(D_array)"]SBBnonsphere,color='green',linestyle="-.",label="SB snowSBBnonsphere")
#axes2.semilogx(D_array,(mDADvD_dict_SBcloudice["m(D_array)"]-mDADvD_dict_MC["m(D_array)"])/mDADvD_dict_MC["m(D_array)"],color='red',label="|SB cloud ice - McSnow|/MC_snow")
#define labels
axes[5].set_xlabel("mass / m")
#axes2.set_xlabel("diameter / m")
axes[5].set_ylabel("terminal velocity / m s-1")
#axes2.set_ylabel("relative diff. of \n SB cloud ice to McSnow")
#set limits
axes[5].set_xlim([1e-12,1e-2])
axes[5].set_ylim([0,2])
#show legend
axes[5].legend()
#axes2.legend()
plt.tight_layout()
'''
dir_save = '/home/mkarrer/Dokumente/plots/'
if not os.path.exists(dir_save): #create direktory if it does not exists
    os.makedirs(dir_save)
out_filestring = "mD_rel"
plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
plt.savefig(dir_save + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])