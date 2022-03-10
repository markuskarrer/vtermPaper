
# coding: utf-8
#import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.pylab as pylab
import os
import subprocess
#import other self-defined functions

#from IPython.core.debugger import Tracer ; Tracer()()


'''
plot different sticking efficiencies as they are used in McSnow (including embedded SB-1D model); .f90 files are all located in MCSNOW/src
'''

n_steps=1001
temp_array_C = np.linspace(-40,0,n_steps)
temp_array_K = temp_array_C+273.15
#initialize different sticking efficiencies e_... (names specify source (paper) and where in the schemes it is used)
e_lin83_snowself_particleparticle = np.zeros(temp_array_C.shape)
e_cotton86_iceself = np.zeros(temp_array_C.shape)
e_pruppklett98_MC = np.zeros(temp_array_C.shape)
e_pruppklett98_MC_mod = np.zeros(temp_array_C.shape)
e_connolly12_MC = np.zeros(temp_array_C.shape)
e_connolly12_0p5_MC = np.zeros(temp_array_C.shape)

for i, (temp_C,temp_K) in enumerate(zip(temp_array_C,temp_array_K)):
    #Lin #mo_2mom_mcrph_processes.f90 SUBROUTINE snow_selfcollection() e_coll & SUBROUTINE particle_particle_collection() e_coll
    e_lin83_snowself_particleparticle[i] = min(np.exp(0.09*(temp_C)),1.0) 
    #Cotton #mo_2mom_mcrph_processes.f90 SUBROUTINE ice_selfcollection() e_coll
    e_cotton86_iceself[i] = min(10**(0.035*(temp_C)-0.7),0.2)
    #Pruppacher & Klett  #mo_colleffi.f90 PURE FUNCTION stick_effi_tempOnly(T) se
    if  temp_C >= 0.0 :
        e_pruppklett98_MC[i] = 1.0
    elif  temp_C >= -4. :
        e_pruppklett98_MC[i] = 0.1
    elif  temp_C >= -6. :
        e_pruppklett98_MC[i] = 0.6
    elif  temp_C >= -9. :
        e_pruppklett98_MC[i] = 0.1
    elif  temp_C >= -12.5 :
        e_pruppklett98_MC[i] = 0.4
    elif  temp_C >= -17 :
        e_pruppklett98_MC[i] = 1.0
    elif  temp_C >= -20 :
        e_pruppklett98_MC[i] = 0.4
    elif  temp_C >= -30 :
        e_pruppklett98_MC[i] = 0.25
    elif  temp_C >= -40 :
        e_pruppklett98_MC[i] = 0.1
    else:
        e_pruppklett98_MC[i] = 0.02

    if temp_C >= -2 :
        e_pruppklett98_MC_mod[i] = 1.0
    elif  temp_C >= -6 :
        e_pruppklett98_MC_mod[i] = 0.19*(temp_C+6)+0.24
    elif  temp_C >= -10 :
        e_pruppklett98_MC_mod[i] = 0.24
    elif  temp_C >= -12.5 :
        e_pruppklett98_MC_mod[i] = -0.304*(temp_C+12.5)+1.0
    elif  temp_C >= -17:
        e_pruppklett98_MC_mod[i] = 1.0
    elif  temp_C >= -20 :
        e_pruppklett98_MC_mod[i] =  0.286666*(temp_C+20)+0.14
    elif  temp_C >= -40 :
        e_pruppklett98_MC_mod[i] = 0.005*(temp_C+40)+0.04
    else:
        e_pruppklett98_MC_mod[i] =  0.04
    #Connolly #mo_colleffi.f90 PURE FUNCTION stick_effi_tempOnly(T) se
    if temp_C >= 0 :
        e_connolly12_MC[i] = 0.14
    elif  temp_C >= -10 :
        e_connolly12_MC[i] = -0.01*(temp_C+10)+0.24
    elif  temp_C >= -15 :
        e_connolly12_MC[i] = -0.08*(temp_C+15)+0.64
    elif  temp_C >= -20 :
        e_connolly12_MC[i] =  0.10*(temp_C+20)+0.14
    elif  temp_C >= -40 :
        e_connolly12_MC[i] = 0.005*(temp_C+40)+0.04
    else:
        e_connolly12_MC[i] =  0.04
    #Connolly*0.5 (for istick=3) #mo_colleffi.f90 PURE FUNCTION stick_effi_tempOnly(T) se
    e_connolly12_0p5_MC[i] = 0.5*e_connolly12_MC[i]
        
#increase font sizes
params = {'legend.fontsize': 'large',
    'figure.figsize': (15, 5),
    'axes.labelsize': 'x-large', #size: Either a relative value of 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large' or an absolute font size, e.g., 12
    'axes.titlesize':'x-large',
    'xtick.labelsize':'x-large',
    'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

#get number of plot automatically
figsize_height = 6.0/2.0 #*(number_of_plots)
fig = plt.figure(figsize=(8.0,figsize_height))
plt.plot(temp_array_C,e_lin83_snowself_particleparticle,label="Lin83_SBsnowself_SBpartpart",color="r")

#plt.plot(temp_array_C,e_cotton86_iceself,label="Cotton86_SBiceself",color="brown")
#plt.plot(temp_array_C,e_pruppklett98_MC,label="PruppKlett98_MC",linestyle="-",color="g")
plt.plot(temp_array_C,e_pruppklett98_MC_mod,label="PruppKlett98_MC_mod",linestyle="--",color="darkgreen")
#plt.plot(temp_array_C,e_connolly12_MC,label="Connolly12_MC",linestyle="--",color="navy")
#plt.plot(temp_array_C,e_connolly12_0p5_MC,label="0.5*Connolly12_MC",linestyle="--",color="deepskyblue")
plt.ylim([0,1.01])
plt.xlabel("temperature / ($^\circ$C)")
plt.ylabel("sticking efficiency")
plt.legend()
#save figure
plt.tight_layout()
dir_save = '/home/mkarrer/Dokumente/plots/'
if not os.path.exists(dir_save): #create direktory if it does not exists
    os.makedirs(dir_save)
out_filestring = "stickeff_few"
plt.savefig(dir_save + out_filestring + '.pdf', dpi=400)
plt.savefig(dir_save + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + dir_save + out_filestring + '.pdf'
#subprocess.Popen(['evince',dir_save + out_filestring + '.pdf'])
