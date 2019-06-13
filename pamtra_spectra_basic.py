'''
this script is a basic example to generate idealized spectra with PAMTRA
'''

import pyPamtra
import numpy as np
import matplotlib.pyplot as plt
import sys
from subprocess import call
#from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug

def run_pamtra(which_mode="12"):
    '''
    #run pamtra with switchin on and off different categories
    INPUT: which_mode: if string contains 1->1st mode active; if contains 2-> 2nd mode active
    OUTPUT: pam_output: dictionary with moments and spectra of pamtra calculations
    '''
    #Initialise pamtra
    pam = pyPamtra.pyPamtra()

    #create some artificial athmosphere
    pam = pyPamtra.importer.createUsStandardProfile(pam,hgt_lev=[100,200])

    #calculate the spectra with PAMTRA
    #see https://pamtra.readthedocs.io/en/latest/descriptorFile.html table on bottom for more infos
    pam.df.addHydrometeor(("mode1", 0.6, -1, -99, 0.0772, 2.2284, 0.3971,1.88, 13, 100, "exp", -99.,    -99.,     -99.,    -99.,  -99,    -99,     "ss-rayleigh-gans",  "corAtlas_1.356_1.356_1281",    -99.))
    pam.df.addHydrometeor(("mode2", 0.6, -1, -99.,0.0772, 2.2284 ,0.3971,1.88, 13, 100, "exp", -99.,    -99.,     -99.,    -99.,  -99,    -99,     "ss-rayleigh-gans",  "corAtlas_1.356_1.356_1281",    -99.))


    pam.set["verbose"] = 0

    #specify hydrometeor bulk properties (this defines your size distribution)
    #number concentration [#/m-3]
    num_conc_mode1 = 1e6
    num_conc_mode2 = 1e3
    #mean mass [kg]
    mean_mass_mode1 = 1e-10
    mean_mass_mode2 = 1e-8
    
    #deactivate modes    
    if not "1" in which_mode:
        num_conc_mode1 = 0
    if not "2" in which_mode:
        num_conc_mode2 = 0


    #specify hydrometeor bulk properties (this defines your size distribution)
    pam.p["hydro_q"][0,0,0,0] = mean_mass_mode1*num_conc_mode1 #mode1
    pam.p["hydro_n"][0,0,0,0] = num_conc_mode1  #mode1
    pam.p["hydro_q"][0,0,0,1] = mean_mass_mode2*num_conc_mode2 #mode2
    pam.p["hydro_n"][0,0,0,1] = num_conc_mode2  #mode2
    #from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug
    pam.p["airturb"][:] = 0
    pam.p["wind_w"][:] =0.0
    pam.nmlSet["passive"] = False  # Activate this for Microwave radiometer
    pam.nmlSet["radar_mode"] = "spectrum"
    #pam.nmlSet["radar_pnoise0"]=-100 #set radar noise to an arbitrary low number
    pam.nmlSet["radar_noise_distance_factor"] = 2.0
    pam.nmlSet["radar_save_noise_corrected_spectra"]=  False
    pam.nmlSet["randomseed"]=  10
    pam.nmlSet['radar_airmotion'] = True
    pam.nmlSet['radar_airmotion_model'] = 'constant'
    pam.nmlSet['radar_airmotion_vmax'] = 0.0
    pam.nmlSet['radar_airmotion_vmin'] = 0.0
    pam.nmlSet["radar_nfft"] = 1024
    pam.nmlSet['radar_aliasing_nyquist_interv'] = 1
    pam.nmlSet['radar_no_Ave'] = 15 #15
    pam.nmlSet["save_psd"] = True
    pam.nmlSet["hydro_adaptive_grid"]=True
    #pam.nmlSet["radar_max_V"] = 3.
    #pam.nmlSet["radar_min_V"] = -3.

    pam.runParallelPamtra([9.6,35.5,95], pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers='auto')

    pam_output = dict() #dictionary of pamtra output with some more meaningful names than in pam.r
    #save the result of the radar moments into pam_output
    #reflectivity
    #comb_dic["Ze_X"] = np.array([pam.r["Ze"][0,0,1,0,0,0]])
    pam_output["Ze_X"] = pam.r["Ze"][0,0,0,0,0,0]
    pam_output["Ze_Ka"] = pam.r["Ze"][0,0,0,1,0,0]
    pam_output["Ze_W"] = pam.r["Ze"][0,0,0,2,0,0]
    pam_output["DWR_X_Ka"] = pam.r["Ze"][0,0,0,0,0,0]-pam.r["Ze"][0,0,0,1,0,0]
    pam_output["DWR_Ka_W"] = pam.r["Ze"][0,0,0,1,0,0]-pam.r["Ze"][0,0,0,2,0,0]


    #higher moments
    pam_output["vDoppler_X"] = pam.r["radar_moments"][0,0,0,0,0,0,0]
    pam_output["vDoppler_Ka"] = pam.r["radar_moments"][0,0,0,1,0,0,0]
    pam_output["vDoppler_W"] = pam.r["radar_moments"][0,0,0,2,0,0,0]
    pam_output["swidth_X"] = pam.r["radar_moments"][0,0,0,0,0,0,0]
    pam_output["swidth_Ka"] = pam.r["radar_moments"][0,0,0,1,0,0,1]
    pam_output["swidth_W"] = pam.r["radar_moments"][0,0,0,2,0,0,2]
    pam_output["skewn_X"] = pam.r["radar_moments"][0,0,0,0,0,0,2]
    pam_output["skewn_Ka"] = pam.r["radar_moments"][0,0,0,1,0,0,2]
    pam_output["skewn_W"] = pam.r["radar_moments"][0,0,0,2,0,0,2]

        
    #full spectra
    pam_output["psd_d"]         = pam.r['psd_d'][0,0,0,0]
    pam_output["psd_n_mode1"]   = np.maximum(0,pam.r['psd_n'][0,0,0,0]) #get rid of -9999 if deactivated
    pam_output["psd_n_mode2"]   = np.maximum(0,pam.r['psd_n'][0,0,0,1]) #get rid of -9999 if deactivated
    pam_output["psd_n_all"]     = pam_output["psd_n_mode1"]+pam_output["psd_n_mode2"]
    #from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug
    pam_output["psd_m_mode1"]   = np.maximum(0,pam.r['psd_mass'][0,0,0,0]) #get rid of -9999 if deactivated #this is the mass at the center of the bin
    pam_output["psd_m_mode2"]   = np.maximum(0,pam.r['psd_mass'][0,0,0,1]) #get rid of -9999 if deactivated
    pam_output["psd_m_all"]     = pam_output["psd_m_mode1"]+pam_output["psd_m_mode2"]

    pam_output["velbins_X"]  = pam.r["radar_vel"][0]
    pam_output["velbins_Ka"] = pam.r["radar_vel"][1]
    pam_output["velbins_W"]  = pam.r["radar_vel"][2]
    pam_output["spectra_X"]  = pam.r["radar_spectra"][0,0,0,0,0]
    pam_output["spectra_Ka"] = pam.r["radar_spectra"][0,0,0,1,0]
    pam_output["spectra_W"]  = pam.r["radar_spectra"][0,0,0,2,0]

    return pam_output

#run pamtra with both categories
pam_output = run_pamtra(which_mode="12")

# Initialise the spectra plot
figsize_height = 8.0/2.0*(2.)
figspec = plt.figure(figsize=(8.0,figsize_height))
#plot the size distributions of both categories
#N(D)
axsizedist = plt.subplot(311)
axsizedist.loglog(pam_output["psd_d"],pam_output["psd_n_all"],color='r',label='all')
axsizedist.loglog(pam_output["psd_d"],pam_output["psd_n_mode1"],color='orange',linestyle='--',marker='o',markersize=6,markevery=3,label='mode1')
axsizedist.loglog(pam_output["psd_d"],pam_output["psd_n_mode2"],color='magenta',linestyle='-.',marker='*',markersize=6,markevery=3,label='mode2')
plt.legend()
plt.xlabel('diameter m')
plt.ylabel('number concentration / m-4')

#m(D)
axsizedist = plt.subplot(312)
axsizedist.loglog(pam_output["psd_d"],pam_output["psd_m_mode1"]*pam_output["psd_n_mode1"]+pam_output["psd_m_mode2"]*pam_output["psd_n_mode2"],color='r',label='all')
axsizedist.loglog(pam_output["psd_d"],pam_output["psd_m_mode1"]*pam_output["psd_n_mode1"],color='orange',linestyle='--',marker='o',markersize=6,markevery=3)
axsizedist.loglog(pam_output["psd_d"],pam_output["psd_m_mode2"]*pam_output["psd_n_mode2"],color='magenta',linestyle='-.',marker='*',markersize=6,markevery=3)
axsizedist.plot(np.nan,np.nan,color='orange',linestyle='--',label='mode1',marker='o')
axsizedist.plot(np.nan,np.nan,color='magenta',linestyle='-.',label='mode2',marker='*')
plt.legend()
plt.xlabel('diameter m')
plt.ylabel('mass concentration / kg m-4')

#plot the Doppler spectra in all three frequencies
axspec = plt.subplot(313)
axspec.plot(pam_output["velbins_X"],pam_output["spectra_X"],label='9.5GHz',color='b',linestyle='-')
axspec.plot(pam_output["velbins_Ka"],pam_output["spectra_Ka"],label='35GHz',color='r',linestyle='-')
axspec.plot(pam_output["velbins_W"],pam_output["spectra_W"],label='95GHz',color='g',linestyle='-')

#save moments in string to add them in plot
moments_string = "Ze: " + str(pam_output["Ze_Ka"]) + "\n" #save all moments in one string
moments_string+="vDoppler_Ka: " + str(pam_output["vDoppler_Ka"]) + "\n"
moments_string+="swidth_Ka: " + str(pam_output["swidth_Ka"]) + "\n"
moments_string+="skewn_Ka: " + str(pam_output["skewn_Ka"]) + "\n"

#run pamtra with the first mode only
pam_output = run_pamtra(which_mode="1")
#plot the spectra in all three frequencies if there is only mode1
axspec.plot(pam_output["velbins_X"],pam_output["spectra_X"],label='__None',color='b',linestyle='--',marker='o',markerfacecolor='k',markeredgecolor='k',markersize=6,markevery=3)

#run pamtra with the first mode only
pam_output = run_pamtra(which_mode="2")
#plot the spectra in all three frequencies if there is only mode2
#from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug
axspec.plot(pam_output["velbins_X"],pam_output["spectra_X"],label='__None',color='b',linestyle='-.',marker='*',markerfacecolor='k',markeredgecolor='k',markersize=6,markevery=3)

#adjust the plot range
axspec.set_ylim([-40,20]); axspec.set_xlim([-0.2,2])

plt.xlabel('Doppler velocity m s-1')
plt.ylabel('spectral reflectivity dB')
#add invisible line for labelling mode1 and mode1
axspec.plot(np.nan,np.nan,color='k',linestyle='--',label='mode1',marker='o')
axspec.plot(np.nan,np.nan,color='k',linestyle='-.',label='mode2',marker='*')



#add string of moments
axspec.text(0,1.0,moments_string,
     horizontalalignment='left',
     verticalalignment='top',
     transform = axspec.transAxes) #put text in upper left corner

#create legend
plt.legend()
plt.tight_layout()
output_file = "/home/mkarrer/Dokumente/pamtrafiles/ideal_spectra/spectra.png"
plt.savefig(output_file,dpi=400)
print "the spectra is at: " + output_file
plt.clf()
sys.exit()