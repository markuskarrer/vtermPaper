import pyPamtra
import numpy as np
import matplotlib.pyplot as plt
import sys
from subprocess import call
#from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug

#functions which project categories of variables to the interval [0,1]
def project_mixrat_to_0to1(value,log_minval,log_maxval):
    number = (np.log10(value)-log_minval)/(log_maxval-log_minval)
    return number
def revert_0to1_to_mixrat(number,log_minval,log_maxval):
    mixrat_val = 10**((log_maxval-log_minval)*number+log_minval)
    return mixrat_val
#def project_numcon_to_0to1(value,log_minval,log_maxval):#ATTENTION: log_minval is not used (and assumed to be 0)
#    number = np.log10(value)/log_maxval #(np.log10(value)-log_minval)/(log_maxval-log_minval)
#    return number
#def revert_0to1_to_numcon(number,log_minval,log_maxval): #ATTENTION: log_minval is not used (and assumed to be 0)
#    mixrat_val = 10**(number*log_maxval)
#    return mixrat_val
def project_linear(value,minval,maxval): 
    number = (value-minval)/(maxval-minval)
    return number
def revert_linear(number,minval,maxval):
    value = number*(maxval-minval)+minval
    return value
def project_linear_symmetricaboutzero(value,minval): #minval and maxval must be the same here (to be symmetric)
    number = (value-minval)/(-2*minval) #curr_line[i_key]/(2*max(-min_for_key[i_key],max_for_key[i_key]))+1
    return number
def revert_linear_symmetricaboutzero(number,minval): #minval and maxval must be the same here (to be symmetric)
    value = number*(-2*minval)+minval
    return value

def calc_spectra_from_twogammadist(a_ms_ice, b_ms_ice,a_ms_snow, b_ms_snow ,fallspeedmodel,qi,qni,qs,qns,mu_ice,gam_ice,mu_snow,gam_snow,turb,first_spectra,comb_dic=dict(),return_also_fullspectra=False):
    if first_spectra:
        comb_dic["qni"] = [qni]; 
        if qni==0:
            comb_dic["meanmass_ice"] = [np.nan]
        else:
            comb_dic["meanmass_ice"] = [qi/qni] #comb_dic["qni"] = [qni]
        if qi==0:
            comb_dic["qs/qi"] = [np.nan]
        else:
            comb_dic["qs/qi"] = [qs/qi]
        comb_dic["qs"] = [qs]
        comb_dic["qns"] = [qns]
        if qns==0:
            comb_dic["meanmass_snow"] = [np.nan]
        else:
            comb_dic["meanmass_snow"] = [qs/qns] #comb_dic["qns"] = [qns]
        comb_dic["mu_ice"] = [mu_ice]
        comb_dic["mu_snow"] = [mu_snow]
        comb_dic["gam_ice"] = [gam_ice]
        comb_dic["gam_snow"] = [gam_snow]
        comb_dic["turb"] = [turb]
    else:
        
        comb_dic["qni"].append(qni); 
        if qni==0:
            comb_dic["meanmass_ice"].append([np.nan])
        else:
            comb_dic["meanmass_ice"].append([qi/qni])
        if qi==0:
             comb_dic["qs/qi"].append([np.nan])
        else:
            comb_dic["qs/qi"].append([qs/qi]); 
            comb_dic["qs"].append([qs])
        if qns==0:
            comb_dic["meanmass_snow"].append([np.nan])
        else:
            comb_dic["meanmass_snow"].append([qs/qns])
        comb_dic["qns"].append(qns);
        comb_dic["mu_ice"].append(mu_ice)
        comb_dic["mu_snow"].append(mu_snow)
        comb_dic["gam_ice"].append(gam_ice)
        comb_dic["gam_snow"].append(gam_snow)
        comb_dic["turb"].append(turb)

    pam = pyPamtra.pyPamtra()

    a_vel_icecosmo5 = 2.60e1; b_vel_icecosmo5 = 0.215790 #ATTENTION 2.77e1 or 2.60d
    a_velD_icecosmo5 = a_vel_icecosmo5*(a_ms_ice)**(b_vel_icecosmo5); b_velD_icecosmo5 = b_ms_ice*b_vel_icecosmo5 #derive constant as a function of D (v=a_velD*diam**v_velD) instead of m
    #print "power law ice (a,b): ", a_velD_icecosmo5,b_vel_icecosmo5
    #from IPython.core.debugger import Tracer ; Tracer()()
    if fallspeedmodel=="powerlaw":
        pam.df.addHydrometeor(("ice", 0.6, -1, -99, a_ms_ice, b_ms_ice,0.684, 2., 13, 100, "mgamma", -99.,    -99.,    mu_ice,   gam_ice,  1e-8,    1e-1,     "ss-rayleigh-gans",  "corPowerLaw_" + "{:.3f}".format(a_velD_icecosmo5)+"_" + "{:.3f}".format(b_velD_icecosmo5),    -99.)) #from Ulis ice
        #ATTENTION: for setting set powerlaw of ice to snow category: uncomment next line
        #pam.df.addHydrometeor(("snow",0.6, -1, -99.,  a_ms_snow, b_ms_snow ,0.3971,1.88,13,100, "mgamma", -99.,    -99.,     mu_snow ,  gam_snow  ,  1e-8,    1e-1,     "ss-rayleigh-gans",  "corPowerLaw_" + "{:.3f}".format(a_velD_icecosmo5)+"_" + "{:.3f}".format(b_velD_icecosmo5),    -99.))
        pam.df.addHydrometeor(("snow",0.6, -1, -99.,  a_ms_snow, b_ms_snow ,0.3971,1.88,13,100, "mgamma", -99.,    -99.,     mu_snow ,  gam_snow  ,  1e-8,    1e-1,     "ss-rayleigh-gans",  "corPowerLaw_5.511054_0.25",    -99.))
    if fallspeedmodel=="powerlawGaussvspread":
        import scipy
        
        #from IPython.core.debugger import Tracer ; Tracer()() 
        #TODO: add more categories according to student-t distribution (assuming Gaussian spread in fall speed)
        pam.df.addHydrometeor(("ice", 0.6, -1, -99, a_ms_ice, b_ms_ice,0.684, 2., 13, 100, "mgamma", -99.,    -99.,    mu_ice,   gam_ice,  1e-8,    1e-1,     "ss-rayleigh-gans",  "corPowerLaw_" + "{:.3f}".format(a_velD_icecosmo5)+"_" + "{:.3f}".format(b_velD_icecosmo5),    -99.)) #from Ulis ice
        pam.df.addHydrometeor(("snow",0.6, -1, -99.,  a_ms_snow, b_ms_snow ,0.3971,1.88,13,100, "mgamma", -99.,    -99.,     mu_snow ,  gam_snow  ,  1e-8,    1e-1,     "ss-rayleigh-gans",  "corPowerLaw_5.511054_0.25",    -99.))
    
    if fallspeedmodel=="Atlas":
        #ATTENTION: for testing set to cvel of Atlas type ice to those of snow: uncomment next line
        #pam.df.addHydrometeor(("ice", 0.6, -1, -99, a_ms_ice, b_ms_ice,0.684, 2., 13, 100, "mgamma", -99.,    -99.,    mu_ice,   gam_ice,  1e-8,    1e-1,     "ss-rayleigh-gans",  "corAtlas_1.872_1.872_3698",    -99.))
        pam.df.addHydrometeor(("ice", 0.6, -1, -99, a_ms_ice, b_ms_ice,0.684, 2., 13, 100, "mgamma", -99.,    -99.,    mu_ice,   gam_ice,  1e-8,    1e-1,     "ss-rayleigh-gans",  "corAtlas_1.872_1.872_932.5",    -99.))
        #a_vel_snowSBBnonsphere = 1.271; b_vel_snowSBBnonsphere = 1.252; c_vel_snowSBBnonsphere = 3.698e3 #fits to Ulis ice (despite the name)
        #ATTENTION: for testing set to cvel of Atlas type snow to those of ice
        #pam.df.addHydrometeor(("snow",0.6, -1, -99.,  a_ms_snow, b_ms_snow ,0.3971,1.88,13,100, "mgamma", -99.,    -99.,     mu_snow ,    gam_snow  ,  1e-8,    1e-1,     "ss-rayleigh-gans",  "corAtlas_1.271_1.252_932.5",    -99.))
        pam.df.addHydrometeor(("snow",0.6, -1, -99.,  a_ms_snow, b_ms_snow ,0.3971,1.88,13,100, "mgamma", -99.,    -99.,     mu_snow ,    gam_snow  ,  1e-8,    1e-1,     "ss-rayleigh-gans",  "corAtlas_1.271_1.252_3698",    -99.))
#iwc_q     0.6      -1        -99.     1.58783  2.56   0.684    2.     13           100       mgamma     -99.    -99.    1.564   0.8547  1.d-8    1.d-0     ss-rayleigh-gans  corPowerLaw_30.606_0.5533    -99.

    pam = pyPamtra.importer.createUsStandardProfile(pam,hgt_lev=[100,200,300])

    pam.set["verbose"] = 0
    pam.p["hydro_q"][0,0,1,0] = qi #1e-4 #cloud ice
    pam.p["hydro_n"][0,0,1,0] = qni #1e5 #cloud ice
    pam.p["hydro_q"][0,0,1,1] = qs #1e-4 #snow
    pam.p["hydro_n"][0,0,1,1] = qns #1e4#snow
    
    pam.p["airturb"][:] = turb
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
    pam.nmlSet["hydro_adaptive_grid"]=False
    #pam.nmlSet["radar_max_V"] = 3.
    #pam.nmlSet["radar_min_V"] = -3.

    pam.runParallelPamtra([9.6,35.5,95], pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers='auto')

    #save the result of the radar moments into comb_dic
    if first_spectra:
        #reflectivity
        comb_dic["Ze_X"] = np.array([pam.r["Ze"][0,0,1,0,0,0]])
        comb_dic["Ze_Ka"] = np.array([pam.r["Ze"][0,0,1,1,0,0]])
        comb_dic["Ze_W"] = np.array([pam.r["Ze"][0,0,1,2,0,0]])
        comb_dic["DWR_X_Ka"] = np.array([pam.r["Ze"][0,0,1,0,0,0]-pam.r["Ze"][0,0,1,1,0,0]])
        comb_dic["DWR_Ka_W"] = np.array([pam.r["Ze"][0,0,1,1,0,0]-pam.r["Ze"][0,0,1,2,0,0]])


        #higher moments
        comb_dic["vDoppler_X"] = np.array([pam.r["radar_moments"][0,0,1,0,0,0,0]])
        comb_dic["vDoppler_Ka"] = np.array([pam.r["radar_moments"][0,0,1,1,0,0,0]])
        comb_dic["vDoppler_W"] = np.array([pam.r["radar_moments"][0,0,1,2,0,0,0]])
        comb_dic["swidth_X"] = np.array([pam.r["radar_moments"][0,0,1,0,0,0,1]])
        comb_dic["swidth_Ka"] = np.array([pam.r["radar_moments"][0,0,1,1,0,0,1]])
        comb_dic["swidth_W"] = np.array([pam.r["radar_moments"][0,0,1,2,0,0,1]])
        comb_dic["skewn_X"] = np.array([pam.r["radar_moments"][0,0,1,0,0,0,2]])
        comb_dic["skewn_Ka"] = np.array([pam.r["radar_moments"][0,0,1,1,0,0,2]])
        comb_dic["skewn_W"] = np.array([pam.r["radar_moments"][0,0,1,2,0,0,2]])
    else:

        #reflectivity
        comb_dic["Ze_X"] = np.append(comb_dic["Ze_X"],pam.r["Ze"][0,0,1,0,0,0])
        comb_dic["Ze_Ka"] = np.append(comb_dic["Ze_Ka"],pam.r["Ze"][0,0,1,1,0,0])
        comb_dic["Ze_W"] = np.append(comb_dic["Ze_W"],pam.r["Ze"][0,0,1,2,0,0])
        comb_dic["DWR_X_Ka"] = np.append(comb_dic["DWR_X_Ka"],pam.r["Ze"][0,0,1,0,0,0]-pam.r["Ze"][0,0,1,1,0,0])
        comb_dic["DWR_Ka_W"] = np.append(comb_dic["DWR_Ka_W"],pam.r["Ze"][0,0,1,1,0,0]-pam.r["Ze"][0,0,1,2,0,0])
        #higher moments
        comb_dic["vDoppler_X"] = np.append(comb_dic["vDoppler_X"],pam.r["radar_moments"][0,0,1,0,0,0,0])  
        comb_dic["vDoppler_Ka"] = np.append(comb_dic["vDoppler_Ka"],pam.r["radar_moments"][0,0,1,1,0,0,0])
        comb_dic["vDoppler_W"] = np.append(comb_dic["vDoppler_W"],pam.r["radar_moments"][0,0,1,2,0,0,0])
        comb_dic["swidth_X"] = np.append(comb_dic["swidth_X"],pam.r["radar_moments"][0,0,1,0,0,0,1])
        comb_dic["swidth_Ka"] = np.append(comb_dic["swidth_Ka"],pam.r["radar_moments"][0,0,1,1,0,0,1])
        comb_dic["swidth_W"] = np.append(comb_dic["swidth_W"],pam.r["radar_moments"][0,0,1,2,0,0,1])
        comb_dic["skewn_X"] = np.append(comb_dic["skewn_X"],pam.r["radar_moments"][0,0,1,0,0,0,2])
        comb_dic["skewn_Ka"] = np.append(comb_dic["skewn_Ka"],pam.r["radar_moments"][0,0,1,1,0,0,2])
        comb_dic["skewn_W"] = np.append(comb_dic["skewn_W"],pam.r["radar_moments"][0,0,1,2,0,0,2])

        
    if return_also_fullspectra: #ATTENTION: this is only applicable for a single spectra (so far)
        #from IPython.core.debugger import Tracer ; Tracer()()
        comb_dic["psd_d"] = pam.r['psd_d'][0,0,1,0] #this is the same for ice and snow because hydro_adaptive_grid is set to False
        comb_dic["psd_n_ice"] =  pam.r['psd_n'][0,0,1,0]
        comb_dic["psd_n_snow"] =  pam.r['psd_n'][0,0,1,1]
        comb_dic["psd_n_all"] = pam.r['psd_n'][0,0,1,0]+pam.r['psd_n'][0,0,1,1]
        comb_dic["psd_m_ice"] =  pam.r['psd_mass'][0,0,1,0] #this is the mass at the center of the bin
        comb_dic["psd_m_snow"] =  pam.r['psd_mass'][0,0,1,1]
        comb_dic["psd_m_all"] = pam.r['psd_mass'][0,0,1,0]+pam.r['psd_mass'][0,0,1,1]

        comb_dic["velbins_X"] = pam.r["radar_vel"][0]
        comb_dic["velbins_Ka"] = pam.r["radar_vel"][1]
        comb_dic["velbins_W"] = pam.r["radar_vel"][2]
        comb_dic["spectra_X"]  = pam.r["radar_spectra"][0,0,1,0,0]
        comb_dic["spectra_Ka"] = pam.r["radar_spectra"][0,0,1,1,0]
        comb_dic["spectra_W"]  = pam.r["radar_spectra"][0,0,1,2,0]

    return comb_dic

def collect_values_of_comb(comb_dic,selected_keys,min_vals,max_vals,values_flag):
    '''
    collect the values which belong to a combination of parameters and save them in an array
    '''
    curr_line = np.array([]) #this is where the values are collected
    for i_key,key in enumerate(selected_keys): 
        curr_line = np.append(curr_line,comb_dic[key][i_comb])
        if values_flag[i_key]==0:
            curr_line[i_key] = project_linear(curr_line[i_key],min_vals[i_key],max_vals[i_key]) #normalize to maximum
        if values_flag[i_key]==1:
            curr_line[i_key] = project_linear_symmetricaboutzero(curr_line[i_key],min_vals[i_key]) #make it symmetric about zero
        #if values_flag[i_key]==2:
        #    curr_line[i_key] = project_numcon_to_0to1(curr_line[i_key],min_vals[i_key],max_vals[i_key])  #logarithmate the values
        if values_flag[i_key]==2: #this is f.e. used for the mixing ratio: mixrat=[10**-6,10**-1] is mapped to [0,1]
            curr_line[i_key] = project_mixrat_to_0to1(curr_line[i_key],min_vals[i_key],max_vals[i_key]) #logarithmate the values and add a constant number to the logarithms to be in the positive space

    return curr_line
create_spectra = False #calculate the spectra based on the input values 
a_ms_ice = 1.58783; b_ms_ice = 2.56
a_ms_snow =  0.038; b_ms_snow = 2.0
for fallspeedmodel in ["Atlas","powerlaw"]:
    print fallspeedmodel

    comb_dic = dict() #initialize dictionary which includes the input parameter (f.e.  parameters of the size distribution) and the forward operated quantities

    #set up figure
    number_of_plots = 1
    figsize_height = 6.0/2.0*(number_of_plots)
    fig	=	plt.figure(figsize=(8.0,figsize_height))#figsize=(4, 4))

    #define values for a highlighted combination (this is the red line)
    plot_highlighted_spectra = True
    qni_highlighted = 1e5
    qs_highlighted = 1e-4
    qns_highlighted = 1e4
    mu_ice_highlighted = 1.564
    gam_ice_highlighted = 0.6 #0.8547
    mu_snow_highlighted = 1.0
    gam_snow_highlighted = 2.0 #1.0
    turb_highlighted = 0.001

    qi=1e-4 #this is held fix because it is just scaling the spectra

    first_spectra = True #
    if create_spectra: #do not run this expensive loop if the data is there already
            for qni in [1e1,1e2,1e3,1e4,1e5]: #,1e5]:
                for qs in [1e-4,1e-5,1e-3]: #,1e-2]: #,1e-3]:
                    for qns in [1e4,1e1,1e2,1e3,1e5]: #,1e4]:
                        for mu_ice in [1.564]: #,1.0,0.8,1.8]: #first is always the default from SB
                            for gam_ice in [0.8547]: #,0.4,0.6,1.0]: #,1.0,1.2]: #first is always the default from SB
                                for mu_snow in [1.0]: #,0.8,1.2,1.4,1.6]: #first is always the default from SB
                                    for gam_snow in [1.0]: #,0.8,1.2,1.4,1.6,1.8]: #first is always the default from SB
                                        for turb in [0.0, 0.1,0.2,0.5,1.0]:
                                            print "qni",qni,"turb",turb
                                            if first_spectra:
                                                comb_dic = calc_spectra_from_twogammadist(a_ms_ice, b_ms_ice,a_ms_snow, b_ms_snow,fallspeedmodel,qi,qni,qs,qns,mu_ice,gam_ice,mu_snow,gam_snow,turb,first_spectra,comb_dic=dict())
                                            else:
                                                comb_dic = calc_spectra_from_twogammadist(a_ms_ice, b_ms_ice,a_ms_snow, b_ms_snow,fallspeedmodel,qi,qni,qs,qns,mu_ice,gam_ice,mu_snow,gam_snow,turb,first_spectra,comb_dic=comb_dic)
                                            first_spectra=False #from now on the dictionaries are present already
                                            

                                                #sys.exit()

            np.save('saved_dicts/two_mgamma_variousparams_' + fallspeedmodel + '.npy', comb_dic) 
                
    if not create_spectra: #load the dictionary if it was not created before
        print 'saved_dicts/two_mgamma_variousparams_' + fallspeedmodel + '.npy'
        comb_dic = np.load('saved_dicts/two_mgamma_variousparams_' + fallspeedmodel + '.npy', comb_dic).item()
    #select list of keys (and define thereby the order) for the spider diagram
    #selected_keys = ["qs","qni","qns","mu_ice","mu_snow","gam_ice","gam_snow","swidth_X","swidth_Ka","swidth_W","skewn_X","skewn_Ka","skewn_W"]
    selected_keys = ["qs/qi","meanmass_ice","meanmass_snow","mu_ice","mu_snow","gam_ice","gam_snow","turb","swidth_X","swidth_Ka","swidth_W","skewn_Ka","vDoppler_Ka","DWR_X_Ka","DWR_Ka_W"] #,"skewn_Ka","skewn_W"]

    values_flag =    [ 2,           2,            2,            0,        0,         0,       0,       0,         0,          0,          0,        1,         0,            0,       0     ]#        1,          1]# 0-> linear and only positive 1-> linear (positive or negative) 2-> logarithmic (only bigger than 1) 3-> logarithmic (only smaller than 1)
    min_vals =       [ -2,         -10 ,           -10,         0,       0,        0,       0,         0,   0,          0,          0,         -5,         0,           0,      0      ]#,       -5,         -5]       #for the logarithmicly scaled variables this is the logarithm
    max_vals =       [ 3,           -5 ,            -5,          2,       2,        2,       2,        2,  0.5,        0.5,        0.5,        5,         2,            10,     10     ]#,        5,          5]       #for the logarithmicly scaled variables this is the logarithm

    len_of_dict = len(comb_dic[selected_keys[0]]) #the variable to get the length is chosen arbitrarily here
    num_of_dictentries = len(selected_keys)
    # Initialise the spider plot
    figspider = plt.figure(figsize=(10.0,10.0))
    axspider = plt.subplot(111, polar=True)

    #set angles for direction of parameters
    angles = [n / float(num_of_dictentries) * 2 * np.pi for n in range(num_of_dictentries)]

    # Draw one axes per variable + add labels labels yet
    parameters = selected_keys


    labelparam = ["" for x in range(len(selected_keys))]
    for i_key in range(0,len(selected_keys)):
        labelparam[i_key] = parameters[i_key] #+ ' ' + str(max_for_key[i_key])
        #print parameters[i_key] + str(max_for_key[i_key])
    #label each angle with the corresponding parameter and maximum value
    plt.xticks(angles, labelparam, color='grey', size=12)

    ###############
    #plot each combination of parameters in the polar plot
    ###############
    for i_comb in range(len_of_dict-1,-1,-1): #loop over each combination 

        curr_line = collect_values_of_comb(comb_dic,selected_keys,min_vals,max_vals,values_flag)

        #plot lines which not fullfill the constraints of the observations in gray and those which fullfill it in blue
        if comb_dic["swidth_Ka"][i_comb]<0.2  and (-1.<comb_dic["skewn_Ka"][i_comb]<1.) and comb_dic["DWR_Ka_W"][i_comb]<4: # and comb_dic["qni"][i_comb]==1e3:
            axspider.plot(angles,curr_line, linewidth=2, linestyle='-',color='blue',alpha=0.8)
        else:
            axspider.plot(angles,curr_line, linewidth=1, linestyle='-',color='black',alpha=0.1)


    #calculate a specific combination and plot also the spectra from that combination
    first_spectra = True #reset it, because we are creating a new dictionary now
    highlighted_dic = calc_spectra_from_twogammadist(a_ms_ice, b_ms_ice,a_ms_snow, b_ms_snow,fallspeedmodel,qi,qni_highlighted,qs_highlighted,qns_highlighted,mu_ice_highlighted,gam_ice_highlighted,mu_snow_highlighted,gam_snow_highlighted,turb_highlighted,first_spectra,comb_dic=dict(),return_also_fullspectra=True)
    highlighted_dic_iceonly = calc_spectra_from_twogammadist(a_ms_ice, b_ms_ice,a_ms_snow, b_ms_snow,fallspeedmodel,qi,qni_highlighted,0,0,mu_ice_highlighted,gam_ice_highlighted,mu_snow_highlighted,gam_snow_highlighted,turb_highlighted,first_spectra,comb_dic=dict(),return_also_fullspectra=True)
    highlighted_dic_snowonly = calc_spectra_from_twogammadist(a_ms_ice, b_ms_ice,a_ms_snow, b_ms_snow,fallspeedmodel,0,0,qs_highlighted,qns_highlighted,mu_ice_highlighted,gam_ice_highlighted,mu_snow_highlighted,gam_snow_highlighted,turb_highlighted,first_spectra,comb_dic=dict(),return_also_fullspectra=True)
    if plot_highlighted_spectra:
        # Initialise the spectra plot
        figsize_height = 8.0/2.0*(2.)
        figspec = plt.figure(figsize=(8.0,figsize_height))
        #plot the size distributions of both categories
        #N(D)
        axsizedist = plt.subplot(311)
        axsizedist.semilogx(highlighted_dic["psd_d"],highlighted_dic["psd_n_all"],color='r',label='all')
        axsizedist.semilogx(highlighted_dic["psd_d"],highlighted_dic["psd_n_ice"],color='orange',linestyle='--',marker='o',markersize=6,markevery=3,label='only ice')
        axsizedist.semilogx(highlighted_dic["psd_d"],highlighted_dic["psd_n_snow"],color='magenta',linestyle='-.',marker='*',markersize=6,markevery=3,label='only snow')
        plt.legend()
        plt.xlabel('diameter m')
        plt.ylabel('number concentration / m-4')
        #m(D)
        axsizedist = plt.subplot(312)
        #axsizedist.semilogx(highlighted_dic["psd_d"],1.58783*highlighted_dic["psd_d"]**2.56*highlighted_dic["psd_n_all"],color='k')
        #m_ice = a_ms_ice*highlighted_dic["psd_d"]**b_ms_ice*highlighted_dic["psd_n_all"]
        #m_snow = a_ms_snow*highlighted_dic["psd_d"]**b_ms_snow*highlighted_dic["psd_n_all"]
        axsizedist.semilogx(highlighted_dic["psd_d"],highlighted_dic["psd_m_ice"]*highlighted_dic["psd_n_ice"]+highlighted_dic["psd_m_snow"]*highlighted_dic["psd_n_snow"],color='r',label='all')
        axsizedist.semilogx(highlighted_dic["psd_d"],highlighted_dic["psd_m_ice"]*highlighted_dic["psd_n_ice"],color='orange',linestyle='--',marker='o',markersize=6,markevery=3)
        axsizedist.semilogx(highlighted_dic["psd_d"],highlighted_dic["psd_m_snow"]*highlighted_dic["psd_n_snow"],color='magenta',linestyle='-.',marker='*',markersize=6,markevery=3)
        axsizedist.plot(np.nan,np.nan,color='orange',linestyle='--',label='only ice',marker='o')
        axsizedist.plot(np.nan,np.nan,color='magenta',linestyle='-.',label='only snow',marker='*')
        plt.legend()
        plt.xlabel('diameter m')
        plt.ylabel('mass concentration / kg m-4')
        #plot the Doppler spectra in all three frequencies
        axspec = plt.subplot(313)
        axspec.plot(highlighted_dic["velbins_X"],highlighted_dic["spectra_X"],label='9.5GHz',color='b',linestyle='-')
        axspec.plot(highlighted_dic["velbins_Ka"],highlighted_dic["spectra_Ka"],label='35GHz',color='r',linestyle='-')
        axspec.plot(highlighted_dic["velbins_W"],highlighted_dic["spectra_W"],label='95GHz',color='g',linestyle='-')
        #plot the spectra in all three frequencies if there is only ice
        axspec.plot(highlighted_dic_iceonly["velbins_X"],highlighted_dic_iceonly["spectra_X"],label='__None',color='b',linestyle='--',marker='o',markerfacecolor='k',markeredgecolor='k',markersize=6,markevery=3)
        #axspec.plot(highlighted_dic_iceonly["velbins_Ka"],highlighted_dic_iceonly["spectra_Ka"],label='__None',color='r',linestyle='--')
        #axspec.plot(highlighted_dic_iceonly["velbins_W"],highlighted_dic_iceonly["spectra_W"],label='__None',color='g',linestyle='--')
        #plot the spectra in all three frequencies if there is only snow
        axspec.plot(highlighted_dic_snowonly["velbins_X"],highlighted_dic_snowonly["spectra_X"],label='__None',color='b',linestyle='-.',marker='*',markerfacecolor='k',markeredgecolor='k',markersize=6,markevery=3)
        #axspec.plot(highlighted_dic_snowonly["velbins_Ka"],highlighted_dic_snowonly["spectra_Ka"],label='__None',color='r',linestyle='-.')
        #axspec.plot(highlighted_dic_snowonly["velbins_W"],highlighted_dic_snowonly["spectra_W"],label='__None',color='g',linestyle='-.')
        
        axspec.set_ylim([-30,20]); axspec.set_xlim([-0.2,2])
        axspec.text(-2.8,20,"Ze: " + str(highlighted_dic["Ze_Ka"])) #+ "\n vDoppler      spectr. width     skewness      kurtosis \n" + str(pam.r["radar_moments"][0,0,1,0,0]))
        plt.xlabel('Doppler velocity m s-1')
        plt.ylabel('spectral reflectivity dB')
        #add invisible line for labelling iceonly and snowonly
        axspec.plot(np.nan,np.nan,color='k',linestyle='--',label='only ice',marker='o')
        axspec.plot(np.nan,np.nan,color='k',linestyle='-.',label='only snow',marker='*')
        #
        axspec.text(0,10,"vDoppler_Ka: " + str(highlighted_dic["vDoppler_Ka"][0]))
        axspec.text(0,6,"swidth_Ka: " + str(highlighted_dic["swidth_Ka"][0]))
        axspec.text(0,2,"skewn_Ka: " + str(highlighted_dic["skewn_Ka"][0]))

        #create legend
        plt.legend()
        plt.tight_layout()
        plt.savefig("/home/mkarrer/Dokumente/pamtrafiles/ideal_spectra/" + fallspeedmodel+ "/" + "spectra_"+ fallspeedmodel + ".png",dpi=400) # + "/" + "qni" + str(qni) + "qs" + str(qs) + "qs" + str(qni)  + "mu_ice" + str(mu_ice) + "gam_ice" + str(gam_ice) + "mu_snow" + str(mu_snow) + "gam_snow" + str(gam_snow) + ".png")
        print "the spectra is at: " + "/home/mkarrer/Dokumente/pamtrafiles/ideal_spectra/" + fallspeedmodel  + "/" + "spectra_"+ fallspeedmodel + ".png"#+ "qni" + str(qni) + "qs" + str(qni) + "qs" + str(qni)  + "mu_ice" + str(mu_ice) + "gam_ice" + str(gam_ice) + "mu_snow" + str(mu_snow) + "gam_snow" + str(gam_snow) + ".png"
        plt.clf()
    ###           
    #plot the highlighted combination in the spider-plot
    ###
    #also plot the ice and snow only
    highlighted_line_iceonly_line = collect_values_of_comb(highlighted_dic_iceonly,selected_keys,min_vals,max_vals,values_flag)
    axspider.plot(angles,highlighted_line_iceonly_line, linewidth=1.5, linestyle='--',color='orange',alpha=0.8)
    highlighted_line_snowonly = collect_values_of_comb(highlighted_dic_snowonly,selected_keys,min_vals,max_vals,values_flag)
    axspider.plot(angles,highlighted_line_snowonly, linewidth=1.5, linestyle='-.',color='magenta',alpha=0.8)
    #no plot the supersition of both categories above
    highlighted_line = collect_values_of_comb(highlighted_dic,selected_keys,min_vals,max_vals,values_flag)
    axspider.plot(angles,highlighted_line, linewidth=1.5, linestyle='-',color='red',alpha=0.8)
    ############
    #mark regions with the same radius colorscale and plot labels for each of them
    ############

    #mark regions with the same radius colorscale (therefore the selected keys must be ordered accordingly)
    change_indices =[0,1,5,10,11,12,13] #define where the data range changes
    angles_change = [None] * (len(change_indices)+1)
    angle1=-angles[change_indices[0]+1]/2; angle2=(angles[change_indices[0]]+angles[change_indices[0]+1])/2; angle3=(angles[change_indices[1]]+angles[change_indices[1]+1])/2; angle4=(angles[change_indices[2]]+angles[change_indices[2]])/2;
    colors_for_colorscale=["firebrick","darksalmon"]*10 #,"firebrick","darksalmon","firebrick","darksalmon"] #["rosybrown","firebrick","darksalmon","sienna","sandybrown"]

    #add filling and labels for each section
    for i_section in range(0,len(change_indices)):
        if i_section==0: #snow mixing ratio
            #calculate the angle at the end of the section
            angles_change[i_section] = (angles[change_indices[i_section]]+ angles[change_indices[i_section+1]])/2.0
            #fill the section
            axspider.fill_between(   np.linspace(-angles[1]/2., (angles_change[i_section]), 100),0,1,alpha=0.2,color=colors_for_colorscale[i_section],linewidth=2)
            #create ticks for each section
            for number in np.arange(0.2,1.1,0.2):
                axspider.text(0,number,'{:.0E}'.format(revert_0to1_to_mixrat(number,min_vals[change_indices[0]],max_vals[change_indices[0]])),color=colors_for_colorscale[i_section])
        if i_section==1: #number concentration
            #calculate the angle at the end of the section
            angles_change[i_section] = (angles[change_indices[i_section]]+ angles[change_indices[i_section+1]])/2.0
            #fill the section        
            axspider.fill_between(  np.linspace(angles_change[i_section], angles_change[i_section-1], 100),0,1,alpha=0.2,color=colors_for_colorscale[i_section],linewidth=2)
            #create ticks for each section
            for number in np.arange(0.2,1.1,0.2):
                axspider.text((angles_change[i_section-1]+angles_change[i_section])/2.0,number,'{:.0E}'.format(revert_0to1_to_mixrat(number,min_vals[change_indices[i_section]],max_vals[change_indices[i_section]])),color=colors_for_colorscale[i_section])
        if i_section==2: #distribution parameters
            #calculate the angle at the end of the section
            angles_change[i_section] = (angles[change_indices[i_section]]+ angles[change_indices[i_section+1]])/2.0
            #fill the section
            axspider.fill_between(np.linspace(angles_change[i_section], angles_change[i_section-1], 100),0,1,alpha=0.2,color=colors_for_colorscale[i_section],linewidth=2)
            #create ticks for each section
            for number in np.arange(0.2,1.1,0.2):
                axspider.text((angles_change[i_section-1]+angles_change[i_section])/2.0,number,'{:1.2f}'.format(revert_linear(number,min_vals[change_indices[i_section]],max_vals[change_indices[i_section]])),color=colors_for_colorscale[i_section])
        if i_section==3: #spectral widths
            #calculate the angle at the end of the section
            angles_change[i_section] = (angles[change_indices[i_section]]+ angles[change_indices[i_section+1]])/2.0
            #fill the section
            axspider.fill_between(np.linspace(angles_change[i_section], angles_change[i_section-1], 100),0,1,alpha=0.2,color=colors_for_colorscale[i_section],linewidth=2)
            #create ticks for each section
            for number in np.arange(0.2,1.1,0.2):
                axspider.text((angles_change[i_section-1]+angles_change[i_section])/2.0,number,'{:1.2f}'.format(revert_linear(number,min_vals[change_indices[i_section]],max_vals[change_indices[i_section]])),color=colors_for_colorscale[i_section])
        if i_section==4: #and (i_section<(len(change_indices)-1)):
            #calculate the angle at the end of the section
            angles_change[i_section] = (angles[change_indices[i_section]]+ angles[change_indices[i_section+1]])/2.0 #+ angles[change_indices[i_section+1]])/2.0
            #fill the section        
            axspider.fill_between(  np.linspace(angles_change[i_section], angles_change[i_section-1], 100),0,1,alpha=0.2,color=colors_for_colorscale[i_section],linewidth=2)
            #create ticks for each section
            for number in np.arange(0.2,1.1,0.2):
                axspider.text((angles_change[i_section-1]+angles_change[i_section])/2.0,number,'{:1.2f}'.format(revert_linear_symmetricaboutzero(number,min_vals[change_indices[i_section]])),color=colors_for_colorscale[i_section])
        if i_section==5: #and (i_section<(len(change_indices)-1)):
            #calculate the angle at the end of the section
            angles_change[i_section] = (angles[change_indices[i_section]]+ angles[change_indices[i_section+1]])/2.0 #+ angles[change_indices[i_section+1]])/2.0
            #fill the section        
            axspider.fill_between(  np.linspace(angles_change[i_section], angles_change[i_section-1], 100),0,1,alpha=0.2,color=colors_for_colorscale[i_section],linewidth=2)
            #create ticks for each section
            for number in np.arange(0.2,1.1,0.2):
                axspider.text((angles_change[i_section-1]+angles_change[i_section])/2.0,number,'{:1.2f}'.format(revert_linear(number,min_vals[change_indices[i_section]],max_vals[change_indices[i_section]])),color=colors_for_colorscale[i_section])
        if i_section==6: #and (i_section<(len(change_indices)-1)):
            #calculate the angle at the end of the section
            angles_change[i_section] = 2.*np.pi-(angles[1]) #+ angles[change_indices[i_section+1]])/2.0 #this is for the last one
            #fill the section        
            axspider.fill_between(  np.linspace(angles_change[i_section], angles_change[i_section-1], 100),0,1,alpha=0.2,color=colors_for_colorscale[i_section],linewidth=2)
            #create ticks for each section
            for number in np.arange(0.2,1.1,0.2):
                axspider.text((angles_change[i_section-1]+angles_change[i_section])/2.0,number,'{:1.2f}'.format(revert_linear(number,min_vals[change_indices[i_section]],max_vals[change_indices[i_section]])),color=colors_for_colorscale[i_section])

        
    #remove the default radius label
    axspider.set_yticklabels([])

    axspider.set_ylim(0,1)
    figspider.savefig("/home/mkarrer/Dokumente/pamtrafiles/ideal_spectra/spider_plot" + fallspeedmodel +".png",dpi=400)