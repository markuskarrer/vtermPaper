'''
plot the percentiles of diameter D and terminal velocity v from McSnow simulations resolved by height
'''

#from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import math
import os
import subprocess
from netCDF4 import Dataset
import traceback
import sys
#self-written functions
from functions import __plotting_functions
from functions import __postprocess_McSnow
from functions import __postprocess_PAMTRA
from functions import __postprocess_SB
from functions import __general_utilities
#read variables passed by shell script
tstep = int(os.environ["tstep"]) #start of the averaging period
tstep_end = int(os.environ["tstep_end"]) #string with 4 numbers and 'min'
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"]
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output
MC_dir = os.environ["MC"]
if "MCtermvel_specifier_onestring" in os.environ.keys():
    MCtermvel_specifier_onestring = os.environ["MCtermvel_specifier_onestring"]
else:
    MCtermvel_specifier_onestring = ''
if "McSnow_geom_specifier_onestring" in  os.environ.keys():
    McSnow_geom_specifier_onestring = os.environ["McSnow_geom_specifier_onestring"]
else:
    McSnow_geom_specifier_onestring = ''
if "separated_by_sensrunsMC" in os.environ.keys():
    separated_by_sensrunsMC = (os.environ["separated_by_sensrunsMC"]=="True") #plot a line for different sensitivity runs
    separated_by_sensruns_onestring= os.environ["model_setup_specifier_onestring"]
else:
    print "separated_by_sensruns not found in environ: set to False"
    separated_by_sensruns = False
if "separated_by_fallspeedsens" in os.environ.keys():
    separated_by_fallspeedsens = (os.environ["separated_by_fallspeedsens"]=="True") #plot a line for different sensitivity runs
else:
    print "separated_by_fallspeedsens not found in environ: set to False"
    separated_by_fallspeedsens = False
#from https://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy
def weighted_quantile(values, quantiles, sample_weight=None, 
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)



#directory of experiments
directory = MC_dir + "/experiments/"

##general settings for the plotting
number_of_plots = 5

#optimize the appearance of the plot (figure size, fonts)
[fig,axes] = __plotting_functions.proper_font_and_fig_size(number_of_plots,legend_fontsize='medium')

#loop over different SB sensitivity runs
#plot different lines for different McSnow geometries
McSnow_geom_list = McSnow_geom_specifier_onestring.split('_')
MCtermvel_list = MCtermvel_specifier_onestring.split('_')
for i_sensMC, sensrun_now_MC in  enumerate(McSnow_geom_list): #loop over different McSnow geometry sensitivity runs
    for i_sensMCfallspeed, sensrun_now_MC_fallspeed in enumerate(MCtermvel_list): #loop over different fall speed models (boehm, KC05,HW10,powerlaw, powerlawSB)
        if separated_by_sensrunsMC:#modify the experiment string
            #from IPython.core.debugger import Tracer ; Tracer()()
            #raw_input(sensrun_now_MC_fallspeed)
            
            if sensrun_now_MC_fallspeed=="powerlawSBdefault": #"fallspeed model" is the old powerlaw
                experiment_splitted = experiment.split('_vt',1) #this results e.g. in [_rm10_rt0_', '3_at2_
                experiment = "_".join([experiment_splitted[0],'vt5',experiment_splitted[1][2:]]) #[1:] cuts the old velocity index plus the following '_'                
            elif sensrun_now_MC_fallspeed=="powerlaw": ##"fallspeed model" is the new powerlaw
                experiment_splitted = experiment.split('_vt',1) #this results e.g. in [_rm10_rt0_', '3_at2_
                experiment = "_".join([experiment_splitted[0],'vt4',experiment_splitted[1][2:]]) #[1:] cuts the old velocity index
                print "powerlaw",experiment
            elif sensrun_now_MC_fallspeed=="Atlas": ##"fallspeed model" is the new powerlaw
                experiment_splitted = experiment.split('_vt',1) #this results e.g. in [_rm10_rt0_', '3_at2_
                experiment = "_".join([experiment_splitted[0],'vt6',experiment_splitted[1][2:]]) #[1:] cuts the old velocity index
                print "powerlaw",experiment
            elif sensrun_now_MC_fallspeed=="boehm": #other sensruns (not modifying the fallspeed model
                #fix the fall speed model again
                experiment_splitted = experiment.split('_vt',1) #this results e.g. in [_rm10_rt0_', '3_at2_
                experiment = "_".join([experiment_splitted[0],'vt3',experiment_splitted[1][2:]]) #[1:] cuts the old velocity index #ATTENTION: vt2 and vt1 not implemented
            McSnow_geom_list_str = MCtermvel_specifier_onestring

        else:
            McSnow_geom_list_str = MCtermvel_specifier_onestring
        #load file
        SP_file = Dataset(directory + experiment + '/mass2fr_' + str(tstep).zfill(4) + '-' + str(tstep_end).zfill(4) + 'min_avtstep_' + str(av_tstep) + '.ncdf',mode='r')
        #create dictionary for all superparticle variables
        SP = dict()

        #if necessary change name of variables
        varlist = SP_file.variables
        #read PAMTRA variables to pamData dictionary
        for var in varlist:#read files and write it with different names in Data
            SP[var] = np.squeeze(SP_file.variables[var])
        #get maximum height for plotting
        model_top = np.nanmax(SP['height'])

        #define height array
        height_bound_vec= np.linspace(0,model_top,251)
        height_cent_vec = (height_bound_vec[:-1]+height_bound_vec[1:])/2
        #initialize arrays of D and v for profiles
        SP_percentiles = dict()
        #quantiles=[0.01,0.1,0.5,0.9,0.99] #
        #quantiles=[0.25,0.5,0.75]
        quantiles=[0.1,0.5,0.9]

        if len(quantiles)==3:
            linestyle_quantiles=["--","-","--"]
            linewidth_quantiles=[0.5,1.5,0.5]
        elif len(quantiles)==5:
            linestyle_quantiles=["-.","--","-","--","-."]
            linewidth_quantiles=[0.5,0.5,1.5,0.5,0.5]
        else:
            print "check quantiles and choose an appropriate linestyle_quantiles vector"; sys.exit(1)
        #print quantiles,linestyle_quantiles; raw_input()

        SP_percentiles["diam"] = np.zeros([len(quantiles),height_bound_vec.shape[0]-1]) #three different percentiles at each height #max. dimension
        SP_percentiles["m_tot"] = np.zeros([len(quantiles),height_bound_vec.shape[0]-1]) #three different percentiles at each height #max. dimension
        SP_percentiles["vt"] = np.zeros([len(quantiles),height_bound_vec.shape[0]-1]) #three different percentiles at each height #terminal velocity
        
        debug_upper_bound=False
        if debug_upper_bound:
            prop="diam"
            in_bin = ((SP['height']>model_top-0.01*model_top))
            prop_upper_boundary = SP[prop][in_bin]
            mult_of_prop_upper_boundary = SP["xi"][in_bin]
            if prop_upper_boundary.shape[0]>1:
                SP_percentiles[prop] =  weighted_quantile(prop_upper_boundary,quantiles,mult_of_prop_upper_boundary) #TODO: that is percentiles of SP not RP!!
            else:
                SP_percentiles[prop] =  np.nan
            print sensrun_now_MC_fallspeed,SP_percentiles[prop]
            raw_input(); continue
            
        #loop over the above defined height vector
        for i_height,height in enumerate(height_bound_vec[:-1]):
            for key in SP_percentiles.keys():
                #select only all SP in the corresponding height
                in_bin = ((SP['height']>height_bound_vec[i_height]) & (SP['height']<height_bound_vec[i_height+1]))
                prop_in_bin = SP[key][in_bin] #[((SP['height']>height_bound_vec[i_height]) & (SP['height']<height_bound_vec[i_height+1]))]
                mult_of_prop_in_bin = SP["xi"][in_bin]*SP["m_tot"][in_bin] #[((SP['height']>height_bound_vec[i_height]) & (SP['height']<height_bound_vec[i_height+1]))]
                #for i_SP_in_bin,in_bin_bool in enumerate(in_bin): #loop over all SP to check whether they are in the right heidht and get their position in the SP dictionary
                #    if in_bin_bool: #if this
                #        for mult_now in mult_of_prop_in_bin:
                #            prop_in_bin = np.append(prop_in_bin,SP[key][i_SP_in_bin])
                if prop_in_bin.shape[0]>1:
                    SP_percentiles[key][:,i_height] =  weighted_quantile(prop_in_bin,quantiles,mult_of_prop_in_bin) #TODO: that is percentiles of SP not RP!!
                else:
                    SP_percentiles[key][:,i_height] =  np.nan
                
            print i_height,height,prop_in_bin.shape, SP_percentiles[key][:,i_height]
            
        #read MCSnows moments profile to see how the percentiles would look like with the size dist prescribed in SB
        filestring_hei2massdens = directory + experiment + "/hei2massdens.dat"
        timestep = tstep/10 #TODO: do not hardcode the 30 minute output interval here
        
        hei2massdens = __postprocess_McSnow.read_hei2massdens(filestring_hei2massdens,timestep=timestep)  
        #create a pseudo-SB-dictionary with MC moments
        twomomMC = dict()
        twomomMC["qs"] = np.expand_dims(hei2massdens["Md"],axis=0)
        twomomMC["qns"] = np.expand_dims(hei2massdens["Nd"],axis=0)
        twomomMC["heights"] = hei2massdens["z"]
        
        #######
        #read twomoment-data
        #######
        #define filestring
        filestring = directory + experiment + "/twomom_d.dat"
        #load netCDF4 file
        twomom_file = Dataset(directory + experiment + '/twomom_d.ncdf',mode='r')
        #create dictionary for all variables from twomom_d.ncdf
        twomom = dict()
        #get relevant timestep
        i_timestep=(tstep/10)-1 #there is no output for t=0min after that there are output steps in 30 minute steps (this could vary)
        #if necessary change name of variables
        varlist = twomom_file.variables
        #read 2mom variables to twomom dictionary
        for var in varlist:#read files and write it with different names in Data

            twomom[var] = np.squeeze(twomom_file.variables[var])
        #initialize arrays for the percentiles with MC moments but SB distribution
        diam = np.logspace(-4,-1,1000)
        m_array = np.logspace(-12,-6,1000)
        N_heights= twomom["qi"].shape[1]
        Nd_snowjplatesnonsphere = np.ones([N_heights,diam.shape[0]])*np.nan
        Md_snowjplatesnonsphere = np.ones([N_heights,diam.shape[0]])*np.nan
        Nm_snowjplatesnonsphere = np.ones([N_heights,m_array.shape[0]])*np.nan
        Mm_snowjplatesnonsphere = np.ones([N_heights,m_array.shape[0]])*np.nan
        
        D_SB_quantiles = np.ones([N_heights,len(quantiles)])*np.nan
        m_SB_quantiles = np.ones([N_heights,len(quantiles)])*np.nan

        #initialize arrays for the percentiles SB runs
        NdMCSB_snowjplatesnonsphere = np.ones([height_bound_vec.shape[0],diam.shape[0]])*np.nan
        MdMCSB_snowjplatesnonsphere = np.ones([height_bound_vec.shape[0],diam.shape[0]])*np.nan
        D_MCSB_quantiles = np.ones([height_bound_vec.shape[0],len(quantiles)])*np.nan
        NmMCSB_snowjplatesnonsphere = np.ones([height_bound_vec.shape[0],diam.shape[0]])*np.nan
        MmMCSB_snowjplatesnonsphere = np.ones([height_bound_vec.shape[0],diam.shape[0]])*np.nan
        m_MCSB_quantiles = np.ones([height_bound_vec.shape[0],len(quantiles)])*np.nan
        
        for category in ["snowjplatesnonsphere"]: #ATTENTION: this so far only works for snow; CHECK the coefficients in __postprocess_SB.init_particles before introducing other categories
            ###calculate the percentiles for the SB run            
            for i_height in range(0,N_heights):
                #get the PSD from SB
                [Nd_snowjplatesnonsphere[i_height,:],Md_snowjplatesnonsphere[i_height,:]] = __postprocess_SB.calc_distribution_from_moments(twomom,category,diam,i_time=i_timestep,i_height=i_height)
                [Nm_snowjplatesnonsphere[i_height,:],Mm_snowjplatesnonsphere[i_height,:]] = __postprocess_SB.calc_fmass_distribution_from_moments(twomom,category,m_array,i_time=i_timestep,i_height=i_height) 
                #i_time=twomom["qi"].shape[0] means we take the last timestep

                #calculate the percentiles from the CDF of the number size distribution
                Md_CDF = np.cumsum(Md_snowjplatesnonsphere[i_height])
                Md_CDF_norm = Md_CDF/Md_CDF[-1]
                Mm_CDF = np.cumsum(Mm_snowjplatesnonsphere[i_height])
                Mm_CDF_norm = Mm_CDF/Mm_CDF[-1]
                
                for i_quantile,quantile_now in enumerate(quantiles):
                    try: #there might be heights with no too low concentration
                        D_SB_quantiles[i_height,i_quantile] = diam[np.where(Md_CDF_norm>quantile_now)[0][0]] #get first diameter exceeding the quantile in the normed CDF
                        m_SB_quantiles[i_height,i_quantile] = m_array[np.where(Mm_CDF_norm>quantile_now)[0][0]] #get first mass exceeding the quantile in the normed CDF
                    except:
                        D_SB_quantiles[i_height,i_quantile] = np.nan
                        m_SB_quantiles[i_height,i_quantile] = np.nan

                #from IPython.core.debugger import Tracer ; Tracer()()
            ###calculate the percentiles for the MC moments with the SB parameters
            for i_height in range(0,height_bound_vec.shape[0]):    
                #get the PSD (N(D)) with the MC moments but the SB size distribution parameters
                [NdMCSB_snowjplatesnonsphere[i_height,:],MdMCSB_snowjplatesnonsphere[i_height,:]] = __postprocess_SB.calc_distribution_from_moments(twomomMC,category,diam,i_time=0,i_height=i_height) #i_time=0 -> the time has already been selected before 
                #get the PSD (N(m)) with the MC moments but the SB size distribution parameters
                [NmMCSB_snowjplatesnonsphere[i_height,:],MmMCSB_snowjplatesnonsphere[i_height,:]] = __postprocess_SB.calc_fmass_distribution_from_moments(twomomMC,category,m_array,i_time=0,i_height=i_height)
                
                #calculate the percentiles from the CDF of the number size distribution
                #for N(D)
                Md_CDF_MCSB = np.cumsum(MdMCSB_snowjplatesnonsphere[i_height])
                Md_CDF_MCSB_norm = Md_CDF_MCSB/Md_CDF_MCSB[-1]
                #for N(m)
                Mm_CDF_MCSB = np.cumsum(MmMCSB_snowjplatesnonsphere[i_height])
                Mm_CDF_MCSB_norm = Mm_CDF_MCSB/Mm_CDF_MCSB[-1]
                #if Mm_CDF_MCSB[-1]>0:
                #    from IPython.core.debugger import Tracer ; Tracer()()

                for i_quantile,quantile_now in enumerate(quantiles):
                    try: #there might be heights with no too low concentration
                        D_MCSB_quantiles[i_height,i_quantile] = diam[np.where(Md_CDF_MCSB_norm>quantile_now)[0][0]] #get first diameter exceeding the quantile in the normed CDF
                        m_MCSB_quantiles[i_height,i_quantile] = m_array[np.where(Mm_CDF_MCSB_norm>quantile_now)[0][0]]
                    except:
                        D_MCSB_quantiles[i_height,i_quantile] = np.nan
                        m_MCSB_quantiles[i_height,i_quantile] = np.nan


        #from IPython.core.debugger import Tracer ; Tracer()()
        

        ######
        ##  plot the percentiles of the simulated diameters
        #####
        i_ax = 0
        for key in ["diam","vt","m_tot"]:

            #axes[i_ax].semilogx(SP_percentiles[key][2,:],height_cent_vec,color=np.array(['b','r','g'])[i_sensMC+i_sensMCfallspeed])
            for i_quantile,quantile_now in enumerate(quantiles):
                axes[i_ax].plot(SP_percentiles[key][i_quantile,:],height_cent_vec,color=np.array(['r','b','g','y'])[0],linestyle=linestyle_quantiles[i_quantile],linewidth=linewidth_quantiles[i_quantile])
                if key=="diam": #TODO: vterm quantiles are not implemented yet for the SB scheme
                    axes[i_ax].plot(D_MCSB_quantiles[:,i_quantile],height_bound_vec,color=np.array(['r','b','g','y'])[1],linestyle=linestyle_quantiles[i_quantile],linewidth=linewidth_quantiles[i_quantile])#percentiles from MCrun with SB PSD coefficients
                    axes[i_ax].plot(D_SB_quantiles[:,i_quantile],twomom["heights"],color=np.array(['r','b','g','y'])[2],linestyle=linestyle_quantiles[i_quantile],linewidth=linewidth_quantiles[i_quantile]) #percentiles from SB runs
                    #from IPython.core.debugger import Tracer ; Tracer()()
                if key=="m_tot": #TODO: vterm quantiles are not implemented yet for the SB scheme
                    axes[i_ax].plot(m_MCSB_quantiles[:,i_quantile],height_bound_vec,color=np.array(['r','b','g','y'])[1],linestyle=linestyle_quantiles[i_quantile],linewidth=linewidth_quantiles[i_quantile])#percentiles from MCrun with SB PSD coefficients
                    axes[i_ax].plot(m_SB_quantiles[:,i_quantile],twomom["heights"],color=np.array(['r','b','g','y'])[2],linestyle=linestyle_quantiles[i_quantile],linewidth=linewidth_quantiles[i_quantile])
                    
                    
                if key=="diam":
                    axes[i_ax].set_xscale("log")
                elif key=="vt":
                    axes[i_ax].set_xscale("linear")
                if key=="m_tot":
                    axes[i_ax].set_xscale("log")
            if key=="diam":
                #make labels
                axes[i_ax].set_xlabel("diameter D / m")

                #change the axis
                axes[i_ax].set_xlim([1e-4,4e-2]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                axes[i_ax].set_ylim([0,model_top]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                axes[i_ax].grid(b=True,which="both",axis="both")
            elif key=="vt":
                #make labels
                axes[i_ax].set_xlabel("velocity v / m  s-1")
                axes[i_ax].set_ylabel("height / m" ) #TODO: plot also the untis of these properties

                #change the axis
                axes[i_ax].set_xlim([0,2.5]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                axes[i_ax].set_ylim([0,model_top]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                axes[i_ax].grid(b=True,which="both",axis="both")
            elif key=="m_tot":
                #make labels
                axes[i_ax].set_xlabel("mass m / kg")
                axes[i_ax].set_ylabel("height / m" ) #TODO: plot also the untis of these properties

                #change the axis
                axes[i_ax].set_xlim([1e-12,1e-6]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                axes[i_ax].set_ylim([0,model_top]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                axes[i_ax].grid(b=True,which="both",axis="both")
        
            i_ax+=1
            
        ####
        #plot the difference between the percentiles
        ####
        #i_ax=2
        for key in ["diam","vt"]:

            #axes[i_ax].semilogx(SP_percentiles[key][2,:],height_cent_vec,color=np.array(['b','r','g'])[i_sensMC+i_sensMCfallspeed])
            axes[i_ax].plot(SP_percentiles[key][2,:]-SP_percentiles[key][0,:],height_cent_vec,color=np.array(['r','b','g','y'])[0],linestyle="-",linewidth=1.5)
            if key=="diam": #TODO: vterm quantiles are not implemented yet for the SB scheme
                axes[i_ax].plot(D_MCSB_quantiles[:,2]-D_MCSB_quantiles[:,0],height_bound_vec,color=np.array(['r','b','g','y'])[1],linestyle="-",linewidth=1.5)
                axes[i_ax].plot(D_SB_quantiles[:,2]-D_SB_quantiles[:,0],twomom["heights"],color=np.array(['r','b','g','y'])[2],linestyle="-",linewidth=1.5)

            if key=="diam":
                axes[i_ax].set_xscale("linear")
            elif key=="vt":
                axes[i_ax].set_xscale("linear")

            if key=="diam":
                #make labels
                axes[i_ax].set_xlabel(r"$\Delta$ diameter D / m")

                #change the axis
                #axes[i_ax].set_xlim([1e-4,4e-2]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                axes[i_ax].set_ylim([0,model_top]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                axes[i_ax].grid(b=True,which="both",axis="both")
            elif key=="vt":
                #make labels
                axes[i_ax].set_xlabel(r"$\Delta$ velocity m s-1 / m")
                axes[i_ax].set_ylabel("height / m" ) #TODO: plot also the untis of these properties

                #change the axis
                #axes[i_ax].set_xlim([0,1.0]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                axes[i_ax].set_ylim([0,model_top]) #np.array([1e-5,1e-4,2.0,2.0,2.0])[i_prop]])
                axes[i_ax].grid(b=True,which="both",axis="both")
        
            i_ax+=1

#add labels and legend
for i_ax,ax_now in enumerate(axes):
    #generate some dummy labels indicating mean and quantiles
    if i_ax<2: #absolute values
        ax_now.plot(np.nan,np.nan,color="k",linestyle="-",linewidth=2.0,label="median")
        ax_now.plot(np.nan,np.nan,color="k",linestyle="--",linewidth=0.5,label="quantiles(mass): " + str(quantiles[0]) + " and " + str(quantiles[-1]))
    else: #differences
        ax_now.plot(np.nan,np.nan,color="k",linestyle="-",linewidth=2.0,label="diff of quantiles: " + str(quantiles[0]) + " and " + str(quantiles[-1]))
        
    #ax_now.plot(np.nan,np.nan,color="r",linestyle="-",linewidth=2.0,label="Boehm")
    ax_now.plot(np.nan,np.nan,color="r",linestyle="-",linewidth=2.0,label="Atlas in MC")
    ax_now.plot(np.nan,np.nan,color="b",linestyle="-",linewidth=2.0,label="SB with MC moments")
    ax_now.plot(np.nan,np.nan,color="g",linestyle="-",linewidth=2.0,label="SB")
    #ax_now.plot(np.nan,np.nan,color="y",linestyle="-",linewidth=2.0,label="power law SB")


    #show legend
    ax_now.legend() 
    ax_now.set_ylabel("height / m" ) #TODO: plot also the untis of these properties  
############
#save figure
############
plt.tight_layout()
output_string_folder = '/home/mkarrer/Dokumente/plots/SBMC_percentiles/' + experiment
if not os.path.exists(output_string_folder): #create direktory if it does not exists
    os.makedirs(output_string_folder)
out_filestring = "Dv_perc_" + McSnow_geom_list_str + separated_by_sensruns_onestring + '_' + testcase
plt.savefig(output_string_folder + out_filestring + '.pdf', dpi=400)
plt.savefig(output_string_folder + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + output_string_folder + out_filestring + '.pdf'
subprocess.Popen(['evince',output_string_folder + out_filestring + '.pdf'])