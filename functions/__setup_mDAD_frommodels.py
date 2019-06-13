'''
functions which setup mass-diameter and area-diameter relationships as used in the microphysic schemes /models (SB, P3, McSnow)
'''
import numpy as np
from numpy import pi, r_
import matplotlib.pyplot as plt
from scipy import optimize
import sys
import __postprocess_McSnow
import __postprocess_SB
import __fallspeed_relations
#from IPython.core.debugger import Tracer ; Tracer()()

def setup_diam_arrays(n_steps = 1001,min_D = 1e-4, max_D = 1e-2):
    '''
    setup a common diameter array
    INPUT:  n_steps: number of diameters in the array
            min_D: minimal diameter (in m)
            max_D: maximum diameter (in m)
    OUTPU:
            log-spaced diameter array
    '''
    
    D_array = np.logspace(np.log10(min_D),np.log10(max_D),n_steps)
    
    return D_array


def get_model_mDADs(model="MC"):
    '''
    get the mass-diameter and area-diameter relationships as assumed from the model
    INPUT: model: model from which the vD is returned
    
    OUTPUT: mDAD_dict: dictionary which contains the diameter and corresponding mass and area arrays for different models
    '''

    mDAD_dict = dict() #initialize dictionary which contains diameter and corresponding mass and area arrays for different models
    
    mDAD_dict["model_conf_name"] = model #name of the model and this specific configuration
    
    #for model_now in modellist:
    if model=="MC":
        
        mDAD_dict["num_piecewise_fit"] = 2 #number of different m-D and A-D fits which are piecewise connected
        
        #get parameter from McSnow in the state of Brdar&Seifert2018
        mth,unr_alf,unr_bet,rhoi,rhol,Dth,unr_sig,unr_gam,sph_sig,sph_gam =  __postprocess_McSnow.return_parameter_mD_AD_rel("1d__xi") #the argument is a workaround to get the standard McSnow settings
        sph_alf=np.pi/6.*rhoi
        sph_bet=3.
        print "used constants (control MC)"
        print "mth,               unr_alf,    unr_bet,rhoi,  rhol,       Dth,     unr_sig,     unr_gam,  sph_sig,    sph_gam"
        print mth,unr_alf,unr_bet,rhoi,rhol,Dth,unr_sig,unr_gam,sph_sig,sph_gam
        
        #copy constants to dictionary and replace name following a general naming convention
        #general names are following this naming convention: if piecewise: m1=a1Db1; m2=a2Db2; ... A1=c1D**d1,... Dth1,Dth2,...mth1,mth2,... (limiter between piecewise linear functions in diameter and mass space)
        for general_name,varlocal in zip(('a1','b1','a2', 'b2','c1','d1','c2','d2','Dth1','mth1'),('sph_alf','sph_bet','unr_alf', 'unr_bet','sph_sig', 'sph_gam','unr_sig', 'unr_gam','Dth','mth')):
            print varlocal,general_name,locals()[varlocal]
            mDAD_dict[general_name] = locals()[varlocal]
    elif model=="MCJagg_old": #parameters which were derived from Jussis aggregate models (using dendrites); found in McSnow code
        
        mDAD_dict["num_piecewise_fit"] = 2 #number of different m-D and A-D fits which are piecewise connected

        
        #get parameter from McSnow in the state of Brdar&Seifert2018
        mth,unr_alf,unr_bet,rhoi,rhol,Dth,unr_sig,unr_gam,sph_sig,sph_gam =  __postprocess_McSnow.return_parameter_mD_AD_rel("1d_xi") #the argument is a workaround to get the standard McSnow settings
        sph_alf=np.pi/6.*rhoi
        sph_bet=3.

        #overwrite variables!!!
        #add other used parameter here
        unr_alf = 0.01243; unr_bet = 2.00000 #found commented in McSnows mo_mass_diam.f90
        unr_sig = 0.05625; unr_gam = 1.81000
        print "##################"
        print "added parameter"
        print "unr_alf,unr_bet,unr_sig,unr_gam"
        print unr_alf,"         ",unr_bet,"             ",unr_sig,"           ",unr_gam
        #calculate Dth and mth for the added coefficients
        Dth=(sph_alf/unr_alf)**(1./(unr_bet-3.)) #D_th_Jaggdent=(sph_alf/unr_alf_Jaggdent)**(1./(unr_bet_Jaggdent-3.))
        mth=sph_alf*Dth**3 #m_th_Jaggdent=sph_alf*D_th_Jaggdent**3
        print "Dth,m_th"
        print Dth,mth
            
        #copy constants to dictionary and replace name following a general naming convention
        #general names are following this naming convention: if piecewise: m1=a1Db1; m2=a2Db2; ... A1=c1D**d1,... Dth1,Dth2,...mth1,mth2,... (limiter between piecewise linear functions in diameter and mass space)
        for general_name,varlocal in zip(('a1','b1','a2', 'b2','c1','d1','c2','d2','Dth1','mth1'),('sph_alf','sph_bet','unr_alf', 'unr_bet','sph_sig', 'sph_gam','unr_sig', 'unr_gam','Dth','mth')):
            print varlocal,general_name,locals()[varlocal]
            mDAD_dict[general_name] = locals()[varlocal]
        
        
    elif model=='SBcloudice':
        
        mDAD_dict["num_piecewise_fit"] = 1 #number of different m-D and A-D fits which are piecewise connected

        
        #get parameters from the SB-categories
        cloud_water,rain,cloud_ice,snow,graupel,hail = __postprocess_SB.init_particles()
        cloud_water,rain,cloud_ice,snow,graupel,hail = __postprocess_SB.convert_Nm_to_ND(cloud_water,rain,cloud_ice,snow,graupel,hail)
        # terminal velocity parameter
        a_vel_icecosmo5 = 2.77e1; b_vel_icecosmo5 = 0.215790

        
        print "##################"
        print "SB parameters: cloud ice"
        
        #copy constants to dictionary and replace name following a general naming convention
        #general names are following this naming convention: if piecewise: m1=a1Db1; m2=a2Db2; ... A1=c1D**d1,... Dth1,Dth2,...mth1,mth2,... (limiter between piecewise linear functions in diameter and mass space)
        for general_name,varlocal in zip(('a1','b1','am1','bm1',"aterm","bterm"),('cloud_ice.a_ms','cloud_ice.b_ms','cloud_ice.a_geo','cloud_ice.b_geo','a_vel_icecosmo5','b_vel_icecosmo5')): 
            #the variables are actually in a class, so we convert then here to normal vvariables
            exec("varnow = " + varlocal)
            print general_name,varnow
            mDAD_dict[general_name] = varnow #locals()[varlocal]
    
    elif model=='SBcloudice_uli':
        
        mDAD_dict["num_piecewise_fit"] = 1 #number of different m-D and A-D fits which are piecewise connected

        
        #get parameters from the SB-categories
        cloud_water,rain,cloud_ice,snow,graupel,hail = __postprocess_SB.init_particles()
        cloud_water,rain,cloud_ice,snow,graupel,hail = __postprocess_SB.convert_Nm_to_ND(cloud_water,rain,cloud_ice,snow,graupel,hail)
        # terminal velocity parameters
        a_vel_ulis_ice = 2.60e1; b_vel_ulis_ice = 0.215790

        print "##################"
        print "SB parameters: cloud ice"
        
        #copy constants to dictionary and replace name following a general naming convention
        #general names are following this naming convention: if piecewise: m1=a1Db1; m2=a2Db2; ... A1=c1D**d1,... Dth1,Dth2,...mth1,mth2,... (limiter between piecewise linear functions in diameter and mass space)
        for general_name,varlocal in zip(('a1','b1','am1','bm1',"aterm","bterm"),('cloud_ice.a_ms','cloud_ice.b_ms','cloud_ice.a_geo','cloud_ice.b_geo','a_vel_ulis_ice','b_vel_ulis_ice')): 
            #the variables are actually in a class, so we convert then here to normal vvariables
            exec("varnow = " + varlocal)
            print varnow
            mDAD_dict[general_name] = varnow #locals()[varlocal]
            
    elif model=='SBcloudice_Atlas':
        
        mDAD_dict["num_piecewise_fit"] = 1 #number of different m-D and A-D fits which are piecewise connected

        
        #get parameters from the SB-categories
        cloud_water,rain,cloud_ice,snow,graupel,hail = __postprocess_SB.init_particles()
        cloud_water,rain,cloud_ice,snow,graupel,hail = __postprocess_SB.convert_Nm_to_ND(cloud_water,rain,cloud_ice,snow,graupel,hail)
        # terminal velocity parameters
        a_vel_icecosmo5nonsphere = 1.860; b_vel_icecosmo5nonsphere = 1.872; c_vel_icecosmo5nonsphere = 9.325e2 #fits to Ulis ice (despite the name) #this is as a function of the molten diameter

        print "##################"
        print "SB parameters: cloud ice"
        
        #copy constants to dictionary and replace name following a general naming convention
        #general names are following this naming convention: if piecewise: m1=a1Db1; m2=a2Db2; ... A1=c1D**d1,... Dth1,Dth2,...mth1,mth2,... (limiter between piecewise linear functions in diameter and mass space)
        for general_name,varlocal in zip(('a1','b1','am1','bm1',"aterm","bterm","cterm"),('cloud_ice.a_ms','cloud_ice.b_ms','cloud_ice.a_geo','cloud_ice.b_geo','a_vel_icecosmo5nonsphere','b_vel_icecosmo5nonsphere','c_vel_icecosmo5nonsphere')): 
            #the variables are actually in a class, so we convert then here to normal vvariables
            exec("varnow = " + varlocal)
            print varnow
            mDAD_dict[general_name] = varnow #locals()[varlocal]
            
    
    
    elif model=='SBsnow':
    
        mDAD_dict["num_piecewise_fit"] = 1 #number of different m-D and A-D fits which are piecewise connected

        
        #get parameters from the SB-categories
        cloud_water,rain,cloud_ice,snow,graupel,hail = __postprocess_SB.init_particles()
        cloud_water,rain,cloud_ice,snow,graupel,hail = __postprocess_SB.convert_Nm_to_ND(cloud_water,rain,cloud_ice,snow,graupel,hail)
        
        print "##################"
        print "SB parameters: snow"
        # terminal velocity parameters
        a_vel_snowSBB = 8.294000; b_vel_snowSBB = 0.125000

        #copy constants to dictionary and replace name following a general naming convention
        #general names are following this naming convention: if piecewise: m1=a1Db1; m2=a2Db2; ... A1=c1D**d1,... Dth1,Dth2,...mth1,mth2,... (limiter between piecewise linear functions in diameter and mass space) #here we have additionally the coefficients as a function of mass and the fall speed coefficients, but no area (as in the scheme itself)
        for general_name,varlocal in zip(('a1','b1','am1','bm1',"aterm","bterm"),('snow.a_ms','snow.b_ms','snow.a_geo','snow.b_geo','a_vel_snowSBB','b_vel_snowSBB')): 
            #the variables are actually in a class, so we convert then here to normal vvariables
            exec("varnow = " + varlocal)
            print varnow
            mDAD_dict[general_name] = varnow #locals()[varlocal]
            
    elif model=='SBsnow_Atlas':
    
        mDAD_dict["num_piecewise_fit"] = 1 #number of different m-D and A-D fits which are piecewise connected

        
        #get parameters from the SB-categories
        cloud_water,rain,cloud_ice,snow,graupel,hail = __postprocess_SB.init_particles()
        cloud_water,rain,cloud_ice,snow,graupel,hail = __postprocess_SB.convert_Nm_to_ND(cloud_water,rain,cloud_ice,snow,graupel,hail)
        # terminal velocity parameters
        a_vel_snowSBBnonsphere = 1.271; b_vel_snowSBBnonsphere = 1.252; c_vel_snowSBBnonsphere = 3.698e3

        print "##################"
        print "SB parameters: snow"
        
        #copy constants to dictionary and replace name following a general naming convention
        #general names are following this naming convention: if piecewise: m1=a1Db1; m2=a2Db2; ... A1=c1D**d1,... Dth1,Dth2,...mth1,mth2,... (limiter between piecewise linear functions in diameter and mass space) #here we have additionally the coefficients as a function of mass and the fall speed coefficients, but no area (as in the scheme itself)
        for general_name,varlocal in zip(('a1','b1','am1','bm1',"aterm","bterm","cterm"),('cloud_ice.a_ms','cloud_ice.b_ms','cloud_ice.a_geo','cloud_ice.b_geo','a_vel_snowSBBnonsphere','b_vel_snowSBBnonsphere','c_vel_snowSBBnonsphere')): 
            #the variables are actually in a class, so we convert then here to normal vvariables
            exec("varnow = " + varlocal)
            print varnow
            mDAD_dict[general_name] = varnow #locals()[varlocal]
    elif model=='P3':
        
        mDAD_dict["num_piecewise_fit"] = 2 #number of different m-D and A-D fits which are piecewise connected

        
        #get parameter from P3 (TAKEN FROM THE P3 LOOKUP-TABLE)
        ##############################################
        ####1:spherical ice m(D)=cs1*D**ds1###########
        ##############################################
        #upper limit: dcrit=6.71e-05=67mu m (constant)
        ##############################################
        dcrit_P3 = 6.71e-05
        mcrit_P3 = 1000. * dcrit_P3**3 
        #set up m-D relationship for solid ice with D < Dcrit
        cs1_P3  = np.pi*1./6.*900.
        ds1_P3  = 3.
        #########################################
        ####2:dense nonspherical m(D)=cs*D**ds###
        #########################################	
        #lower limit: dcrit, upper limit:dcrits (if we consider rimed##
        #########################################
        # Brown and Francis (1995)
        ds_P3	=	1.9
        # cs = 0.01855 # original (pre v2.3), based on assumption of Dmax
        cs_P3	=	0.0121 # scaled value based on assumtion of Dmean from Hogan et al. 2012, JAMC

        ################
        #A-D-parameters
        ################
        ####1:spherical ice: 	 A(D)=aas1*D**bas1
        aas1_P3=np.pi/4.
        bas1_P3=2

        ###2:dense nonspherical: A(D)=aas2*D**bas2
        bas2_P3=1.88
        aas2_P3=0.2285*100.**bas2_P3/(100.**2)
        
        print "##################"
        print "P3 parameters: snow"
        
        #copy constants to dictionary and replace name following a general naming convention
        #general names are following this naming convention: if piecewise: m1=a1Db1; m2=a2Db2; ... A1=c1D**d1,... Dth1,Dth2,...mth1,mth2,... (limiter between piecewise linear functions in diameter and mass space)
        for general_name,varlocal in zip(('a1','b1','a2', 'b2','c1','d1','c2','d2','Dth1','mth1'),('cs1_P3','ds1_P3','cs_P3', 'ds_P3','aas1_P3', 'bas1_P3','aas2_P3', 'bas2_P3','dcrit_P3','mcrit_P3')):
            print varlocal,general_name,locals()[varlocal]
            mDAD_dict[general_name] = locals()[varlocal]
            
    elif model=='GSFCcloudice':
    
        mDAD_dict["num_piecewise_fit"] = 1 #number of different m-D and A-D fits which are piecewise connected

         
        print "##################"
        print "GSFC (Godard scheme) parameters: cloud ice"
        # terminal velocity parameters
        a_vel = 1.30493; b_vel = 0.11 #ATTENTION: check these parameters

        #copy constants to dictionary and replace name following a general naming convention
        #general names are following this naming convention: if piecewise: m1=a1Db1; m2=a2Db2; ... A1=c1D**d1,... Dth1,Dth2,...mth1,mth2,... (limiter between piecewise linear functions in diameter and mass space) #here we have additionally the coefficients as a function of mass and the fall speed coefficients, but no area (as in the scheme itself)
        for general_name,varlocal in zip(('aterm','bterm'),('a_vel','b_vel')): 
            #the variables are actually in a class, so we convert then here to normal vvariables
            exec("varnow = " + varlocal)
            print varnow
            mDAD_dict[general_name] = varnow
            
    elif model=='GSFCsnow':
    
        mDAD_dict["num_piecewise_fit"] = 1 #number of different m-D and A-D fits which are piecewise connected

         
        print "##################"
        print "GSFC (Godard scheme) parameters: snow"
        # terminal velocity parameters
        a_vel = 1.30493; b_vel = 0.11 #ATTENTION: check these parameters

        #copy constants to dictionary and replace name following a general naming convention
        #general names are following this naming convention: if piecewise: m1=a1Db1; m2=a2Db2; ... A1=c1D**d1,... Dth1,Dth2,...mth1,mth2,... (limiter between piecewise linear functions in diameter and mass space) #here we have additionally the coefficients as a function of mass and the fall speed coefficients, but no area (as in the scheme itself)
        for general_name,varlocal in zip(('aterm','bterm'),('a_vel','b_vel')): 
            #the variables are actually in a class, so we convert then here to normal vvariables
            exec("varnow = " + varlocal)
            print varnow
            mDAD_dict[general_name] = varnow

    elif model=='morr2mom_cloudice':
    
        mDAD_dict["num_piecewise_fit"] = 1 #number of different m-D and A-D fits which are piecewise connected
         
        print "##################"
        print "Morrison2mom parameters: ice"
        # terminal velocity parameters
        a_vel = 700.; b_vel = 1.0

        #copy constants to dictionary and replace name following a general naming convention
        #general names are following this naming convention: if piecewise: m1=a1Db1; m2=a2Db2; ... A1=c1D**d1,... Dth1,Dth2,...mth1,mth2,... (limiter between piecewise linear functions in diameter and mass space) #here we have additionally the coefficients as a function of mass and the fall speed coefficients, but no area (as in the scheme itself)
        for general_name,varlocal in zip(('aterm','bterm'),('a_vel','b_vel')): 
            #the variables are actually in a class, so we convert then here to normal vvariables
            exec("varnow = " + varlocal)
            print varnow
            mDAD_dict[general_name] = varnow #locals()[varlocal]

    elif model=='morr2mom_snow':
    
        mDAD_dict["num_piecewise_fit"] = 1 #number of different m-D and A-D fits which are piecewise connected

         
        print "##################"
        print "Morrison2mom parameters: snow"
        # terminal velocity parameters
        a_vel = 11.72; b_vel = 0.41

        #copy constants to dictionary and replace name following a general naming convention
        #general names are following this naming convention: if piecewise: m1=a1Db1; m2=a2Db2; ... A1=c1D**d1,... Dth1,Dth2,...mth1,mth2,... (limiter between piecewise linear functions in diameter and mass space) #here we have additionally the coefficients as a function of mass and the fall speed coefficients, but no area (as in the scheme itself)
        for general_name,varlocal in zip(('aterm','bterm'),('a_vel','b_vel')): 
            #the variables are actually in a class, so we convert then here to normal vvariables
            exec("varnow = " + varlocal)
            print varnow
            mDAD_dict[general_name] = varnow #locals()[varlocal]
            
    elif model=='thompson_cloudice':
    
        mDAD_dict["num_piecewise_fit"] = 1 #number of different m-D and A-D fits which are piecewise connected
         
        print "##################"
        print "Thompson parameters: ice"
        # terminal velocity parameters
        a_vel = 1847.5; b_vel = 1.0

        #copy constants to dictionary and replace name following a general naming convention
        #general names are following this naming convention: if piecewise: m1=a1Db1; m2=a2Db2; ... A1=c1D**d1,... Dth1,Dth2,...mth1,mth2,... (limiter between piecewise linear functions in diameter and mass space) #here we have additionally the coefficients as a function of mass and the fall speed coefficients, but no area (as in the scheme itself)
        for general_name,varlocal in zip(('aterm','bterm'),('a_vel','b_vel')): 
            #the variables are actually in a class, so we convert then here to normal vvariables
            exec("varnow = " + varlocal)
            print varnow
            mDAD_dict[general_name] = varnow #locals()[varlocal]

    elif model=='thompson_snow':
    
        mDAD_dict["num_piecewise_fit"] = 1 #number of different m-D and A-D fits which are piecewise connected

         
        print "##################"
        print "Thompson parameters: snow"
        # terminal velocity parameters
        a_velthom = 40.0; b_velthom = 0.55; c_velthom = 100

        #copy constants to dictionary and replace name following a general naming convention
        #general names are following this naming convention: if piecewise: m1=a1Db1; m2=a2Db2; ... A1=c1D**d1,... Dth1,Dth2,...mth1,mth2,... (limiter between piecewise linear functions in diameter and mass space) #here we have additionally the coefficients as a function of mass and the fall speed coefficients, but no area (as in the scheme itself)
        for general_name,varlocal in zip(('atermthom','btermthom', 'ctermthom'),('a_velthom','b_velthom','c_velthom')): 
            #the variables are actually in a class, so we convert then here to normal vvariables
            exec("varnow = " + varlocal)
            print varnow
            mDAD_dict[general_name] = varnow #locals()[varlocal]
    
    elif model=="Jussis_needles":
        mDAD_dict["num_piecewise_fit"] = 1 #number of different m-D and A-D fits which are piecewise connected

        
        #theoretical parameters from Jussis needles (monomers)
        #m-D
        a=4.15e-3; b=1.87
        #A-D
        c=3.18e-3; d=1.44
        
        print "##################"
        print "Jussis needles parameters: "
        
        #copy constants to dictionary and replace name following a general naming convention
        #general names are following this naming convention: if piecewise: m1=a1Db1; m2=a2Db2; ... A1=c1D**d1,... Dth1,Dth2,...mth1,mth2,... (limiter between piecewise linear functions in diameter and mass space) #here we have additionally the coefficients as a function of mass and the fall speed coefficients, but no area (as in the scheme itself)
        for general_name,varlocal in zip(('a1','b1','c1','d1'),('a','b','c','d')): 
            #the variables are actually in a class, so we convert then here to normal vvariables
            exec("varnow = " + varlocal)
            print varnow
            mDAD_dict[general_name] = varnow #locals()[varlocal]
        
    return mDAD_dict
        
def calc_area_mass_vterm_arrays(diam_array,mDADvD_dict):
    '''
    calculate the area, mass and terminal velocity from the constants
    INPUT:  diam_array: array of diameters
            mDADvD_dict: dictionary which contains the model name (and its configuration), coefficients and the arrays
    OUTPU:  mDADvD_dict: same dictionary as input , but with the added arrays of area, mass and velocity
    '''
    ###
    #calculate the mass, area and terminal speed corresponding to the diam_array
    ###

    #seperate piecewise fitted power-laws from full-range fits
    if mDADvD_dict["num_piecewise_fit"]==1:

        if "c1" in mDADvD_dict.keys(): #we have explicit information about the A-D relation and do not need a v-D fits
            mDADvD_dict["m(D_array)"] = mDADvD_dict["a1"]*diam_array**mDADvD_dict["b1"]
            mDADvD_dict["A(D_array)"] = mDADvD_dict["c1"]*diam_array**mDADvD_dict["d1"]
            
            #now: calculate terminal velocity from the D,m,A-arrays
            for fallspeedmodel in ['HW10','mitch_heym','bohm','KC05']:#do this for all fall speed models in general
                mDADvD_dict["v_" + fallspeedmodel + "(D_array)"] = __fallspeed_relations.calc_vterm(fallspeedmodel,mDADvD_dict["m(D_array)"],diam_array,mDADvD_dict["A(D_array)"])
        elif "cterm" in mDADvD_dict.keys(): #then we have an Atlas-type
            D_array_molten = np.logspace(-12,-1,1001) #span a wide range of molten diameters here
            m_from_moltenD = 1000.*np.pi/6.*D_array_molten**3

            mDADvD_dict["D_max_from_moltenD"] = mDADvD_dict["am1"]*m_from_moltenD**mDADvD_dict["bm1"] #calculate maximum dimension from molten diameter
            mDADvD_dict["v(D_array)"] = mDADvD_dict["aterm"]-mDADvD_dict["bterm"]*np.exp(-mDADvD_dict["cterm"]*D_array_molten) #these are formulated as a melted diameter!!
            #from IPython.core.debugger import Tracer ; Tracer()()

        elif "a1" in mDADvD_dict.keys() and "aterm" in mDADvD_dict.keys(): #power-law which has to be converted from Nm to ND space
            mDADvD_dict["m(D_array)"] = mDADvD_dict["a1"]*diam_array**mDADvD_dict["b1"]
            mDADvD_dict["v(D_array)"] = mDADvD_dict["aterm"]*(mDADvD_dict["a1"]*diam_array**mDADvD_dict["b1"])**mDADvD_dict["bterm"] #aterm and bterm are formulated as a function of mass
        elif "aterm" in mDADvD_dict.keys(): #just the powerlaw in v-D formulation
            
            mDADvD_dict["v(D_array)"] = mDADvD_dict["aterm"]*diam_array**mDADvD_dict["bterm"] #aterm and bterm are formulated as a function of mass
        elif "atermthom": #the Thompson scheme has a special formulation for snow (and rain)
            mDADvD_dict["v(D_array)"] = mDADvD_dict["atermthom"]*diam_array**mDADvD_dict["btermthom"]*np.exp(-mDADvD_dict["ctermthom"]*diam_array) #aterm and bterm are formulated as a function of mass

            

                
    elif mDADvD_dict["num_piecewise_fit"]==2:
        #initialize arrays
        mDADvD_dict["m(D_array)"] = np.zeros(diam_array.shape)
        mDADvD_dict["A(D_array)"] = np.zeros(diam_array.shape)
        mDADvD_dict["v(D_array)"] = np.zeros(diam_array.shape)

        
        for i, D in enumerate(diam_array):
            if D<=mDADvD_dict["Dth1"]:
                mDADvD_dict["m(D_array)"][i] = mDADvD_dict["a1"]*D**mDADvD_dict["b1"]
                mDADvD_dict["A(D_array)"][i] = mDADvD_dict["c1"]*D**mDADvD_dict["d1"]

            elif D>mDADvD_dict["Dth1"]:
                mDADvD_dict["m(D_array)"][i] = mDADvD_dict["a2"]*D**mDADvD_dict["b2"]
                mDADvD_dict["A(D_array)"][i] = mDADvD_dict["c2"]*D**mDADvD_dict["d2"]
                
        #now: calculate terminal velocity from the D,m,A-arrays
        for fallspeedmodel in ['HW10','mitch_heym','bohm','KC05']:#do this for all fall speed models in general
            mDADvD_dict["v_" + fallspeedmodel + "(D_array)"] = __fallspeed_relations.calc_vterm(fallspeedmodel,mDADvD_dict["m(D_array)"],diam_array,mDADvD_dict["A(D_array)"])
        

    return mDADvD_dict
    