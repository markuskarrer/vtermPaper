#import modulus
import numpy as np
import sys
from scipy.special import gamma
#from IPython.core.debugger import Tracer ; Tracer()()

# define functions to postprocess output from the SB-scheme

def calc_Dmean(twomom,category):
    '''
    this scripts calculates the mean diameter (based on the two moments) for the ICON-SB scheme as done in mo_2mom_mcrph_processes within ICON
    INPUT:  twomom: dictionary with (McSnows) two-moment variables
            category: hydrometeor category for which D_mean should be calculated
    '''
    ''' adapted from mo_2mom_mcrph_processes
    ! mean mass with limiters, Eq. (94) of SB2006
        ELEMENTAL FUNCTION particle_meanmass(this,q,n) RESULT(xmean)
        CLASS(particle), INTENT(in) :: this
        REAL(wp),        INTENT(in) :: q,n
        REAL(wp)                    :: xmean
        REAL(wp), PARAMETER         :: eps = 1e-20_wp

        xmean = MIN(MAX(q/(n+eps),this%x_min),this%x_max)
        END FUNCTION particle_meanmass
        ! mass-diameter relation, power law, Eq. (32) of SB2006
        ELEMENTAL FUNCTION particle_diameter(this,x) RESULT(D)
        CLASS(particle), INTENT(in) :: this
        REAL(wp),        INTENT(in) :: x
        REAL(wp)                    :: D

        D = this%a_geo * exp(this%b_geo*log(x))    ! v = a_geo * x**b_geo
        END FUNCTION particle_diameter
    '''
    #get hyd_type depending on category
    index_in_cat_array = ['c','i','r','s','g','h'].index(category)
    hyd_type = ['cloud_nue1mue1','rainSBB','snow_cosmo5','ice_cosmo5','graupelhail_cosmo5','hail_cosmo5'][index_in_cat_array] #ATTENTION: change this when something changes in SB
    
    if hyd_type=='cloud_nue1mue1':
        a_geo = 1.24e-1   #..a_geo..Koeff.
        b_geo = 0.333333   #..b_geo..Koeff. Geometrie = 1/3

        x_max = 2.60e-10
        x_min = 4.20e-15
    if hyd_type=='rainSBB':
        a_geo = 1.24e-1   #..a_geo..Koeff.
        b_geo = 0.333333   #..b_geo..Koeff. Geometrie = 1/3

        x_max = 3.00e-6
        x_min = 2.60e-10
    if hyd_type=='snowSBB':
        a_geo = 5.130000   #..a_geo..Koeff. Geometrie, x = 0.038*D**2
        b_geo = 0.500000   #..b_geo..Koeff. Geometrie = 1/2

        x_max = 2.00e-5
        x_min = 1.00e-10
    if hyd_type=='snow_cosmo5': #not used in HDCP2 version of SB2006
        a_geo = 2.400000   #..a_geo..Koeff. Geometrie, 
        b_geo = 0.455000   #..b_geo..Koeff. Geometrie 
        x_max = 2.00e-5
        x_min = 1.00e-10

    if hyd_type=='ice_cosmo5':
        a_geo = 0.835000   #..a_geo..Koeff. Geometrie, 
        b_geo = 0.390000   #..b_geo..Koeff. Geometrie 
        x_max = 1.00e-5
        x_min = 1.00e-12

    if hyd_type=='graupelhail_cosmo5':
        a_geo = 1.42e-01   #..a_geo..Koeff. Geometrie, 
        b_geo = 0.314000   #..b_geo..Koeff. Geometrie 
        x_max = 5.00e-4
        x_min = 1.00e-09

    if hyd_type=='hail_cosmo5':
        a_geo = 0.1366   #..a_geo..Koeff. Geometrie, 
        b_geo = 0.333333   #..b_geo..Koeff. Geometrie 
        x_max = 5.00e-4
        x_min = 2.60e-9
    #copy moments to q and N    
    q = twomom['q' + category]
    N = twomom['qn' + category]
    eps = 1e-20
    #calculate mean mass first
    xmean = np.where(q>0,np.minimum(np.maximum(q/(N+eps),x_min),x_max),np.nan)

    #calculate mean mass
    twomom["D_mean_" + category] = a_geo * np.exp(b_geo*np.log(xmean))    # v = a_geo * x**b_geo

    return twomom["D_mean_" + category]

def init_particles():
    #define particle class
    class particle(object):
            def __init__(self, **kwargs):
                    self.__dict__.update(kwargs)
    #define objects for particle types
    cloud_water = particle(nu_SB	=  1.,
                            mu_SB       = 1.0,
                            a_geo       =  0.124,
                            b_geo       =  0.33333,
                            xmax       = 2.60e-10, #& !..x_max..maximale Teilchenmasse D=80e-6m
                            xmin       = 4.20e-15, #& !..x_min..minimale Teilchenmasse D=2.e-6m
                            mixrat_var = 'qc',
                            numcon_var = 'qnc')
    rain	 = particle(nu_SB		=  0.0,
                            mu_SB        = 0.33333,
                            a_geo      =  0.124,
                            b_geo      =  0.33333,
                            xmax      = 3.00e-06,  #& !..x_max
                            xmin      = 2.60e-10,  #& !..x_min
                            mixrat_var = 'qr',
                            numcon_var = 'qnr')
    cloud_ice = particle(nu_SB	        =  0.,
                            mu_SB        = 0.3333,
                            a_geo      =  0.835,
                            b_geo      =  0.39,
                            xmax       =  1.00e-05, #& !..x_max..maximale Teilchenmasse D=???e-2m
                            xmin       =  1.00e-12, #& !..x_min..minimale Teilchenmasse D=200e-6m
                            mixrat_var = 'qi',
                            numcon_var = 'qni')
    snow      = particle(nu_SB	=  0,
                            mu_SB = 0.5,
                            a_geo  =  5.13,
                            b_geo  =  0.5,
                            xmax   =  2.00e-05, #& !..x_max..maximale Teilchenmasse
                            xmin   =  1.00e-10, #& !..x_min..minimale Teilchenmasse
                            mixrat_var = 'qs', #if set to 'qi': ATTENTION: set for testing 'qs',
                            numcon_var = 'qns') 
    graupel = particle(nu_SB	        =  1.0,
                    mu_SB        = 0.33333,
                    a_geo      =  0.142,
                    b_geo      =  0.314,
                    xmax       =  5.00e-04, #& #!..x_max..maximale Teilchenmasse
                    xmin       =  1.00e-09, #& !..x_min..minimale Teilchenmasse
                    mixrat_var = 'qg',
                    numcon_var = 'qng')
    hail = particle(nu_SB	        =  1.0,
                    mu_SB        =  0.33333,
                    a_geo      =  0.1366,
                    b_geo      =  0.3333333,
                    xmax       =  5.00e-04, #& !..x_max..maximale Teilchenmasse
                    xmin       =  2.60e-9,  #& !..x_min..minimale Teilchenmasse
                    mixrat_var = 'qh',
                    numcon_var = 'qnh')
    return cloud_water,rain,cloud_ice,snow,graupel,hail

def convert_Nm_to_ND(cloud_water,rain,cloud_ice,snow,graupel,hail):

    #convert from N(m) to N(D) space
    for n_cat,curr_cat in enumerate((cloud_water,rain,cloud_ice,snow,graupel,hail)):	
        curr_cat.a_ms = (1./curr_cat.a_geo)**(1./curr_cat.b_geo)
        curr_cat.b_ms = 1./curr_cat.b_geo
        curr_cat.mu = curr_cat.b_ms*curr_cat.nu_SB+curr_cat.b_ms-1
        curr_cat.gam = curr_cat.b_ms*curr_cat.mu_SB
        curr_cat.Dmax = (curr_cat.xmax/curr_cat.a_ms)**(1./curr_cat.b_ms)
        curr_cat.Dmin = (curr_cat.xmin/curr_cat.a_ms)**(1./curr_cat.b_ms)
    
    return cloud_water,rain,cloud_ice,snow,graupel,hail

def convert_ND_to_Nm_from_coeff(a_mD,b_mD,a_vD,b_vD):
    '''
    this function converts coefficients (m-D, A-D) from N(D) space to N(m) space (e.g. to use it in the Sb-scheme)
    INPUT:
        a_mD: a in mass-diameter relationship (m=a*D**b)
        b_mD: b in mass-diameter relationship (m=a*D**b)
        a_vD: a in velocity-diameter relationship (v=a*D**b)
        b_vD: b in velocity-diameter relationship (v=a*D**b)
    OUTPUT:
        a_mm: a in diameter-mass relationship (D=a*m**b)
        b_mm: b in diameter-mass relationship (D=a*m**b)
        a_vm: a in velocity-mass relationship (v=a*m**b)
        b_vm: b in velocity-mass relationship (v=a*m**b)
    '''

    #from IPython.core.debugger import Tracer ; Tracer()()

    #convert coefficients from N(D) to N(m) space
    b_mm = 1./b_mD
    a_mm = (1./a_mD)**(1./b_mD)
    a_vm = a_vD*(1./a_mD)**(b_vD/b_mD)
    b_vm = b_vD/b_mD
    
    return a_mm,b_mm,a_vm,b_vm

def calc_distribution_from_moments(twomom,category,d_ds,i_time=0,i_height=249):
    '''
    calculate the normalized number concentration (as a function of diameter) corresponding to moments of the SB06 categories
    INPUT:  twomom: dictionary containing the moments
            category: name of the category used (f.e. cloud_ice, snow,... (see below))
            d_ds: diameter array at which the number concentration should be evaluated
            i_time: timestep index of the entries in the twomom dictionary which should be analyzed
            i_height: height index of the entries in the twomom dictionary which should be analyzed
    '''
    cloud_water,rain,cloud_ice,snow,graupel,hail = init_particles() #get all parameters from the SB-categories
    
    [cloud_water,rain,cloud_ice,snow,graupel,hail] = convert_Nm_to_ND(cloud_water,rain,cloud_ice,snow,graupel,hail)

        
    #select the category
    if category=="icecosmo5":
        curr_cat = cloud_ice #if set to "
    elif category=="snowSBB":
        curr_cat = snow
    else:
        print category + " currently not implemented in __postprocess_SB.calc_distribution_from_moments"
    
    #for debugging: output the calculated converted parameters
    #print curr_cat.mixrat_var 
    #print 'a_ms,b_ms,mu,gam,Dmin,Dmax'
    #print curr_cat.a_ms,curr_cat.b_ms,curr_cat.mu,curr_cat.gam,curr_cat.Dmin,curr_cat.Dmax  

    ###
    #calculate the normalized number concentration 
    ###
    diam	=	d_ds
    #initialize and calculate bin width (del_diam)
    del_diam = np.zeros(diam.shape)
    del_diam[0:(diam.shape[0]-1)]=	diam[1:]-diam[0:(diam.shape[0]-1)]

    #copy the mass density and the number concentration to PAMTRA conventions
    q_h  =  twomom[curr_cat.mixrat_var][i_time,i_height]
    n_tot =  twomom[curr_cat.numcon_var][i_time,i_height]

    #calculate the distribution based on PAMTRA code
    #taken from PAMTRA make_dist_params.f90	
    work2 = gamma((curr_cat.mu + curr_cat.b_ms + 1.0) / curr_cat.gam)
    work3 = gamma((curr_cat.mu + 1.0) / curr_cat.gam)
    lam	=	(curr_cat.a_ms / q_h * n_tot * work2 / work3)**(curr_cat.gam / curr_cat.b_ms)
    N_0 = curr_cat.gam * n_tot / work3 * lam**((curr_cat.mu + 1.0) / curr_cat.gam)
    N_D	=	N_0*diam**curr_cat.mu*np.exp(-lam*diam**curr_cat.gam) #/ del_diam		#normalized number concentrations with *del_diam-> not normalized
    #apply diameter (mass) limits #or better dont because they just limit the mean mass not truncate the spectrum
    #N_D[np.where( diam < curr_cat.Dmin)] = 0.0 #N_D[diam<curr_cat.Dmin or
    #N_D[np.where( diam > curr_cat.Dmax)] = 0.0 #N_D[diam<curr_cat.Dmin or
    M_D = N_D * curr_cat.a_ms*diam**curr_cat.b_ms
    return N_D,M_D

def calc_fmass_distribution_from_moments(twomom,category,m_ds,i_time=0,i_height=249):
    '''
    calculate the normalized number concentration (as a function of mass) corresponding to moments of the SB06 categories
    INPUT:  twomom: dictionary containing the moments
            category: name of the category used (f.e. cloud_ice, snow,... (see below))
            d_ds: diameter array at which the number concentration should be evaluated
            i_time: timestep index of the entries in the twomom dictionary which should be analyzed
            i_height: height index of the entries in the twomom dictionary which should be analyzed
    '''
    cloud_water,rain,cloud_ice,snow,graupel,hail = init_particles() #get all parameters from the SB-categories

    #select the category
    if category=="icecosmo5":
        curr_cat = cloud_ice #if set to "
    elif category=="snowSBB":
        curr_cat = snow
    else:
        print category + " currently not implemented in __postprocess_SB.calc_distribution_from_moments"
    ###
    #calculate the normalized number concentration 
    ###
    #mass_array	=	m_ds
    #initialize and calculate bin width (del_diam)
    del_mass = np.zeros(m_ds.shape)
    del_mass[0:(m_ds.shape[0]-1)]=	m_ds[1:]-m_ds[0:(m_ds.shape[0]-1)]

    #copy the mass density and the number concentration to PAMTRA conventions
    q_h  =  twomom[curr_cat.mixrat_var][i_time,i_height]
    n_tot =  twomom[curr_cat.numcon_var][i_time,i_height]
    
    #calucate the A and lambda parameter
    x_mean = q_h/n_tot
    lam = (gamma((curr_cat.nu_SB+1)/curr_cat.mu_SB)/gamma((curr_cat.nu_SB+2)/curr_cat.mu_SB)*x_mean)**(-curr_cat.mu_SB) #based on Seifert&Beheng 2006
    A = curr_cat.mu_SB*n_tot/gamma((curr_cat.nu_SB+1)/curr_cat.mu_SB)*lam**((curr_cat.nu_SB+1)/curr_cat.mu_SB) #based on Seifert&Beheng 2006
    f_m = A*m_ds**curr_cat.nu_SB*np.exp(-lam*m_ds**curr_cat.mu_SB)
    #apply mass limits
    #f_m[np.where( m_ds < curr_cat.xmin)] = 0.0 #N_D[diam<curr_cat.Dmin or
    #f_m[np.where( m_ds > curr_cat.xmax)] = 0.0 #N_D[diam<curr_cat.Dmin or
    return f_m