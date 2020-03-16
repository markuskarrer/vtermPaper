import numpy as np
import sys
'''
stores the coefficients of the m-D & v-D parameterizations from the vterm study in a a class
'''
from IPython.core.debugger import Tracer ; debug=Tracer()

class particle(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def init_class(monotype):
    if monotype=="plate":
        return particle(a_ms = 0.075588,
                        b_ms = 2.218742,
                        A_vAtlas = 1.366326,
                        B_vAtlas = 1.390639,
                        C_vAtlas = 1285.591040,
                        a_vpower = 12.114902,
                        b_vpower = 0.5090708
                        )
    elif monotype=="dendrite":
        return particle(a_ms = 0.026505, 
                        b_ms = 2.217536,
                        A_vAtlas = 1.034505, #1.034505&1.067134&1388.834943
                        B_vAtlas = 1.067134,
                        C_vAtlas = 1388.834943,
                        a_vpower = 24.348411, #24.348411&0.69751846
                        b_vpower = 0.69751846
                        )
    elif monotype=="mixcolumndend":
        return particle(a_ms = 0.016676, #0.016676&1.940710
                        b_ms = 1.940710,
                        A_vAtlas = 1.108004, #1.108004&1.111288&2358.022852
                        B_vAtlas = 1.111288,
                        C_vAtlas = 2358.022852,
                        a_vpower = 8.567136, #8.567136&0.3932750694
                        b_vpower = 0.3932750694
                        )
    elif monotype=="column":
        return particle(a_ms = 0.045801, #0.045801&2.073687 
                        b_ms = 2.073687,
                        A_vAtlas = 1.614318, #1.614318&1.654455&1612.423152
                        B_vAtlas = 1.654455,
                        C_vAtlas = 1612.423152,
                        a_vpower = 22.800083, #22.800083&0.52146112999152605
                        b_vpower = 0.5214611
                        )
    elif monotype=="powerlawOLD":
        return particle(a_geo=5.130000, #
                        b_geo=0.500000,
                        nu_SB=0.000000,
                        mu_SB=0.500000,
                        a_vpower =8.294000, 
                        b_vpower = 0.125000
                        )
    else:
        print monotype + " is not defined in vtermstudy_coeffs.init_class"
        sys.exit()

def add_size_dist_param_ND(particle_class,mu,gam):
    '''
    INPUT: 
        particle_class: class containing (at least) the m-D  parameter
        mu: paramter in N(D)=N0*D**mu*exp(-lam*D**gam)
        gam: paramter in N(D)=N0*D**mu*exp(-lam*D**gam)
    '''
    
    
    particle_class.mu=mu
    particle_class.gam=gam


    return particle_class

def complement_mDDM(particle_class):
    '''
    calculate m-D or D-m coefficients if missing
    '''

    if "a_ms" in particle_class.__dict__.keys() and not "a_geo" in particle_class.__dict__.keys(): #then we need to fill the D-m parameters
        
        #calculate D(m) coefficients
        particle_class.a_geo = (1./particle_class.a_ms)**(1./particle_class.b_ms)
        particle_class.b_geo = 1./particle_class.b_ms
    elif "a_geo" in particle_class.__dict__.keys() and not "a_ms" in particle_class.__dict__.keys(): #then we need to fill the D-m parameters
        
        #calculate m(D) coefficients
        particle_class.a_ms = (1./particle_class.a_geo)**(1./particle_class.b_geo)
        particle_class.b_ms = 1./particle_class.b_geo

    return particle_class
def complement_distparam(particle_class):
    '''
    calculate missing size distribution parameter
    '''
    if "mu" in particle_class.__dict__.keys() and not "nu_SB" in particle_class.__dict__.keys(): #then we need to fill the N(m) parameters

        particle_class.nu_SB =  particle_class.b_geo*particle_class.mu+particle_class.b_geo-1
        particle_class.mu_SB = particle_class.b_geo*particle_class.gam
    elif "nu_SB" in particle_class.__dict__.keys() and not "mu" in particle_class.__dict__.keys(): #then we need to fill the N(m) parameters
        particle_class.mu =  particle_class.b_ms*particle_class.nu_SB+particle_class.b_ms-1
        particle_class.gam = particle_class.b_ms*particle_class.mu_SB
    
    return particle_class

def return_class(monotype="plate",mu=None,gam=None,nu_SB=None,mu_SB=None):
    '''
    main part
    '''
    
    particle_class = init_class(monotype)
    
    particle_class = add_size_dist_param_ND(particle_class,mu=mu,gam=gam)

    particle_class = complement_mDDM(particle_class)
    
    particle_class = complement_distparam(particle_class)
    
    return particle_class

if __name__ == "__main__":
    #for testing
    for monotype in ["plate","dendrite","column","mixcolumndend"]:
        particle_class = return_class(monotype=monotype,mu=0.0,gam=1.0)

        print monotype,particle_class.__dict__
        print "\n",monotype,particle_class.A_vAtlas,particle_class.B_vAtlas,particle_class.C_vAtlas,"\n " 
