

def vterm_hw10_modified_turbcorr(mtot,diam,area,rho=1.287,eta=1.717696e-5):
    '''
    returns terminal velocity according to Heymsfield and Westbrook (2010) as implemented in McSnow but with turbulence correction from Boehm (1992)
    '''
    #define constants
    do_i = 8.0
    co_i = 0.35
    grav  = 9.81 #from mo_atmo_types()
    # modified Best number eq. on p. 2478
    Ar = area * 4.0/np.pi
    Xbest = rho* 8.0*mtot*grav*diam/(eta**2*np.pi*np.sqrt(Ar))
    #print "Xbest",Xbest
    #additional turbulence correction (from B?hm)
    X_0 = 2.8e6
    psi = (1.0+1.6*(Xbest/X_0)**2)/(1.0+(Xbest/X_0)**2)
    #psi = 1 #uncomment this line to get back to original HW10
    Xbest = Xbest / psi

    # Re-X eq. on p. 2478
    c1 = 4.0 / ( do_i**2 * np.sqrt(co_i) )
    c2 = 0.25 * do_i**2
    bracket = np.sqrt(1.0 + c1*np.sqrt(Xbest)) - 1.0
    #turbulence correction from MH05
    a0=1.0e-5
    b0=1.0
    Re = c2*bracket**2 #- a0*Xbest**b0 

    vt = eta * Re / (rho * diam)
    
    return vt
