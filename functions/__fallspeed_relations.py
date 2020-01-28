# coding: utf-8
#import packages
import numpy as np
import os 
#import other self-defined functions
import __postprocess_McSnow
import __postprocess_SB
#from IPython.core.debugger import Tracer ; Tracer()()

'''
this code calculates different fallspeed models (f.e. taken from McSnow)
e.g. calc_vterm("bohm",mass,diam,area) 
mass: numpy-array containing masses
diam: numpy-array containing maximum dimension
area: numpy-array containing projected area 
ALL IN SI-UNITS
'''
#wrapper around different fall speed models

def calc_vterm(velocity_model,mass,diam,area,as_ratio=1.0,turb_corr="all"):
    #turb_corr is only for mitch_heym
    if velocity_model=='HW10':
        vterm = vterm_hw10(mass,diam,area,rho=1.287,eta=1.717696e-5)
    elif velocity_model=='KC05':
        vterm = vterm_kc05(mass,diam,area,rho=1.287,eta=1.717696e-5) 
    elif velocity_model=='bohm89':
        vterm = vterm_bohm89(mass,diam,area,rho=1.287,eta=1.717696e-5)
    elif velocity_model=='bohm':
        vterm = vterm_bohm(mass,diam,area,np.ones_like(mass),rho=1.287,eta=1.717696e-5)
    elif velocity_model=='bohmar':
        vterm = vterm_bohm(mass,diam,area,as_ratio,rho=1.287,eta=1.717696e-5)
    elif velocity_model=='mitch_heym': #the turbulence correction produces a jump at 0.5mm
        vterm = vterm_mitch_heym(mass,diam,area,rho=1.287,eta=1.717696e-5,turb_corr=turb_corr)
    else:
        print velocity_model + " not implemented (in functions/__fallspeed_relations.py"
        
    return vterm

'''#copied here from McSnow mo_mass2diam.f90
!
! terminal fall velocity 
!   from Heymsfield and Westbrook (2010)
!
FUNCTION vterm_hw10( atmo, sp ) RESULT( vt )
  TYPE(t_atmo), INTENT(in)    :: atmo
  TYPE(t_sp),   INTENT(inout) :: sp    ! super-droplet
  REAL(mp)                    :: vt

  REAL(wp), PARAMETER :: do_i = 8.0_wp
  REAL(wp), PARAMETER :: co_i = 0.35_wp

  REAL(wp) :: mtot, Re, c1, c2, Xbest, bracket, Ar

  mtot = sp%m_r + sp%m_i + sp%m_w

  ! modified Best number eq. on p. 2478
  Ar = sp%p_area * z4_pi
  Xbest = atmo%rho * 8._wp * mtot * grav * sp%d / (atmo%eta**2 * pi * SQRT(Ar))

  ! Re-X eq. on p. 2478
  c1 = 4.0_wp / ( do_i**2 * SQRT(co_i) )
  c2 = 0.25_wp * do_i**2
  bracket = SQRT(1.0 + c1*SQRT(Xbest)) - 1.0_wp
  Re = c2*bracket**2

  vt = atmo%eta * Re / (atmo%rho * sp%d)
END FUNCTION vterm_hw10
'''

def vterm_hw10(mtot,diam,area,rho=1.287,eta=1.717696e-5):
    '''
    returns terminal velocity according to Heymsfield and Westbrook (2010) as implemented in McSnow
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
    #psi = (1.0+1.6*(Xbest/X_0)**2)/(1.0+(Xbest/X_0)**2)
    psi = 1 #uncomment this line to get back to original HW10
    Xbest = Xbest / psi

    # Re-X eq. on p. 2478
    c1 = 4.0 / ( do_i**2 * np.sqrt(co_i) )
    #c1 = 0.292**0.5/6 #as in Boehm
    #print "c1",c1
    c2 = 0.25 * do_i**2
    bracket = np.sqrt(1.0 + c1*np.sqrt(Xbest)) - 1.0
    #turbulence correction from MH05
    a0=1.0e-5
    b0=1.0
    Re = c2*bracket**2 #- a0*Xbest**b0 
    #print "Re1",c2*bracket**2,"Re2",a0*Xbest**b0,"rel (Re-Re1)/Re1",(a0*Xbest**b0)/(c2*bracket**2)

    vt = eta * Re / (rho * diam)
    
    return vt
#copied here from McSnow mo_mass2diam.f90
'''
FUNCTION X2Cd_kc05rough( Xbest ) RESULT( Cd )
  REAL(wp), INTENT(in)  :: Xbest
  REAL(wp)              :: Cd

  REAL(wp), PARAMETER :: do_i = 5.83_wp
  REAL(wp), PARAMETER :: co_i = 0.6_wp
  REAL(wp), PARAMETER :: Ct = 1.6_wp
  REAL(wp), PARAMETER :: X0_i = .35714285714285714285e-6 !1.0_wp/2.8e6_wp
  ! derived constants
  REAL(wp), PARAMETER :: c1 = 4.0 / ( do_i**2 * SQRT(co_i) )
  REAL(wp), PARAMETER :: c2 = 0.25_wp * do_i**2

  REAL(wp) :: bracket, psi, Re

  ! Re-X eq. (2.5)
  bracket = SQRT(1.0 + c1*SQRT(Xbest)) - 1.0
  ! turbulent Reynold's number, eq (3.3)
  psi = (1+(Xbest*X0_i)**2) / (1+Ct*(Xbest*X0_i)**2)
  Re  = c2*bracket**2 ! * SQRT(psi) ! TODO remove psi in Re?
  ! eq. (2.1) from KC05 with (3.2)
  Cd = co_i * (1._wp + do_i/SQRT(Re))**2 / psi
END FUNCTION X2Cd_kc05rough
'''
def X2Cd_kc05rough( Xbest ):
    Ct = 1.6
    #original  
    do_i = 5.83
    co_i = 0.6
    '''# test from HW #with that vterm gets even faster
    do_i = 8.0
    co_i = 0.35
    #'''    
    X0_i = .35714285714285714285e-6
    c1 = 4.0 / ( do_i**2 * np.sqrt(co_i) )
    c2 = 0.25 * do_i**2
    
    #! Re-X eq. (2.5)
    bracket = np.sqrt(1.0 + c1*np.sqrt(Xbest)) - 1.0
    #! turbulent Reynold's number, eq (3.3)
    psi = (1+(Xbest*X0_i)**2) / (1+Ct*(Xbest*X0_i)**2)
    Re  = c2*bracket**2 #! * np.sqrt(psi) #! TODO remove psi in Re?
    #! eq. (2.1) from KC05 with (3.2)
    Cd = co_i * (1. + do_i/np.sqrt(Re))**2 / psi
    
    return Cd
#copied here from McSnow mo_mass2diam.f90
'''
!
! terminal fall velocity of ice from Khvorostyanov and Curry (2005/2002)
! optional parameter to switch to smooth surface parameters
! 
! Note: Rassmussen and Heymsfield JAS 44, 2754pp, 1987 interpolate to Cd 0.6 at 
! Re > 2.5e4 on wet particles. This is not done here to avoid a transition to pure water drops
!
FUNCTION vterm_kc05( atmo, sp, smooth ) RESULT( vt )
  TYPE(t_atmo), INTENT(in)    :: atmo  ! background atmosphere
  TYPE(t_sp),   INTENT(inout) :: sp    ! super-droplet
  LOGICAL,      INTENT(in), OPTIONAL :: smooth
  REAL(mp)                    :: vt

  REAL(wp) :: mtot, Xbest, Cd, Vb, Fb

  mtot = sp%m_r + sp%m_i + sp%m_w

  ! Best number eq. (2.4b) with buoyancy
  Vb = sp%m_i*rhoii + MAX(sp%v_r, sp%m_r*rhoii + sp%m_w*rholi)
  Fb = atmo%rho * Vb * grav
  Xbest = 2._wp * ABS(mtot*grav-Fb) * atmo%rho * sp%d**2 / (sp%p_area * atmo%eta**2)

  IF( PRESENT(smooth) ) THEN
    IF( smooth ) THEN
      Cd  = X2Cd_kc05smooth( Xbest )
    ELSE
      Cd  = X2Cd_kc05rough( Xbest )
    ENDIF
  ELSE
      Cd  = X2Cd_kc05rough( Xbest )
  ENDIF
  vt = SQRT( 2*ABS(mtot*grav - Fb)/(atmo%rho * sp%p_area * Cd) )
END FUNCTION vterm_kc05
'''
def vterm_kc05(mtot,diam,area,rho=1.287,eta=1.717696e-5):
    '''
    returns terminal velocity according to Khvorostyanov and Curry (2005/2002) as implemented in McSnow
    '''
    rhoi = 919. #from mo_atmo_types()
    rhoii = 1./rhoi #from mo_atmo_types()
    grav  = 9.81 #from mo_atmo_types()

    #! Best number eq. (2.4b) with buoyancy
    Vb = mtot*rhoii #TODO: this function can not handle riming yet+ MAX(sp%v_r, sp%m_r*rhoii + sp%m_w*rholi)
    Fb = rho * Vb * grav
    Xbest = 2. * abs(mtot*grav-Fb) * rho * diam**2 / (area * eta**2)
    #Xbest = rho* 8.0*mtot*grav/(eta**2*np.pi*(area * 4.0/np.pi / diam**2)**(1./4)) #HW10 formulation, but exponent as in Bohm
    #Xbest = rho* 8.0*mtot*grav/(eta**2*np.pi*(area * 4.0/np.pi / diam**2)**(1./2)) #HW10


    Cd  = X2Cd_kc05rough( Xbest )

    vt = np.sqrt( 2*abs(mtot*grav - Fb)/(rho * area * Cd) )
    return vt

#copied here from McSnow mo_mass2diam.f90
'''
!
! test, compare boehm to KC fall velocities
!
REAL(wp) FUNCTION vterm_bohm( atmo, sp )
  TYPE(t_atmo), INTENT(in) :: atmo
  TYPE(t_sp),   INTENT(in) :: sp

  REAL(wp) :: alpha, X, mtot, q, C_DP, k, gama_big, C_DP_prim, beta, gama_small, C_DO, N_Re0, N_Re
  !REAL(wp), PARAMETER :: X_0 = 2.8e6_wp
  REAL(wp), PARAMETER :: lambda = 4.7_wp/1000.0_wp
  REAL(wp) :: X_0

  mtot = sp%m_w + sp%m_r + sp%m_i
  q = sp%p_area / (pi/4.0_wp * sp%d**2)
  alpha = 1.0_wp
  IF( sp%m_r + sp%m_i < TINY(sp%m_i) ) THEN ! pure water
    alpha = exp(-sp%d/lambda) + (1.0_wp - exp(-sp%d/lambda))/(1.0_wp+sp%d/lambda)
    !alpha = MIN(1.0_wp,1.05_wp-0.655_wp*sp%d*100.0_wp)
    X_0 = 6.7e6_wp
  ELSE ! rough ice
    alpha = 1.0_wp
    X_0 = 2.8e6_wp
  END IF

  ! 99 (7) or I (17)/(18)
  X = 8.0_wp*mtot*grav*atmo%rho/(pi*(atmo%eta**2)*MAX(alpha,1.0_wp)*MAX(q**(1.0_wp/4.0_wp),q))

  ! 99 (9) attention sqrt(alpha) should be just alpha, see I (11)/(12)/(13)/(14)
  k = min(max(0.82_wp+0.18_wp*alpha,0.85_wp),0.37_wp+0.63_wp/alpha,1.33_wp/(max(log(alpha),0.0_wp)+1.19_wp))
  ! 99 (11) or I (11)
  gama_big = max(1.0_wp, min(1.98_wp,3.76_wp-8.41_wp*alpha+9.18*alpha**2-3.53*alpha**3))
  ! 99 (10) or I (7)/(8)
  C_DP = max(0.292_wp*k*gama_big,0.492_wp-0.2_wp/sqrt(alpha))
  C_DP = max(1.0_wp,q*(1.46_wp*q-0.46_wp))*C_DP
alias umount_MC='fusermount -u /home/mkarrer/Dokumente/McSnow'
  ! I (23) turbulence correction
  C_DP_prim = C_DP*(1.0_wp+1.6_wp*(X/X_0)**2)/(1.0_wp+(X/X_0)**2)

  ! I (21)
  beta = SQRT(1.0_wp+C_DP_prim/6.0_wp/k*SQRT(X/C_DP_prim))-1
  ! I (20)
  N_Re0 = 6.0_wp*k/C_DP_prim*beta**2


  ! low Reynolds number
  ! I (16)
  C_DO = 4.5_wp*k**2*MAX(alpha,1.0_wp)
  ! I (26)
  gama_small = (C_DO - C_DP)/4.0_wp/C_DP
  ! I (27)
  N_Re  = N_Re0*(1.0_wp + (2.0_wp*beta*EXP(-beta*gama_small))/((2.0_wp+beta)*(1.0_wp+beta)) )

  ! I (22)
  vterm_bohm = N_Re*atmo%eta/sp%d/atmo%rho
END FUNCTION vterm_bohm
'''
def vterm_bohm(mtot,diam,area,as_ratio,rho=1.287,eta=1.717696e-5):
    '''
    returns the fall speed with assumptions from B?hm (1992)
    '''
    #!REAL(wp), PARAMETER :: X_0 = 2.8e6_wp
    #REAL(wp), PARAMETER :: lambda = 4.7_wp/1000.0_wp
    #REAL(wp) :: X_0

    q = area / (np.pi/4.0 * diam**2)

    
    alpha = np.array(as_ratio) #1.0
    #print a
    X_0 = 2.8e6
    grav  = 9.81 #from mo_atmo_types()

    #! 99 (7) or I (17)/(18)
    X = 8.0*mtot*grav*rho/(np.pi*(eta**2)*np.maximum(alpha,np.ones_like(alpha)*1.0)*np.maximum(q**(1.0/4.0),q)) #reduced to 8.0*mtot*grav*rho/(np.pi*(eta**2)*q**(1/4) = 8.0*mtot*grav*rho/(np.pi*(eta**2)*(area / (np.pi/4.0 * diam**2))**(1/4) for alpha=1 and q<1 (which is usually the case)
    #print "q",q,"np.maximum(q**(1.0/4.0),q))",np.maximum(q**(1.0/4.0),q)
    #! 99 (9) attention sqrt(alpha) should be just alpha, see I (11)/(12)/(13)/(14)
    k = np.minimum(np.maximum(0.82+0.18*alpha,np.ones_like(alpha)*0.85),0.37+0.63/alpha,1.33/(np.maximum(np.log(alpha),np.ones_like(alpha)*0.0)+1.19)) #k is 1 for alpha=1
    #print "k",k
    #! 99 (11) or I (11)
    gama_big = np.maximum(np.ones_like(alpha)*1.0, np.minimum(np.ones_like(alpha)*1.98,3.76-8.41*alpha+9.18*alpha**2-3.53*alpha**3)) #1 for alpha=1
    #print "gama_big",gama_big
    #! 99 (10) or I (7)/(8)
    C_DP = np.maximum(0.292*k*gama_big,0.492-0.2/np.sqrt(alpha)) #0.292 for alpha=1
    C_DP = np.maximum(1.0,q*(1.46*q-0.46))*C_DP #0.292 for alpha=1
    #print "C_DP",C_DP
    #! I (23) turbulence correction
    C_DP_prim = C_DP*(1.0+1.6*(X/X_0)**2)/(1.0+(X/X_0)**2) #0.292 for small particles; larger for bigger particles 
    #print "C_DP_prim",C_DP_prim
    #! I (21)
    beta = np.sqrt(1.0+C_DP_prim/6.0/k*np.sqrt(X/C_DP_prim))-1
    #! I (20)
    N_Re0 = 6.0*k/C_DP_prim*beta**2


    #! low Reynolds number
    #! I (16)
    C_DO = 4.5*k**2*np.maximum(alpha,np.ones_like(alpha)*1.0)
    #! I (26)
    gama_small = (C_DO - C_DP)/4.0/C_DP
    #! I (27)
    N_Re  = N_Re0*(1.0 + (2.0*beta*np.exp(-beta*gama_small))/((2.0+beta)*(1.0+beta)) )

    #! I (22)
    vterm_bohm = N_Re*eta/diam/rho
    return vterm_bohm


def vterm_bohm89(mass,diam, area, rho=1.287, eta=1.717696e-5):
    grav = 9.81
    rho_ice = 917.0
    q = area / (np.pi/4.0 * diam**2)

    X = 8.0*mass*grav*rho/(np.pi*(eta**2)*q**0.25)

    Re = 8.5*((1.0+0.1519*X**0.5)**0.5-1.0)**2
    vterm_bohm = Re*eta/diam/rho

    return vterm_bohm


'''
#copied here from the P3 look-up table creator create_p3_lookupTable_1.f90
! assume 600 hPa, 253 K for p and T for fallspeed calcs (for reference air density) #TODO: do we have to take this into account??
 g   = 9.861                           ! gravity
# p   = 60000.                          ! air pressure (pa)
# t   = 253.15                          ! temp (K)
# rho = p/(287.15*t)                    ! air density (kg m-3)
# mu  = 1.496E-6*t**1.5/(t+120.)/rho    ! viscosity of air
# dv  = 8.794E-5*t**1.81/p              ! diffusivity of water vapor in air
# dt  = 10.                             ! time step for collection (s)

! parameters for surface roughness of ice particle used for fallspeed
! see mitchell and heymsfield 2005
 del0 = 5.83
 c0   = 0.6
 c1   = 4./(del0**2*c0**0.5)
 c2   = del0**2/4.

! correction for turbulence
!            if (d1.lt.500.e-6) then
          a0 = 0.
          b0 = 0.
!            else
!               a0=1.7e-3
!               b0=0.8
!            end if

! fall speed for ice
! Best number
          xx = 2.*cs1*g*rho*d1**(ds1+2.-bas1)/(aas1*(mu*rho)**2)

! drag terms
          b1 = c1*xx**0.5/(2.*((1.+c1*xx**0.5)**0.5-1.)*(1.+c1*xx**0.5)**0.5)-a0*b0*xx** &
               b0/(c2*((1.+c1*xx**0.5)**0.5-1.)**2)

          a1 = (c2*((1.+c1*xx**0.5)**0.5-1.)**2-a0*xx**b0)/xx**b1

! velocity in terms of drag terms
          fall2(jj) = a1*mu**(1.-2.*b1)*(2.*cs1*g/(rho*aas1))**b1*d1**(b1*(ds1-bas1+2.)-1.)
'''
def vterm_mitch_heym(mtot,diam,area,rho=1.287,eta=1.717696e-5,turb_corr="all"):
    '''
    returns fall speed as assumed by Mitchel and Heymsfield (2005)
    INPUT: 
        mtot: array of masses
        diam: array of diameters
        area: array of projected areas
    (OPTIONAL)
        rho: air density
        eta: air viscosity
        turb_corr: apply turbulence correction
    '''
    #some constants (as assumed in P3)
    g   = 9.861                           # gravity
    
    # parameters for surface roughness of ice particle used for fallspeed
    #! see mitchell and heymsfield 2005
    del0 = 5.83
    c0   = 0.6
    c1   = 4./(del0**2*c0**0.5)
    c2   = del0**2/4.
    
    if turb_corr=="small": #turbulence correction only for particles >500mum
        '''
        if False: #TODO: apply correction only for larger part of arrays # (diam < 500e-6): #no turbulence correction for small particles
            a0 = 0.
            b0 = 0.
        else:
            a0=1.7e-3
            b0=0.8
        '''    
        #a0 and b0 should be arrays here
        a0 = np.zeros_like(diam); b0 = np.zeros_like(diam)
        a0[diam>500e-6]=1.7e-3
        b0[diam>500e-6]=0.8
    elif turb_corr=="all": #turbulence correction for all
        a0=1.7e-3*np.ones_like(diam)
        b0=0.8*np.ones_like(diam)
    # Best number
    X = 2.*mtot* g*rho*diam**2 / (area*(eta*rho)**2) #restructured from: 2.*cs1*g*rho*d1**(ds1+2.-bas1)/(aas1*(mu*rho)**2)# cs1*d1**ds1=mass; d1=diam; aas1*d1**bas1=area
    #drag terms
    b1 = c1*X**0.5/(2.*((1.+c1*X**0.5)**0.5-1.)*(1.+c1*X**0.5)**0.5)-a0*b0*X**b0/(c2*((1.+c1*X**0.5)**0.5-1.)**2) #eq. 7
    a1 = (c2*((1.+c1*X**0.5)**0.5-1.)**2-a0*X**b0)/X**b1 #eq. 6 
    #velocity in terms of drag terms
    vterm_mitch_heym = a1*eta**(1.-2.*b1)*(2*g/rho)**b1*(mtot/area)**b1*diam**(2*b1-1)#(2.*cs1*g/(rho*aas1))**b1*d1**(b1*(ds1-bas1+2.)-1.) #a1*eta**(1.-2.*b1)*(2.*cs1*g/(rho*aas1))**b1*d1**(b1*(ds1-bas1+2.)-1.) #eq. 11 and 12 inserted in eq 10
    
    return vterm_mitch_heym
    
