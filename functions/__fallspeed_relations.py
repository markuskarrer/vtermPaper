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
'''
#wrapper around different fall speed models

def calc_vterm(velocity_model,mass,diam,area):
    if velocity_model=='HW10':
        vterm = vterm_hw10(mass,diam,area,rho=1.287,eta=1.717696e-5)
    elif velocity_model=='KC05':
        vterm = vterm_kc05(mass,diam,area,rho=1.287,eta=1.717696e-5) 
    elif velocity_model=='bohm':
        vterm = vterm_bohm(mass,diam,area,rho=1.287,eta=1.717696e-5)         
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
    
    # Re-X eq. on p. 2478
    c1 = 4.0 / ( do_i**2 * np.sqrt(co_i) )
    c2 = 0.25 * do_i**2
    bracket = np.sqrt(1.0 + c1*np.sqrt(Xbest)) - 1.0
    Re = c2*bracket**2

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
    
    do_i = 5.83
    co_i = 0.6
    Ct = 1.6
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
def vterm_bohm(mtot,diam,area,rho=1.287,eta=1.717696e-5):
    #!REAL(wp), PARAMETER :: X_0 = 2.8e6_wp
    #REAL(wp), PARAMETER :: lambda = 4.7_wp/1000.0_wp
    #REAL(wp) :: X_0

    q = area / (np.pi/4.0 * diam**2)


    alpha = 1.0
    X_0 = 2.8e6
    grav  = 9.81 #from mo_atmo_types()

    #! 99 (7) or I (17)/(18)
    X = 8.0*mtot*grav*rho/(np.pi*(eta**2)*max(alpha,1.0)*np.maximum(q**(1.0/4.0),q))

    #! 99 (9) attention sqrt(alpha) should be just alpha, see I (11)/(12)/(13)/(14)
    k = min(max(0.82+0.18*alpha,0.85),0.37+0.63/alpha,1.33/(max(np.log(alpha),0.0)+1.19))
    #! 99 (11) or I (11)
    gama_big = max(1.0, min(1.98,3.76-8.41*alpha+9.18*alpha**2-3.53*alpha**3))
    #! 99 (10) or I (7)/(8)
    C_DP = np.maximum(0.292*k*gama_big,0.492-0.2/np.sqrt(alpha))
    C_DP = np.maximum(1.0,q*(1.46*q-0.46))*C_DP
    #! I (23) turbulence correction
    C_DP_prim = C_DP*(1.0+1.6*(X/X_0)**2)/(1.0+(X/X_0)**2)

    #! I (21)
    beta = np.sqrt(1.0+C_DP_prim/6.0/k*np.sqrt(X/C_DP_prim))-1
    #! I (20)
    N_Re0 = 6.0*k/C_DP_prim*beta**2


    #! low Reynolds number
    #! I (16)
    C_DO = 4.5*k**2*max(alpha,1.0)
    #! I (26)
    gama_small = (C_DO - C_DP)/4.0/C_DP
    #! I (27)
    N_Re  = N_Re0*(1.0 + (2.0*beta*np.exp(-beta*gama_small))/((2.0+beta)*(1.0+beta)) )

    #! I (22)
    vterm_bohm = N_Re*eta/diam/rho
    return vterm_bohm
