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