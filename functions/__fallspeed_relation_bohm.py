# coding: utf-8
#import packages
import numpy as np
import os 

'''
this code calculates terminal velocity according to Boehm (1992) (taken from McSnow)
e.g. vterm_bohm(np.array([1e-6]),np.array([1e-2]),np.array([1e-5]),np.array([1.0]))
                 mass              max. dimension    proj. area       as_ratio          rho=1.287,eta=1.717696e-5)

'''

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

	X_0 = 2.8e6
	grav  = 9.81 #from mo_atmo_types()

	#! 99 (7) or I (17)/(18)
	X = 8.0*mtot*grav*rho/(np.pi*(eta**2)*np.maximum(alpha,np.ones_like(alpha)*1.0)*np.maximum(q**(1.0/4.0),q)) #reduced to 8.0*mtot*grav*rho/(np.pi*(eta**2)*q**(1/4) = 8.0*mtot*grav*rho/(np.pi*(eta**2)*(area / (np.pi/4.0 * diam**2))**(1/4) for alpha=1 and q<1 (which is usually the case)

	#! 99 (9) attention sqrt(alpha) should be just alpha, see I (11)/(12)/(13)/(14)
	k = np.minimum(np.maximum(0.82+0.18*alpha,np.ones_like(alpha)*0.85),0.37+0.63/alpha,1.33/(np.maximum(np.log(alpha),np.ones_like(alpha)*0.0)+1.19)) #k is 1 for alpha=1

	#! 99 (11) or I (11)
	gama_big = np.maximum(np.ones_like(alpha)*1.0, np.minimum(np.ones_like(alpha)*1.98,3.76-8.41*alpha+9.18*alpha**2-3.53*alpha**3)) #1 for alpha=1

	#! 99 (10) or I (7)/(8)
	C_DP = np.maximum(0.292*k*gama_big,0.492-0.2/np.sqrt(alpha)) #0.292 for alpha=1
	C_DP = np.maximum(1.0,q*(1.46*q-0.46))*C_DP #0.292 for alpha=1

	#! I (23) turbulence correction
	C_DP_prim = C_DP*(1.0+1.6*(X/X_0)**2)/(1.0+(X/X_0)**2) #0.292 for small particles; larger for bigger particles 

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

