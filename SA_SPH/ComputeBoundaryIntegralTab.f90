!cfile ComputeBoundaryIntegralTab.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ComputeBoundaryIntegralTab
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!
!************************************************************************************
! Module purpose : Module to compute local coordinates x, y, z of a grid of points,
!                  regularly distributed on the semisphere
!                  z<0 (radius = 2h), whose centre is the origin O of local axis.
!
! Calling routine: Gest_Trans
!
! Called routines: 
!
!************************************************************************************
!
subroutine ComputeBoundaryIntegralTab

!Computes local coordinates x, y, z of a grid of points, regularly distributed
!on the semisphere z<0 (radius = 2h), whose centre is the origin O of local axis.
!The semisphere will be superposed to the influence sphere of the
!generic particle near a plane boundary face, and oriented in such a way that
!the axis x, y, z coincide with the face local axes r, s, n.
!In the first three columns of the array BoundaryIntegralTab() the coordinates
!x, y, z, of each point is stored; in the forth column then relative d_alpha
!(portion of solid angle relative to the point, necessary fo integrations)
!is stored
!BITcols = 4
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Scalars ..
integer(4) :: i, j, points
double precision ::  delta_fi, delta_teta, fi, teta, xp, yp, zp, delta_alpha, &
                     delta_fimz, duesendelta_fimz
!
!.. Executable Statements ..
!
  delta_fi = FI_INTERVAL / FI_STEPS
  delta_fimz = half * delta_fi
  duesendelta_fimz = two * Sin(delta_fimz)
  delta_teta = TETA_INTERVAL / TETA_STEPS
!
  points = 0
  fi = -half * delta_fi
  Do i = 1, FI_STEPS
    fi = fi + delta_fi
    teta = -half * delta_teta
    zp = -doubleh * Cos(fi)
    Do j = 1, TETA_STEPS
      teta = teta + delta_teta
      xp = doubleh * Sin(fi) * Cos(teta)
      yp = doubleh * Sin(fi) * Sin(teta)
      delta_alpha = delta_teta * duesendelta_fimz * Sin(fi)
      points = points + 1
      BoundIntegralTab(points, 1) = xp
      BoundIntegralTab(points, 2) = yp
      BoundIntegralTab(points, 3) = zp
      BoundIntegralTab(points, 4) = delta_alpha
!.. coseni direttori del raggio OP di componenti (xp,yp,zp)
      BoundIntegralTab(points, 5) = xp / doubleh
      BoundIntegralTab(points, 6) = yp / doubleh
      BoundIntegralTab(points, 7) = zp / doubleh
    end do
  end do
!
return
end subroutine ComputeBoundaryIntegralTab
!---split

