!cfile J2Wro2.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : J2Wro2
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
! Module purpose : Module to compute the boundary volume integrals IntWdV only
!
! Calling routine: BoundaryVolumeIntegrals2D
!
! Called routines: 
!
!************************************************************************************
!
double precision function J2Wro2 (ro)
!
!.. assign modules
use GLOBAL_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
double precision,parameter :: a1 = 0.00833333333333333d0 != 1 / 120
double precision,parameter :: a2 = 0.26666666666666667d0 != 8 / 30
!
!.. Formal Arguments ..
double precision,intent(IN) :: ro
!
!.. Local Scalars ..
double precision :: ro2, ro3
!
!.. Executable Statements ..
!
  ro2 = ro * ro
  ro3 = ro2 * ro
!
  if (zero <= ro .and. ro < one) then
    J2Wro2 = KERNELCONST2D * (0.25d0 - (a1 * (40.0d0 - 36.0d0 * ro2 + 15.0d0 * ro3) * ro3))
  else if (one <= ro .and. ro < two) then
    J2Wro2 = KERNELCONST2D * (a2 - (a1 * (80.0d0 - 90.0d0 * ro + 36.0d0 * ro2 - 5.0d0 * ro3) * ro3))
  else
    J2Wro2 = zero
  end if
!
return
End function J2Wro2
!---split

