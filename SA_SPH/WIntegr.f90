!cfile WIntegr.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : WIntegr
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
! Module purpose : Module for Computing the definite integral
!
! Calling routine: ComputeVolumeIntegral_WdV2D
!
! Called routines: 
!
!************************************************************************************
!
double precision Function WIntegr (ri, h) 

!Computes the definite integral
!2h
!S W(r,h)rdr
!ri
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
double precision,parameter :: a   =  0.5d0
double precision,parameter :: a2  = -0.375d0
double precision,parameter :: a3  =  0.15d0
double precision,parameter :: aa  = -0.125d0
double precision,parameter :: aa1 =  0.05d0
double precision,parameter :: b   =  0.35d0
!
!.. Formal Arguments ..
double precision,intent(IN) :: ri
double precision,intent(IN) :: h
!
!.. Local Scalars ..
double precision :: WIntegr1, WIntegr2, ro, ro2, ro3, q1, q2, q4
!
!.. Executable Statements ..
!

!KERNELCONST2D = 0.454728408833987d0     =10/(7*pigreco) definita in AdM_USER_TYPE

  ro = ri / h

  if ( zero <= ro .AND. ro < one ) then
    ro2 = ro * ro
    ro3 = ro2 * ro
    WIntegr1 = (a + a2 * ro2 + a3 * ro3) * ro2
    WIntegr = KERNELCONST2D * (b - WIntegr1)
  else if ( one <= ro .AND. ro < two ) then
    q1 = two - ro
    q2 = q1 * q1
    q4 = q2 * q2
    WIntegr2 = (aa + aa1 * q1) * q4
    WIntegr = -KERNELCONST2D * WIntegr2
  Else
    WIntegr = zero
  end if

return
End Function WIntegr
!---split

