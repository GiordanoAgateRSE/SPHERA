!cfile VelLaw.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : VelLaw
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
! Module purpose : Module for velocity law calculation
!
! Calling routine: Loop_Irre_2D, Loop_Irre3D, SetParticles
!
! Called routines: 
!
!************************************************************************************
!
subroutine VelLaw (vlaw,vel,np)
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
  implicit none

!
!.. Formal Arguments ..
double precision,intent(IN), dimension(0:3,MAXPOINTSVLAW) :: vlaw
double precision,intent(OUT),dimension(3)                 :: vel
integer(4),      intent(IN)                               :: np
!
!.. Local Scalars ..
double precision :: fra
integer(4)       :: n
!
!.. Executable Statements ..
!
 if ( np <= 1 ) return
 do n = 2,np
   if ( tempo > vlaw(0,n) ) cycle
   fra = ( tempo - vlaw(0,n-1) ) / ( vlaw(0,n) - vlaw(0,n-1) )   
   vel(1:3) = vlaw(1:3,n-1) +  ( vlaw(1:3,n) - vlaw(1:3,n-1) ) * fra
   return
 end do
 vel(1:3) = vlaw(1:3,np)

return
end subroutine VelLaw
!---split

