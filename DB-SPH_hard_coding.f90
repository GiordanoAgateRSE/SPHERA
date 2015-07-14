!AA503 the whole subroutine
!cfile DB-SPH_hard_coding.f90
!************************************************************************************
!                          S P H E R A 5.0.3 
!
!                Smoothed Particle Hydrodinamics Code
!
!************************************************************************************
!
! File name     : sloshing_tank_control_points
!
! Creation      : Amicarelli-Agate, 12apr13
!
!************************************************************************************
! Module purpose : hard coding to define the control points for the sloshing tank 
!
! Calling routine: calc_pelo
!
! Called routines: ParticleCellNumber
!
!************************************************************************************

subroutine sloshing_tank_control_points

! Modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE

! Declarations ..
implicit none
double precision :: Amp,t0,T_per 
integer(4) :: i,j

!.. External Routines ..
integer(4),external :: ParticleCellNumber

!Initializations
T_per = 1.875  !1.3 !1.875     
Amp = 0.032       
t0 = 1.

!Statements
do i = 1,nlines
!.. loop sui punti della linea
   do j = control_lines(i)%Icont(1),control_lines(i)%Icont(2) 
         control_points(j)%coord(1) = 0.05 -Amp + Amp * cos(2.*PIGRECO*(tempo-t0)/T_per)
      control_points(j)%cella = ParticleCellNumber(control_points(j)%coord(:))
   end do
end do

return
end subroutine sloshing_tank_control_points
!---split
