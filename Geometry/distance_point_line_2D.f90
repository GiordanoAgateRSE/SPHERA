!AA501b the whole subroutine
!AA504 sub
!cfile distance_point_line_2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : distance_point_line_2D
!
! Creation : Amicarelli, 24Oct12
!
!************************************************************************************
! Module purpose : Computation of the distance between a point and a plane
!
! Calling routine: point_inout_polygone
!
! Called routines: /
!
!************************************************************************************

subroutine distance_point_line_2D(P0,P1_line,P2_line,dis,normal)

! Declarations
 implicit none
 double precision,intent(IN) :: P0(2),P1_line(2),P2_line(2)
 double precision,intent(INOUT) :: dis
 double precision,intent(INOUT) :: normal(2)
 double precision               :: a,b,c
 
! Statements
 a = (P2_line(2)-P1_line(2))
 b = -(P2_line(1)-P1_line(1)) 
 c = -(a*P1_line(1)+b*P1_line(2))
 normal(1) = a/dsqrt(a**2+b**2)
 normal(2) = b/dsqrt(a**2+b**2)
 dis = (a*P0(1)+b*P0(2)+c)/dsqrt(a**2+b**2)
 
return
end subroutine distance_point_line_2D
!---split

