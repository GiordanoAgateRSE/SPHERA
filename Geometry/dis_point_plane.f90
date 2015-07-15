!AA501b the whole subroutine
!AA504 sub
!cfile dis_point_plane.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : dis_point_plane
!
! Creation : Amicarelli-Agate, 25Jul12
!
!************************************************************************************
! Module purpose : Computation of the distance between a point and a plane
!AA601 sub
! Calling routine: RHS_body_dynamics,DBSPH_inlet_outlet
!
! Called routines: /
!
!************************************************************************************

subroutine dis_point_plane(P0,P1_plane,P2_plane,P3_plane,dis,normal)

! Declarations
 implicit none
 double precision,intent(IN) :: P0(3),P1_plane(3),P2_plane(3),P3_plane(3) 
 double precision,intent(INOUT) :: dis
 double precision,intent(INOUT) :: normal(3)
 double precision               :: a,b,c,d
 
! Statements
 a = (P2_plane(2)-P1_plane(2)) * (P3_plane(3)-P1_plane(3)) - (P3_plane(2)-P1_plane(2)) * (P2_plane(3)-P1_plane(3))
 b = -((P2_plane(1)-P1_plane(1)) * (P3_plane(3)-P1_plane(3)) - (P3_plane(1)-P1_plane(1)) * (P2_plane(3)-P1_plane(3)))
 c = (P2_plane(1)-P1_plane(1)) * (P3_plane(2)-P1_plane(2)) - (P3_plane(1)-P1_plane(1)) * (P2_plane(2)-P1_plane(2))
 d = -(a*P1_plane(1)+b*P1_plane(2)+c*P1_plane(3))
 normal(1) = a/dsqrt(a**2+b**2+c**2)
 normal(2) = b/dsqrt(a**2+b**2+c**2)
 normal(3) = c/dsqrt(a**2+b**2+c**2)
 dis = abs((a*P0(1)+b*P0(2)+c*P0(3)+d)/dsqrt(a**2+b**2+c**2))
 
return
end subroutine dis_point_plane
!---split

