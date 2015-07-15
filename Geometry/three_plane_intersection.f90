!AA504 the whole subroutine, startig from distance_point_line_3D of AA501b
!cfile three_plane_intersection.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : three_plane_intersection
!
! Versions:
! 01   Amicarelli   08Apr14   (v5.04) New subroutine from adaptation of "distance_point_line_3D" (from v5.01, Body dynamics)
!
!************************************************************************************
! Module purpose : Computation of the intersection between 3 planes
!
! Calling routine: distance_point_line_3D,line_plane_intersection
!
! Called routines: Matrix_Inversion_3x3 
!
!************************************************************************************

subroutine three_plane_intersection(a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,test_point,intersection_pl)

! Declarations
 implicit none
 double precision,intent(in)    :: a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3
 double precision,intent(inout) :: intersection_pl(3)
 integer(4), intent(inout)         :: test_point
 integer(4)                     :: test
 double precision               :: b_vec(3)
 double precision               :: matrix(3,3),inverted(3,3)
 
!Statements 
!Solving the 3x3 linear system 
 matrix(1,1) = a1
 matrix(1,2) = b1
 matrix(1,3) = c1
 matrix(2,1) = a2
 matrix(2,2) = b2
 matrix(2,3) = c2
 matrix(3,1) = a3
 matrix(3,2) = b3
 matrix(3,3) = c3 
 call Matrix_Inversion_3x3(matrix,inverted,test)
 if (test==1) then
 test_point = 1
 b_vec(1) = -d1
 b_vec(2) = -d2
 b_vec(3) = -d3
 intersection_pl(1) = dot_product(inverted(1,:),b_vec) 
 intersection_pl(2) = dot_product(inverted(2,:),b_vec) 
 intersection_pl(3) = dot_product(inverted(3,:),b_vec) 
 else
 test_point = 0 
 endif
 
 return
end subroutine three_plane_intersection
!---split

