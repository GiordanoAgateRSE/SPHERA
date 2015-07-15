!AA601 the whole subroutine
!cfile Granular_flows.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : area_triangle
!
! Creation        : Amicarelli A., 31July14
!
!************************************************************************************
! Module purpose : Computation of the area of a generic triangle, 
!                  provided the coordinates of its vertices
!
! Calling routine: area_quadrilateral,Import_ply_surface_meshes
!
! Called routines: Vector_Product
!
!************************************************************************************

subroutine area_triangle(P1,P2,P3,area,normal)

! Declarations
 implicit none
 double precision,intent(IN)    :: P1(3),P2(3),P3(3)
 double precision,intent(OUT)   :: area
 double precision,intent(OUT)   :: normal(3)
 double precision               :: vec_1(3),vec_2(3)

 !Statements 
 vec_1(:) = P2(:)-P1(:)
 vec_2(:) = P3(:)-P1(:) 
 call Vector_Product(vec_1,vec_2,normal,3)
 area = 0.5d0 * dsqrt(dot_product(normal,normal))
 normal(:) = normal(:)/(2.d0*area)

 return
end subroutine area_triangle
!---split

