!AA601 the whole subroutine
!cfile Granular_flows.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : area_quadrilateral
!
! Versions: 
! 01 Amicarelli 26Jan15 (creation, AA601, DBSPH-input)
!
!************************************************************************************
! Module purpose : Computation of the area of a generic quadrilateral from the coordinates of its vertices
!
! Calling routine: ReadSectionFlowRate
!
! Called routines: area_triangle
!
!************************************************************************************

subroutine area_quadrilateral(P1,P2,P3,P4,area)

! Declarations
 implicit none
 double precision,intent(IN)    :: P1(3),P2(3),P3(3),P4(3)
 double precision,intent(INOUT) :: area
 double precision               :: area_triangle_1,area_triangle_2
!AA601
 double precision               :: aux_normal(3)
  
!Statements 
!Area of the triangle P1,P2,P3: 0.5*vector_product(vec_a1,vec_b), vec_a1=(P2-P1),vec_b=(P3-P1)
 call area_triangle(P1,P2,P3,area_triangle_1,aux_normal)
!Area of the trinagle P1,P4,P3: 0.5*vector_product(vec_a2,vec_b), vec_a2=(P4-P1),vec_b=(P3-P1)
 call area_triangle(P1,P4,P3,area_triangle_2,aux_normal)
!Area of the quadrilateral: sum of the areas of the trinagles
 area = area_triangle_1 + area_triangle_2
 return
end subroutine area_quadrilateral
!---split

