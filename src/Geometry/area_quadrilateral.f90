!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-) 
!      
!     
!   
!      
!  

! This file is part of SPHERA.
!  
!  
!  
!  
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
!  
!  
!  
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: area_quadrilateral 
! Description: Computation of the area of a generic quadrilateral from the coordinates of its vertices.  
!----------------------------------------------------------------------------------------------------------------------------------

subroutine area_quadrilateral(P1,P2,P3,P4,area)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(IN) :: P1(3),P2(3),P3(3),P4(3)
double precision,intent(INOUT) :: area
double precision :: area_triangle_1,area_triangle_2
double precision :: aux_normal(3)
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
! Area of the triangle P1,P2,P3: 0.5*vector_product(vec_a1,vec_b), 
! vec_a1=(P2-P1),vec_b=(P3-P1)
call area_triangle(P1,P2,P3,area_triangle_1,aux_normal)
! Area of the trinagle P1,P4,P3: 0.5*vector_product(vec_a2,vec_b),
! vec_a2=(P4-P1),vec_b=(P3-P1)
call area_triangle(P1,P4,P3,area_triangle_2,aux_normal)
! Area of the quadrilateral: sum of the areas of the trinagles
area = area_triangle_1 + area_triangle_2
!------------------------
! Deallocations
!------------------------
return
end subroutine area_quadrilateral

