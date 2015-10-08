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
! Program unit: three_plane_intersection   
! Description: Computation of the intersection of 3 planes.         
!----------------------------------------------------------------------------------------------------------------------------------

subroutine three_plane_intersection(a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,       &
                                    test_point,intersection_pl)
!------------------------
! Modules
!------------------------ 
implicit none
double precision,intent(in) :: a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3
integer(4),intent(inout) :: test_point
double precision,intent(inout) :: intersection_pl(3)
integer(4) :: test
double precision :: b_vec(3)
double precision :: matrix(3,3),inverted(3,3)
!------------------------
! Declarations
!------------------------
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
   b_vec(1) = - d1
   b_vec(2) = - d2
   b_vec(3) = - d3
   intersection_pl(1) = dot_product(inverted(1,:),b_vec) 
   intersection_pl(2) = dot_product(inverted(2,:),b_vec) 
   intersection_pl(3) = dot_product(inverted(3,:),b_vec) 
   else
      test_point = 0 
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine three_plane_intersection

