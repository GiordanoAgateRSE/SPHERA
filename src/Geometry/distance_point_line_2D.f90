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
! Program unit: distance_point_line_2D  
! Description: Computation of the distance between a point and a plane.     
!----------------------------------------------------------------------------------------------------------------------------------

subroutine distance_point_line_2D(P0,P1_line,P2_line,dis,normal)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(IN) :: P0(2),P1_line(2),P2_line(2)
double precision,intent(INOUT) :: dis
double precision,intent(INOUT) :: normal(2)
double precision :: a,b,c
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
a = (P2_line(2) - P1_line(2))
b = -(P2_line(1) - P1_line(1)) 
c = -(a * P1_line(1) + b * P1_line(2))
normal(1) = a / dsqrt(a ** 2 + b ** 2)
normal(2) = b / dsqrt(a ** 2 + b ** 2)
dis = (a * P0(1) + b * P0(2) + c) / dsqrt(a ** 2 + b ** 2)
!------------------------
! Deallocations
!------------------------
return
end subroutine distance_point_line_2D

