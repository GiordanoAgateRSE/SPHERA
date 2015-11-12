!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-Ricerca di Sistema)



! SPHERA authors and email contact are provided on SPHERA documentation.

! This file is part of SPHERA v.8.0.
! SPHERA v.8.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with SPHERA. If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: area_triangle  
! Description: Computation of the area of a generic triangle, provided the coordinates of its vertices.   
!----------------------------------------------------------------------------------------------------------------------------------

subroutine area_triangle(P1,P2,P3,area,normal)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(IN) :: P1(3),P2(3),P3(3)
double precision,intent(OUT) :: area
double precision,intent(OUT) :: normal(3)
double precision :: vec_1(3),vec_2(3)
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
vec_1(:) = P2(:)-P1(:)
vec_2(:) = P3(:)-P1(:) 
!------------------------
! Statements
!------------------------
call Vector_Product(vec_1,vec_2,normal,3)
area = 0.5d0 * dsqrt(dot_product(normal,normal))
normal(:) = normal(:) / (2.d0 * area)
!------------------------
! Deallocations
!------------------------
return
end subroutine area_triangle

