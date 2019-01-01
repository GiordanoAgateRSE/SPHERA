!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2019 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.9.0.0
! SPHERA v.9.0.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with SPHERA. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: area_pentagon
! Description: Computation of the area of a generic pentagon from the
!              coordinates of its vertices.
!-------------------------------------------------------------------------------
subroutine area_pentagon(P1,P2,P3,P4,P5,area)
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(IN) :: P1(3),P2(3),P3(3),P4(3),P5(3)
double precision,intent(INOUT) :: area
double precision :: aux_scalar
double precision :: aux_normal(3)
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine area_triangle(P1,P2,P3,area,normal)
      implicit none
      double precision,intent(IN) :: P1(3),P2(3),P3(3)
      double precision,intent(OUT) :: area
      double precision,intent(OUT) :: normal(3)
   end subroutine area_triangle
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
area = 0.d0
!------------------------
! Statements
!------------------------
! Area of the triangle P1,P2,P3
call area_triangle(P1,P2,P3,aux_scalar,aux_normal)
! Update of the area of the pentagon
area = area + aux_scalar
! Area of the triangle P1,P3,P4
call area_triangle(P1,P3,P4,aux_scalar,aux_normal)
! Update of the area of the pentagon
area = area + aux_scalar
! Area of the triangle P1,P4,P5
call area_triangle(P1,P4,P5,aux_scalar,aux_normal)
! Update of the area of the pentagon
area = area + aux_scalar
!------------------------
! Deallocations
!------------------------
return
end subroutine area_pentagon

