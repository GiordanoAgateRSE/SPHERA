!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
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
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: dis_point_plane  
! Description: Computation of the distance between a point and a plane.    
!-------------------------------------------------------------------------------
subroutine dis_point_plane(P0,P1_plane,P2_plane,P3_plane,dis,normal)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(IN) :: P0(3),P1_plane(3),P2_plane(3),P3_plane(3) 
double precision,intent(INOUT) :: dis
double precision,intent(INOUT) :: normal(3)
double precision :: a,b,c,d
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
a = (P2_plane(2) - P1_plane(2)) * (P3_plane(3) - P1_plane(3)) -                &
    (P3_plane(2) - P1_plane(2)) * (P2_plane(3) - P1_plane(3))
b = - ((P2_plane(1) - P1_plane(1)) * (P3_plane(3) - P1_plane(3)) -             &
    (P3_plane(1) - P1_plane(1)) * (P2_plane(3) - P1_plane(3)))
c = (P2_plane(1) - P1_plane(1)) * (P3_plane(2) - P1_plane(2)) -                &
    (P3_plane(1) - P1_plane(1)) * (P2_plane(2) - P1_plane(2))
d = - (a * P1_plane(1) + b * P1_plane(2) + c * P1_plane(3))
normal(1) = a / dsqrt(a ** 2 + b ** 2 + c ** 2)
normal(2) = b / dsqrt(a ** 2 + b ** 2 + c ** 2)
normal(3) = c / dsqrt(a ** 2 + b ** 2 + c ** 2)
dis = abs((a * P0(1) + b * P0(2) + c * P0(3) + d) /                            &
      dsqrt(a ** 2 + b ** 2 + c ** 2))
!------------------------
! Deallocations
!------------------------
return
end subroutine dis_point_plane

