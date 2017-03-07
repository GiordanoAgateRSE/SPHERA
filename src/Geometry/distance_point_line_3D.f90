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
! Program unit: distance_point_line_3D  
! Description: Computation of the distance between a point and a line in 3D.      
!-------------------------------------------------------------------------------
subroutine distance_point_line_3D(P0,P1_line,P2_line,dis)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(IN) :: P0(3),P1_line(3),P2_line(3) 
double precision,intent(INOUT) :: dis
integer(4) :: test_point
double precision :: a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3
double precision :: P3_plane1(3),P3_plane2(3),intersection_pl(3),aux_vec(3)
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
! Fictitious points to define the fictitious planes whose intersection is the 
! given line
P3_plane1(1) = P1_line(1) - 999.d0
P3_plane1(2) = P1_line(2)
P3_plane1(3) = P1_line(3)
P3_plane2(1) = P1_line(1) 
P3_plane2(2) = P1_line(2) - 999.d0
if (P2_line(3)==P2_line(3)) then
   P3_plane2(3) = P1_line(3) - 999.d0
   else
   P3_plane2(3) = P1_line(3)
endif 
! Finding the coefficients for the Cartesian equation of the planes
! Plane 1: first auxiliary plane passing for the line 
! (containing the points P1_line, P2_line and P3_plane1)
a1 = (P2_line(2) - P1_line(2)) * (P3_plane1(3) - P1_line(3)) -                 &
     (P3_plane1(2) - P1_line(2)) * (P2_line(3) - P1_line(3))
b1 = -((P2_line(1) - P1_line(1)) * (P3_plane1(3) - P1_line(3)) -               &
     (P3_plane1(1) - P1_line(1)) * (P2_line(3) - P1_line(3)))
c1 = (P2_line(1) - P1_line(1)) * (P3_plane1(2) - P1_line(2)) -                 &
     (P3_plane1(1) - P1_line(1)) * (P2_line(2) - P1_line(2))
d1 = - (a1 * P1_line(1) + b1 * P1_line(2) + c1 * P1_line(3))
! Plane 2: second auxiliary plane passing for the line (containing the points 
! P1_line, P2_line and P3_plane2)
a2 = (P2_line(2) - P1_line(2)) * (P3_plane2(3) - P1_line(3)) - (P3_plane2(2) - &
     P1_line(2)) * (P2_line(3) - P1_line(3))
b2 = - ((P2_line(1) - P1_line(1)) * (P3_plane2(3) - P1_line(3)) -              &
     (P3_plane2(1) - P1_line(1)) * (P2_line(3) - P1_line(3)))
c2 = (P2_line(1) - P1_line(1)) * (P3_plane2(2) - P1_line(2)) - (P3_plane2(1) - &
     P1_line(1)) * (P2_line(2) - P1_line(2))
d2 = - (a2 * P1_line(1) + b2 * P1_line(2) + c2 * P1_line(3))
! Plane 3: plane passing for P0 and perpendicular to P1-P2 
a3 = P2_line(1) - P1_line(1) 
b3 = P2_line(2) - P1_line(2) 
c3 = P2_line(3) - P1_line(3) 
d3 = - (a3 * P0(1) + b3 * P0(2) + c3 * P0(3)) 
! Intersection line-plane (intersection_pl): intersection of 3 planes
call three_plane_intersection(a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,test_point,  &
                              intersection_pl)
! distance point-intersection
aux_vec(:) = intersection_pl(:)-P0(:)
dis = sqrt(dot_product(aux_vec,aux_vec))
!------------------------
! Deallocations
!------------------------
return
end subroutine distance_point_line_3D

