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
! Program unit: line_plane_intersection  
! Description: Computation of the intersection point, if unique, between a line and a plane. 
!----------------------------------------------------------------------------------------------------------------------------------

subroutine line_plane_intersection(P1_line,P2_line,P1_plane3,P2_plane3,        &
                                   P3_plane3,test_point,intersection_pl)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(IN) :: P1_line(3),P2_line(3),P1_plane3(3),P2_plane3(3)
double precision,intent(IN) :: P3_plane3(3) 
integer(4),intent(inout) :: test_point
double precision,intent(INOUT) :: intersection_pl(3)
double precision :: a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3
double precision :: P3_plane1(3),P3_plane2(3),b_vec(3)
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
! Fictitious points to define the fictitious planes, whose intersection is the 
! given line.
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
! Plane 1: first auxiliary plane passing for the line (containing the points 
! P1_line, P2_line and P3_plane1)
a1 = (P2_line(2) - P1_line(2)) * (P3_plane1(3) - P1_line(3)) -                 &
     (P3_plane1(2) - P1_line(2)) * (P2_line(3) - P1_line(3))
b1 = - ((P2_line(1) - P1_line(1)) * (P3_plane1(3) - P1_line(3)) -              &
     (P3_plane1(1) - P1_line(1)) * (P2_line(3) - P1_line(3)))
c1 = (P2_line(1) - P1_line(1)) * (P3_plane1(2) - P1_line(2)) - (P3_plane1(1) - &
     P1_line(1)) * (P2_line(2) - P1_line(2))
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
! Plane 3: plane passing for P1_plane3, P2_plane3 and P3_plane3
a3 = (P2_plane3(2) - P1_plane3(2)) * (P3_plane3(3) - P1_plane3(3)) -           &
     (P3_plane3(2) - P1_plane3(2)) * (P2_plane3(3) - P1_plane3(3))
b3 = - ((P2_plane3(1) - P1_plane3(1)) * (P3_plane3(3) - P1_plane3(3)) -        &
     (P3_plane3(1) - P1_plane3(1)) * (P2_plane3(3) - P1_plane3(3)))
c3 = (P2_plane3(1) - P1_plane3(1)) * (P3_plane3(2) - P1_plane3(2)) -           &
     (P3_plane3(1) - P1_plane3(1)) * (P2_plane3(2) - P1_plane3(2))
d3 = - (a3 * P1_plane3(1) + b3 * P1_plane3(2) + c3 * P1_plane3(3))
! Intersection line-plane (intersection_pl): intersection of 3 planes
call three_plane_intersection(a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,test_point,  &
                              intersection_pl)
!------------------------
! Deallocations
!------------------------
return
end subroutine line_plane_intersection

