!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2021 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: point_inout_convex_non_degenerate_polygon
! Description: Test to evaluate if a point lies inside or strictly outside a 
!              polygon. A point is internal to the polygon if its distances from
!              the lines passing for the polygon sides (no matter about the 
!              number of sides, but they must be taken in either a clockwise or 
!              an anti-clockwise order), have all the same sign of a generic 
!              polygon vertex not belonging to the selected side -a null 
!              distance is always a positive test for internal points-). The 
!              maximum number of polygon sides is now equal to 6 (triangles, 
!              quadrilaterals, pentagons and hexagons can be treated). Polygons
!              must be convex and non-degenerate (a n-side polygon should have n
!              vertices, not more).
!-------------------------------------------------------------------------------
subroutine point_inout_convex_non_degenerate_polygon(point,n_sides,point_pol_1,&
                                                     point_pol_2,point_pol_3,  &
                                                     point_pol_4,point_pol_5,  &
                                                     point_pol_6,test)
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: n_sides
double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
double precision,intent(in) :: point_pol_3(2),point_pol_4(2),point_pol_5(2)
double precision,intent(in) :: point_pol_6(2)
integer(4),intent(inout) :: test
double precision :: dis1,dis2
double precision :: normal(2)
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
dis2 = -9999999.d0
!------------------------
! Statements
!------------------------
call distance_point_line_2D(point,point_pol_1,point_pol_2,dis1,normal)
if (dis1/=0.d0) call distance_point_line_2D(point_pol_3,point_pol_1,           &
                                            point_pol_2,dis2,normal)
if (dsign(dis1,dis2)/=(dis1)) then
   test = 0
   else
      call distance_point_line_2D(point,point_pol_2,point_pol_3,dis1,normal)
      if (dis1/=0.d0) call distance_point_line_2D(point_pol_1,point_pol_2,     &
                                                  point_pol_3,dis2,normal)
      if (dsign(dis1,dis2)/=(dis1)) then
         test = 0
         elseif (n_sides==3) then
            call distance_point_line_2D(point,point_pol_3,point_pol_1,         &
                                        dis1,normal)
            if (dis1/=0.d0) call distance_point_line_2D(point_pol_2,           &
                                                        point_pol_3,           &
                                                        point_pol_1,dis2,      &
                                                        normal) 
            if (dsign(dis1,dis2)/=(dis1)) then
               test = 0
               else
                  test = 1
            endif
            else
               call distance_point_line_2D(point,point_pol_3,point_pol_4,      &
                                           dis1,normal)
               if (dis1/=0.d0) call distance_point_line_2D(point_pol_1,        &
                                                           point_pol_3,        &
                                                           point_pol_4,        &
                                                           dis2,normal) 
               if (dsign(dis1,dis2)/=(dis1)) then
                  test = 0
                  elseif (n_sides==4) then
                     call distance_point_line_2D(point,point_pol_4,            &
                                                 point_pol_1,dis1,normal)
                     if (dis1/=0.d0) call distance_point_line_2D(              &
                                        point_pol_2,point_pol_4,               &
                                        point_pol_1,dis2,normal) 
                     if (dsign(dis1,dis2)/=(dis1)) then
                        test = 0
                        else
                           test = 1
                     endif
                     else
                        call distance_point_line_2D(point,point_pol_4,         &
                                                    point_pol_5,dis1,normal)
                        if (dis1/=0.d0) call distance_point_line_2D            &
                           (point_pol_1,point_pol_4,point_pol_5,dis2,normal) 
                        if (dsign(dis1,dis2)/=(dis1)) then
                           test = 0
                           elseif (n_sides==5) then
                              call distance_point_line_2D(point,point_pol_5,   &
                                                 point_pol_1,dis1,normal)
                              if (dis1/=0.d0) call distance_point_line_2D      &
                                 (point_pol_2,point_pol_5,point_pol_1,dis2,    &
                                 normal) 
                              if (dsign(dis1,dis2)/=(dis1)) then
                                 test = 0
                                 else
                                    test = 1
                              endif
                              else
                                 call distance_point_line_2D(point,point_pol_5,&
                                                    point_pol_6,dis1,normal)
                                 if (dis1/=0.d0) call distance_point_line_2D   &
                                    (point_pol_1,point_pol_5,point_pol_6,dis2, &
                                    normal) 
                                 if (dsign(dis1,dis2)/=(dis1)) then
                                    test = 0
                                    else
                                       call distance_point_line_2D             &
                                          (point,point_pol_6,point_pol_1,dis1, &
                                          normal)
                                       if (dis1/=0.d0) call                    &
                                          distance_point_line_2D(point_pol_2,  &
                                          point_pol_6,point_pol_1,dis2,normal) 
                                       if (dsign(dis1,dis2)/=(dis1)) then
                                          test = 0
                                          else
                                             test = 1
                                       endif 
                                 endif
                        endif
               endif
      endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine point_inout_convex_non_degenerate_polygon
