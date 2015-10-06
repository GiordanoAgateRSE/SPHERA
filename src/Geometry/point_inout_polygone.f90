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
! Program unit: point_inout_polygone   
! Description: Test to evaluate if a point lies inside or strictly outside a polygone (a triangle or a quadrilateral).      
!----------------------------------------------------------------------------------------------------------------------------------

subroutine point_inout_polygone(point,n_sides,point_pol_1,point_pol_2,         &
                                point_pol_3,point_pol_4,test)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: n_sides
double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
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
         else
            if (n_sides==3) then
               call distance_point_line_2D(point,point_pol_3,point_pol_1,      &
                                           dis1,normal)
               if (dis1/=0.d0) call distance_point_line_2D(point_pol_2,        &
                                                           point_pol_3,        &
                                                           point_pol_1,dis2,   &
                                                           normal) 
               if (dsign(dis1,dis2)/=(dis1)) then
                  test = 0
                  else
                     test = 1
               endif
               else
                  call distance_point_line_2D(point,point_pol_3,point_pol_4,   &
                                              dis1,normal)
                  if (dis1/=0.d0) call distance_point_line_2D(point_pol_1,     &
                                                              point_pol_3,     &
                                                              point_pol_4,     &
                                                              dis2,normal) 
                  if (dsign(dis1,dis2)/=(dis1)) then
                     test = 0
                     else
                        call distance_point_line_2D(point,point_pol_4,         &
                                                    point_pol_1,dis1,normal)
                        if (dis1/=0.d0) call distance_point_line_2D(           &
                                           point_pol_2,point_pol_4,            &
                                           point_pol_1,dis2,normal) 
                        if (dsign(dis1,dis2)/=(dis1)) then
                           test = 0
                           else
                              test = 1
                        endif
                  endif
            endif
      endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine point_inout_polygone

