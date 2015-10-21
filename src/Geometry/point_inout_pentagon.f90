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
! Program unit: point_inout_pentagon
! Description: Test to evaluate if a point lies inside or strictly outside a generic pentagon. The pentagon is 
!              partitioned into 3 triangles (P1P2P5,P2P3P5,P3P4P5). A point is internal to the pentagon if it is internal to one of
!              its triangles.   
!----------------------------------------------------------------------------------------------------------------------------------
subroutine point_inout_pentagon(point,point_pol_1,point_pol_2,point_pol_3,     &
                                point_pol_4,point_pol_5,test)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
double precision,intent(in) :: point_pol_3(2),point_pol_4(2),point_pol_5(2)
integer(4),intent(inout) :: test
integer(4) :: test1,test2,test3
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine point_inout_convex_non_degenerate_polygon(point,n_sides,         &
                                                        point_pol_1,           &
                                                        point_pol_2,           &
                                                        point_pol_3,           &
                                                        point_pol_4,           &
                                                        point_pol_5,           &
                                                        point_pol_6,test)
      implicit none
      integer(4),intent(in) :: n_sides
      double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
      double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
      double precision,intent(in) :: point_pol_5(2),point_pol_6(2)
      integer(4),intent(inout) :: test
      double precision :: dis1,dis2
      double precision :: normal(2)
   end subroutine point_inout_convex_non_degenerate_polygon
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
test1 = 0
test2 = 0
test3 = 0
!------------------------
! Statements
!------------------------
call point_inout_convex_non_degenerate_polygon(point,3,point_pol_1,            &
                                               point_pol_2,point_pol_5,        &
                                               point_pol_5,point_pol_5,        &
                                               point_pol_5,test1)
test = test1
if (test/=1) then 
   call point_inout_convex_non_degenerate_polygon(point,3,point_pol_2,         &
                                                  point_pol_3,point_pol_5,     &
                                                  point_pol_5,point_pol_5,     &
                                                  point_pol_5,test2)
   test = test2
endif
if (test/=1) then
   call point_inout_convex_non_degenerate_polygon(point,3,point_pol_3,         &
                                                  point_pol_4,point_pol_5,     &
                                                  point_pol_5,point_pol_5,     &
                                                  point_pol_5,test3)
   test = test3
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine point_inout_pentagon

