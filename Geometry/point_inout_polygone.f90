!AA501b the whole subroutine
!AA504 sub
!cfile Body_dynamics.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : point_inout_polygone
!
!AA504 sub start
! Versions: 
! 01  Amicarelli      24Oct12       (creation) 
! 02  Amicarelli      08apr14       (v5.04) removed minor errors      
!AA504 sub end
!
!************************************************************************************
! Module purpose : test to evaluate if a point lies inside or strictly outside a polygone
!                  (triangle or quadrilateral)
!
!AA504sub
!AA601sub
! Calling routine: RHS_body_dynamics,GeneratePart,DBSPH_inlet_outlet
!
! Called routines: distance_point_line_2D
!
!************************************************************************************

subroutine point_inout_polygone(point,n_sides,point_pol_1,point_pol_2,point_pol_3,point_pol_4,test)

! Declarations
 implicit none
 integer(4),intent(in) :: n_sides
 double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2),point_pol_3(2),point_pol_4(2)
 integer(4),intent(inout) :: test
 double precision :: dis1,dis2
 double precision :: normal(2)
 
! Statements
  call distance_point_line_2D(point,point_pol_1,point_pol_2,dis1,normal)
!AA504 rm line 
  if (dis1.ne.0.) call distance_point_line_2D(point_pol_3,point_pol_1,point_pol_2,dis2,normal)
  if (dsign(dis1,dis2).ne.(dis1)) then
     test = 0
     else
        call distance_point_line_2D(point,point_pol_2,point_pol_3,dis1,normal)
!AA504 rm line
        if (dis1.ne.0.) call distance_point_line_2D(point_pol_1,point_pol_2,point_pol_3,dis2,normal)
        if (dsign(dis1,dis2).ne.(dis1)) then
           test = 0
           else
              if (n_sides == 3) then
                 call distance_point_line_2D(point,point_pol_3,point_pol_1,dis1,normal)
                 if (dis1.ne.0.) call distance_point_line_2D(point_pol_2,point_pol_3,point_pol_1,dis2,normal) 
                 if (dsign(dis1,dis2).ne.(dis1)) then
                    test = 0
                    else
                    test = 1
                 endif
                 else
                    call distance_point_line_2D(point,point_pol_3,point_pol_4,dis1,normal)
                    if (dis1.ne.0.) call distance_point_line_2D(point_pol_1,point_pol_3,point_pol_4,dis2,normal) 
                    if (dsign(dis1,dis2).ne.(dis1)) then
                       test = 0
                       else
                          call distance_point_line_2D(point,point_pol_4,point_pol_1,dis1,normal)
                          if (dis1.ne.0.) call distance_point_line_2D(point_pol_2,point_pol_4,point_pol_1,dis2,normal) 
                          if (dsign(dis1,dis2).ne.(dis1)) then
                             test = 0
                             else
                                test = 1
                          endif
                    endif
              endif
        endif
  endif
        
return
end subroutine point_inout_polygone
!---split

