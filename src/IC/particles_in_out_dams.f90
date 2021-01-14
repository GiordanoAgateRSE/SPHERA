!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2020 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: particles_in_out_dams
! Description: To test if the SPH particles (preliminarily designed at this 
!              stage) lie out of (test_dam=0) or within (test_dam=1) the 
!              volume of a dam (initial conditions).
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine particles_in_out_dams(fluid_zone,pos,test_dam)
!------------------------
! Modules
!------------------------
use Dynamic_allocation_module
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: fluid_zone
double precision, intent(in) :: pos(3)
integer(4),intent(out) :: test_dam
integer(4) :: test_xy,test_face,i_face,test_xy_2
double precision :: aux_scal
double precision :: aux_vec(3)
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
!------------------------
! Statements
!------------------------
! Check the presence of a dam zone
test_dam = 0
if (Partz(fluid_zone)%dam_zone_ID>0) then
! Test if the point lies inside the plan projection of the dam zone
   call point_inout_convex_non_degenerate_polygon(pos(1:2),                    &
      Partz(fluid_zone)%dam_zone_n_vertices,                                   &
      Partz(fluid_zone)%dam_zone_vertices(1,1:2),                              &
      Partz(fluid_zone)%dam_zone_vertices(2,1:2),                              &
      Partz(fluid_zone)%dam_zone_vertices(3,1:2),                              &
      Partz(fluid_zone)%dam_zone_vertices(4,1:2),                              &
      Partz(fluid_zone)%dam_zone_vertices(4,1:2),                              &
      Partz(fluid_zone)%dam_zone_vertices(4,1:2),test_xy)
   if (test_xy==1) then
      test_face = 0
! Loop on the faces of the dam in the dam zone
!$omp parallel do default(none)                                                &
!$omp shared(NumFacce,test_dam,Tratto,BoundaryFace,Partz,fluid_zone,pos)       &
!$omp shared(test_face)                                                        &
!$omp private(i_face,aux_vec,aux_scal,test_xy_2)
      do i_face=1,NumFacce
         if (test_face==1) cycle
! Check if the face belongs to the dam_zone boundary and in particular to its 
! top
if ((Tratto(BoundaryFace(i_face)%stretch)%zone==Partz(fluid_zone)%dam_zone_ID).and.    &
            (BoundaryFace(i_face)%T(3,3)<0.d0)) then 
! Test if the particle horizontal coordinates lie inside the horizontal 
! projection of the face
            call point_inout_convex_non_degenerate_polygon(pos(1:2),           &
               BoundaryFace(i_face)%nodes,BoundaryFace(i_face)%Node(1)%GX(1:2),&
               BoundaryFace(i_face)%Node(2)%GX(1:2),                           &
               BoundaryFace(i_face)%Node(3)%GX(1:2),                           &
               BoundaryFace(i_face)%Node(4)%GX(1:2),                           &
               BoundaryFace(i_face)%Node(4)%GX(1:2),                           &
               BoundaryFace(i_face)%Node(4)%GX(1:2),test_xy_2)
            if (test_xy_2==1) then
! Test if the particle position is invisible (test_dam=1) or visible 
! (test_dam=0) to the dam. Visibility to the dam means invisibility to the 
! current dam top face. 
               aux_vec(:) = pos(:) - BoundaryFace(i_face)%Node(1)%GX(:)
               aux_scal = dot_product(BoundaryFace(i_face)%T(:,3),aux_vec)
! Note: even if a particle simultanously belongs to 2 faces test_dam will 
! provide the same results: no matter the face order
!$omp critical (GeneratePart_cs)
               test_face = 1
               if (aux_scal<0.d0) then 
                  test_dam=0
                  else
                     test_dam = 1
               endif
!$omp end critical (GeneratePart_cs)
            endif   
         endif
      enddo
!$omp end parallel do                      
   endif   
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine particles_in_out_dams
#endif
