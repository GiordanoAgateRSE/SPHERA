!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: particle_position_extrusion
! Description: Particle positions extruded from DEM-DTM or any 3D solid bottom.
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine particle_position_extrusion(i_vertex,aux_factor,ii,jj,kk,DEM_zone,  &
   zmax,test_z,pos)
!------------------------
! Modules
!------------------------
use Dynamic_allocation_module
use Hybrid_allocation_module
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: i_vertex,aux_factor,ii,jj,kk,DEM_zone
double precision,intent(in) :: zmax
integer(4),intent(inout) :: test_z
double precision,intent(inout) :: pos(3)
integer(4) :: test_face,i_face,test_xy
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
! Setting coordinates. One avoids particles to be located on the very diagonals 
! of topography (otherwise some particles could be not removed even if below 
! the topography; it is unclear why this rare eventuality would happen, but no 
! problem remains anymore).
pos(1) = Vertice(1,i_vertex) - aux_factor / 2.d0 * Domain%dx + (ii - 1 +       &
         0.501d0) * Domain%dx
pos(2) = Vertice(2,i_vertex) - aux_factor / 2.d0 * Domain%dx + (jj - 1 + 0.5d0)&
         * Domain%dx
pos(3) = (zmax - Domain%dx / 2.d0) - (kk - 1) * Domain%dx
! Test if the particle is below the reservoir
test_z = 0
test_face = 0
! Loop over boundaries
!$omp parallel do default(none)                                                &
!$omp shared(NumFacce,Tratto,BoundaryFace,DEM_zone,test_face,test_z,pos)       &
!$omp private(i_face,aux_vec,aux_scal,test_xy)
do i_face=1,NumFacce
   if (test_face==1) cycle
   if (Tratto(BoundaryFace(i_face)%stretch)%zone==DEM_zone) then  
! Test if the point lies inside the plan projection of the face     
      call point_inout_convex_non_degenerate_polygon(                          &
         pos(1:2),BoundaryFace(i_face)%nodes,                                  &
         BoundaryFace(i_face)%Node(1)%GX(1:2),                                 &
         BoundaryFace(i_face)%Node(2)%GX(1:2),                                 &
         BoundaryFace(i_face)%Node(3)%GX(1:2),                                 &
         BoundaryFace(i_face)%Node(4)%GX(1:2),                                 &
         BoundaryFace(i_face)%Node(4)%GX(1:2),                                 &
         BoundaryFace(i_face)%Node(4)%GX(1:2),test_xy)
      if (test_xy==1) then
! No need for a critical section to update test_face: only one face/process 
! is interested (no conflict); in particular cases there could be a conflict, 
! but no face has a preference and test_face cannot come back to zero 
! (no actual problem)
         test_face = 1
! Test if the point is invisible (test_z=1, below the topography) or visible 
! (test_z=0) to the face
         aux_vec(:) = pos(:) - BoundaryFace(i_face)%Node(1)%GX(:)
         aux_scal = dot_product(BoundaryFace(i_face)%T(:,3),aux_vec)
! No need for a critical section to update test_z: in particular cases there 
! could be a conflict, but the "if" condition would always provide the same 
! result
         if (aux_scal<=0.d0) test_z=1
      endif
   endif
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine particle_position_extrusion
#endif
