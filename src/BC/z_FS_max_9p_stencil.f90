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
! Program unit: z_FS_max_9p_stencil
! Description: To assess the maximum free surface height in the 9-point stencil 
! around the current vertex around which the extrusion takes place. 
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine z_FS_max_9p_stencil(i_zone,i_vertex,z_FS_max)
!------------------------
! Modules
!------------------------
use Dynamic_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: i_zone,i_vertex
double precision,intent(inout) :: z_FS_max
integer(4) :: j_vertex,GridColumn
double precision :: distance_hor
double precision :: pos(3)
integer(4),external :: ParticleCellNumber
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
! Reference minimum height in the presence of dry bottom
z_FS_max = Partz(i_zone)%BC_zmax_vertices(i_vertex,3)
! In order to prevent the unefficient introduction of an array to store the 
! bottom vertices per cell and keep the computation fater and lighter, this 
! loop only concerns the zone vertices. Under this frame, mass penetration can 
! occur at the edges of the "zmax" zone, where the free surface grows with the 
! distance from the zone. However, the open sections associated to the "zmax" 
! zone can properly avoid this possible shortcoming in terms of input file 
! configuration.
!$omp parallel do default(none)                                                &
!$omp shared(Partz,i_zone,i_vertex,Grid,z_FS_max,Z_fluid_step)                 &
!$omp private(j_vertex,distance_hor,pos,GridColumn)
do j_vertex=1,size(Partz(i_zone)%BC_zmax_vertices,1)
   distance_hor = dsqrt((Partz(i_zone)%BC_zmax_vertices(i_vertex,1) -          &
                  Partz(i_zone)%BC_zmax_vertices(j_vertex,1)) ** 2 +           &
                  (Partz(i_zone)%BC_zmax_vertices(i_vertex,2) -                &
                  Partz(i_zone)%BC_zmax_vertices(j_vertex,2)) ** 2)
   if (distance_hor<=(1.5*Partz(i_zone)%dx_CartTopog)) then
! The test vertex is one of the 9 close vertices of the current vertex 
      pos(1) = Partz(i_zone)%BC_zmax_vertices(j_vertex,1)
      pos(2) = Partz(i_zone)%BC_zmax_vertices(j_vertex,2)
      pos(3) = Grid%extr(3,1) + 1.d-7
      GridColumn = ParticleCellNumber(pos)
!$omp critical (omp_z_FS_max_9p_stencil_cs)
! "Z_fluid_step" refers to the revious time step. This time integration 
! truncation error cannot be major issue.
      z_FS_max = max(z_FS_max,(Z_fluid_step(GridColumn,1)))
!$omp end critical (omp_z_FS_max_9p_stencil_cs)
   endif
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine z_FS_max_9p_stencil
#endif
