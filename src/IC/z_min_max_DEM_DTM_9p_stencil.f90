!-------------------------------------------------------------------------------
! SPHERA v.10.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.10.0.0
! SPHERA v.10.0.0 is free software: you can redistribute it and/or modify
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
! Program unit: z_min_max_DEM_DTM_9p_stencil
! Description: To assess the minimum/maximum height at the 9 vertices of the 
!              DEM-DTM (or any 3D solid bottom) around the given DEM-DTM vertex 
!              (included). The minimum local bottom height is useful for the 
!              fluid body extrusion from topography (IC). The maximum local 
!              bottom height is useful for the BCs of the "zmax" zones 
!              (where the fluid height is imposed).
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine z_min_max_DEM_DTM_9p_stencil(min_flag,i_zone,i_vertex,z_aux)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
logical,intent(in) :: min_flag
integer(4),intent(in) :: i_zone,i_vertex
double precision,intent(inout) :: z_aux
integer(4) :: j_vertex
double precision :: distance_hor
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
if (min_flag.eqv..true.) then
   z_aux = max_positive_number
   else
      z_aux = max_negative_number
endif
! Loop over the sub-selection of vertices of the associated bottom
do j_vertex=Partz(i_zone)%ID_first_vertex_sel,Partz(i_zone)%ID_last_vertex_sel
   distance_hor = dsqrt((Vertice(1,i_vertex) - Vertice(1,j_vertex)) ** 2 +     &
                  (Vertice(2,i_vertex) - Vertice(2,j_vertex)) ** 2)
   if (distance_hor<=(1.5*Partz(i_zone)%dx_CartTopog)) then
      if (min_flag.eqv..true.) then
         z_aux = min(z_aux,Vertice(3,j_vertex))
         else
            z_aux = max(z_aux,Vertice(3,j_vertex))
      endif
   endif
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine z_min_max_DEM_DTM_9p_stencil
#endif
