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
! Program unit: main_wall_info
! Description: To compute/recall the number of vertices of the main wall (i.e., 
!              the solid boundary used for extrusion of water bodies from 
!              topography). In case of no extrusion, the first "fixed" (wall) 
!              boundary/zone is selected and set up (with all its vertices) for 
!              assessing the 2D synthetic quantities (e.g., maximum water 
!              depth, maximum specific flow rate).
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine main_wall_info(n_vertices_main_wall)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(out) :: n_vertices_main_wall
logical :: aux_logical
integer(4) :: aux_integer,i_zone
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
aux_logical = .false.
aux_integer = 0
n_vertices_main_wall = 0
!------------------------
! Statements
!------------------------
do i_zone=1,NPartZone
   if (Partz(i_zone)%IC_source_type==2) then
! A main wall had already been set up
      aux_logical = .true.
      n_vertices_main_wall = Partz(i_zone)%ID_last_vertex -                    &
                             Partz(i_zone)%ID_first_vertex + 1
      exit
   endif
enddo
if (aux_logical.eqv..false.) then
   do i_zone=1,NPartZone
      if (Partz(i_zone)%tipo=="fixe") then
         Partz(i_zone)%ID_first_vertex = aux_integer + 1
         Partz(i_zone)%ID_last_vertex = Partz(i_zone)%ID_first_vertex +        &
                                        Tratto(i_zone)%numvertices - 1
         n_vertices_main_wall = Tratto(i_zone)%numvertices
         exit
      endif
! Temporary count of all the vertices of the previous zones
      aux_integer = aux_integer + Tratto(i_zone)%numvertices
   enddo
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine main_wall_info
#endif
