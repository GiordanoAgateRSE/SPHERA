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
! Program unit: main_wall_info
! Description: To count the number of vertices of the main wall (i.e., 
!              the solid ("fixe") boundary with the minimum zone ID).
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine main_wall_info(n_vertices_main_wall)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Dynamic_allocation_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(out) :: n_vertices_main_wall
integer(4) :: i_zone
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
n_vertices_main_wall = 0
!------------------------
! Statements
!------------------------
do i_zone=1,NPartZone
   if (Partz(i_zone)%tipo=="fixe") then
      n_vertices_main_wall = Partz(i_zone)%ID_last_vertex_sel -                &
                             Partz(i_zone)%ID_first_vertex_sel + 1
      write(ulog,'(2a,i4,a5,i11)') 'Program unit "main_wall_info": i_zone, ',  &
         'Partz(i_zone)%tipo, n_vertices_main_wall: ',                         &
         i_zone,Partz(i_zone)%tipo,n_vertices_main_wall
      exit
   endif
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine main_wall_info
#endif
