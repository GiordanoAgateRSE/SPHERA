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
! Program unit: IC_CAE_bodies
! Description: Initial Conditions for the CAE-made solid bodies      
!-------------------------------------------------------------------------------
#if (defined SPACE_3D) && (defined SOLID_BODIES)
subroutine IC_CAE_bodies
!------------------------
! Modules
!------------------------
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_file_module
use Memory_I_O_interface_module
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: io_stat,n_vtu_grids,ier,i_vtu_grid,n_vtu_cells,aux_scal
character(100) :: file_name,array_name,aux_char,aux_char_2
type(body_particle),dimension(:),allocatable :: aux_bp_arr
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine vtu_file_reading(i_vtu_grid,file_name,n_vtu_cells)
      implicit none
      integer(4),intent(in) :: i_vtu_grid
      character(100),intent(in) :: file_name
      integer(4),intent(inout) :: n_vtu_cells
   end subroutine vtu_file_reading
   logical function ReadCheck(IoErr,Ier,Nrighe,ainp,listadati,ninp,ulog)
      implicit none
      integer(4) :: IoErr,Ier,Nrighe,ninp,ulog
      character(*) :: ainp,listadati
   end function ReadCheck
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
n_bodies_CAE = 0
n_body_part_CAE = 0
!------------------------
! Statements
!------------------------
! Reading the ".vtu" files: start
! Counting of the number of «.vtu» files
call system ("find ./input/10_body_dynamics/vtu/*.vtu | wc -l > test_vtu.txt")
file_name = "test_vtu.txt"
call open_close_file(.true.,max_file_unit_booked+1,file_name)
read(max_file_unit_booked+1,*,iostat=io_stat) n_vtu_grids
if (.not.ReadCheck(io_stat,ier,1,file_name,"n_vtu_grids",                      &
   max_file_unit_booked+1,ulog)) then
   write(uerr,*) 'Error in reading the file ',file_name,'. The execution ',    &
      'stops here (program unit "IC_CAE_bodies").'
   stop
endif
call open_close_file(.false.,max_file_unit_booked+1,file_name)
call system ("rm -f test_vtu.txt")
! Number of CAE-made bodies
n_bodies_CAE = n_vtu_grids
! Exit the program unit in the absence of CAE-made solid bodies
if (n_vtu_grids==0) return
! Allocation of the ".vtu" grids
array_name = "vtu_grids"
call allocate_de_vtu_grid_r1(.true.,vtu_grids,n_vtu_grids,array_name,          &
   ulog_flag=.true.)
call system("ls input/10_body_dynamics/vtu/*.vtu > test_vtu.txt")
file_name = "test_vtu.txt"
call open_close_file(.true.,max_file_unit_booked+1,file_name)
do i_vtu_grid=1,n_vtu_grids
! Loop over all the ".vtu" file names 
! Read the file name
! "(a27)" is the substring format for the file path 
   read(max_file_unit_booked+1,'(a27,a)',iostat=io_stat) aux_char,aux_char_2
   if (.not.ReadCheck(io_stat,ier,1,file_name,file_name,                       &
      max_file_unit_booked+1,ulog)) then
      write(uerr,*) "Error in reading the file ",file_name,". The execution ", &
         'stops here (program unit "IC_CAE_bodies").'
      stop
   endif
! Read the ".vtu" file
   file_name = "input/10_body_dynamics/vtu/" // trim(adjustl(aux_char_2))
   n_vtu_cells = 0
   call vtu_file_reading(i_vtu_grid,file_name,n_vtu_cells)
! Counting update for the number of CAE-made body particles
   n_body_part_CAE = n_body_part_CAE + n_vtu_cells
! Number of body particles of the current body
   body_arr(n_bodies-n_bodies_CAE+i_vtu_grid)%npart = n_vtu_cells
enddo
! Close and remove the temporary files with the ".vtu" file names
file_name = "test_vtu.txt"
call open_close_file(.false.,max_file_unit_booked+1,file_name)
call system("rm -f test_vtu.txt")
! Reading the ".vtu" files: end
! Counting update for the number of body particles
n_body_part = n_body_part + n_body_part_CAE
! Reallocation of the array of body particles: start
! Allocation of the auxiliary array of body particles
array_name = "aux_bp_arr"
aux_scal = n_body_part - n_body_part_CAE
call allocate_de_BodPar_r1(.true.,aux_bp_arr,aux_scal,array_name,              &
   ulog_flag=.true.)
aux_bp_arr(1:aux_scal) = bp_arr(1:aux_scal)
! Temporary deallocation of the array of body particles
array_name = "bp_arr"
call allocate_de_BodPar_r1(.false.,bp_arr,array_name=array_name,               &
   ulog_flag=.true.)
! Reallocation of the array of body particles
array_name = "bp_arr"
call allocate_de_BodPar_r1(.true.,bp_arr,n_body_part,array_name,               &
   ulog_flag=.true.)
bp_arr(1:aux_scal) = aux_bp_arr(1:aux_scal) 
! Reallocation of the array of body particles: end
! Assignment of the variables of the CAE-made SPH body particles
call IC_CAE_body_particles
!------------------------
! Deallocations
!------------------------
do i_vtu_grid=1,n_vtu_grids
! Deallocation of the array of the ".vtu" cell types
   array_name = "vtu_cells_points_IDs"
   call allocate_de_int4_r2(.false.,vtu_grids(i_vtu_grid)%cells%points_IDs,    &
      array_name=array_name,ulog_flag=.true.)
! Deallocation of the array of the ".vtu" cell "surface"
   array_name = "vtu_cells_surface"
   call allocate_de_log_r2(.false.,vtu_grids(i_vtu_grid)%cells%surface,        &
      array_name=array_name,ulog_flag=.true.)
! Deallocation of the array of the ".vtu" point positions
   array_name = "vtu_points_vertices"
   call allocate_de_vertex_r1(.false.,vtu_grids(i_vtu_grid)%points%vertex,     &
      array_name=array_name,ulog_flag=.true.)
! Deallocation of the array of the ".vtu" point "surface"
   array_name = "vtu_points_surface"
   call allocate_de_log_r1(.false.,vtu_grids(i_vtu_grid)%points%surface,       &
      array_name=array_name,ulog_flag=.true.)
enddo
! Deallocation of the ".vtu" grids
array_name = "vtu_grids"
call allocate_de_vtu_grid_r1(.false.,vtu_grids,array_name=array_name,          &
   ulog_flag=.true.)
! Deallocation of the auxiliary array of body particles
array_name = "aux_bp_arr"
call allocate_de_BodPar_r1(.false.,aux_bp_arr,array_name=array_name,           &
   ulog_flag=.true.)
return
end subroutine IC_CAE_bodies
#endif
