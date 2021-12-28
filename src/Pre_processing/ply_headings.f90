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
! Program unit: ply_headings                               
! Description: Reading the headings of a generic ".ply" file to obtain the 
!              number of vertices and faces
!-------------------------------------------------------------------------------
subroutine ply_headings(I_O_unit,file_name,n_vertices,n_faces)
!------------------------
! Modules
!------------------------
use I_O_file_module
use Memory_I_O_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: I_O_unit
character(100),intent(in) :: file_name
integer(4),intent(out) :: n_vertices,n_faces
integer(4) :: io_stat,ier
character(100) :: aux_char,aux_char_2
logical,external :: ReadCheck
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
call open_close_file(.true.,I_O_unit,file_name)
read(I_O_unit,"(3/)",IOSTAT=io_stat)
if (.not.ReadCheck(io_stat,ier,1,file_name,'".ply" headings',I_O_unit,ulog))   &
   then
   write(uerr,*) "Error in reading the file ",file_name,". The execution ",    &
      "stops here."
   stop
endif
read(I_O_unit,*,IOSTAT=io_stat) aux_char,aux_char_2,n_vertices
if (.not.ReadCheck(io_stat,ier,5,file_name,'".ply" headings',I_O_unit,ulog))   &
   then
   write(uerr,*) "Error in reading the file ",file_name,". The execution ",    &
      "stops here."
   stop
endif
read(I_O_unit,"(2/)",IOSTAT=io_stat)
if (.not.ReadCheck(io_stat,ier,6,file_name,'".ply" headings',I_O_unit,ulog))   &
   then
   write(uerr,*) "Error in reading the file ",file_name,". The execution ",    &
      "stops here."
   stop
endif
read(I_O_unit,*,IOSTAT=io_stat) aux_char,aux_char_2,n_faces
if (.not.ReadCheck(io_stat,ier,9,file_name,'".ply" headings',I_O_unit,ulog))   &
   then
   write(uerr,*) "Error in reading the file ",file_name,". The execution ",    &
      "stops here."
   stop
endif
read(I_O_unit,"(1/)",IOSTAT=io_stat)
if (.not.ReadCheck(io_stat,ier,10,file_name,'".ply" headings',I_O_unit,ulog))  &
   then
   write(uerr,*) "Error in reading the file ",file_name,". The execution ",    &
      "stops here."
   stop
endif
call open_close_file(.false.,I_O_unit,file_name)
!------------------------
! Deallocations
!------------------------
return
end subroutine ply_headings
