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
! Program unit: check_max_file_unit_ID
! Description: Check on the machine/OS-dependent number of existing file units
!-------------------------------------------------------------------------------
subroutine check_max_file_unit_ID(max_file_unit_ID,file_unit_requested,        &
   file_name)
!------------------------
! Modules
!------------------------
use I_O_file_module
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(inout) :: max_file_unit_ID
integer(4),intent(in) :: file_unit_requested
character(100),intent(in) :: file_name
integer(4) :: io_stat,ier
character(100) :: file_name_2
logical,external :: ReadCheck
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine open_close_file(open_flag,I_O_unit,file_name)
      implicit none
      character(100),intent(inout) :: file_name
      logical,intent(in) :: open_flag
      integer(4),intent(in) :: I_O_unit
   end subroutine open_close_file
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
write(ulog,'(a,i9)') "File unit or maximum file unit requested: ",             &
   file_unit_requested
call system ("ulimit -n > ulimit_n.txt")
file_name_2 = "ulimit_n.txt"
call open_close_file(.true.,max_file_unit_booked+1,file_name_2)
read(max_file_unit_booked+1,*,iostat=io_stat) max_file_unit_ID
if (.not.ReadCheck(io_stat,ier,1,"ulimit_n.txt","max_file_unit_ID",            &
   max_file_unit_booked+1,ulog)) then
   write(uerr,*) "Error in reading the file ",file_name_2,". The execution ",  &
      "stops here."
   stop
endif
call open_close_file(.false.,max_file_unit_booked+1,file_name_2)
call system ("rm -f ulimit_n.txt")
write(ulog,'(a,i9)') "Machine/OS-dependent maximum number of file units: ",    &
   max_file_unit_ID
if (max_file_unit_ID<file_unit_requested) then
   write(ulog,'(4a)') "The file unit requested for ",trim(adjustl(file_name)), &
      "is larger than the machine/OS-dependent maximum number of file units: ",&
      "a further check is mandatory."
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine check_max_file_unit_ID
