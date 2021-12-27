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
! Program unit: open_close_file
! Description: File opening or closing.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Further Copyright acknowledgments
! This program unit represents a modification of the program unit 
! open_close_file of Grid Interpolator v.2.0 (RSE SpA). Please refer to 
! the git log file for further details.
!-------------------------------------------------------------------------------
subroutine open_close_file(open_flag,I_O_unit,file_name)
!------------------------
! Modules
!------------------------
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
logical,intent(in) :: open_flag
integer(4),intent(in) :: I_O_unit
character(100),intent(in) :: file_name
integer(4) :: open_stat
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
if (open_flag.eqv..true.) then
   open(I_O_unit,file=trim(file_name),IOSTAT=open_stat)
   if (open_stat/=0) then
      write(uerr,*) "Error in opening the file ",trim(file_name),              &
         ". The program stops."
      stop
   endif
   else
      close(I_O_unit)
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine open_close_file
