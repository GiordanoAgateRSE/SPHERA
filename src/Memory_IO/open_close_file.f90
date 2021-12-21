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
! On the SPHERA git commit associated with the first appearance of this file in 
! SPHERA (RSE SpA): this file is copied and pasted from Grid Interpolator v.2.0 
! (RSE SpA); the distribution of this file under the GNU-GPL license is 
! authorized by the Copyright owner of Grid Interpolator (RSE SpA).
!-------------------------------------------------------------------------------
subroutine open_close_file(open_flag,file_unit,file_name,uerr)
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
character(100),intent(inout) :: file_name
logical,intent(in) :: open_flag
integer(4),intent(in) :: file_unit,uerr
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
   open(file_unit,file=trim(file_name),IOSTAT=open_stat)
   if (open_stat/=0) then
      write(uerr,*) "Error in opening the file ",trim(file_name),              &
         ". The program stops."
      stop
   endif
   else
      close(file_unit)
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine open_close_file
