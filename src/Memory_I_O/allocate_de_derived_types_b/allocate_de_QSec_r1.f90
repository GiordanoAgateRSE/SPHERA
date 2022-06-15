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
! Program unit: allocate_de_QSec_r1
! Description: Allocation/Deallocation of a generic allocatable array of 
!              derived type tyQ_section_array and range (number of 
!              dimensions) 1
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine allocate_de_QSec_r1(allocation_flag,array,extent_1,array_name,      &
   ulog_flag)
!------------------------
! Modules
!------------------------
use I_O_file_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
type (tyQ_section_array),dimension(:),allocatable,intent(inout) :: array
logical,intent(in) :: allocation_flag,ulog_flag
integer(4),intent(in),optional :: extent_1
character(100),intent(in) :: array_name
integer(4) :: alloc_stat
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
if (allocation_flag.eqv..true.) then
   if(.not.allocated(array)) then
      allocate(array(extent_1),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(uerr,*) "Allocation of ",trim(adjustl(array_name)),             &
            " failed; the execution stops here."
         stop
         else
            if (ulog_flag.eqv..true.) then
!$omp critical (omp_Memory_I_O_ulog)
               write(ulog,*) "Allocation of ",trim(adjustl(array_name)),       &
                  " completed."
!$omp end critical (omp_Memory_I_O_ulog)
            endif
      endif
   endif
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
!------------------------
! Deallocations
!------------------------
   else
      if(allocated(array)) then
         deallocate(array,STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(uerr,*) "Deallocation of ",trim(adjustl(array_name)),        &
               " failed; the execution stops here."
            stop
            else
               if (ulog_flag.eqv..true.) then
!$omp critical (omp_Memory_I_O_ulog)
                  write(ulog,*) "Deallocation of ",trim(adjustl(array_name)),  &
                     " completed."
!$omp end critical (omp_Memory_I_O_ulog)
               endif
         endif
      endif
endif
return
end subroutine allocate_de_QSec_r1
#endif