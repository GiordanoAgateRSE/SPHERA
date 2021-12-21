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
! Program unit: allocate_de_ch100_r1
! Description: Allocation/deallocation of a generic allocatable array of type 
!              character(100) and range (number of dimensions) 1
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Further Copyright acknowledgments
! On the SPHERA git commit associated with the first appearance of this file in 
! SPHERA (RSE SpA): this file is copied and pasted from Grid Interpolator v.2.0 
! (RSE SpA); the distribution of this file under the GNU-GPL license is 
! authorized by the Copyright owner of Grid Interpolator (RSE SpA).
!-------------------------------------------------------------------------------
subroutine allocate_de_ch100_r1(allocate_flag,array,extent_1,uerr,array_name)
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
character(100),dimension(:),allocatable,intent(inout) :: array
logical,intent(in) :: allocate_flag
integer(4),intent(in) :: extent_1
! error unit
integer(4),intent(in) :: uerr
character(100),intent(in) :: array_name
integer(4) :: alloc_stat
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
if((allocate_flag.eqv..true.).and.(.not.allocated(array))) then
   allocate(array(extent_1),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) "Allocation of ",trim(array_name),                         &
         " failed; the execution stops here."
      stop
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
if((allocate_flag.eqv..false.).and.(allocated(array))) then
   deallocate(array,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) "Deallocation of ",trim(array_name),                       &
         " failed; the execution stops here."
      stop
   endif
endif
return
end subroutine allocate_de_ch100_r1
