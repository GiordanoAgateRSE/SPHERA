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
! Program unit: time_elapsed_IC
! Description: To assess the elapsed time at the end of the initial conditions    
!-------------------------------------------------------------------------------
subroutine time_elapsed_IC
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: machine_Julian_day,machine_hour,machine_minute,machine_second
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
if (exetype=="linux") then
   if (input_any_t%tmax>1.d-15) then
      call system("date +%j%H%M%S>date_pre_iterations.txt")
      open(unit_time_elapsed,file='date_pre_iterations.txt',status="unknown",  &
         form="formatted")
      read(unit_time_elapsed,'(i3,i2,i2,i2)') machine_Julian_day,machine_hour, &
         machine_minute,machine_second
      close(unit_time_elapsed)
      Domain%t_pre_iter = machine_Julian_day * 24 * 60 * 60 + machine_hour *   &
                          60 * 60 + machine_minute * 60 + machine_second
   endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine time_elapsed_IC
