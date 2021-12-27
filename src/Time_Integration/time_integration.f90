!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2021 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)

! SPHERA authors and email contact are provided on SPHERA documentation.

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
! Program unit: time_integration                                          
! Description: Explicit Runge-Kutta time integration schemes.
!-------------------------------------------------------------------------------
subroutine time_integration  
!------------------------
! Modules
!------------------------
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use Memory_I_O_module
!------------------------
! Declarations
!------------------------
implicit none
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
! Explicit Runge-Kutta time integration schemes (Euler-RK1-, Heun-RK2)
call start_and_stop(2,17)
select case (Domain%RKscheme)
   case (1) 
      call Euler
   case (2) 
      if (Domain%time_stage==1) then
         call Euler
         else
            call Heun
      endif
endselect
call start_and_stop(3,17)
! Equation of State 
call start_and_stop(2,13)
call CalcPre  
call start_and_stop(3,13)
!------------------------
! Deallocations
!------------------------
return
end subroutine time_integration
