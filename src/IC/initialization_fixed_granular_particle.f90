!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-Ricerca di Sistema)



! SPHERA authors and email contact are provided on SPHERA documentation.

! This file is part of SPHERA v.8.0.
! SPHERA v.8.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with SPHERA. If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: initialization_fixed_granular_particle     
! Description: To initialize the most of the fixed SPH mixture particles (bed-load transport).              
!----------------------------------------------------------------------------------------------------------------------------------

subroutine initialization_fixed_granular_particle(npi)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module 
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: npi
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
pg(npi)%Beta_slope = -999.d0
pg(npi)%Gamma_slope = -999.d0 
pg(npi)%u_star = 0.d0
pg(npi)%C_L = 0.d0
pg(npi)%C_D = 0.d0
pg(npi)%k_BetaGamma = -999.d0
pg(npi)%tau_tauc = 0.0d0
!------------------------
! Statements
!------------------------
!------------------------
! Deallocations
!------------------------
return
end subroutine initialization_fixed_granular_particle

