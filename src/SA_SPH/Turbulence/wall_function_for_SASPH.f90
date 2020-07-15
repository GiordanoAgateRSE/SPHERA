!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2020 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: wall_function_for_SASPH
! Description: Assessment of the mixing-length turbulent viscosity and the slip 
!              coefficient only for the SASPH boundary terms of the momentum 
!              equation, depending on the wall function of the Surface Neutral 
!              Boundary Layer (SNBL) for rough walls. The slip coefficient also 
!              affects the partial smoothing of the velocity field.
!-------------------------------------------------------------------------------
subroutine wall_function_for_SASPH(r_0w,u_t_0,d_50,rho_0,slip_coefficient_0w,  &
   mu_T_0w)
!------------------------
! Modules
!------------------------
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
! Tangential velocity of the computational particle
double precision,intent(in) :: u_t_0
! Mean diameter of the wall roughness elements around the wall frontier
double precision,intent(in) :: d_50
! Density of the computational particle
double precision,intent(in) :: rho_0
! Distance between the computational particle and the neighbouring wall element
double precision,intent(in) :: r_0w
! Slip coefficient for the current particle-frontier interaction
double precision,intent(out) :: slip_coefficient_0w
! Mixing-length turbulent viscosity for the current particle-frontier 
! interaction
double precision,intent(out) :: mu_T_0w
! Roughness length
double precision :: z_0
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
z_0 = d_50 / 30.d0
if (r_0w<(z_0*Nepero_number)) then
   slip_coefficient_0w = 1.d0
   mu_T_0w = rho_0 * (k_v ** 2) * u_t_0 * z_0 * Nepero_number
   else
      slip_coefficient_0w = 1.d0 / log(r_0w / z_0)
      mu_T_0w = rho_0 * (k_v ** 2) * u_t_0 * r_0w / log(r_0w / z_0)
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine wall_function_for_SASPH
