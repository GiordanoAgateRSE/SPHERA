!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
subroutine wall_function_for_SASPH(u_t_0,d_50,r_0w,slip_coefficient_0w,ni_T_0w)
!------------------------
! Modules
!------------------------
use Hybrid_allocation_module
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
! Tangential velocity of the computational particle
double precision,intent(in) :: u_t_0
! Mean diameter of the wall roughness elements around the wall frontier
double precision,intent(in) :: d_50
! Distance between the computational particle and the neighbouring wall element
double precision,intent(in) :: r_0w
! Slip coefficient for the current particle-frontier interaction
double precision,intent(out) :: slip_coefficient_0w
! Mixing-length turbulent kinematic viscosity for the current particle-frontier 
! interaction
double precision,intent(out) :: ni_T_0w
! Roughness length
double precision :: z0
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
slip_coefficient_0w = 0.d0
! Turbulent viscosity is initialized to a non-null value just to avoid the 
! product "0*0" when assessing the shear stress out of the wall-function depth.
ni_T_0w = 1.d-15
!------------------------
! Statements
!------------------------
if (r_0w<=(0.5d0*Domain%dx)) then
! Only the particle layer close to the boundary is selected
   z0 = d_50 / 10.d0
   if (r_0w<(z0*Nepero_number)) then
! Formulation with underestimation of the wall shear stress, but keeping the 
! slip coefficient non-larger than the unity.
      slip_coefficient_0w = 1.d0
      ni_T_0w = (k_v ** 2) * u_t_0 * z0 * Nepero_number
      else
         slip_coefficient_0w = 1.d0 / log(r_0w / z0)
         ni_T_0w = (k_v ** 2) * u_t_0 * r_0w / log(r_0w / z0)
   endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine wall_function_for_SASPH
