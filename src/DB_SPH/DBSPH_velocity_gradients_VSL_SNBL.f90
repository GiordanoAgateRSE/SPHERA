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
! Program unit: DBSPH_velocity_gradients_VSL_SNBL
! Description: Computation of the velocity gradients in the Viscous Sub-Layer of   
!              the Surface Neutral Boundary Layer. The gradients are used in the  
!              DB-SPH BC shear viscosity term (DB-SPH-NS). For wall elements, 
!              the numerator and the denominator (wall element Shepard 
!              coefficient without contributions from semi-particles) are   
!              updated independently. Their ratio is computed in 
!              "DBSPH_BC_shear_viscosity_term": here are summed their 
!              contributions. To compute the kinematic viscosity of the 
!              semi-particles (before Shepard correction). Contributions to the 
!              discrete Shepard coefficient of wall elements depending on fluid  
!              particles (not on semi-particles).        
!-------------------------------------------------------------------------------
subroutine DBSPH_velocity_gradients_VSL_SNBL(i_0,i_a,npartint)
!------------------------
! Modules
!------------------------
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: i_0,i_a,npartint
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
! Velocity gradient in VSL for fluid particles
grad_vel_VSL_fw(:,npartint) = (pg(i_0)%vel(:)-pg_w(i_a)%vel(:)) /              &
                              dsqrt(dot_product(rag_fw(:,npartint),            &
                              rag_fw(:,npartint)))  
! Contributions to the numerator of the velocity gradient in VSL for surface 
! elements
pg_w(i_a)%grad_vel_VSL_times_mu(:) = pg_w(i_a)%grad_vel_VSL_times_mu(:) +      &
                                     pg(i_0)%mu * grad_vel_VSL_fw(:,npartint)  & 
                                     * kernel_fw(1,npartint) * pg(i_0)%mass /  &
                                     pg(i_0)%dens
! Contributions to the discrete Shepard coefficient of wall elements depending 
! on fluid particles (not on semi-particles)
pg_w(i_a)%sigma = pg_w(i_a)%sigma + kernel_fw(1,npartint) * pg(i_0)%mass /     &
                  pg(i_0)%dens
! Contributions to the numerator of the kinematic viscosity of the 
! semi-particles
pg_w(i_a)%kin_visc_semi_part = pg_w(i_a)%kin_visc_semi_part + pg(i_0)%kin_visc &
                               * kernel_fw(1,npartint) * pg(i_0)%mass /        &
                               pg(i_0)%dens
!------------------------
! Deallocations
!------------------------
return
end subroutine DBSPH_velocity_gradients_VSL_SNBL
