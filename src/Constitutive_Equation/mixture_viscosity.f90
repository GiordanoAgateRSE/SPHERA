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
! Program unit: mixture_viscosity 
! Description: To compute the frictional viscosity for the bed-load transport 
!              layer, according to the formulation of Schaeffer (1987) and the 
!              KTGF in the packing limit.       
!----------------------------------------------------------------------------------------------------------------------------------
subroutine mixture_viscosity
!------------------------
! Modules
!------------------------ 
use Static_allocation_module 
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: npi,i_cell,i_aux,i_grid,j_grid,k_grid
double precision :: z_int,K_0 
integer(4),external :: ParticleCellNumber,CellIndices 
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
!$omp parallel do default(none) shared(nag,Med,pg)                             &
!$omp shared(Granular_flows_options,Domain,ind_interfaces,tempo,K_0)           &                                                        &
!$omp private(npi,i_cell,i_aux,i_grid,j_grid,k_grid,z_int)
do npi=1,nag
   if ((Med(pg(npi)%imed)%tipo=="granular").and.(pg(npi)%state=="flu")) then
      i_cell = ParticleCellNumber(pg(npi)%coord) 
      i_aux = CellIndices(i_cell,i_grid,j_grid,k_grid)
! Vertical effective stress  
      if (ind_interfaces(i_grid,j_grid,3).ne.0) then                  
         z_int = pg(ind_interfaces(i_grid,j_grid,3))%coord(3)
         pg(npi)%sigma_prime_m = ( - Domain%grav(3)) *                         &
            (Med(Granular_flows_options%ID_granular)%den0 -                    &
            Med(Granular_flows_options%ID_main_fluid)%den0) *                  &
            (z_int - pg(npi)%coord(3))
         if (ind_interfaces(i_grid,j_grid,1)==0) then
            pg(npi)%sigma_prime_m = ( - Domain%grav(3)) *                      &
               Med(Granular_flows_options%ID_granular)%den0 * (z_int -         &
               pg(npi)%coord(3))
         endif  
            else
! This case (probably) never occurs
               pg(npi)%sigma_prime_m = pg(npi)%pres *                          &
                  (Med(Granular_flows_options%ID_granular)%den0 -              &
                  Med(Granular_flows_options%ID_main_fluid)%den0)              &
                  / Med(Granular_flows_options%ID_granular)%den0             
      endif 
      if (pg(npi)%sigma_prime_m<0.d0) pg(npi)%sigma_prime_m = 0.d0
      if (Granular_flows_options%t_liq>=0.000001) then
         if ((tempo>=Granular_flows_options%t_q0).and.                         &
            (tempo<=Granular_flows_options%t_liq).and.                         &
            (ind_interfaces(i_grid,j_grid,1)/=0)) then
            pg(npi)%sigma_prime_m = pg(npi)%sigma_prime_m * (1.d0 - (tempo -   &
                                  Granular_flows_options%t_q0) /               &
                                  Granular_flows_options%t_liq) 
         endif
      endif
! Coefficient of lateral earth pressure at rest
      K_0 = 1.d0 - dsin(Med(Granular_flows_options%ID_granular)%phi)
! Mean of the effective stresses
      pg(npi)%sigma_prime_m = pg(npi)%sigma_prime_m * (2.d0 * K_0 + 1.d0) / 3.d0
! secinv=sqrt(I2(e_ij) (secinv is the sqrt of the second inveriant of the 
! strain-rate tensor)
! Frictional viscosity is the mixture viscosity in the bed-load transport layer
      if (pg(npi)%secinv>1.d-9) then 
         pg(npi)%mu = pg(npi)%sigma_prime_m * dsin(                            &
                      Med(Granular_flows_options%ID_granular)%phi) / (2.d0 *   &
                      pg(npi)%secinv)
         else
            pg(npi)%mu = 0.d0       
      endif
! To save computational time in the transition zone of elastic-platic regime
      if (pg(npi)%mu>Med(Granular_flows_options%ID_granular)%mumx)             &
         pg(npi)%mu = Med(Granular_flows_options%ID_granular)%mumx
! Kinematic viscosity is updated
      pg(npi)%visc = pg(npi)%mu / pg(npi)%dens
   endif
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine mixture_viscosity

