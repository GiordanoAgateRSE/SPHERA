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
!-------------------------------------------------------------------------------
! Program unit: mixture_viscosity 
! Description: To compute the frictional viscosity for the bed-load transport 
!              layer, according to the formulation of Schaeffer (1987) and the 
!              KTGF in the packing limit.       
!-------------------------------------------------------------------------------
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
double precision :: mu_main_fluid,eps_fluid_blt,p_fluid_blt_top,gamma_fluid
double precision :: z_blt_top_fluid,alfa_TBT
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
mu_main_fluid = Med(Granular_flows_options%ID_main_fluid)%visc *               &
                Med(Granular_flows_options%ID_main_fluid)%den0
! Volume fraction of the fluid phase in the bed-load transport layer
eps_fluid_blt = 1.d0 -                                                         &
   Med(Granular_flows_options%ID_granular)%gran_vol_frac_max
! Specific weight of the fluid phase
gamma_fluid = Med(Granular_flows_options%ID_main_fluid)%den0 * GI
!------------------------
! Statements
!------------------------
!$omp parallel do default(none) shared(nag,Med,pg)                             &
!$omp shared(Granular_flows_options,ind_interfaces,tempo,mu_main_fluid)        &
!$omp shared(eps_fluid_blt,gamma_fluid)                                        &
!$omp private(npi,i_cell,i_aux,i_grid,j_grid,k_grid,p_fluid_blt_top)           &
!$omp private(z_blt_top_fluid,alfa_TBT)
do npi=1,nag
   if ((Med(pg(npi)%imed)%tipo=="granular").and.(pg(npi)%state=="flu")) then
      i_cell = ParticleCellNumber(pg(npi)%coord) 
      i_aux = CellIndices(i_cell,i_grid,j_grid,k_grid)
! Fluid pressure (simplifying assumption: 1D filtration with piezometric lines  
!    in the bed-load transport layer parallel to the local 3D slope of the 
!    bed-load transport layer top)
      if (ind_interfaces(i_grid,j_grid,1)==0) then
         pg(npi)%pres_fluid = 0.d0
         else
            p_fluid_blt_top = pg(ind_interfaces(i_grid,j_grid,2))%pres
            z_blt_top_fluid = pg(ind_interfaces(i_grid,j_grid,2))%coord(3)
! "alfa_TBT" means “Topographic angle at the Bed-load transport layer Top” 
! between the local interface normal (computed on the mixture particle, 
! which locally defines this interface) and the vertical. The formulation for 
! the fluid pressure only depends on this angle, no matter abut the flow 
! direction.
            alfa_TBT = dacos(                                                  &
               pg(ind_interfaces(i_grid,j_grid,3))%normal_int_mixture_top(3))
! The limiter of PI/2 simply allows assigning null fluid pressure in case of 
! "anti-gravity" slopes (very rare and very local in case of bed-load transport)
            alfa_TBT = max(alfa_TBT,PIGRECO/2.d0) 
            pg(npi)%pres_fluid = p_fluid_blt_top + gamma_fluid *               &
                                 (z_blt_top_fluid - pg(npi)%coord(3)) *        &
                                 ((dcos(alfa_TBT)) ** 2)
      endif
! Mean of the effective stresses as the difference between the total pressure 
! and the fluid pressure
! The mean effective stress cannot be negative
! In the presence of a negative fluid pressure, this is zeroed
      pg(npi)%sigma_prime_m = max(pg(npi)%pres - max(pg(npi)%pres_fluid,0.d0), &
                              0.d0) 
! Liquefaction model
      if (Granular_flows_options%t_liq>=0.000001) then
         if ((tempo>=Granular_flows_options%t_q0).and.                         &
            (tempo<=Granular_flows_options%t_liq).and.                         &
            (ind_interfaces(i_grid,j_grid,1)/=0)) then
            pg(npi)%sigma_prime_m = pg(npi)%sigma_prime_m * (1.d0 - (tempo -   &
                                  Granular_flows_options%t_q0) /               &
                                  Granular_flows_options%t_liq) 
         endif
      endif
! secinv=sqrt(I2(e_ij) (secinv is the sqrt of the second inveriant of the 
! strain-rate tensor)
! Frictional viscosity is the mixture viscosity in the bed-load transport layer
      if (pg(npi)%secinv>1.d-9) then 
         pg(npi)%mu = pg(npi)%sigma_prime_m * dsin(                            &
                      Med(Granular_flows_options%ID_granular)%phi) / (2.d0 *   &
                      pg(npi)%secinv) + mu_main_fluid * eps_fluid_blt                 
         else
! Fictitious value representative of perfect uniform flow conditions
            pg(npi)%mu = 0.d0       
      endif
! To save computational time in the transition zone of elastic-platic regime
      if (pg(npi)%mu>Med(Granular_flows_options%ID_granular)%mumx) then
         pg(npi)%mu = Med(Granular_flows_options%ID_granular)%mumx
         pg(npi)%vel(:) = 0.d0
         if (Granular_flows_options%erosion_flag/=1) pg(npi)%state = "sol"
      endif
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

