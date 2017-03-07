!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
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
!-------------------------------------------------------------------------------
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
integer(4) :: npi,i_cell,i_aux,i_grid,j_grid,k_grid,alloc_stat
double precision :: mu_main_fluid,eps_fluid_blt,p_fluid_blt_top,gamma_fluid
! "alfa_TBT" means “Topographic angle at the Bed-load transport layer Top” 
! between the local interface normal (computed on the mixture particle, 
! which locally defines this interface) and the vertical. The formulation for 
! the fluid pressure depends on this angle, which approximately represents
! the underground flow direction (Eulerian saturation schemes).
double precision :: alfa_TBT,z_blt_top_fluid,z_soil_bottom,z_inf_sat,z_sat_top
! "alfa_WT" means “Water Table slope (computed on the mixture particle, 
! which locally defines this interface). The formulation for the fluid pressure
! depends on this angle, which approximately represents the underground flow 
! direction (Lagrangian saturation scheme).
double precision :: alfa_WT
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
! Specific weight of the fluid phase
gamma_fluid = Med(Granular_flows_options%ID_main_fluid)%den0 * GI
!------------------------
! Statements
!------------------------
!$omp parallel do default(none) shared(nag,Med,pg)                             &
!$omp shared(Granular_flows_options,ind_interfaces,simulation_time)            &
!$omp shared(eps_fluid_blt,gamma_fluid,mu_main_fluid,nout)                     &
!$omp private(npi,i_cell,i_aux,i_grid,j_grid,k_grid,p_fluid_blt_top)           &
!$omp private(z_blt_top_fluid,alfa_TBT,alfa_WT,z_soil_bottom,z_inf_sat)        &
!$omp private(z_sat_top)
do npi=1,nag
   if ((Med(pg(npi)%imed)%tipo=="granular").and.(pg(npi)%state=="flu")) then
      i_cell = ParticleCellNumber(pg(npi)%coord) 
      i_aux = CellIndices(i_cell,i_grid,j_grid,k_grid)
! Volume fraction of the fluid phase in the bed-load transport layer
      eps_fluid_blt = 1.d0 - Med(pg(npi)%imed)%gran_vol_frac_max
      select case (Granular_flows_options%saturation_scheme)
         case(0) 
! Dry soil
            pg(npi)%pres_fluid = 0.d0         
         case(1)
! Saturation scheme for uniform fully saturated soil
! Fluid pressure
! Phreatic zone (simplifying assumption: 1D filtration with piezometric lines  
! in the bed-load transport layer parallel to the local 3D slope of the 
! bed-load transport layer top)
            if (ind_interfaces(i_grid,j_grid,2)/=0) then
               p_fluid_blt_top = pg(ind_interfaces(i_grid,j_grid,2))%pres
               z_blt_top_fluid = pg(ind_interfaces(i_grid,j_grid,2))%coord(3)
               else
                  p_fluid_blt_top = 0.d0
                  z_blt_top_fluid = pg(ind_interfaces(i_grid,j_grid,3))%coord(3)
            endif
            alfa_TBT = dacos(                                                  &
               pg(ind_interfaces(i_grid,j_grid,3))%normal_int_mixture_top(3))
! The limiter of PI/2 simply allows assigning null fluid pressure in case of 
! "anti-gravity" slopes (very rare and very local in case of bed-load transport)
            alfa_TBT = min(alfa_TBT,PIGRECO/2.d0) 
            pg(npi)%pres_fluid = p_fluid_blt_top + gamma_fluid *               &
                                 (z_blt_top_fluid - pg(npi)%coord(3)) *        &
                                 ((dcos(alfa_TBT)) ** 2)      
         case(2)
! Saturation scheme depending on t_minimum_saturation and t_maximum_sateration
            if (Granular_flows_options%saturation_conditions(i_grid,j_grid)==1)&
               then
! Fluid pressure
! Phreatic zone (simplifying assumption: 1D filtration with piezometric lines  
! in the bed-load transport layer parallel to the local 3D slope of the 
! bed-load transport layer top)
               p_fluid_blt_top = pg(ind_interfaces(i_grid,j_grid,2))%pres
               z_blt_top_fluid = pg(ind_interfaces(i_grid,j_grid,2))%coord(3)
               alfa_TBT = dacos(                                               &
                  pg(ind_interfaces(i_grid,j_grid,3))%normal_int_mixture_top(3))
! The limiter of PI/2 simply allows assigning null fluid pressure in case of 
! "anti-gravity" slopes (very rare and very local in case of bed-load transport)
               alfa_TBT = min(alfa_TBT,PIGRECO/2.d0) 
               pg(npi)%pres_fluid = p_fluid_blt_top + gamma_fluid *            &
                                    (z_blt_top_fluid - pg(npi)%coord(3)) *     &
                                    ((dcos(alfa_TBT)) ** 2)
               elseif (Granular_flows_options%saturation_conditions(i_grid,    &
                  j_grid)==2) then
! Infiltration zone
                  z_blt_top_fluid = pg(ind_interfaces(i_grid,j_grid,2))%coord(3)
                  z_soil_bottom = pg(ind_interfaces(i_grid,j_grid,5))%coord(3)
                  z_inf_sat = z_blt_top_fluid - (simulation_time -             &
                              Granular_flows_options%time_minimum_saturation) /&
                              (Granular_flows_options%time_maximum_saturation -&
                              Granular_flows_options%time_minimum_saturation) *&
                              (z_blt_top_fluid - z_soil_bottom)
                  if (pg(npi)%coord(3)>=z_inf_sat) then
                     p_fluid_blt_top = pg(ind_interfaces(i_grid,j_grid,2))%pres
                     pg(npi)%pres_fluid = p_fluid_blt_top * (1.d0 -            &
                                          (z_blt_top_fluid - pg(npi)%coord(3)) &
                                          / (z_blt_top_fluid - z_soil_bottom))
                     else
                        pg(npi)%pres_fluid = 0.d0
                  endif
                  else
! Dry soil
                     if (Granular_flows_options%saturation_conditions(i_grid,  &
                        j_grid)==3) pg(npi)%pres_fluid = 0.d0   
            endif
         case(3)
! Lagrangian scheme for saturation conditions under the hypothesis of 
! stratified flows. Mixture particòes are either fully saturated or dry. 
! Fluid pressure
! Phreatic zone (simplifying assumption: 1D filtration with piezometric lines  
! in the bed-load transport layer parallel to the local 3D slope of the 
! bed-load transport layer top)
            pg(npi)%pres_fluid = 0.d0
            if ((ind_interfaces(i_grid,j_grid,6)/=0).and.(                     &
               Med(pg(npi)%imed)%saturated_medium_flag.eqv..true.)) then
               z_sat_top = pg(ind_interfaces(i_grid,j_grid,6))%coord(3)
               if (z_sat_top>pg(npi)%coord(3)) then
                  alfa_WT = dacos(                                             &
                   pg(ind_interfaces(i_grid,j_grid,3))%normal_int_sat_top(3))
! The limiter of PI/2 simply allows assigning null fluid pressure in case of 
! "anti-gravity" slopes (very rare and very local in case of bed-load transport)
                  alfa_WT = min(alfa_WT,PIGRECO/2.d0) 
                  pg(npi)%pres_fluid = gamma_fluid * (z_sat_top -              &
                                       pg(npi)%coord(3)) * ((dcos(alfa_WT)) ** &
                                       2)
               endif
            endif
         case default 
            write(nout,*) 'The saturation scheme chosen in the input file is ',&
               'wrong. The program terminates here. '
            stop
      endselect
! Mean of the effective stresses as the difference between the total pressure 
! and the fluid pressure
! The mean effective stress cannot be negative
! In the presence of a negative fluid pressure, this is zeroed
      pg(npi)%sigma_prime_m = max(pg(npi)%pres - max(pg(npi)%pres_fluid,0.d0), &
                              0.d0) 
! Liquefaction model
      if (Granular_flows_options%t_liq>=0.000001) then
         if ((simulation_time>=Granular_flows_options%t_q0).and.               &
            (simulation_time<=Granular_flows_options%t_liq).and.               &
            (pg(npi)%pres_fluid>0.d0)) then
            pg(npi)%sigma_prime_m = pg(npi)%sigma_prime_m * (1.d0 -            &
                                    (simulation_time -                         &
                                    Granular_flows_options%t_q0) /             &
                                    Granular_flows_options%t_liq) 
         endif
      endif
! secinv=sqrt(I2(e_ij)) (secinv is the sqrt of the second inveriant of the 
! strain-rate tensor)
! Mixture viscosity in the bed-load transport layer
      if (pg(npi)%secinv>1.d-12) then
         pg(npi)%mu = pg(npi)%sigma_prime_m * dsin(Med(pg(npi)%imed)%phi) /    &
                      (2.d0 * pg(npi)%secinv) + mu_main_fluid * eps_fluid_blt
         else
! This condition simply avoids a division by zero (elastic strain rate regime)
            pg(npi)%mu = Med(pg(npi)%imed)%mumx
      endif
! No matter about the presence/absence of an erosion criterion, the particles 
! in the transition zone of elastic-plastic regime are held fixed. This is 
! consistent with negible strain rates and saves computational time in the 
! transition zone in the elasto-platic strain rate regime.
      if (pg(npi)%mu>=Med(pg(npi)%imed)%mumx) then
         pg(npi)%mu = Med(pg(npi)%imed)%mumx
         pg(npi)%vel(:) = 0.d0
         elseif(pg(npi)%mu>Med(pg(npi)%imed)%limiting_viscosity) then
            pg(npi)%mu = Med(pg(npi)%imed)%limiting_viscosity
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

