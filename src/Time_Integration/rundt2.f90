!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)

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
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: rundt2                                          
! Description: Time step computation according to stability constraints (CFL 
!              condition, visosity stability criterion, interface diffusion 
!              criterion -not recommended-). 
!              Plus, a special treatment for Monaghan artificial viscosity term 
!              and management of low-velocity SPH mixture particles for bed-load
!              transport phenomena.
!-------------------------------------------------------------------------------
subroutine rundt2
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
use I_O_file_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: j,npi,ii
double precision :: dtmin,dt_CFL,dt_dif,dt_vis,diffmax,U,celiq,visc_Mon_j
double precision :: dt_Mon_j,dt_Mon
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
dtmin = 1.0d+30                                              
diffmax = zero
dt_Mon = max_positive_number 
!------------------------
! Statements
!------------------------
if (indarrayFlu==0) then
! In case there is no fluid (or mobile) particle
   dt = 0.d0
   else
! Time step depends on 3 stability criteria:
! 1) CFL criterion: dt_CFL=CFL*(2h/(c+U))
! 2) viscous term criterion: dt_vis=vsc_coeff*(rho*h**2/(0.5*mu)) 
! 3) interface diffusion criterion: dt_diff=(h**2/2*teta)
      do ii=1,indarrayFlu
         npi = Array_Flu(ii)
         if ((Granular_flows_options%ID_erosion_criterion==1).and.             &
            (Med(pg(npi)%imed)%tipo=="granular")) then
! Redundant and safety check
            if (pg(npi)%state=="sol") cycle
            if ((pg(npi)%coord(1)<Granular_flows_options%x_min_dt).or.         &
               (pg(npi)%coord(1)>Granular_flows_options%x_max_dt).or.          &
               (pg(npi)%coord(2)<Granular_flows_options%y_min_dt).or.          &
               (pg(npi)%coord(2)>Granular_flows_options%y_max_dt).or.          &
               (pg(npi)%coord(3)<Granular_flows_options%z_min_dt).or.          &
               (pg(npi)%coord(3)>Granular_flows_options%z_max_dt)) then
               cycle
            endif
         endif
         U = sqrt(pg(npi)%vel(1) ** 2 + pg(npi)%vel(2) ** 2 + pg(npi)%vel(3)   &
             ** 2)
         dt_CFL = Domain%CFL * 2.d0 * Domain%h / (Med(pg(npi)%imed)%celerita + &
                  U)
         dt_vis = Domain%vsc_coeff * pg(npi)%dens * squareh / (half *          &
                  pg(npi)%mu)
         dtmin = min(dtmin,dt_CFL,dt_vis)
         celiq = Med(pg(npi)%imed)%eps / pg(npi)%dens
         if (celiq>=zero) pg(npi)%csound = dsqrt(celiq)
      enddo
      do j=1,NMedium
         diffmax = max(diffmax,Med(j)%codif)
         if (dt_alfa_Mon.eqv..true.) then
            visc_Mon_j = Med(j)%alfaMon * Med(j)%celerita * Domain%h /         &
               Med(j)%den0
            dt_Mon_j = squareh / (half * visc_Mon_j)
            dt_Mon = min(dt_Mon,dt_Mon_j)
         endif
      enddo
      if (diffusione) then
         dt_dif = half * squareh / (diffmax+0.000000001d0)
         dtmin = min(dtmin,dt_dif)
      endif
      if (dt_alfa_Mon.eqv..true.) dtmin = min(dtmin,dt_Mon)
      dt = (one - pesodt) * dtmin + pesodt * dt_average
endif
dt_average = (dt_average * (on_going_time_step - 1) + dt) / on_going_time_step
!------------------------
! Deallocations
!------------------------
return
end subroutine rundt2

