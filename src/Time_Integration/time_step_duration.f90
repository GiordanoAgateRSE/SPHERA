!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2020 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: time_step_duration                                 
! Description: Computation of the time step duration (dt) according to 
!              stability constraints (CFL condition, viscosity term stability 
!              criterion). Plus, a special treatment for Monaghan artificial 
!              viscosity term and management of low-velocity SPH mixture 
!              particles for bed-load transport phenomena.
!-------------------------------------------------------------------------------
subroutine time_step_duration
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
double precision :: dtmin,dt_CFL,dt_vis,U,celiq,visc_Mon_j
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
dtmin = max_positive_number                                  
dt_Mon = max_positive_number
dt_vis = max_positive_number
!------------------------
! Statements
!------------------------
if (indarrayFlu==0) then
! In case there is no fluid (or mobile) particle
   dt = 0.d0
   else
! Time step depends on 2 stability criteria:
! 1) CFL criterion: dt_CFL<=CFL*(2h/(c+U))
! 2) viscous term criterion: dt_vis<=vsc_coeff*(rho*h**2/(0.5*mu))
! Loop over the media
      do ii=1,NMedium
         if ((Med(ii)%tipo=="liquid  ").and.(Med(ii)%kin_visc>1.d-24)) then
            dt_vis = Domain%vsc_coeff * squareh / (half * Med(ii)%kin_visc)
            dtmin = min(dtmin,dt_vis)
         endif
      enddo
! Loop over the particles
      do ii=1,indarrayFlu
         npi = Array_Flu(ii)
         if (Med(pg(npi)%imed)%tipo=="granular") then
! Redundant and safety check
            if ((pg(npi)%state=="sol").or.(pg(npi)%mu==Med(pg(npi)%imed)%mumx))&
               cycle
            if ((pg(npi)%coord(1)<Granular_flows_options%x_min_dt).or.         &
               (pg(npi)%coord(1)>Granular_flows_options%x_max_dt).or.          &
               (pg(npi)%coord(2)<Granular_flows_options%y_min_dt).or.          &
               (pg(npi)%coord(2)>Granular_flows_options%y_max_dt).or.          &
               (pg(npi)%coord(3)<Granular_flows_options%z_min_dt).or.          &
               (pg(npi)%coord(3)>Granular_flows_options%z_max_dt)) then
               cycle
            endif
            if (pg(npi)%mu>1.d-24) dt_vis = Domain%vsc_coeff * pg(npi)%dens *  &
                                            squareh / (half * pg(npi)%mu)
         endif
         U = sqrt(pg(npi)%vel(1) ** 2 + pg(npi)%vel(2) ** 2 + pg(npi)%vel(3)   &
             ** 2)
         dt_CFL = Domain%CFL * 2.d0 * Domain%h / (Med(pg(npi)%imed)%celerita + &
                  U)
         dtmin = min(dtmin,dt_CFL,dt_vis)
         celiq = Med(pg(npi)%imed)%eps / pg(npi)%dens
         if (celiq>=zero) pg(npi)%csound = dsqrt(celiq)
      enddo
      do j=1,NMedium
         if (dt_alfa_Mon.eqv..true.) then
            visc_Mon_j = Med(j)%alfaMon * Med(j)%celerita * Domain%h /         &
               Med(j)%den0
            dt_Mon_j = squareh / (half * visc_Mon_j)
            dt_Mon = min(dt_Mon,dt_Mon_j)
         endif
      enddo
      if (dt_alfa_Mon.eqv..true.) dtmin = min(dtmin,dt_Mon)
      if (on_going_time_step==0) then
         dt = dtmin
         else
            dt = (one - pesodt) * dtmin + pesodt * dt_average
      endif
endif
if (on_going_time_step==0) then
   dt_average = dt
   else
      dt_average = (dt_average * (on_going_time_step - 1) + dt) /              &
                   on_going_time_step
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine time_step_duration
