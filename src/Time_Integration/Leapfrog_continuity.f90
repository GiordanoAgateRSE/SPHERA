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
! Program unit: Leapfrog_continuity
! Description: Leapfrog time integration scheme for: the continuity equation; 
!              the volume balance equation. Mass update.
!-------------------------------------------------------------------------------
subroutine Leapfrog_continuity
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Dynamic_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: npi,ii
double precision :: rho_old,dvolume_dt
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
!$omp parallel do default(none)                                                &
!$omp shared(pg,dt,indarrayFlu,Array_Flu,Domain,Med,input_any_t)               &
!$omp private(npi,ii,rho_old,dvolume_dt)
do ii=1,indarrayFlu
   npi = Array_Flu(ii)
   if ((pg(npi)%cella<=0).or.(pg(npi)%vel_type/="std")) cycle
   rho_old = pg(npi)%dens
   if (Domain%tipo=="bsph") pg(npi)%dden = pg(npi)%dden / pg(npi)%uni
! Continuity equation
   if (Domain%tipo=="semi") pg(npi)%dens = pg(npi)%dens + dt * pg(npi)%dden
   if (input_any_t%density_thresholds==1) then        
      if (pg(npi)%dens<(0.9d0*Med(pg(npi)%imed)%den0)) then
         pg(npi)%dens = 0.9d0*Med(pg(npi)%imed)%den0
      endif
      if (pg(npi)%dens>(1.1d0*Med(pg(npi)%imed)%den0)) then
         pg(npi)%dens = 1.1d0 * Med(pg(npi)%imed)%den0
      endif
   endif
! Mass equation and volume update
   if (Domain%tipo=="semi") then
      if (input_any_t%ALE3) then
         dvolume_dt = -pg(npi)%volume / rho_old * (pg(npi)%dden -              &
                      pg(npi)%dden_ALE12)
         pg(npi)%dmass_dt = rho_old * dvolume_dt + pg(npi)%volume * pg(npi)%dden
         pg(npi)%mass = pg(npi)%mass + pg(npi)%dmass_dt * dt
         pg(npi)%volume = pg(npi)%mass / pg(npi)%dens
      endif
   endif
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine Leapfrog_continuity
