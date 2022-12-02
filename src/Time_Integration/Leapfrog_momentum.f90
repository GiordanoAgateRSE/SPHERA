!-------------------------------------------------------------------------------
! SPHERA v.10.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.10.0.0
! SPHERA v.10.0.0 is free software: you can redistribute it and/or modify
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
! Program unit: Leapfrog_momentum
! Description: Leapfrog time integration scheme (momentum equation). Detection 
!              of the frozen-mass particles (ALE3).
!-------------------------------------------------------------------------------
subroutine Leapfrog_momentum(dt_previous_step,dtvel)
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
double precision,intent(in) :: dt_previous_step
double precision,intent(inout) :: dtvel
integer(4) :: ii,npi,aux_int
! Minimum and maximum value for the particle mass (ALE3)
double precision :: ALE3_min_p_mass,ALE3_max_p_mass
double precision,dimension(3) :: dmom_dt
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
! Initializations for the selection of the frozen-mass particles
#ifdef SPACE_3D
aux_int = 3
#elif defined SPACE_2D
aux_int = 2
#endif
ALE3_min_p_mass = (Domain%dx ** aux_int) * Med(1)%den0 * (1.d0 -               &
                  input_any_t%max_delta_mass / 100.d0) 
ALE3_max_p_mass = (Domain%dx ** aux_int) * Med(1)%den0 * (1.d0 +               &
                  input_any_t%max_delta_mass / 100.d0)
!------------------------
! Statements
!------------------------
! dt computation
dtvel = half * (dt + dt_previous_step) 
!$omp parallel do default(none)                                                &
!$omp shared(pg,dtvel,indarrayFlu,Array_Flu,input_any_t,mass_frozen_count)     &
!$omp shared(ALE3_min_p_mass,ALE3_max_p_mass)                                  &
!$omp private(npi,ii,dmom_dt)
do ii=1,indarrayFlu
   npi = Array_Flu(ii)
#ifdef SOLID_BODIES
   if (pg(npi)%cella==-2) cycle
#endif
   if (input_any_t%ALE3) then
      dmom_dt(1:3) = pg(npi)%mass * pg(npi)%acc(1:3) + pg(npi)%vel(1:3) *      &
                     pg(npi)%dmass_dt
      pg(npi)%mom(1:3) = pg(npi)%mom(1:3) + dtvel * dmom_dt(1:3)
      pg(npi)%vel(1:3) = pg(npi)%mom(1:3) / pg(npi)%mass
      else
         pg(npi)%vel(:) = pg(npi)%vel(:) + dtvel * pg(npi)%acc(:)
   endif
! kodvel=0: the particle is internal to the domain. 
! kodvel=1: the particle has a critical flux condition. The vertical 
! velocity component is assigned.
   if (pg(npi)%kodvel==1) then
      pg(npi)%vel(3) = pg(npi)%velass(3)
      if (input_any_t%ALE3) then
         pg(npi)%mom(3) = pg(npi)%vel(3) * pg(npi)%mass
      endif
! kodvel=2: the particle has an assigned normal velocity or source 
! condition. All the velocity components are assigned.
      elseif (pg(npi)%kodvel==2) then
         pg(npi)%vel(1:3) = pg(npi)%velass(1:3)
         if (input_any_t%ALE3) then
            pg(npi)%mom(1:3) = pg(npi)%vel(1:3) * pg(npi)%mass
         endif
   endif
! Selection of the frozen-mass particles: start
   if ((input_any_t%max_delta_mass>-1.d-21).and.(input_any_t%ALE3)) then
      if ((.not.(pg(npi)%mass_frozen)).and.((pg(npi)%mass<ALE3_min_p_mass).or. &
         (pg(npi)%mass>ALE3_max_p_mass))) then
         pg(npi)%mass_frozen = .true.
!$omp critical (omp_mass_frozen_count)
         mass_frozen_count = mass_frozen_count + 1
!$omp end critical (omp_mass_frozen_count)
      endif
   endif
! Selection of the frozen-mass particles: end
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine Leapfrog_momentum
