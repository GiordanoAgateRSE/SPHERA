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
! Program unit: Leapfrog_trajectories
! Description: Leapfrog time integration scheme (fluid particle trajectories)
!-------------------------------------------------------------------------------
subroutine Leapfrog_trajectories
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: npi
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
!$omp shared(nag,pg,dt,input_any_t)                                            &
!$omp private(npi)
! Loop over the active particles
do npi=1,nag
   if (pg(npi)%cella==0) cycle
! To save the old coordinates
   pg(npi)%CoordOld(:) = pg(npi)%coord(:)
   if (pg(npi)%vel_type/="std") then
! If the motion type is not "std", velocities are assigned
      pg(npi)%coord(:) = pg(npi)%coord(:) + dt * pg(npi)%vstart(:)
      else
! Otherwise, the partial smoothed velocity field is integrated in time
         pg(npi)%coord(:) = pg(npi)%coord(:) + dt * pg(npi)%vel(:)  
   endif
! ALE velocity increment (it was already integrated in the acceleration)
   if (input_any_t%ALE3) then
      pg(npi)%dvel_ALE1(1:3) = pg(npi)%dvel_ALE1(1:3) * dt
      pg(npi)%dvel_ALE3(1:3) = pg(npi)%dvel_ALE3(1:3) * dt
   endif
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine Leapfrog_trajectories
