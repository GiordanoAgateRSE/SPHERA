!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2020 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: velocity_smoothing_2
! Description: Partial smoothing of the fluid velocity field (part 2)
!-------------------------------------------------------------------------------
subroutine velocity_smoothing_2
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
double precision :: TetaV1
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
!$omp shared(pg,Med,Domain,dt,indarrayFlu,Array_Flu,input_any_t)               &
!$omp private(npi,ii,TetaV1)
! Loop over all the active particles
do ii=1,indarrayFlu
   npi = Array_Flu(ii)
   if (input_any_t%TetaV>1.d-9) then
! TetaV depending on the time step
      TetaV1 = input_any_t%TetaV * Med(pg(npi)%imed)%Celerita * dt / Domain%h
      if (pg(npi)%kodvel==0) then              
! The particle is inside the domain and far from boundaries
         pg(npi)%var(:) = pg(npi)%vel(:) + TetaV1 * pg(npi)%var(:)     
         pg(npi)%vel(:) = pg(npi)%var(:)                                        
         else
! The particle is close to a "source", "level" or "normal velocity boundary
! (kodvel=1 or kodvel=2): the final velocity is kept unmodified
            pg(npi)%var(:) = pg(npi)%vel(:)                                
      endif
      else
         pg(npi)%var(:) = pg(npi)%vel(:)
   endif
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine velocity_smoothing_2
