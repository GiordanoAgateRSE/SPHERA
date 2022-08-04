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
! Program unit: KTGF_update
! Description: To update some quantities associated with dense granular flows 
!              (Kinetic Theory of Granular Flow under the packing limit, 
!              Amicarelli et al., 2017, IJCFD)
!-------------------------------------------------------------------------------
subroutine KTGF_update
!------------------------
! Modules
!------------------------
use Hybrid_allocation_module
use Static_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: npi,ncel,aux,igridi,jgridi,kgridi
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
if (Granular_flows_options%KTGF_config==1) then
!$omp parallel do default(none)                                                &
!$omp shared(pg,nag)                                                           &
!$omp private(npi,ncel)
   do npi=1,nag
      pg(npi)%vel_old(:) = pg(npi)%vel(:)
      pg(npi)%normal_int_old(:) = pg(npi)%normal_int(:)
      call initialization_fixed_granular_particle(npi)             
   enddo
!$omp end parallel do
endif
!$omp parallel do default(none)                                                &
!$omp shared(pg,nag)                                                           &
!$omp private(npi) 
do npi=1,nag
   call Shields(npi)
enddo
!$omp end parallel do
if (Granular_flows_options%KTGF_config==1) then
! Initializing viscosity for fixed particles
!$omp parallel do default(none)                                                &
!$omp shared(pg,nag,Granular_flows_options,Med)                                &
!$omp private(npi,ncel,aux,igridi,jgridi,kgridi)
   do npi=1,nag
      ncel = ParticleCellNumber(pg(npi)%coord)
      aux = CellIndices(ncel,igridi,jgridi,kgridi)
      if (pg(npi)%state=="sol") then
         pg(npi)%mu = Med(pg(npi)%imed)%mumx
         pg(npi)%kin_visc = pg(npi)%mu / pg(npi)%dens
      endif
   enddo
!$omp end parallel do
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine KTGF_update
