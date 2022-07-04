!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2021 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: EoS_barotropic_linear 
! Description: Barotropic linear Equation of State (EoS): pressure estimation      
!-------------------------------------------------------------------------------
subroutine EoS_barotropic_linear
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: npi
double precision :: rhorif,c2
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
if (on_going_time_step>0) then  
! Loop over all the particles
!$omp parallel do default(none)                                                & 
!$omp shared(nag,pg,Domain,Med)                                                &
!$omp private(npi,rhorif,c2)
   do npi=1,nag
! Skip the outgone particles and the particles with velocity type different 
! from "std"
      if ((pg(npi)%cella==0).or.(pg(npi)%vel_type/="std")) cycle
      rhorif = Med(pg(npi)%imed)%den0
      c2 = Med(pg(npi)%imed)%eps / rhorif
      pg(npi)%pres = c2 * (pg(npi)%dens - rhorif) + Domain%prif
   enddo
!$omp end parallel do
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine EoS_barotropic_linear
