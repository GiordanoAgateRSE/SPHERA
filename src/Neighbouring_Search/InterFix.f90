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
! Program unit: InterFix              
! Description: 
!-------------------------------------------------------------------------------
subroutine InterFix(npi,appo,unity)
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
! Local maximum number of particles within the kernel support 
integer(4),parameter :: local_d = 500
integer(4),intent(in) :: npi
double precision,intent(inout) :: unity
double precision,intent(inout),dimension(3) :: appo
integer(4) :: npj,contj,npartint   
double precision :: rhoj,amassj,pesoj
double precision,dimension(3) :: pesogradj
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
unity = zero
appo(:) = zero
!------------------------
! Statements
!------------------------
do contj=1,nPartIntorno(npi)
   npartint = (npi - 1) * NMAXPARTJ + contj
   npj = PartIntorno(npartint)
   if (pg(npj)%vel_type=="std") cycle    
   rhoj = pg(npj)%dens
   amassj = pg(npj)%mass
   pesoj = amassj * Partkernel(4,npartint) / rhoj
   pesogradj(1:3) = amassj * rag(1:3,npartint) * PartKernel(1,npartint) / rhoj
   unity = unity + pesoj 
   appo(:) = appo(:) + pesogradj(:)
enddo
appo(:) = -appo(:)
!------------------------
! Deallocations
!------------------------
return
end subroutine InterFix
