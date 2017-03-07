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
! Program unit: diffumorris
! Description:  
!-------------------------------------------------------------------------------
subroutine diffumorris (npi,npj,npartint,dervol,factdiff,rvw)
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
integer(4) :: npi,npj,npartint
double precision :: dervol,factdiff,rvw,anuitilde,rhotilde,amassj,coei,coej
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
pg(npj)%coefdif = med(pg(npj)%imed)%codif
pg(npi)%coefdif = med(pg(npi)%imed)%codif
coej = pg(npj)%coefdif
coei = pg(npi)%coefdif
amassj = pg(npj)%mass
anuitilde = 4.0d0 * (coei * coej) / (coei + coej + 0.00001d0)
rhotilde  = pg(npj)%dens
if ((pg(npi)%visc==med(2)%numx).or.(pg(npj)%visc==med(2)%numx)) then
   anuitilde = zero
end if
if ((pg(npj)%vel_type/="std")) then 
   rhotilde = pg(npi)%dens
   amassj = pg(npi)%mass
   anuitilde = zero
end if
factdiff = amassj * anuitilde / rhotilde
rvw = - dervol * PartKernel(2,npartint) * (rag(1,npartint) * rag(1,npartint)   &
      + rag(2,npartint) * rag(2,npartint) + rag(3,npartint) * rag(3,npartint)) 
!------------------------
! Deallocations
!------------------------
return
end subroutine diffumorris

