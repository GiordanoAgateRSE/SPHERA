!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-Ricerca di Sistema)



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
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: viscomon 
! Description: Monaghan (2005) artificial viscosity term. It is also active for separating particles. Volume viscosity term 
!              is neglected in the momentum equation.     
!----------------------------------------------------------------------------------------------------------------------------------

subroutine viscomon (npi,npj,npartint,dervel,rvwalfa,rvwbeta)
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
integer(4), intent(IN)  :: npi,npj,npartint
double precision,intent(IN) :: dervel(3)
double precision,intent(OUT) :: rvwalfa(3),rvwbeta(3)
double precision :: celtilde,rhotilde,amassj,vrij,TermMon
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
vrij = - dervel(1) * rag(1,npartint) - dervel(2) * rag(2,npartint) - dervel(3) &
       * rag(3,npartint)
! AA504 removed part: Monaghan's term is always active even if particle are 
! separating 
if (pg(npj)%vel_type/="std") then
   amassj   = pg(npi)%mass
   rhotilde = two * pg(npi)%dens
   celtilde = Med(pg(npi)%imed)%celerita + Med(pg(npi)%imed)%celerita
   else
      amassj = pg(npj)%mass
      rhotilde = pg(npi)%dens + pg(npj)%dens
      if (esplosione) then  
         celtilde = pg(npi)%Csound + pg(npj)%Csound
         else
         celtilde = Med(pg(npi)%imed)%celerita + Med(pg(npj)%imed)%celerita
      end if
end if
TermMon = Med(pg(npi)%imed)%alfaMon * celtilde * Domain%h / rhotilde
! Volume viscosity term (compressible flows)
! can be neglected and may causes several problems when activated
rvwalfa(1:3) = amassj * TermMon * vrij * rag(1:3,npartint) *                   &
               PartKernel(2,npartint)
! beta_Monaghan is not validated
rvwbeta(1:3) = zero
!------------------------
! Deallocations
!------------------------
return
end subroutine viscomon

