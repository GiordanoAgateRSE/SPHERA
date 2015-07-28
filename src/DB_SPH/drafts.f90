!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-; SPHERA has been authored for RSE SpA by 
!    Andrea Amicarelli, Antonio Di Monaco, Sauro Manenti, Elia Bon, Daria Gatti, Giordano Agate, Stefano Falappi, 
!    Barbara Flamini, Roberto Guandalini, David Zuccalà).
! Main numerical developments of SPHERA: 
!    Amicarelli et al. (2015,CAF), Amicarelli et al. (2013,IJNME), Manenti et al. (2012,JHE), Di Monaco et al. (2011,EACFM). 
! Email contact: andrea.amicarelli@rse-web.it

! This file is part of SPHERA.
! SPHERA is free software: you can redistribute it and/or modify
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
! Program unit: viscomorris_wall_elements
! Description: Wall element contributions to Morris' viscosity term.              
!----------------------------------------------------------------------------------------------------------------------------------

subroutine viscomorris_wall_elements(npi,npj,npartint,dervel,rvw)
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
integer(4),intent(IN) :: npi,npj,npartint
double precision,intent(IN)  :: dervel(3)
double precision,intent(OUT) :: rvw(3)
double precision :: rhoWw,rhotilde,anuitilde,factivis,dis2
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
rhoWw = pg_w(npj)%dens * pg_w(npj)%weight * kernel_fw(1,npartint) 
rhotilde  = (pg(npi)%visc * pg(npi)%dens + pg(npi)%visc * pg_w(npj)%dens +     &
            0.001d0)
anuitilde = 4.0d0 * (pg(npi)%visc * pg(npi)%visc)    
factivis = rhoWw * anuitilde / rhotilde
dis2 = (rag_fw(1,npartint) * rag_fw(1,npartint) + rag_fw(2,npartint) *         &
       rag_fw(2,npartint) + rag_fw(3,npartint) * rag_fw(3,npartint))
rvw(:) = factivis * dervel(:) * (rag_fw(:,npartint) * pg_w(npj)%normal(:)) /   &
         dis2
!------------------------
! Deallocations
!------------------------
return
end subroutine viscomorris_wall_elements

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: viscomon_wall_elements
! Description: Wall element contributions for Monaghan artificial viscosity term.               
!----------------------------------------------------------------------------------------------------------------------------------

subroutine viscomon_wall_elements(npi,npj,npartint,dervel,rvwalfa,rvwbeta)
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
integer(4),intent(IN) :: npi,npj,npartint
double precision,intent(IN) :: dervel(3)
double precision,intent(OUT) :: rvwalfa(3), rvwbeta(3)
double precision :: rhoWw,rhotilde,celtilde,vrij,TermMon,dis2
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
vrij = - dervel(1) * rag_fw(1,npartint) - dervel(2) * rag_fw(2,npartint) -     &
       dervel(3) * rag_fw(3,npartint)
dis2 = (rag_fw(1,npartint) * rag_fw(1,npartint) + rag_fw(2,npartint) *         &
       rag_fw(2,npartint) + rag_fw(3,npartint) * rag_fw(3,npartint))
if (vrij>zero) then
   rvwalfa = zero
   rvwbeta = zero
   else
      rhotilde = pg(npi)%dens + pg_w(npj)%dens
      celtilde = 2.d0 * Med(pg(npi)%imed)%celerita 
      rhoWw = pg_w(npj)%dens * pg_w(npj)%weight * kernel_fw(1,npartint) 
      TermMon = Med(pg(npi)%imed)%alfaMon * celtilde * Domain%h / rhotilde
      rvwalfa(1:3) = rhoWw * TermMon * vrij * pg_w(npj)%normal(:) / dis2
      rvwbeta(1:3) = zero
end if
!------------------------
! Deallocations
!------------------------
return
end subroutine viscomon_wall_elements

