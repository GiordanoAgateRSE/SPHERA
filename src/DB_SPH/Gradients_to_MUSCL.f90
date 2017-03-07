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
! Program unit: Gradients_to_MUSCL
! Description: 0th-order consistency estimation of velocity and density 
!              gradients for the MUSCL reconstruction (to feed the Partial 
!              Linearized Riemann Solver; Amicarelli et al., 2013, IJNME).                
!-------------------------------------------------------------------------------
subroutine Gradients_to_MUSCL
!------------------------
! Modules
!------------------------ 
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: npi,npj,npartint,contj
double precision :: volume,Ww_Shep
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
!$omp private(npi,contj,npj,npartint,volume)                                   &
!$omp shared(nag,pg,NMAXPARTJ,rag,nPartIntorno,Partintorno,PartKernel)         &
!$omp shared(nPartIntorno_fw,DBSPH,NMedium)
! Loop over the fluid computational particles
do npi=1,nag 
   if (nPartIntorno_fw(npi)>0) then   
! Initializing the gradients 
      pg(npi)%drho = zero
      pg(npi)%dvel = zero
! Loop over the fluid neighbouring particles
      do contj=1,nPartIntorno(npi)   
         npartint = (npi - 1) * NMAXPARTJ + contj
         npj = PartIntorno(npartint)
         if (pg(npi)%imed==pg(npj)%imed) then
            volume = pg(npj)%mass / pg(npj)%dens 
! Computation of the density gradient: fluid particle contributions
            pg(npi)%drho(:) = pg(npi)%drho(:) + rag(:,npartint) * (pg(npj)%dens&
                              - pg(npi)%dens) * Partkernel(1,npartint)  *      &
                              volume
! Computation of the velocity gradient: fluid particle contributions
            pg(npi)%dvel(1,:) = pg(npi)%dvel(1,:) + rag(:,npartint) *          &
                                (pg(npj)%var(1) - pg(npi)%var(1)) *            &
                                Partkernel(1,npartint) * volume
            pg(npi)%dvel(2,:) = pg(npi)%dvel(2,:) + rag(:,npartint) *          &
                                (pg(npj)%var(2) - pg(npi)%var(2)) *            &
                                Partkernel(1,npartint) * volume
            pg(npi)%dvel(3,:) = pg(npi)%dvel(3,:) + rag(:,npartint) *          &
                                (pg(npj)%var(3) - pg(npi)%var(3)) *            &
                                Partkernel(1,npartint) * volume
         endif
      enddo
      if ((DBSPH%MUSCL_boundary_flag.eqv.(.true.)).and.(NMedium==1)) then
         call Gradients_to_MUSCL_boundary(npi)
      endif
   endif
enddo 
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine Gradients_to_MUSCL

