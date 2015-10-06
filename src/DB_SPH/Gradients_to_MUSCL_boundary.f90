!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-) 
!      
!     
!   
!      
!  

! This file is part of SPHERA.
!  
!  
!  
!  
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
!  
!  
!  
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: Gradients_to_MUSCL_boundary
! Description: Estimation of the boundary terms for the MUSCL reconstruction scheme (DB-SPH), in case they are required in input.                 
!----------------------------------------------------------------------------------------------------------------------------------

subroutine Gradients_to_MUSCL_boundary(npi)
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
integer(4),intent(in) :: npi
integer(4) :: contj,npartint,npj
double precision :: vol_Shep,Ww_Shep
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
! Loop over the wall neighbouring elements 
! (parallelization makes no sense for this loop)
do contj=1,nPartIntorno_fw(npi)   
   npartint = (npi - 1) * NMAXPARTJ + contj
   npj = PartIntorno_fw(npartint)
! Boundary term contributions for the density and velocity gradients
   Ww_Shep = pg_w(npj)%weight * kernel_fw(1,npartint) / pg(npi)%sigma
! Computation of the density gradient: wall contributions
   pg(npi)%drho(:) = pg(npi)%drho(:) +  (pg_w(npj)%dens-pg(npi)%dens) *        &
                     pg_w(npj)%normal(:) * Ww_Shep
! Computation of the velocity gradient: wall contributions
   pg(npi)%dvel(1,:) = pg(npi)%dvel(1,:) + (pg_w(npj)%vel(1) - pg(npi)%var(1)) &
                       * pg_w(npj)%normal(:) * Ww_Shep
   pg(npi)%dvel(2,:) = pg(npi)%dvel(2,:) + (pg_w(npj)%vel(2) - pg(npi)%var(2)) &
                       * pg_w(npj)%normal(:) * Ww_Shep
   pg(npi)%dvel(3,:) = pg(npi)%dvel(3,:) + (pg_w(npj)%vel(3) - pg(npi)%var(3)) & 
                       * pg_w(npj)%normal(:) * Ww_Shep
! Semi-particles
   vol_Shep = pg_w(npj)%mass / pg_w(npj)%dens / pg(npi)%sigma
! Computation of the density gradient: wall contributions
   pg(npi)%drho(:) = pg(npi)%drho(:) + rag_fw(:,npartint) * (pg_w(npj)%dens -  &
                     pg(npi)%dens) * kernel_fw(2,npartint) * vol_Shep
! Computation of the velocity gradient: wall contributions
   pg(npi)%dvel(1,:) = pg(npi)%dvel(1,:) + rag_fw(:,npartint) *                &
                       (pg_w(npj)%vel(1) - pg(npi)%var(1)) *                   &
                       kernel_fw(2,npartint) * vol_Shep
   pg(npi)%dvel(2,:) = pg(npi)%dvel(2,:) + rag_fw(:,npartint) *                &
                       (pg_w(npj)%vel(2) - pg(npi)%var(2)) *                   &
                       kernel_fw(2,npartint) * vol_Shep
   pg(npi)%dvel(3,:) = pg(npi)%dvel(3,:) + rag_fw(:,npartint) *                &
                       (pg_w(npj)%vel(3) - pg(npi)%var(3)) *                   &
                       kernel_fw(2,npartint) * vol_Shep
end do   
!------------------------
! Deallocations
!------------------------
return
end subroutine Gradients_to_MUSCL_boundary

