!AA601 sub the whole subroutine
!cfile BC_wall_elements.f90
!************************************************************************************
!                          S P H E R A 5.0.3 
!
!                Smoothed Particle Hydrodinamics Code
!
!************************************************************************************
!
! File name     : Gradients_to_MUSCL_boundary
!
! Creation      : Amicarelli A., 30Nov11 (hard-coding)
! Updates       : Amicarelli A., 26Jan15 (integration with no hard-coding)
!
!************************************************************************************
! Module purpose : Boundary terms for the MUSCL reconstruction (DBSPH), if required in input
!
! Calling subroutines: Gradients_to_MUSCL
! Called subroutines:  /  
!
!************************************************************************************

subroutine Gradients_to_MUSCL_boundary(npi)

! Modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE

! Declarations ..
implicit none
integer(4), intent(in) :: npi
integer(4) :: contj,npartint,npj
double precision :: vol_Shep,Ww_Shep

! Statements
! Loop over the wall neighbouring elements: start (parallelization has no sense for this loop)
do contj=1,nPartIntorno_fw(npi)   
   npartint = (npi-1)* NMAXPARTJ + contj
   npj = PartIntorno_fw(npartint)
! Boundary term contributions for the density and velocity gradients
   Ww_Shep = pg_w(npj)%weight * kernel_fw(1,npartint) / pg(npi)%sigma
! Computation of the density gradient: wall contributions
   pg(npi)%drho(:) = pg(npi)%drho(:) +  (pg_w(npj)%dens-pg(npi)%dens) * pg_w(npj)%normal(:) * Ww_Shep
! Computation of the velocity gradient: wall contributions
   pg(npi)%dvel(1,:) = pg(npi)%dvel(1,:) + (pg_w(npj)%vel(1) - pg(npi)%var(1)) * pg_w(npj)%normal(:) * Ww_Shep
   pg(npi)%dvel(2,:) = pg(npi)%dvel(2,:) + (pg_w(npj)%vel(2) - pg(npi)%var(2)) * pg_w(npj)%normal(:) * Ww_Shep
   pg(npi)%dvel(3,:) = pg(npi)%dvel(3,:) + (pg_w(npj)%vel(3) - pg(npi)%var(3)) * pg_w(npj)%normal(:) * Ww_Shep
! Semi-particles
   vol_Shep = pg_w(npj)%mass / pg_w(npj)%dens / pg(npi)%sigma
! Computation of the density gradient: wall contributions
   pg(npi)%drho(:) = pg(npi)%drho(:) + rag_fw(:,npartint) * (pg_w(npj)%dens-pg(npi)%dens) &
                     * kernel_fw(2,npartint)  * vol_Shep
! Computation of the velocity gradient: wall contributions
   pg(npi)%dvel(1,:) = pg(npi)%dvel(1,:) + rag_fw(:,npartint) * (pg_w(npj)%vel(1) - pg(npi)%var(1)) * &
                       kernel_fw(2,npartint) * vol_Shep
   pg(npi)%dvel(2,:) = pg(npi)%dvel(2,:) + rag_fw(:,npartint) * (pg_w(npj)%vel(2) - pg(npi)%var(2)) * &
                       kernel_fw(2,npartint) * vol_Shep
   pg(npi)%dvel(3,:) = pg(npi)%dvel(3,:) + rag_fw(:,npartint) * (pg_w(npj)%vel(3) - pg(npi)%var(3)) * &
                       kernel_fw(2,npartint) * vol_Shep
end do   
! Loop over the wall neighbouring elements: end

return
end subroutine Gradients_to_MUSCL_boundary
!---split

