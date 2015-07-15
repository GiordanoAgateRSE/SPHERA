!cfile Gradients_to_MUSCL.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : Gradients_to_MUSCL
!
! Last updating : November 30, 2011
!
! Creation : Amicarelli/Agate  30 Nov 2011
!
!************************************************************************************
! Module purpose : Pseudo consistent estimation, with boundary contributions and Shepard kernel,
!                  of the velocity and density gradients for the MUSCL 
!                  reconstruction (to feed the partial Riemann solver). 
!
!AA501
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!AA601 sub
! Called subroutines: DBSPH_inlet_outlet
!
!************************************************************************************
!
subroutine Gradients_to_MUSCL
!
!.. assign modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE
!
! Implicit Declarations
implicit none
!
! Local Parameters
integer(4)       :: npi,npj,npartint,contj
double precision :: vol_Shep,Ww_Shep
!
! Executable Statements
!AA601 sub
!$omp parallel do  &
!$omp default(none) &
!$omp private(npi,contj,npj,npartint,vol_Shep) &
!$omp shared (nag,pg,NMAXPARTJ,rag,nPartIntorno,Partintorno,PartKernel,nPartIntorno_fw,DBSPH)
 do npi = 1, nag ! loop over the fluid computational particles: start
!
!AA501
    if (nPartIntorno_fw(npi) > 0) then   
!    
! Initializing the gradients 
    pg(npi)%drho = zero
    pg(npi)%dvel = zero
!
    do contj = 1, nPartIntorno(npi)   !loop over the fluid neighbouring particles: start
       npartint = (npi-1)* NMAXPARTJ + contj
       npj = PartIntorno(npartint)
!AA601
       if (pg(npi)%imed==pg(npj)%imed) then
!AA406!!!test
          vol_Shep = pg(npj)%mass / pg(npj)%dens 
!       vol_Shep = pg(npj)%mass / pg(npj)%dens / pg(npi)%sigma
! Computation of the density gradient: fluid particle contributions
          pg(npi)%drho(:) = pg(npi)%drho(:) + rag(:,npartint) * (pg(npj)%dens-pg(npi)%dens) &
                            * Partkernel(1,npartint)  * vol_Shep
! Computation of the velocity gradient: fluid particle contributions
          pg(npi)%dvel(1,:) = pg(npi)%dvel(1,:) + rag(:,npartint) * (pg(npj)%var(1) - pg(npi)%var(1)) * &
                              Partkernel(1,npartint) * vol_Shep
          pg(npi)%dvel(2,:) = pg(npi)%dvel(2,:) + rag(:,npartint) * (pg(npj)%var(2) - pg(npi)%var(2)) * &
                              Partkernel(1,npartint) * vol_Shep
          pg(npi)%dvel(3,:) = pg(npi)%dvel(3,:) + rag(:,npartint) * (pg(npj)%var(3) - pg(npi)%var(3)) * &
                              Partkernel(1,npartint) * vol_Shep
!AA601 
       endif
    end do !loop over the fluid neighbouring particles: end
!AA601 sub
    if (DBSPH%MUSCL_boundary_flag==1) call Gradients_to_MUSCL_boundary(npi)
!AA501
    endif
! 
 end do ! loop over the fluid computational particles: end
!$omp end parallel do
!
return
!
return
end subroutine Gradients_to_MUSCL
!---split


