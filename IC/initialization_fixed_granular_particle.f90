!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: initialization_fixed_granular_particle
!
! Versions: 
! 01 Amicarelli 08Apr14 (creation)
!
!************************************************************************************
! Module purpose : (v5.04) To initialize the most of the fixed SPH granular particles 
!
! Calling routines: Shields,SetParticleParameters,Loop_Irre_2D,Loop_Irre_3D
!
! Called subroutines: / 
!
!************************************************************************************

 subroutine initialization_fixed_granular_particle(npi)

! Assigning modules
 use GLOBAL_MODULE 
 use AdM_USER_TYPE
 use ALLOC_MODULE

! Declarations
 implicit none
 integer(4),intent(in)  :: npi
 
! Statements
 pg(npi)%Beta_slope = -999.d0
 pg(npi)%Gamma_slope = -999.d0 
 pg(npi)%u_star = 0.d0
 pg(npi)%C_L = 0.d0
 pg(npi)%C_D = 0.d0
 pg(npi)%k_BetaGamma = -999.d0
 pg(npi)%tau_tauc = 0.0d0
 
 return
 end subroutine initialization_fixed_granular_particle
!---split

