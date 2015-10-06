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
! Program unit: initialization_fixed_granular_particle     
! Description: To initialize the most of the fixed SPH mixture particles (bed-load transport).              
!----------------------------------------------------------------------------------------------------------------------------------

subroutine initialization_fixed_granular_particle(npi)
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
integer(4),intent(in) :: npi
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
pg(npi)%Beta_slope = -999.d0
pg(npi)%Gamma_slope = -999.d0 
pg(npi)%u_star = 0.d0
pg(npi)%C_L = 0.d0
pg(npi)%C_D = 0.d0
pg(npi)%k_BetaGamma = -999.d0
pg(npi)%tau_tauc = 0.0d0
!------------------------
! Statements
!------------------------
!------------------------
! Deallocations
!------------------------
return
end subroutine initialization_fixed_granular_particle

