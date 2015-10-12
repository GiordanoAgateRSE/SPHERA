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
! Program unit: viscomorris
! Description: To compute the volume inter-particle contributions to the shear 
! viscosity term in the momentum equation (Morris, 1997, JCP). Both interactions
! "fluid particle - fluid particle" and "fluid particle - semi-particle" (DBSPH)
! are considered.   
!----------------------------------------------------------------------------------------------------------------------------------
subroutine viscomorris(npi,npj,npartint,mass_comput_part,dens_comput_part,     &
kin_visc_comput_part,mass_neighbour,dens_neighbour,kin_visc_neighbour,         &
kernel_der,vel_type,rel_dis,dervel,rvw)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: npi,npj,npartint
double precision,intent(in) :: mass_comput_part,dens_comput_part
double precision,intent(in) :: kin_visc_comput_part,mass_neighbour
double precision,intent(in) :: dens_neighbour,kin_visc_neighbour,kernel_der
double precision,intent(in) :: rel_dis(3),dervel(3)
character(3),intent(in) :: vel_type
double precision,intent(out) :: rvw(3)
double precision :: amassj,rhotilde,anuitilde,factivis
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
if (vel_type/="std") then 
   amassj = mass_comput_part
   rhotilde = dens_comput_part
   anuitilde = two * kin_visc_comput_part
   else
      amassj = mass_neighbour
      rhotilde  = (kin_visc_comput_part * dens_comput_part + kin_visc_neighbour&
                  * dens_neighbour + 0.001d0)
      anuitilde = 4.0d0 * (kin_visc_comput_part * kin_visc_neighbour)      
endif
factivis = amassj * anuitilde / rhotilde
rvw(1:3) = factivis * (-dervel(1:3) * kernel_der * dot_product(rel_dis,rel_dis))
!------------------------
! Deallocations
!------------------------
return
end subroutine viscomorris

