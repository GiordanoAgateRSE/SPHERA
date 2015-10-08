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
! Program unit: InterFix              
! Description:      
!----------------------------------------------------------------------------------------------------------------------------------

subroutine InterFix(npi,appo,unity)
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
integer(4),parameter :: local_d = 500 ! Local maximum number of particles 
                                      ! within the kernel support 
integer(4),intent(IN) :: npi
double precision,intent(INOUT) :: unity
double precision,intent(INOUT),dimension(3) :: appo
integer(4) :: npj,contj,npartint   
double precision :: rhoj,amassj,pesoj
double precision,dimension(3) :: pesogradj
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
unity = zero
appo(:) = zero
!------------------------
! Statements
!------------------------
do contj=1,nPartIntorno(npi)
   npartint = (npi - 1)* NMAXPARTJ + contj
   npj = PartIntorno(npartint)
   if ( pg(npj)%vel_type=="std") cycle    
   rhoj = pg(npj)%dens
   amassj = pg(npj)%mass
   pesoj = amassj * Partkernel(4,npartint) / rhoj
   pesogradj(1:3) = amassj * rag(1:3,npartint) * PartKernel(1,npartint) / rhoj
   unity = unity + pesoj  
   appo(:) = appo(:) + pesogradj(:)  
enddo
appo(:) = -appo(:)
!------------------------
! Deallocations
!------------------------
return
end subroutine InterFix 

