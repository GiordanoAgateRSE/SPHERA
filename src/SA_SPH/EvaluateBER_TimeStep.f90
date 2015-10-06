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
! Program unit: EvaluateBER_TimeStep                                      
! Description:    
!              (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine EvaluateBER_TimeStep
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
integer(4) :: imed
double precision :: DTlocal
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
DTminBER = 1.0d30
!------------------------
! Statements
!------------------------
do imed=1,NMedium
   DTlocal = Domain%h / Med(imed)%celerita
   if (DTlocal<DTminBER) DTminBER = DTlocal
enddo
DTminBER = sqrttwo * DTminBER
!------------------------
! Deallocations
!------------------------
return
end subroutine EvaluateBER_TimeStep

