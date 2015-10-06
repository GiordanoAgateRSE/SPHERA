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
! Program unit: BoundaryMassForceMatrix3D                                
! Description: Generation of the generalised boundary mass force matrix RN, on the base of the cosine matrix 
!              T and the parameter Fi.
!              (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine BoundaryMassForceMatrix3D(T,RMF,Fi)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(INOUT),dimension(SPACEDIM) :: Fi
double precision,intent(INOUT),dimension(SPACEDIM,SPACEDIM) :: T,RMF
integer(4) :: i
double precision,dimension(SPACEDIM,SPACEDIM) :: Diag,FiR,TTR
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
Diag = zero
!------------------------
! Statements
!------------------------
do i=1,SPACEDIM
   Diag(i,i) = Fi(i)
enddo
call MatrixTransposition(T,TTR,SPACEDIM,SPACEDIM)
call MatrixProduct(Diag,TTR,FiR,SPACEDIM,SPACEDIM,SPACEDIM)
call MatrixProduct(T,FiR,RMF,SPACEDIM,SPACEDIM,SPACEDIM)
!------------------------
! Deallocations
!------------------------
return
end subroutine BoundaryMassForceMatrix3D

