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
! Program unit: BoundaryReflectionMatrix2D                                
! Description: Generation of the generalised reflection matrix R, based on the cosine matrix T and the parameters PsiS and PsiN.
!              (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine BoundaryReflectionMatrix2D(T,R,PsiS,PsiN)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
double precision :: PsiS,PsiN 
double precision,dimension(1:SPACEDIM,1:SPACEDIM) :: T,R
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
R(1,1) = PsiS * T(1,1) * T(1,1) + PsiN * T(3,1) * T(3,1) 
R(1,3) = (PsiS - PsiN) * T(1,1) * T(3,1)
R(3,1) = R(1,3)
R(3,3) = PsiS * T(3,1) * T(3,1) + PsiN * T(1,1) * T(1,1)
!------------------------
! Deallocations
!------------------------
return
end subroutine BoundaryReflectionMatrix2D

