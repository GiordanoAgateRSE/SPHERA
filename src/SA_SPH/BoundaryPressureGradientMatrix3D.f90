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
! Program unit: BoundaryPressureGradientMatrix3D                                
! Description: To generate the pressure gradient matrix RRP, based on the cosine matrix T and the parameter vector Psi.
!              (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine BoundaryPressureGradientMatrix3D(T,RGP,Psi)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(INOUT),dimension(SPACEDIM)          :: Psi
double precision,intent(INOUT),dimension(SPACEDIM,SPACEDIM) :: T,RGP
Integer(4) :: i
double precision,dimension(SPACEDIM,SPACEDIM) :: Diag,PsiR,TTR
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
   Diag(i,i) = Psi(i)
enddo
call MatrixTransposition(T,TTR,SPACEDIM,SPACEDIM)
call MatrixProduct(Diag,TTR,PsiR,SPACEDIM,SPACEDIM,SPACEDIM)
call MatrixProduct(T,PsiR,RGP,SPACEDIM,SPACEDIM,SPACEDIM)
!------------------------
! Deallocations
!------------------------
return
end subroutine BoundaryPressureGradientMatrix3D

