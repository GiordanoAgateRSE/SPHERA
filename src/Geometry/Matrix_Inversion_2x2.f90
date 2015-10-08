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
! Program unit: Matrix_Inversion_2x2  
! Description: Computation of the inverse (inv) of a provided 2x2 matrix (mat).   
!----------------------------------------------------------------------------------------------------------------------------------

subroutine Matrix_Inversion_2x2(mat,inv)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(IN) :: mat(2,2) 
double precision,intent(INOUT) :: inv(2,2)
double precision :: det
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
det = mat(1,1) * mat(2,2) - mat(2,1) * mat(1,2)
inv(1,1) = mat(2,2)
inv(1,2) = -mat(2,1)
inv(2,1) = -mat(1,2)
inv(2,2) = mat(1,1)
inv(:,:) = inv(:,:) / det
!------------------------
! Deallocations
!------------------------
return
end subroutine Matrix_Inversion_2x2

