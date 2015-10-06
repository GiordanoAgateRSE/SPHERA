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
! Program unit: Matrix_Inversion_3x3  
! Description: Computation of the inverse (inv) of a provided 3x3 matrix (mat).   
!----------------------------------------------------------------------------------------------------------------------------------

subroutine Matrix_Inversion_3x3(mat,inv,test)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(IN) :: mat(3,3) 
double precision,intent(INOUT) :: inv(3,3)
integer(4),intent(inout) :: test  
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
det = mat(1,1) * mat(2,2) * mat(3,3) + mat(2,1) * mat(3,2) * mat(1,3) +        &
      mat(3,1) * mat(1,2) * mat(2,3) - mat(1,1) * mat(3,2) * mat(2,3) -        &
      mat(3,1) * mat(2,2) * mat(1,3) - mat(2,1) * mat(1,2) * mat(3,3)
if (det/=0.d0) then
   test = 1
   inv(1,1) = mat(2,2) * mat(3,3) - mat(2,3) * mat(3,2)
   inv(1,2) = mat(1,3) * mat(3,2) - mat(1,2) * mat(3,3)
   inv(1,3) = mat(1,2) * mat(2,3) - mat(1,3) * mat(2,2)   
   inv(2,1) = mat(2,3) * mat(3,1) - mat(2,1) * mat(3,3)
   inv(2,2) = mat(1,1) * mat(3,3) - mat(1,3) * mat(3,1)
   inv(2,3) = mat(1,3) * mat(2,1) - mat(1,1) * mat(2,3)
   inv(3,1) = mat(2,1) * mat(3,2) - mat(2,2) * mat(3,1)
   inv(3,2) = mat(1,2) * mat(3,1) - mat(1,1) * mat(3,2)
   inv(3,3) = mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1)
   inv(:,:) = inv(:,:) / det
   else
      test = 0
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine Matrix_Inversion_3x3

