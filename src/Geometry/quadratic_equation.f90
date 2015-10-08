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
! Program unit: quadratic_equation   
! Description: To solve a quadratic equation.       
!----------------------------------------------------------------------------------------------------------------------------------

subroutine quadratic_equation(a,b,c,n_roots,root1,root2)
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
double precision, intent(in) :: a,b,c
integer(4),intent(out) :: n_roots
double precision, intent(out) :: root1,root2
double precision :: discriminant
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
root1=-999.d0
root2=-999.d0
!------------------------
! Statements
!------------------------
if (a/=0.d0) then
   discriminant = b ** 2 - 4.0d0 * a * c
   if (discriminant>0.0d0) then
      n_roots = 2
      root1 = ( - b - dsqrt(discriminant)) / (2.0d0 * a)
      root2 = ( - b + dsqrt(discriminant)) / (2.0d0 * a)
      else
      if (discriminant==0.0d0) then
         n_roots = 1
         root1 = - b / (2.0d0 * a)
         else
            n_roots = 0
      endif      
   endif
   else
      if (b/=0.d0) then
         n_roots = 1
         root1 = - c / b
      endif   
endif       
!------------------------
! Deallocations
!------------------------
return
end subroutine quadratic_equation

