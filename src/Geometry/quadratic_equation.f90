!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-; SPHERA has been authored for RSE SpA by 
!    Andrea Amicarelli, Antonio Di Monaco, Sauro Manenti, Elia Bon, Daria Gatti, Giordano Agate, Stefano Falappi, 
!    Barbara Flamini, Roberto Guandalini, David Zuccalà).
! Main numerical developments of SPHERA: 
!    Amicarelli et al. (2015,CAF), Amicarelli et al. (2013,IJNME), Manenti et al. (2012,JHE), Di Monaco et al. (2011,EACFM). 
! Email contact: andrea.amicarelli@rse-web.it

! This file is part of SPHERA.
! SPHERA is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with SPHERA. If not, see <http://www.gnu.org/licenses/>.
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

