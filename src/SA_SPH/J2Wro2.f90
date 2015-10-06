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
! Program unit: J2Wro2                                         
! Description: 
!               (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

double precision function J2Wro2(ro)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
double precision,parameter :: a1 = 0.00833333333333333d0 != 1 / 120
double precision,parameter :: a2 = 0.26666666666666667d0 != 8 / 30
double precision,intent(IN) :: ro
double precision :: ro2,ro3
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
ro2 = ro * ro
ro3 = ro2 * ro
!------------------------
! Statements
!------------------------
if ((zero<=ro).and.(ro<one)) then
   J2Wro2 = KERNELCONST2D * (0.25d0 - (a1 * (40.0d0 - 36.0d0 * ro2 + 15.0d0 *  &
      ro3) * ro3))
   elseif ((one<=ro).and.(ro<two)) then
      J2Wro2 = KERNELCONST2D * (a2 - (a1 * (80.0d0 - 90.0d0 * ro + 36.0d0 *    &
         ro2 - 5.0d0 * ro3) * ro3))
      else
         J2Wro2 = zero
endif
!------------------------
! Deallocations
!------------------------
return
end function J2Wro2

