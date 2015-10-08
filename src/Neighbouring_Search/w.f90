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
! Program unit: w              
! Description: kernel function 
!----------------------------------------------------------------------------------------------------------------------------------

double precision function w(r,h,coef)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
double precision,parameter :: a1 = 0.666666667d0
double precision :: r,h,s,q,dms,coef
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
s = r / h
if (s<=1.0d0) then
   q = a1 + s * s * (s * 0.5d0 - 1.0d0)
   elseif (s>=2.0d0) then
      q = 0.0d0
      else 
         dms = 2.0d0 - s
         q = dms * dms * dms / 6.0d0
endif
w = q * coef
!------------------------
! Deallocations
!------------------------
return
end function w

