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
! Program unit: NumberSectionPoints
! Description: 
!----------------------------------------------------------------------------------------------------------------------------------

integer(4) function NumberSectionPoints (values,opt)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
implicit none
!------------------------
! Declarations
!------------------------
double precision,dimension(3,2) :: values
character(1) :: opt
integer(4) :: n
integer(4),dimension(3) :: Nmesh
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
! creazione mesh di lato dd
Nmesh = 1
!------------------------
! Statements
!------------------------
do n=1,SPACEDIM  
   if ((n==1).AND.(opt=="x")) cycle
   if ((n==2).AND.(opt=="y")) cycle
   if ((n==3).AND.(opt=="z")) cycle
   Nmesh(n) = nint((values(n,2) - values(n,1)) / Domain%dd)
end do
NumberSectionPoints = Product(Nmesh)
!------------------------
! Deallocations
!------------------------
return
end function NumberSectionPoints

