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
! Program unit: IsPointInternal  
! Description: Checking wheather a point with local normal coordinates csi() is internal to a given face, whose code  
!              is fk (=1 triangle, =2 parallelogram).
!----------------------------------------------------------------------------------------------------------------------------------

Logical Function IsPointInternal(fk,csi)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: i,fk
double precision, dimension(1:SPACEDIM) :: csi
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
IsPointInternal = .FALSE.
!------------------------
! Statements
!------------------------
if (fk==1) then            
! Triangle
   do i=1,3
     if (csi(i)<zero) return
     if (csi(i)>one) return
   end do
   IsPointInternal = .TRUE.
   else if (fk==2) then
! Quadrilateral 
      do i=1,2
         if (csi(i)<zero) return
         if (csi(i)>one) return
      end do
      IsPointInternal = .TRUE.
end if
!------------------------
! Deallocations
!------------------------
return
end Function IsPointInternal

