!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-Ricerca di Sistema)



! SPHERA authors and email contact are provided on SPHERA documentation.

! This file is part of SPHERA v.8.0.
! SPHERA v.8.0 is free software: you can redistribute it and/or modify
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

