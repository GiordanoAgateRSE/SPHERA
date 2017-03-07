!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
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
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: IsPointInternal  
! Description: Checking wheather a point with local normal coordinates csi(1:3)
!              is internal to a given face, whose type code is fk (=1 triangle,
!              =2 quadrilateral), or not. This procedure is based on the 
!              subroutine "LocalNormalCoordinatesGiven". They will be replaced
!              by the subroutine "point_inout_convex_non_degenerate_polygon"
!              (already available in SPHERA), which is more effective and can be
!              applied to any polygon (it already works for triangles and
!              quadrilaterals).
!-------------------------------------------------------------------------------
logical function IsPointInternal(fk,csi)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: fk
double precision,intent(in) :: csi(SPACEDIM)
integer(4) :: i
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
IsPointInternal = .false.
!------------------------
! Statements
!------------------------
if (fk==1) then            
! Triangle
   do i=1,3
     if (csi(i)<zero) return
     if (csi(i)>one) return
   enddo
   IsPointInternal = .true.
   elseif (fk==2) then
! Quadrilateral 
      do i=1,2
         if (csi(i)<zero) return
         if (csi(i)>one) return
      enddo
      IsPointInternal = .true.
endif
!------------------------
! Deallocations
!------------------------
return
end function IsPointInternal

