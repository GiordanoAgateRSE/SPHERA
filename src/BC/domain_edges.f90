!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2021 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.9.0.0
! SPHERA v.9.0.0 is free software: you can redistribute it and/or modify
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
! Program unit: domain_edges
! Description: to find the domain edges (x, y and z maximum and minimum values).
!-------------------------------------------------------------------------------
subroutine domain_edges(MinOfMin,Xmin,Xmax)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
double precision,dimension(SPACEDIM),intent(inout) :: MinOfMin
double precision,dimension(SPACEDIM,NPartZone),intent(inout) :: Xmin,Xmax
integer(4) :: Nz,ii
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
! Loop over the zones in order to set the initial minimum and maximum 
! coordinates of each zone and of the numerical domain
do Nz=1,NPartZone
! To search for the maximum and minimum "partzone" coordinates
   Partz(Nz)%coordMM(1:SPACEDIM,1) = max_positive_number
   Partz(Nz)%coordMM(1:SPACEDIM,2) = max_negative_number
   Xmin(1:SPACEDIM,Nz) = max_positive_number
   Xmax(1:SPACEDIM,Nz) = max_negative_number
! To search for the minimum and maximum coordinates of the current zone 
#ifdef SPACE_3D
      call FindFrame(Xmin,Xmax,Nz)
#elif defined SPACE_2D
         call FindLine(Xmin,Xmax,Nz)
#endif
! To evaluate the minimum and maximum coordinates of the zone
   do ii=1,SPACEDIM
      if (Xmin(ii,Nz)<Partz(Nz)%coordMM(ii,1)) Partz(Nz)%coordMM(ii,1) =       &
         Xmin(ii,Nz)
      if (Xmax(ii,Nz)>Partz(Nz)%coordMM(ii,2)) Partz(Nz)%coordMM(ii,2) =       &
         Xmax(ii,Nz)
! To evaluate the overall minimum among the zones
      MinOfMin(ii) = min(MinOfMin(ii),Xmin(ii,Nz))
   enddo
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine domain_edges
