!-------------------------------------------------------------------------------
! SPHERA v.10.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.10.0.0
! SPHERA v.10.0.0 is free software: you can redistribute it and/or modify
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
! Program unit: cell_indices2pos
! Description: To return the position of the barycentre of the background cell 
!              whose indices are provided as input. So far, it is used only for 
!              CLC.     
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine cell_indices2pos(ix,iy,iz,pos)
!------------------------
! Modules
!------------------------
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: ix,iy,iz
double precision,intent(out) :: pos(3)
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
pos(1) = Grid%extr(1,1) + (dfloat(ix) - 0.5d0) * Grid%dcd(1)
pos(2) = Grid%extr(2,1) + (dfloat(iy) - 0.5d0) * Grid%dcd(2)
pos(3) = Grid%extr(3,1) + (dfloat(iz) - 0.5d0) * Grid%dcd(3)
!------------------------
! Deallocations
!------------------------
return
end subroutine cell_indices2pos
#endif
