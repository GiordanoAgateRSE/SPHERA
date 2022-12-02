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
! Program unit: CellIndices            
! Description: To return the indices (i,j,k) of the cell (nc) in a 3D domain
!              with ni*nj*nk cells and the cell index in the 2D horizontal grid.
!-------------------------------------------------------------------------------
integer(4) function CellIndices(nc,i,j,k)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: nc
integer(4),intent(out) :: i,j,k
integer(4) :: ncij,nucellsij,ni,nj 
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
ni = Grid%ncd(1)
nj = Grid%ncd(2)
nucellsij = ni * nj
! Grid index in the z-direction
k = int((nc - 1) / nucellsij) + 1
ncij = nc - nucellsij * (k - 1)
! Grid index in the y-direction
j = int((ncij - 1) / ni) + 1
! Grid index in the x-direction
i = ncij - ni * (j - 1)
! Cell index in the horizontal projection of the background grid (2D 
! horizontal grid)
CellIndices = ncij
!------------------------
! Deallocations
!------------------------
end function CellIndices
