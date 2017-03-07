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
! Program unit: CellNumber             
! Description: To return the ID of the cell of indices (i,j,k).  
!-------------------------------------------------------------------------------
integer(4) function CellNumber(i,j,k)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(IN) :: i,j,k
integer(4) :: ni,nj,nk
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
ni = Grid%ncd(1)
nj = Grid%ncd(2)
nk = Grid%ncd(3)
!------------------------
! Statements
!------------------------
if ((i<1).or.(i>ni).or.(j<1).or.(j>nj).or.(k<1).or.(k>nk)) then
! .. the cell is outside the grid limits
   CellNumber = 0  
   else
! .. return the cell number
      CellNumber = ((k - 1) * nj + (j - 1)) * ni + i
end if
!------------------------
! Deallocations
!------------------------
end function CellNumber

