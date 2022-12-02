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
! Program unit: ParticleCellNumber              
! Description: To return the ID of the grid cell where particle np is located.
!-------------------------------------------------------------------------------
integer(4) function ParticleCellNumber(coord)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(in) :: coord(3)
integer(4) :: i,j,k
double precision :: xp,yp,zp
integer(4),external :: CellNumber
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
xp = coord(1) - Grid%extr(1,1)
yp = coord(2) - Grid%extr(2,1)
zp = coord(3) - Grid%extr(3,1)
i = ceiling(xp / Grid%dcd(1))
k = ceiling(zp / Grid%dcd(3))
#ifdef SPACE_3D
   j = ceiling(yp / Grid%dcd(2))
#elif defined SPACE_2D
      j = 1
#endif
ParticleCellNumber = CellNumber(i,j,k)
!------------------------
! Deallocations
!------------------------
return
end function ParticleCellNumber
