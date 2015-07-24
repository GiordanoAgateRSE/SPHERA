!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-; SPHERA has been authored for RSE SpA by 
!    Andrea Amicarelli, Antonio Di Monaco, Sauro Manenti, Elia Bon, Daria Gatti, Giordano Agate, Stefano Falappi, 
!    Barbara Flamini, Roberto Guandalini, David Zuccalà).
! Main numerical developments of SPHERA: 
!    Amicarelli et al. (2015,CAF), Amicarelli et al. (2013,IJNME), Manenti et al. (2012,JHE), Di Monaco et al. (2011,EACFM). 
! Email contact: andrea.amicarelli@rse-web.it

! This file is part of SPHERA.
! SPHERA is free software: you can redistribute it and/or modify
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
! Program unit: ParticleCellNumber              
! Description: To return the ID of the grid cell where particle np is located. If particle is outside of the grid, it returns -1 .      
!----------------------------------------------------------------------------------------------------------------------------------

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
integer(4) :: i, j, k
double precision :: xp, yp, zp
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
if (ncord==3) then
   j = ceiling(yp / Grid%dcd(2))
   else
      j = 1
end if
ParticleCellNumber = CellNumber(i, j, k)
!------------------------
! Deallocations
!------------------------
return
end Function ParticleCellNumber

