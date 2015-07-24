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
! Program unit: CreaGrid             
! Description: To create the background positioning grid.    
!----------------------------------------------------------------------------------------------------------------------------------

subroutine CreaGrid
!------------------------
! Modules
!------------------------ 
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: ier
double precision :: epsi
character(len=lencard) :: nomsub = "CreaGrid"
double precision, dimension(3) :: dextr
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
epsi = 0.01d0 * Domain%dd
!------------------------
! Statements
!------------------------
! To set the vertices of the background grid, based on the domain 
! dimensions and h (a tolerance epsi refers to the coordinates 
! for numerical reasons)
Grid%extr(:,1) = Domain%coord(:,1) - doubleh - epsi
Grid%extr(:,2) = Domain%coord(:,2) + doubleh + epsi
! To assess the dimensions of the grid along all the directions
dextr(:) = Grid%extr(:,2) - Grid%extr(:,1)
! To assess the number of grid cells along all the directions, assuming that 
! the cell is a cube having the side length equal to 2h 
Grid%ncd(:) = NINT(dextr(:) / doubleh)
! To assess the actual cell dimensions in all the directions, transforming 
! the cubic cell into a real hesaedric regular cell
Grid%dcd(:) = dextr(:) / Grid%ncd(:)
! In 2D, the number of cells in the y direction is forced to be 1
if (ncord==2) Grid%ncd(2) = 1
! To assess the maximum number of cells covering the grid 
! domain (a parallelepiped)
Grid%nmax = Grid%ncd(1) * Grid%ncd(2) * Grid%ncd(3)
write (nout,'(1x,a)') " "
write (nout,'(1x,a,3i8)') " Number of grid in x, y, z directions : ",          &
   Grid%ncd(1),Grid%ncd(2),Grid%ncd(3)
write (nout,'(1x,a,i10)') " Number of total grid : ",Grid%nmax
write (nout,'(1x,a)') " "
! Allocation of a 2D matrix to detect free surface (erosion criterion)
allocate (ind_interfaces(Grid%ncd(1),Grid%ncd(2),4), stat = ier)
if (ier/=0) then
   write (nout,'(1x,a,i2)')                                                    &
      "    Array ind_interfaces not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
   else
      write (nout,'(1x,a)') "    Array ind_interfaces successfully allocated "
end if
!------------------------
! Deallocations
!------------------------
return
end subroutine CreaGrid

