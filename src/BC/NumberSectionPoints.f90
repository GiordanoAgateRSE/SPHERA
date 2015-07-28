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
! Program unit: NumberSectionPoints
! Description: 
!----------------------------------------------------------------------------------------------------------------------------------

integer(4) function NumberSectionPoints (values,opt)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
implicit none
!------------------------
! Declarations
!------------------------
double precision,dimension(3,2) :: values
character(1) :: opt
integer(4) :: n
integer(4),dimension(3) :: Nmesh
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
! creazione mesh di lato dd
Nmesh = 1
!------------------------
! Statements
!------------------------
do n=1,SPACEDIM  
   if ((n==1).AND.(opt=="x")) cycle
   if ((n==2).AND.(opt=="y")) cycle
   if ((n==3).AND.(opt=="z")) cycle
   Nmesh(n) = nint((values(n,2) - values(n,1)) / Domain%dd)
end do
NumberSectionPoints = Product(Nmesh)
!------------------------
! Deallocations
!------------------------
return
end function NumberSectionPoints

