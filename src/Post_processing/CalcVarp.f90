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
! Program unit: CalcVarp              
! Description: To calculate physical quantities at a monitoring point.     
!-------------------------------------------------------------------------------
subroutine CalcVarp 
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
integer(4) :: i,ii,jj,kk,pointcellnumber
double precision xp,yp,zp
type (TyCtlPoint) :: pglocal
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
if (Npointst<1) return
do i=1, Npointst
   pglocal%coord(:) = control_points(i)%coord(:)
   pglocal%pres = zero
   pglocal%dens = zero
   pglocal%vel(:) = zero
   pglocal%uni = zero
   xp = pglocal%coord(1) - Grid%extr(1,1)
   yp = pglocal%coord(2) - Grid%extr(2,1)
   zp = pglocal%coord(3) - Grid%extr(3,1)
   ii = ceiling(xp / Grid%dcd(1))
   jj = ceiling(yp / Grid%dcd(2))
   kk = ceiling(zp / Grid%dcd(3)) 
   pointcellnumber = CellNumber(ii, jj, kk)
   pglocal%cella = pointcellnumber
   call GetVarPart (pglocal)
   control_points(i)%pres = pglocal%pres
   control_points(i)%dens = pglocal%dens
   control_points(i)%vel(:) = pglocal%vel(:)
   control_points(i)%uni = pglocal%uni
   control_points(i)%cella  = pglocal%cella
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine CalcVarp

