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
if (npointst<1) return
!$omp parallel do default(none)                                                &
!$omp shared(npointst,control_points,Grid)                                     &
!$omp private(i,pglocal,xp,yp,zp,ii,jj,kk,pointcellnumber)
do i=1,npointst
   pglocal%coord(1:3) = control_points(i)%coord(1:3)
   pglocal%pres = 0.d0
   pglocal%dens = 0.d0
   pglocal%vel(1:3) = 0.d0
   pglocal%sigma_fp = 0.d0
#ifdef SOLID_BODIES
   pglocal%sigma_fp_bp = 0.d0
   pglocal%sigma_fp_sbp = 0.d0
#endif
   xp = pglocal%coord(1) - Grid%extr(1,1)
   yp = pglocal%coord(2) - Grid%extr(2,1)
   zp = pglocal%coord(3) - Grid%extr(3,1)
   ii = ceiling(xp / Grid%dcd(1))
   jj = ceiling(yp / Grid%dcd(2))
   kk = ceiling(zp / Grid%dcd(3)) 
   pointcellnumber = CellNumber(ii,jj,kk)
   pglocal%cella = pointcellnumber
   pglocal%B_ren_fp(1:3,1:3) = 0.d0
#ifdef SOLID_BODIES
   pglocal%B_ren_fp_bp(1:3,1:3) = 0.d0
   pglocal%B_ren_fp_sbp(1:3,1:3) = 0.d0
#endif
   call SPH_approximations_at_monitors(pglocal)
   control_points(i)%pres = pglocal%pres
   control_points(i)%dens = pglocal%dens
   control_points(i)%vel(1:3) = pglocal%vel(1:3)
   control_points(i)%sigma_fp = pglocal%sigma_fp
#ifdef SOLID_BODIES
   control_points(i)%sigma_fp_bp = pglocal%sigma_fp_bp
   control_points(i)%sigma_fp_sbp = pglocal%sigma_fp_sbp
#endif
   control_points(i)%B_ren_fp(1:3,1:3) = pglocal%B_ren_fp(1:3,1:3)
#ifdef SOLID_BODIES
   control_points(i)%B_ren_fp_bp(1:3,1:3) = pglocal%B_ren_fp_bp(1:3,1:3)
   control_points(i)%B_ren_fp_sbp(1:3,1:3) = pglocal%B_ren_fp_sbp(1:3,1:3)
#endif
   control_points(i)%cella  = pglocal%cella
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine CalcVarp
