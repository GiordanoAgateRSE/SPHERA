!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: InterpolateBoundaryIntegrals2D                                       
! Description:  Interpolation in table "BoundIntegralTab(:,:)", defined in
!               module "SA_SPH_module", the values in columns "Colmn(nc), nc=1,
!               Ncols" corresponding to the input value "x" to be interpolated,
!               in turn, in column 0. It returns:
!                  Func(nc), nc=1, Ncols : values interpolated in columns 
!                                          Col(nc), nc=1, Ncols
!               (Di Monaco et al., 2011, EACFM)
!-------------------------------------------------------------------------------
#ifdef SPACE_2D
subroutine InterpolateBoundaryIntegrals2D(Ncols,Colmn,xx,Func)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use SA_SPH_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: Ncols
integer(4),dimension(1:NUMCOLS_BIT),intent(in) :: Colmn
double precision,intent(inout) :: xx
double precision,dimension(1:NUMCOLS_BIT),intent(out) :: Func
integer(4) :: nc,ii,jj
double precision :: xi,fi,fip1
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
if (xx<BoundIntegralTab2D(1,0)) xx = BoundIntegralTab2D(1,0)
ii = int((xx - BoundIntegralTab2D(1,0)) / DELTAX_BIT) + 1
xi = BoundIntegralTab2D(ii,0)
do nc=1,Ncols
   jj = Colmn(nc)
   fi = BoundIntegralTab2D(ii,jj)
   fip1 = BoundIntegralTab2D(ii+1,jj)
   Func(nc) = fi + (fip1 - fi) * (xx - xi) / DELTAX_BIT
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine InterpolateBoundaryIntegrals2D
#endif
