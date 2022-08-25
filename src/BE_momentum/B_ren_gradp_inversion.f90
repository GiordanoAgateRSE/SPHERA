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
! Program unit: B_ren_gradp_inversion
! Description: Inversion of the renormalization matrix of an input fluid 
!              particle for the approximation of the pressure-gradient term
!-------------------------------------------------------------------------------
subroutine B_ren_gradp_inversion(npi)
!------------------------
! Modules
!------------------------
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: npi
integer(4) :: test
double precision :: abs_det_thresh
double precision,dimension(3,3) :: aux_mat
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine Matrix_Inversion_3x3(mat,inv,abs_det_thresh,test)
      implicit none
      double precision,intent(in) :: abs_det_thresh
      double precision,intent(in) :: mat(3,3)
      double precision,intent(inout) :: inv(3,3)
      integer(4),intent(inout) :: test
   end subroutine Matrix_Inversion_3x3
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
abs_det_thresh = 1.d-3
!------------------------
! Statements
!------------------------
! Matrix inversion to get the renormalization matrix
call Matrix_Inversion_3x3(pg(npi)%B_ren_gradp,aux_mat,abs_det_thresh,test)
if (test==1) then
   pg(npi)%B_ren_gradp(1:3,1:3) = aux_mat(1:3,1:3)
   pg(npi)%B_ren_gradp_stat = 1
   else
      pg(npi)%B_ren_gradp(1:3,1:3) = 0.d0
      pg(npi)%B_ren_gradp(1,1) = -1.d0
      pg(npi)%B_ren_gradp(2,2) = -1.d0
      pg(npi)%B_ren_gradp(3,3) = -1.d0
      pg(npi)%B_ren_gradp_stat = 0
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine B_ren_gradp_inversion
