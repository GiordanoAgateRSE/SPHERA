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
! Program unit: Matrix_Inversion_3x3  
! Description: Computation of the inverse (inv) of a provided 3x3 matrix (mat).
!              It is also called in 2D.
!-------------------------------------------------------------------------------
subroutine Matrix_Inversion_3x3(mat,inv,abs_det_thresh,test)
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
! abs_det_thresh: absolute value of the minimum threshold below which no 
! inversion applies 
double precision,intent(in) :: abs_det_thresh
double precision,intent(in) :: mat(3,3)
double precision,intent(inout) :: inv(3,3)
integer(4),intent(inout) :: test
double precision :: det,aux_scal
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
det = mat(1,1) * mat(2,2) * mat(3,3) + mat(2,1) * mat(3,2) * mat(1,3) +        &
      mat(3,1) * mat(1,2) * mat(2,3) - mat(1,1) * mat(3,2) * mat(2,3) -        &
      mat(3,1) * mat(2,2) * mat(1,3) - mat(2,1) * mat(1,2) * mat(3,3)
aux_scal = dabs(det)
if (aux_scal>=abs_det_thresh) then
   test = 1
   inv(1,1) = mat(2,2) * mat(3,3) - mat(2,3) * mat(3,2)
   inv(1,2) = mat(1,3) * mat(3,2) - mat(1,2) * mat(3,3)
   inv(1,3) = mat(1,2) * mat(2,3) - mat(1,3) * mat(2,2)   
   inv(2,1) = mat(2,3) * mat(3,1) - mat(2,1) * mat(3,3)
   inv(2,2) = mat(1,1) * mat(3,3) - mat(1,3) * mat(3,1)
   inv(2,3) = mat(1,3) * mat(2,1) - mat(1,1) * mat(2,3)
   inv(3,1) = mat(2,1) * mat(3,2) - mat(2,2) * mat(3,1)
   inv(3,2) = mat(1,2) * mat(3,1) - mat(1,1) * mat(3,2)
   inv(3,3) = mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1)
   inv(:,:) = inv(:,:) / det
   else
      test = 0
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine Matrix_Inversion_3x3
