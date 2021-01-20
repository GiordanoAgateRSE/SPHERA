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
! Program unit: vector_rotation_axis_angle
! Description: Provided 2 vectors, this subroutine computes the rotation axis 
!              and the rotation angle which allow rotating from the unit vector 
!              aligned with the first vector to the unit vector aligned with 
!              the second vector.
!-------------------------------------------------------------------------------
subroutine vector_rotation_axis_angle(vector1,vector2,n_R,teta_R)
!------------------------
! Modules
!------------------------
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(in) :: vector1(3),vector2(3)
double precision,intent(out) :: teta_R
double precision,intent(out) :: n_R(3)
double precision :: abs_vector1,abs_vector2,abs_cross_product,sin_teta_R
double precision :: cos_teta_R
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
abs_vector1 = dsqrt(dot_product(vector1(:),vector1(:)))
abs_vector2 = dsqrt(dot_product(vector2(:),vector2(:)))
call Vector_Product(vector1(:),vector2(:),n_R(:),3)
abs_cross_product = dsqrt(dot_product(n_R,n_R))
if (abs_cross_product>1.d-9) then
   n_R(:) = n_R(:) / abs_cross_product
   if ((abs_vector1<1.d-9).or.(abs_vector2<1.d-9)) then
      write(*,*) 'Subroutine "vector_rotation_axis_angle". The on-going body ',&
         'is made up of only one solid particle. The simulation stops here. '
      stop
      else
         sin_teta_R = abs_cross_product / abs_vector1 / abs_vector2
         cos_teta_R = dot_product(vector1(:),vector2(:)) / abs_vector1 /       &
                      abs_vector2
         if (cos_teta_R>1.d-9) then
            teta_R = datan(sin_teta_R / cos_teta_R)
            elseif (cos_teta_R<-1.d-9) then
               teta_R = datan(sin_teta_R / cos_teta_R) + PIGRECO
               else
                  if (sin_teta_R>1.d-9) then
                     teta_R = PIGRECO / 2.d0
                     elseif (sin_teta_R<-1.d-9) then
                        teta_R = - PIGRECO / 2.d0
                  endif
         endif
   endif         
   else
      n_R(1) = 0.d0
      n_R(2) = 0.d0
      n_R(3) = 1.d0
      teta_R = 0.d0
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine vector_rotation_axis_angle
