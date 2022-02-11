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
! Program unit: vector_rotation_Euler_angles
! Description: 3D rotation of a given vector, provided the vector of Euler's 
!              angles (3D). It is als ocalled in 2D.
!-------------------------------------------------------------------------------
#ifdef SOLID_BODIES
subroutine vector_rotation_Euler_angles(vector,angle)
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(in) :: angle(3) 
double precision,intent(inout) :: vector(3)
double precision :: vec_temp(3)
double precision :: mat1_temp(3,3),mat2_temp(3,3),mat3_temp(3,3),mat4_temp(3,3)
double precision :: cos_dir(3,3)
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
mat1_temp(1,1) = 1.d0
mat1_temp(1,2) = 0.d0
mat1_temp(1,3) = 0.d0
mat1_temp(2,1) = 0.d0
mat1_temp(2,2) = dcos(angle(1))
mat1_temp(2,3) = - dsin(angle(1))
mat1_temp(3,1) = 0.d0
mat1_temp(3,2) = dsin(angle(1))
mat1_temp(3,3) = dcos(angle(1))
mat2_temp(1,1) = dcos(angle(2))
mat2_temp(1,2) = 0.d0
mat2_temp(1,3) = dsin(angle(2))
mat2_temp(2,1) = 0.d0
mat2_temp(2,2) = 1.d0
mat2_temp(2,3) = 0.d0
mat2_temp(3,1) = - dsin(angle(2)) 
mat2_temp(3,2) = 0.d0
mat2_temp(3,3) = dcos(angle(2))
mat3_temp(1,1) = dcos(angle(3))
mat3_temp(1,2) = - dsin(angle(3))
mat3_temp(1,3) = 0.d0
mat3_temp(2,1) = dsin(angle(3))
mat3_temp(2,2) = dcos(angle(3))
mat3_temp(2,3) = 0.d0
mat3_temp(3,1) = 0.d0 
mat3_temp(3,2) = 0.d0
mat3_temp(3,3) = 1.d0    
call MatrixProduct(mat1_temp,mat2_temp,mat4_temp,3,3,3)
call MatrixProduct(mat4_temp,mat3_temp,cos_dir,3,3,3)
vec_temp(:) = vector(:)
vector(1) = dot_product(cos_dir(1,:),vec_temp)
vector(2) = dot_product(cos_dir(2,:),vec_temp)
vector(3) = dot_product(cos_dir(3,:),vec_temp)
!------------------------
! Deallocations
!------------------------
return
end subroutine vector_rotation_Euler_angles
#endif
