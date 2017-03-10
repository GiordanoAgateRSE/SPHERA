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
! Program unit: vector_rotation_Rodrigues
! Description: 3D rotation of a given vector, provided the rotation axis and 
!              the rotation angle, based on Rodrigues formula.           
!-------------------------------------------------------------------------------
subroutine vector_rotation_Rodrigues(n_R,teta_R,vector)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(in) :: teta_R
double precision,intent(in) :: n_R(3)
double precision,intent(inout) :: vector(3)
double precision :: aux_vector(3)
double precision :: R_Rodrigues(3,3)
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
aux_vector(:) = vector(:)
!------------------------
! Statements
!------------------------
R_Rodrigues(1,1) = (n_R(1) ** 2) + (1.d0 - (n_R(1) ** 2)) * dcos(teta_R)
R_Rodrigues(1,2) = n_R(1) * n_R(2) * (1.d0 - dcos(teta_R)) - n_R(3) *          &
                   dsin(teta_R)
R_Rodrigues(1,3) = n_R(1) * n_R(3) * (1.d0 - dcos(teta_R)) + n_R(2) *          &
                   dsin(teta_R)
R_Rodrigues(2,1) = n_R(1) * n_R(2) * (1.d0 - dcos(teta_R)) + n_R(3) *          &
                   dsin(teta_R)
R_Rodrigues(2,2) = (n_R(2) ** 2) + (1.d0 - (n_R(2) ** 2)) * dcos(teta_R)
R_Rodrigues(2,3) = n_R(2) * n_R(3) * (1.d0 - dcos(teta_R)) - n_R(1) *          &
                   dsin(teta_R)
R_Rodrigues(3,1) = n_R(1) * n_R(3) * (1.d0 - dcos(teta_R)) - n_R(2) *          &
                   dsin(teta_R)
R_Rodrigues(3,2) = n_R(2) * n_R(3) * (1.d0 - dcos(teta_R)) + n_R(1) *          &
                   dsin(teta_R)
R_Rodrigues(3,3) = (n_R(3) ** 2) + (1.d0 - (n_R(3) ** 2)) * dcos(teta_R)
vector(1) = dot_product(R_Rodrigues(1,:),aux_vector)
vector(2) = dot_product(R_Rodrigues(2,:),aux_vector)
vector(3) = dot_product(R_Rodrigues(3,:),aux_vector)
!------------------------
! Deallocations
!------------------------
return
end subroutine vector_rotation_Rodrigues

