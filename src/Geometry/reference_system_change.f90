!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-Ricerca di Sistema)



! SPHERA authors and email contact are provided on SPHERA documentation.

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
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: reference_system_change   
! Description: Transformation of coordinates, expressed in a new reference system.        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine reference_system_change(pos_old_ref,new_origin_old_ref,             &
                                   new_cos_dir_old_ref,pos_new_ref)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: i,aux_int
double precision,intent(IN) :: pos_old_ref(3),new_origin_old_ref(3) 
double precision,intent(IN) :: new_cos_dir_old_ref(3,3)
double precision,intent(INOUT) :: pos_new_ref(3)
double precision :: aux_vec(3)
double precision :: inv_new_cos_dir_old_ref(3,3)
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
aux_vec(:) = pos_old_ref(:) - new_origin_old_ref(:)
!------------------------
! Statements
!------------------------
call Matrix_Inversion_3x3(new_cos_dir_old_ref,inv_new_cos_dir_old_ref,aux_int)
do i=1,3
   pos_new_ref(i) = dot_product(inv_new_cos_dir_old_ref(i,:),aux_vec)
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine reference_system_change

