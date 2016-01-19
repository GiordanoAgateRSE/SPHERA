!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-)



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
! Program unit: DBSPH_kinematics
! Description: Imposing input kinematics for the DB-SPH elements (linear interpolation of input data).             
!----------------------------------------------------------------------------------------------------------------------------------
subroutine DBSPH_kinematics
!------------------------
! Modules
!------------------------ 
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: i,j,i_a
double precision :: rel_pos(3)
double precision :: vel_aux(DBSPH%surface_mesh_files,3)
double precision :: omega_aux(DBSPH%surface_mesh_files,3
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine Vector_Product(uu,VV,ww,SPACEDIM)
     implicit none
     integer(4),intent(IN) :: SPACEDIM
     double precision,intent(IN),dimension(SPACEDIM) :: uu,VV
     double precision,intent(INOUT),dimension(SPACEDIM) :: ww
   end subroutine Vector_Product
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
! Loop over the kinematics records to provide a linear interpolation of the 
! imposed velocity to the first surface element.
if (DBSPH%n_w>0) then
!$omp parallel do default(none)                                                &
!$omp shared(DBSPH,tempo)                                                      &
!$omp private(i,j,vel_aux,omega_aux)
   do i=1,DBSPH%surface_mesh_files
      do j=1,DBSPH%n_kinematics_records(i)
         if (DBSPH%kinematics(i,j,1)>=tempo) then
            if (DBSPH%kinematics(i,j,1)==tempo) then
               vel_aux(i,:) = DBSPH%kinematics(i,j,2:4)
               omega_aux(i,:) = DBSPH%kinematics(i,j,5:7)
               else
                  vel_aux(i,:) = DBSPH%kinematics(i,j-1,2:4) +                 &
                               (DBSPH%kinematics(i,j,2:4) -                    &
                               DBSPH%kinematics(i,j-1,2:4)) /                  &
                               (DBSPH%kinematics(i,j,1) -                      &
                               DBSPH%kinematics(i,j-1,1)) * (tempo -           &
                               DBSPH%kinematics(i,j-1,1))
                  omega_aux(i,:) = DBSPH%kinematics(i,j-1,5:7) +               &
                                 (DBSPH%kinematics(i,j,5:7) -                  &
                                 DBSPH%kinematics(i,j-1,5:7)) /                &
                                 (DBSPH%kinematics(i,j,1) -                    &
                                 DBSPH%kinematics(i,j-1,1)) * (tempo -         &
                                 DBSPH%kinematics(i,j-1,1))
            endif
            exit
         endif
      enddo
   enddo
!$omp parallel do default(none)                                                &
!$omp shared(DBSPH,pg_w,i)                                                     &
!$omp private(i_a,rel_pos,vel_aux,omega_aux)
   do i_a=1,DBSPH%n_w
      rel_pos(:) = pg_w(i_a)%coord(:) -                                        &
                DBSPH%rotation_centre(pg_w(i_a)%surface_mesh_file_ID,:)
      call Vector_Product(omega_aux(pg_w(i_a)%surface_mesh_file_ID,1:3),       &
         rel_pos,pg_w(i_a)%vel,3)
      pg_w(i_a)%vel(:) = pg_w(i_a)%vel(:) +                                    &
                         vel_aux(pg_w(i_a)%surface_mesh_file_ID,:)
   enddo
!$omp end parallel do
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine DBSPH_kinematics

