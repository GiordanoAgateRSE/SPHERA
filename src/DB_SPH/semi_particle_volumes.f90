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
! Program unit: semi_particle_volumes
! Description: To compute the semi-particle shape coefficients and volumes.                  
!----------------------------------------------------------------------------------------------------------------------------------

subroutine semi_particle_volumes
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
logical :: aux_false_hyp 
integer(4) :: j,i,ID_P_0_iso,ID_P_b_iso,dealloc_stat,aux_adjacent_faces
integer(4) :: n_adj_faces_max
double precision :: aux_scalar,aux_scalar_2,k_parameter,alfa,alfa_sum
double precision,dimension(3) :: aux_vec
double precision,dimension(4,3) :: aux_face1,aux_face2
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine adjacent_faces_isolated_points(face1,face2,ID_face1_iso,         &
                                             ID_face2_iso,false_hyp)
      implicit none
      double precision, dimension(4,3), intent(in) :: face1,face2
      integer(4),intent(out) :: ID_face1_iso,ID_face2_iso
      logical,intent(out) :: false_hyp 
   end subroutine adjacent_faces_isolated_points
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
if (ncord==3) then
   n_adj_faces_max=3
   else
      n_adj_faces_max=2
endif
!------------------------
! Statements
!------------------------
! Loop over the DBSPH surface elements
!$omp parallel do default(none)                                                &
!$omp shared(DBSPH,pg_w,ncord,Domain,Med,pg,nout,n_adj_faces_max)              &
!$omp private(aux_scalar,i,j,k_parameter,alfa,alfa_sum,aux_adjacent_faces)     &
!$omp private(aux_face1,aux_face2,aux_false_hyp,aux_scalar_2,aux_vec)          &
!$omp private(ID_P_b_iso,ID_P_0_iso)
do i=1,DBSPH%n_w
   alfa_sum = 0.d0
   aux_adjacent_faces = 0
! Loop over the adjacent faces 
   do j=1,n_adj_faces_max
      if (pg_w(i)%adjacent_faces(j)==0) cycle
! Provided 2 adjacent triangular faces, find the 2 vertices not in common, 
! one per face
      aux_face1(1,:) =                                                         &
         DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(1))%pos(:)
      aux_face1(2,:) =                                                         &
         DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(2))%pos(:)
      aux_face1(3,:) =                                                         &
         DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(3))%pos(:)
      aux_face2(1,:) =                                                         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(pg_w(i)%adjacent_faces(j))%vert_list(1))%pos(:)
      aux_face2(2,:) =                                                         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(pg_w(i)%adjacent_faces(j))%vert_list(2))%pos(:)
      aux_face2(3,:) =                                                         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(pg_w(i)%adjacent_faces(j))%vert_list(3))%pos(:)
      if (ncord==2) then
         aux_face1(4,:) =                                                      &
         DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(4))%pos(:)
         aux_face2(4,:) =                                                      &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(pg_w(i)%adjacent_faces(j))%vert_list(4))%pos(:)
      endif
      aux_adjacent_faces = aux_adjacent_faces + 1          
      call adjacent_faces_isolated_points(aux_face1,aux_face2,ID_P_0_iso,      &
                                          ID_P_b_iso,aux_false_hyp)
      if (aux_false_hyp.eqv..true.) then
         write(nout,*) "Error! Two faces are not adjacent, but they should be "&
                       ,"(subroutine semi_particle_volume), Face1: ",i,        &
                       " Face2: ",pg_w(i)%adjacent_faces(j),                   &
                       "; the run terminates here."
         stop
      endif 
! Compute the spread angle between the normal vectors
      aux_scalar = dot_product(pg_w(i)%normal,                                 &
                               pg_w(pg_w(i)%adjacent_faces(j))%normal)
      aux_vec(:) = aux_face2(ID_P_b_iso,:) - aux_face1(ID_P_0_iso,:) 
      aux_scalar_2 = dot_product(pg_w(i)%normal,aux_vec)
      if (aux_scalar_2>=0.000001d0) then
          alfa = PIGRECO + dacos(aux_scalar) 
          else if (aux_scalar_2<=-0.000001d0) then 
             alfa = PIGRECO - dacos(aux_scalar)
             else
                alfa = PIGRECO 
      endif       
! Update alfa_summation (algebric sum)
      alfa_sum = alfa_sum + alfa
   end do
! Compute k_d (shape coefficient) and semi-particle volume (area in 2D) and mass
   if (alfa_sum>(aux_adjacent_faces*PIGRECO*0.5d0)) then
      pg_w(i)%k_d = alfa_sum/(PIGRECO*0.5d0*aux_adjacent_faces)-1.d0
      else
         pg_w(i)%k_d = 0.d0
   endif
   pg_w(i)%volume = pg_w(i)%k_d * DBSPH%k_w * pg_w(i)%weight * Domain%dd /     &
                    DBSPH%dx_dxw   
   pg_w(i)%mass = pg_w(i)%volume * Med(1)%den0  
end do
!$omp end parallel do
if (allocated(DBSPH%surf_mesh%vertices)) then
   deallocate(DBSPH%surf_mesh%vertices,STAT=dealloc_stat) 
   if (dealloc_stat/=0) then
      write(nout,*) 'Deallocation of DBSPH%surf_mesh%vertices in ',            &
                    'DBSPH_find_close_faces failed; the program terminates here'
! Stop the main program
      stop 
   endif   
endif 
if (allocated(DBSPH%surf_mesh%faces)) then
   deallocate(DBSPH%surf_mesh%faces,STAT=dealloc_stat) 
   if (dealloc_stat/=0) then
      write(nout,*) 'Deallocation of DBSPH%surf_mesh%faces in ',               &
                    'DBSPH_find_close_faces failed; the program terminates here'
! Stop the main program
      stop 
   endif   
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine semi_particle_volumes

