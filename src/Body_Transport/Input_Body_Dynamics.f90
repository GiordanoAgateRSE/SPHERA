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
! Program unit: Input_Body_Dynamics
! Description: Input management for body transport in fluid flows.      
!-------------------------------------------------------------------------------
#ifdef SOLID_BODIES
subroutine Input_Body_Dynamics
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: i,nbi,npi,j,k,nei,erased_part,alloc_stat,aux_integer
double precision :: mod_normal,aux_umax
integer(4) :: vec_temp(3)
double precision :: vec2_temp(3)
double precision, dimension(:), allocatable  ::  aux_mass
type (body_particle), dimension(:), allocatable :: aux_bp_arr 
integer(4),external :: ParticleCellNumber
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine vector_rotation_Rodrigues(n_R,teta_R,vector)
      implicit none
      double precision,intent(in) :: teta_R
      double precision,intent(in) :: n_R(3)
      double precision,intent(inout) :: vector(3)
   end subroutine vector_rotation_Rodrigues
end interface
!------------------------
! Allocations
!------------------------
allocate(aux_mass(n_bodies))
!------------------------
! Initializations
!------------------------
n_body_part = 0
aux_mass = 0.d0
erased_part = 0
n_surf_body_part = 0
aux_integer = 0
!------------------------
! Statements
!------------------------
! Loop over the transported bodies
!$omp parallel do default(none)                                                &
!$omp shared(n_bodies,body_arr,Domain,n_body_part,dx_dxbodies)                 &
!$omp private(i)
do i=1,n_bodies
   body_arr(i)%npart = 0
   do j=1,body_arr(i)%n_elem
! Particle length scales for the element
#ifdef SPACE_3D
      body_arr(i)%elem(j)%dx(:) = body_arr(i)%elem(j)%L_geom(:) /              &
         (int(body_arr(i)%elem(j)%L_geom(:) / (Domain%dx/dx_dxbodies)))
#elif defined SPACE_2D
      body_arr(i)%elem(j)%dx(1) = body_arr(i)%elem(j)%L_geom(1) /              &
         (int(body_arr(i)%elem(j)%L_geom(1) / (Domain%dx/dx_dxbodies)))
      body_arr(i)%elem(j)%dx(2) = 1.d0
      body_arr(i)%elem(j)%dx(3) = body_arr(i)%elem(j)%L_geom(3) /              &
         (int(body_arr(i)%elem(j)%L_geom(3) / (Domain%dx/dx_dxbodies)))
#endif
! Number of body particles of the element
#ifdef SPACE_3D
         body_arr(i)%elem(j)%npart = int(body_arr(i)%elem(j)%L_geom(1) /       &
            body_arr(i)%elem(j)%dx(1)) * int(body_arr(i)%elem(j)%L_geom(2) /   &
            body_arr(i)%elem(j)%dx(2)) * int(body_arr(i)%elem(j)%L_geom(3) /   &
            body_arr(i)%elem(j)%dx(3))
#elif defined SPACE_2D
            body_arr(i)%elem(j)%npart = int(body_arr(i)%elem(j)%L_geom(1) /    &
               body_arr(i)%elem(j)%dx(1)) * int(body_arr(i)%elem(j)%L_geom(3)  &
               / body_arr(i)%elem(j)%dx(3))
#endif
! Update of the number of body particles of the body
      body_arr(i)%npart = body_arr(i)%npart + body_arr(i)%elem(j)%npart         
   end do
! Initializing umax (maximum particle velocity of the body)
   body_arr(i)%umax = 0.d0 
   body_arr(i)%pmax = 0.d0
enddo
!$omp end parallel do
! Incrementing the total number of particles 
do i=1,n_bodies
   n_body_part = n_body_part + body_arr(i)%npart      
end do
! Managing body particles  
if (.not.allocated(bp_arr)) then
   allocate(bp_arr(n_body_part),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(ulog,*)                                                            &
         'Allocation of bp_arr in Input_Body_Dynamics failed;',                &
         ' the program terminates here.'
      stop ! Stop the main program
      else
         write(ulog,*)                                                        &
            "Allocation of bp_arr in Input_Body_Dynamics successfully completed"
   endif
endif
npi = 1  
! Loop over the transported bodies (not to be parallelized)
do nbi=1,n_bodies
! Loop over the body elements  
   do nei=1,body_arr(nbi)%n_elem
      vec_temp(:) = int(body_arr(nbi)%elem(nei)%L_geom(:) /                    &
                    body_arr(nbi)%elem(nei)%dx(:))     
#ifdef SPACE_2D 
         vec_temp(2) = 1     
#endif
      do i=1,vec_temp(1)
         do j=1,vec_temp(2)
            do k=1,vec_temp(3)
! Corresponding Body 
               bp_arr(npi)%body = nbi
! Initial relative positions (body reference system: temporary)
               bp_arr(npi)%rel_pos(1) = - body_arr(nbi)%elem(nei)%L_geom(1) /  &
                  2.d0 + body_arr(nbi)%elem(nei)%dx(1) / 2.d0 + (i - 1) *      &
                  body_arr(nbi)%elem(nei)%dx(1)
#ifdef SPACE_3D                  
                  bp_arr(npi)%rel_pos(2) = - body_arr(nbi)%elem(nei)%L_geom(2) &
                     / 2.d0 + body_arr(nbi)%elem(nei)%dx(2) / 2.d0 + (j - 1) * &
                     body_arr(nbi)%elem(nei)%dx(2)  
#elif defined SPACE_2D
                     bp_arr(npi)%rel_pos(2) = 0.d0
#endif
               bp_arr(npi)%rel_pos(3) = - body_arr(nbi)%elem(nei)%L_geom(3) /  &
                  2.d0 + body_arr(nbi)%elem(nei)%dx(3) / 2.d0 + (k - 1) *      &
                  body_arr(nbi)%elem(nei)%dx(3) 
! Initial relative normal (with respect to the body axis: temporary) and 
! area(/length in 2D) only for particles on the body surface
               bp_arr(npi)%area = 0.d0  
               bp_arr(npi)%normal = 0.d0
               if ((i==1).and.(body_arr(nbi)%elem(nei)%normal_act(5)==1)) then
                  bp_arr(npi)%normal(1) = 1.d0
#ifdef SPACE_3D
                     bp_arr(npi)%area = bp_arr(npi)%area +                     &
                                        body_arr(nbi)%elem(nei)%dx(2) *        &
                                        body_arr(nbi)%elem(nei)%dx(3)
#elif defined SPACE_2D
                        bp_arr(npi)%area = bp_arr(npi)%area +                  &
                                           body_arr(nbi)%elem(nei)%dx(3)
#endif
               endif
               if ((i==vec_temp(1)).and.                                       &
                  (body_arr(nbi)%elem(nei)%normal_act(6)==1)) then
                  bp_arr(npi)%normal(1) = -1.d0
#ifdef SPACE_3D
                     bp_arr(npi)%area = bp_arr(npi)%area +                     &
                                        body_arr(nbi)%elem(nei)%dx(2) *        &
                                        body_arr(nbi)%elem(nei)%dx(3)
#elif defined SPACE_2D
                        bp_arr(npi)%area = bp_arr(npi)%area +                  &
                                           body_arr(nbi)%elem(nei)%dx(3)
#endif
               endif
#ifdef SPACE_3D
                  if ((j==1).and.(body_arr(nbi)%elem(nei)%normal_act(3)==1))   &
                     then
                     bp_arr(npi)%normal(2) = 1.d0
                     bp_arr(npi)%area = bp_arr(npi)%area +                     &
                                        body_arr(nbi)%elem(nei)%dx(1) *        &
                                        body_arr(nbi)%elem(nei)%dx(3)
                  endif
                  if ((j==vec_temp(2)).and.                                    &
                     (body_arr(nbi)%elem(nei)%normal_act(4)==1)) then
                     bp_arr(npi)%normal(2) = -1.d0
                     bp_arr(npi)%area = bp_arr(npi)%area +                     &
                                        body_arr(nbi)%elem(nei)%dx(1) *        &
                                        body_arr(nbi)%elem(nei)%dx(3) 
                  endif
#endif
               if ((k==1).and.(body_arr(nbi)%elem(nei)%normal_act(1)==1)) then
                  bp_arr(npi)%normal(3) = 1.d0
#ifdef SPACE_3D                  
                     bp_arr(npi)%area = bp_arr(npi)%area +                     &
                                        body_arr(nbi)%elem(nei)%dx(1) *        &
                                        body_arr(nbi)%elem(nei)%dx(2)
#elif defined SPACE_2D
                        bp_arr(npi)%area = bp_arr(npi)%area +                  &
                                           body_arr(nbi)%elem(nei)%dx(1)
#endif
               endif
               if ((k==vec_temp(3)).and.                                       &
                  (body_arr(nbi)%elem(nei)%normal_act(2)==1)) then
                  bp_arr(npi)%normal(3) = -1.d0
#ifdef SPACE_3D
                     bp_arr(npi)%area = bp_arr(npi)%area +                     &
                                        body_arr(nbi)%elem(nei)%dx(1) *        &
                                        body_arr(nbi)%elem(nei)%dx(2)
#elif defined SPACE_2D
                        bp_arr(npi)%area = bp_arr(npi)%area +                  &
                                           body_arr(nbi)%elem(nei)%dx(1)
#endif
               endif
               mod_normal =                                                    &
                  dsqrt(dot_product(bp_arr(npi)%normal,bp_arr(npi)%normal))
               if (mod_normal>1.d0) bp_arr(npi)%normal(:) =                    &
                  bp_arr(npi)%normal(:) / mod_normal
! Preliminary value of the relative positions (rotation of rel_pos with respect
! to the centre of mass of the element); "relative" means with respect to the 
! element origin, in the world/global reference system 
! (not in the element reference system)
               call vector_rotation_Rodrigues(                                 &
                  body_arr(nbi)%elem(nei)%n_R_IO(:),                           &
                  body_arr(nbi)%elem(nei)%teta_R_IO,bp_arr(npi)%rel_pos(:))
#ifdef SPACE_2D
                  bp_arr(npi)%rel_pos(2) = 0.d0
#endif
! Preliminary value of the absolute positions after
! translation of rel_pos of the distance "corresponding element - 
! reference system origin")
               bp_arr(npi)%pos(:) = body_arr(nbi)%elem(nei)%x_CM(:) +          &
                                    bp_arr(npi)%rel_pos(:)
! Body particle masses
! Deactivating particle masses (before any rotation around the centre of mass 
! of the body)
               if (((bp_arr(npi)%pos(1)>=body_arr(nbi)%elem(nei)%mass_deact(1))&
                  .or.                                                         &
                  (bp_arr(npi)%pos(1)<=body_arr(nbi)%elem(nei)%mass_deact(2))  &
                  .or.                                                         &
                  (bp_arr(npi)%pos(2)>=body_arr(nbi)%elem(nei)%mass_deact(3))  &
                  .or.                                                         &
                  (bp_arr(npi)%pos(2)<=body_arr(nbi)%elem(nei)%mass_deact(4))  &
                  .or.                                                         &
                  (bp_arr(npi)%pos(3)>=body_arr(nbi)%elem(nei)%mass_deact(5))  &
                  .or.                                                         &
                  (bp_arr(npi)%pos(3)<=body_arr(nbi)%elem(nei)%mass_deact(6))) &
                  .and.(mod_normal==0.d0)) then
                  bp_arr(npi)%mass = 0.d0
                  else
                     bp_arr(npi)%mass = body_arr(nbi)%mass / body_arr(nbi)%npart
                     aux_mass(nbi) = aux_mass(nbi) + bp_arr(npi)%mass
               endif
! Since here, the relative position refers to the centre of mass of the body
               bp_arr(npi)%rel_pos(:) = bp_arr(npi)%pos(:) -                   &
                                        body_arr(nbi)%x_CM(:)              
! Initial normal (after rotation around the centre of mass of the element)
               call vector_rotation_Rodrigues(                                 &
                  body_arr(nbi)%elem(nei)%n_R_IO(:),                           &
                  body_arr(nbi)%elem(nei)%teta_R_IO,bp_arr(npi)%normal(:))
#ifdef SPACE_2D
                  bp_arr(npi)%normal(2) = 0.d0
#endif
! Rotations around the centre of rotation provided in input
               bp_arr(npi)%rel_pos(:) = bp_arr(npi)%pos(:) -                   &
                                        body_arr(nbi)%x_rotC(:)
               call vector_rotation_Rodrigues(body_arr(nbi)%n_R_IO(:),         &
               body_arr(nbi)%teta_R_IO,bp_arr(npi)%rel_pos(:))
#ifdef SPACE_2D
                  bp_arr(npi)%rel_pos(2) = 0.d0
#endif
               bp_arr(npi)%pos(:) = bp_arr(npi)%rel_pos(:) +                   &
                                    body_arr(nbi)%x_rotC(:) 
               bp_arr(npi)%rel_pos(:) = bp_arr(npi)%pos(:) -                   &
                  body_arr(nbi)%x_CM(:)
               call vector_rotation_Rodrigues(body_arr(nbi)%n_R_IO(:),         &
               body_arr(nbi)%teta_R_IO,bp_arr(npi)%normal(:))
#ifdef SPACE_2D
                  bp_arr(npi)%normal(2) = 0.d0
#endif
! Initial cell
               bp_arr(npi)%cell =  ParticleCellNumber(bp_arr(npi)%pos)
! Velocity
               call Vector_Product(body_arr(nbi)%omega,bp_arr(npi)%rel_pos,    &
                  vec2_temp,3)
               bp_arr(npi)%vel(:) = body_arr(nbi)%u_CM(:) + vec2_temp(:)
#ifdef SPACE_2D
                  bp_arr(npi)%vel(2) = 0.d0
#endif                   
! Update of umax
               aux_umax = dsqrt(dot_product(bp_arr(npi)%vel,bp_arr(npi)%vel))        
               body_arr(bp_arr(npi)%body)%umax =                               &
                  max(body_arr(bp_arr(npi)%body)%umax,aux_umax)       
! Formal initialization of pressure and acceleration
               bp_arr(npi)%pres = 0.d0              
               bp_arr(npi)%acc(:) = 0.d0
               bp_arr(npi)%vel_mir(:) = 0.d0
               npi = npi + 1                                                                     
            enddo
         enddo
      enddo
   enddo
enddo
! Cancel useless particles
allocate(aux_bp_arr(n_body_part))
do i=1,n_body_part
   if (bp_arr(i)%mass>0.d0) then
      aux_bp_arr(i-erased_part) = bp_arr(i)
      else
         erased_part = erased_part + 1
         body_arr(bp_arr(i)%body)%npart = body_arr(bp_arr(i)%body)%npart -1
   endif
enddo
n_body_part = n_body_part - erased_part 
deallocate(bp_arr)
allocate(bp_arr(n_body_part))
!$omp parallel do default(none) private(i) shared(n_body_part,bp_arr,aux_bp_arr)
do i=1,n_body_part
   bp_arr(i) = aux_bp_arr(i) 
enddo
!$omp end parallel do  
! Counting surface body particles
do i=1,n_body_part
   if (bp_arr(i)%area>0.d0) n_surf_body_part = n_surf_body_part + 1  
enddo
! Allocating the array of the surface body particles
if (.not.allocated(surf_body_part)) then
   allocate(surf_body_part(n_surf_body_part),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(ulog,*)                                                            &
         'Allocation of surf_body_part in Input_Body_Dynamics failed;',        &
         ' the program terminates here.'
      stop ! Stop the main program
      else
         write(ulog,*)                                                        &
            'Allocation of surf_body_part in Input_Body_Dynamics ',            &
            ' successfully completed.'
   endif
endif
surf_body_part = 0
! Listing surface body particles
j = 1
do i=1,n_body_part
   if (bp_arr(i)%area>0.d0) then
      surf_body_part(j) = i
      j = j + 1
   endif
enddo
! Redistributing lost mass 
!$omp parallel do default(none) private(i) shared(n_body_part,bp_arr,aux_mass) &
!$omp shared(body_arr)
do i=1,n_body_part
   if (bp_arr(i)%mass>0.d0) bp_arr(i)%mass = bp_arr(i)%mass -                  &
      (aux_mass(bp_arr(i)%body) - body_arr(bp_arr(i)%body)%mass) / n_body_part
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
deallocate(aux_bp_arr)
do i=1,n_bodies
! To initialize rel_pos_part1_t0
   body_arr(i)%rel_pos_part1_t0(:) = bp_arr(aux_integer+1)%rel_pos(:)
   aux_integer = aux_integer + body_arr(i)%npart
   deallocate(body_arr(i)%elem)
enddo
deallocate(aux_mass)
return
end subroutine Input_Body_Dynamics
#endif
