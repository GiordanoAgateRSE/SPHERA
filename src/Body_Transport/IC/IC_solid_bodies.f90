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
! Program unit: IC_solid_bodies
! Description: Initial Conditions for solid bodies. Both handmade and CAE-made 
!              bodies are considered.
!-------------------------------------------------------------------------------
#ifdef SOLID_BODIES
subroutine IC_solid_bodies
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_file_module
use Memory_I_O_interface_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: i,nbi,npi,j,k,nei,erased_part,aux_integer,ib,i_bp
double precision :: mod_normal,aux_umax
integer(4) :: vec_temp(3)
double precision :: vec2_temp(3)
character(100) :: array_name
double precision, dimension(:), allocatable  ::  aux_mass
type(body_particle),dimension(:),allocatable :: aux_bp_arr 
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
   subroutine Vector_Product(uu,VV,ww,SPACEDIM)
      implicit none
      integer(4),intent(in) :: SPACEDIM
      double precision,intent(in),dimension(SPACEDIM) :: uu,VV
      double precision,intent(inout),dimension(SPACEDIM) :: ww
   end subroutine Vector_Product
end interface
!------------------------
! Allocations
!------------------------
! Auxiliary temporary array for the body masses
array_name = "aux_mass"
call allocate_de_dp_r1(.true.,aux_mass,n_bodies,array_name,ulog_flag=.true.)
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
!$omp shared(n_bodies,body_arr,Domain,n_body_part,dx_dxbodies,Grid)            &
!$omp private(ib)
do ib=1,n_bodies
   body_arr(ib)%npart = 0
   do j=1,body_arr(ib)%n_elem
! Particle length scales for the element
#ifdef SPACE_3D
      body_arr(ib)%elem(j)%dx(:) = body_arr(ib)%elem(j)%L_geom(:) /            &
         (int(body_arr(ib)%elem(j)%L_geom(:) / (Domain%dx/dx_dxbodies)))
#elif defined SPACE_2D
      body_arr(ib)%elem(j)%dx(1) = body_arr(ib)%elem(j)%L_geom(1) /            &
         (int(body_arr(ib)%elem(j)%L_geom(1) / (Domain%dx/dx_dxbodies)))
      body_arr(ib)%elem(j)%dx(2) = 1.d0
      body_arr(ib)%elem(j)%dx(3) = body_arr(ib)%elem(j)%L_geom(3) /            &
         (int(body_arr(ib)%elem(j)%L_geom(3) / (Domain%dx/dx_dxbodies)))
#endif
! Number of body particles of the element
#ifdef SPACE_3D
         body_arr(ib)%elem(j)%npart = int(body_arr(ib)%elem(j)%L_geom(1) /     &
            body_arr(ib)%elem(j)%dx(1)) * int(body_arr(ib)%elem(j)%L_geom(2) / &
            body_arr(ib)%elem(j)%dx(2)) * int(body_arr(ib)%elem(j)%L_geom(3) / &
            body_arr(ib)%elem(j)%dx(3))
#elif defined SPACE_2D
            body_arr(ib)%elem(j)%npart = int(body_arr(ib)%elem(j)%L_geom(1) /  &
               body_arr(ib)%elem(j)%dx(1)) * int(body_arr(ib)%elem(j)%L_geom(3)&
               / body_arr(ib)%elem(j)%dx(3))
#endif
! Update the number of body particles of the body
      body_arr(ib)%npart = body_arr(ib)%npart + body_arr(ib)%elem(j)%npart
   enddo
! Initializations
   body_arr(ib)%umax = 0.d0
   body_arr(ib)%pmax = 0.d0
   body_arr(ib)%volume = 0.d0
   body_arr(ib)%area = 0.d0
   body_arr(ib)%int_ndA(1:3) = 0.d0
   body_arr(ib)%alfa(1:3) = 0.d0
   body_arr(ib)%Force(1:3) = 0.d0
   body_arr(ib)%Moment(1:3) = 0.d0
enddo
!$omp end parallel do
! Incrementing the total number of particles 
do ib=1,n_bodies
   n_body_part = n_body_part + body_arr(ib)%npart
enddo
! Managing body particles
array_name = "bp_arr"
call allocate_de_BodPar_r1(.true.,bp_arr,n_body_part,array_name,               &
   ulog_flag=.true.)
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
! Particle volume
#ifdef SPACE_3D
               bp_arr(npi)%volume = (Domain%dx / dx_dxbodies) ** 3
#elif defined SPACE_2D
               bp_arr(npi)%volume = (Domain%dx / dx_dxbodies) ** 2
#endif
! Update the body volume
               body_arr(nbi)%volume = body_arr(nbi)%volume + bp_arr(npi)%volume
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
               bp_arr(npi)%normal(1:3) = 0.d0
               bp_arr(npi)%surface = .false.
               if ((i==1).and.(body_arr(nbi)%elem(nei)%normal_act(5)==1)) then
                  bp_arr(npi)%surface = .true.
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
                  bp_arr(npi)%surface = .true.
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
                     bp_arr(npi)%surface = .true.
                     bp_arr(npi)%normal(2) = 1.d0
                     bp_arr(npi)%area = bp_arr(npi)%area +                     &
                                        body_arr(nbi)%elem(nei)%dx(1) *        &
                                        body_arr(nbi)%elem(nei)%dx(3)
                  endif
                  if ((j==vec_temp(2)).and.                                    &
                     (body_arr(nbi)%elem(nei)%normal_act(4)==1)) then
                     bp_arr(npi)%surface = .true.
                     bp_arr(npi)%normal(2) = -1.d0
                     bp_arr(npi)%area = bp_arr(npi)%area +                     &
                                        body_arr(nbi)%elem(nei)%dx(1) *        &
                                        body_arr(nbi)%elem(nei)%dx(3) 
                  endif
#endif
               if ((k==1).and.(body_arr(nbi)%elem(nei)%normal_act(1)==1)) then
                  bp_arr(npi)%surface = .true.
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
                  bp_arr(npi)%surface = .true.
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
! Update the body area
               body_arr(nbi)%area = body_arr(nbi)%area + bp_arr(npi)%area
! Body particle normal
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
               if (Grid%dcd(1)>1.d-21) then
! The "if construct" is needed just for a formal assignation in case of restart
                  bp_arr(npi)%cell = ParticleCellNumber(bp_arr(npi)%pos)
               endif
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
! Update the integral of the normal over the external body surface
               body_arr(nbi)%int_ndA(1:3) = body_arr(nbi)%int_ndA(1:3) +       &
                  bp_arr(npi)%area * bp_arr(npi)%normal(1:3)
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
! Allocate the auxiliary array of body particles 
array_name = "aux_bp_arr"
call allocate_de_BodPar_r1(.true.,aux_bp_arr,n_body_part,array_name,           &
   ulog_flag=.true.)
! Cancel useless particles
do i_bp=1,n_body_part
   if (bp_arr(i_bp)%mass>0.d0) then
      aux_bp_arr(i_bp-erased_part) = bp_arr(i_bp)
      else
         erased_part = erased_part + 1
         body_arr(bp_arr(i_bp)%body)%npart = body_arr(bp_arr(i_bp)%body)%npart &
                                             - 1
   endif
enddo
n_body_part = n_body_part - erased_part
! Managing body particles
array_name = "bp_arr"
call allocate_de_BodPar_r1(.false.,bp_arr,array_name=array_name,               &
   ulog_flag=.true.)
call allocate_de_BodPar_r1(.true.,bp_arr,n_body_part,array_name,               &
   ulog_flag=.true.)
!$omp parallel do default(none)                                                &
!$omp shared(n_body_part,bp_arr,aux_bp_arr)                                    &
!$omp private(i_bp) 
do i_bp=1,n_body_part
   bp_arr(i_bp) = aux_bp_arr(i_bp) 
enddo
!$omp end parallel do
! Redistributing lost mass
!$omp parallel do default(none)                                                &
!$omp shared(n_body_part,bp_arr,aux_mass,body_arr)                             &
!$omp private(i_bp)
do i_bp=1,n_body_part
   if (bp_arr(i_bp)%mass>0.d0) bp_arr(i_bp)%mass = bp_arr(i_bp)%mass -         &
      (aux_mass(bp_arr(i_bp)%body) - body_arr(bp_arr(i_bp)%body)%mass) /       &
      n_body_part
enddo
!$omp end parallel do
#ifdef SPACE_3D
! Read the CAE-made bodies and assign ICs
call IC_CAE_bodies
#endif
! Counting surface body particles
do i_bp=1,n_body_part
   if (bp_arr(i_bp)%surface) n_surf_body_part = n_surf_body_part + 1  
enddo
! Allocating and initialization of the array of the surface body particles
array_name = "surf_body_part"
call allocate_de_int4_r1(.true.,surf_body_part,n_surf_body_part,array_name,    &
   ulog_flag=.true.)
surf_body_part(1:n_surf_body_part) = 0
! Listing surface body particles
j = 1
do i_bp=1,n_body_part
   if (bp_arr(i_bp)%surface) then
      surf_body_part(j) = i_bp
      j = j + 1
   endif
enddo
do ib=1,n_bodies
! To initialize rel_pos_part1_t0
   body_arr(ib)%rel_pos_part1_t0(1:3) = bp_arr(aux_integer+1)%rel_pos(1:3)
   aux_integer = aux_integer + body_arr(ib)%npart
   if (.not.body_arr(ib)%CAE) then
      array_name = "body_arr(ib)%elem"
      call allocate_de_Bod_elem_r1(.false.,body_arr(ib)%elem,                  &
         array_name=array_name,ulog_flag=.true.)
   endif
enddo
!------------------------
! Deallocations
!------------------------
! Deallocation of the auxiliary array of body particles
array_name = "aux_bp_arr"
call allocate_de_BodPar_r1(.false.,aux_bp_arr,array_name=array_name,           &
   ulog_flag=.true.)
array_name = "aux_mass"
call allocate_de_dp_r1(.false.,aux_mass,array_name=array_name,ulog_flag=.true.)
return
end subroutine IC_solid_bodies
#endif
