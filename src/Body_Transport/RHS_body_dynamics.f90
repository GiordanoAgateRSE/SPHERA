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
! Program unit: RHS_body_dynamics  
! Description:  To estimate the RHS of the body dynamics equations (Amicarelli 
!               et al.,2015,CAF).     
!-------------------------------------------------------------------------------
subroutine RHS_body_dynamics(dtvel)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(in) :: dtvel
integer(4) :: npartint,i,j,npi,npj,Ncb,Nfzn,aux,nbi,npk,k,nbj,nbk
integer(4) :: n_interactions,aux2,aux3,test,aux_locx_min,aux_locx_max,aux_int
! AA!!! test
integer(4) :: aux_scal_test
double precision :: c2,k_masses,r_per,r_par,temp_dden,temp_acc,alfa_boun
double precision :: aux_impact_vel,aux4,pres_mir,aux_scalar,aux_scalar_2
double precision :: friction_limiter
double precision :: f_pres(3),temp(3),r_par_vec(3),f_coll_bp_bp(3)          
double precision :: f_coll_bp_boun(3),dvar(3),pos_aux(3),normal_plane(3)
double precision :: u_rel(3),x_rel(3),aux_acc(3),aux_vec(3),aux_vec2(3)
double precision :: loc_pos(3),aux_locx_vert(3),sliding_friction_dir(3)
double precision :: sliding_friction(3)
! AA!!! test
double precision :: aux_mat(3,3)
integer(4),dimension(:),allocatable :: inter_front
double precision,dimension(:),allocatable :: interface_sliding_vel_max
double precision,dimension(:,:),allocatable :: Force_mag_sum,r_per_min
double precision,dimension(:,:),allocatable :: aux_gravity,sum_normal_plane
double precision,dimension(:,:,:),allocatable :: Force,Moment
character(255) :: file_name_test
double precision, external :: Gamma_boun
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine point_inout_convex_non_degenerate_polygon(point,n_sides,         &
                                                        point_pol_1,           &
                                                        point_pol_2,           &
                                                        point_pol_3,           &
                                                        point_pol_4,           &
                                                        point_pol_5,           &
                                                        point_pol_6,test)
      implicit none
      integer(4),intent(in) :: n_sides
      double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
      double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
      double precision,intent(in) :: point_pol_5(2),point_pol_6(2)
      integer(4),intent(inout) :: test
      double precision :: dis1,dis2
      double precision :: normal(2)
   end subroutine point_inout_convex_non_degenerate_polygon
end interface
!------------------------
! Allocations
!------------------------
if (ncord==2) aux = n_bodies+NumBSides
if (ncord==3) aux = n_bodies+NumFacce
allocate(Force(n_bodies,aux,3))
allocate(Moment(n_bodies,aux,3))
allocate(Force_mag_sum(n_bodies,aux))
allocate(r_per_min(n_bodies,aux))
allocate(aux_gravity(n_bodies,3))
allocate(sum_normal_plane(n_bodies,3))
! AA!!! test
allocate(inter_front(n_bodies))
allocate(interface_sliding_vel_max(n_bodies))
!------------------------
! Initializations
!------------------------
Force = 0.d0
Moment = 0.d0
Force_mag_sum = 0.d0
r_per_min = 1000000.d0
aux2 = 0
! AA!!! test start
inter_front(:) = 0
aux_scal_test = 0
! AA!!! test end
sum_normal_plane(:,:) = 0.d0
interface_sliding_vel_max(:) = 0.d0
!------------------------
! Statements
!------------------------
! Updating pressure of the body particles
call body_pressure_mirror
! Contributions to fluid dynamics momentum (discretized semi-analytic approach: 
! a mirror particle technique)
do npi=1,n_body_part
   do j=1,nPartIntorno_bp_f(npi)
      npartint = (npi - 1) * NMAXPARTJ + j
      npj = PartIntorno_bp_f(npartint)
      temp_acc = (pg(npj)%pres + bp_arr(npi)%pres) / (pg(npj)%dens *           &
                 pg(npj)%dens)
      pg(npj)%acc(:) = pg(npj)%acc(:) + ( - pg(npj)%mass / (dx_dxbodies **     &
                       ncord) * temp_acc * ( - rag_bp_f(:,npartint)) *         &
                       KerDer_bp_f_Gal(npartint))    
   enddo
enddo
! Loop over the transported bodies (gravity forces and Ic initialization)
!$omp parallel do default(none) private(i)                                     &
!$omp shared(n_bodies,body_arr,Domain,it_start,on_going_time_step,ncord)       &
!$omp shared(aux_gravity,simulation_time,time_max_no_body_gravity_force)
do i=1,n_bodies
   if (body_arr(i)%imposed_kinematics==0) then
      if (simulation_time>time_max_no_body_gravity_force) then
         aux_gravity(i,:) = body_arr(i)%mass * Domain%grav(:)
         else
            aux_gravity(i,:) = 0.d0
      endif
      body_arr(i)%Force = 0.d0
      body_arr(i)%Moment = 0.d0
      if ((body_arr(i)%Ic_imposed==0).and.((ncord==3).or.                      &
         ((it_start+1)==on_going_time_step)) ) then
         body_arr(i)%Ic = 0.d0
      endif
   endif
enddo
!$omp end parallel do
! Loop over the body particles  
do npi=1,n_body_part
   if (body_arr(bp_arr(npi)%body)%imposed_kinematics==0) then
! Only body particles at the body surface     
      if (bp_arr(npi)%area>0.d0) then  
! Computation of the ID of the surface body particles
         aux2 = aux2+1
! Fluid pressure forces
         f_pres(:) = bp_arr(npi)%pres * bp_arr(npi)%area * bp_arr(npi)%normal(:)
         body_arr(bp_arr(npi)%body)%Force(:) =                                 &
            body_arr(bp_arr(npi)%body)%Force(:) + f_pres(:)
         call Vector_Product(bp_arr(npi)%rel_pos(:),f_pres(:),temp(:),3)
         body_arr(bp_arr(npi)%body)%Moment(:) =                                &
            body_arr(bp_arr(npi)%body)%Moment(:) + temp(:)
! Loop over the neighbouring body particles, belonging to other bodies 
! (inter-body impacts, partial contributions) 
         do j=1,nPartIntorno_bp_bp(aux2)
            npartint = (aux2 - 1) * NMAXPARTJ + j
            npj = PartIntorno_bp_bp(npartint)
            r_per = abs(rag_bp_bp(1,npartint) * bp_arr(npj)%normal(1) +        &
                    rag_bp_bp(2,npartint) * bp_arr(npj)%normal(2) +            &
                    rag_bp_bp(3,npartint) * bp_arr(npj)%normal(3))
! (note: rag=-r) 
            r_par_vec(:) = - rag_bp_bp(:,npartint) - r_per *                   &
                           bp_arr(npj)%normal(:)  
            r_par = dsqrt(r_par_vec(1) * r_par_vec(1) + r_par_vec(2) *         &
                    r_par_vec(2) + r_par_vec(3) * r_par_vec(3))
            if ((r_per>0.d0).and.((r_par/(Domain%dx/dx_dxbodies))<=1.d0)) then
               u_rel(:) = bp_arr(npi)%vel(:) - bp_arr(npj)%vel(:) 
               x_rel(:) = - rag_bp_bp(:,npartint) /                            &
                 dsqrt(dot_product(rag_bp_bp(:,npartint),rag_bp_bp(:,npartint)))
               aux_impact_vel = dot_product(u_rel,x_rel)
               if ((impact_vel(aux2,bp_arr(npj)%body)==0.d0) .or.              &
                  (impact_vel(aux2,bp_arr(npj)%body)<aux_impact_vel)) then
                  k = 0
                  nbi = 1
                  do while (nbi.ne.bp_arr(npi)%body) 
                     k = k + body_arr(nbi)%npart
                     nbi = nbi + 1
                  enddo
                  do npk=(k+1),(k+body_arr(bp_arr(npi)%body)%npart)
                     do aux3=1,n_surf_body_part
                        if (surf_body_part(aux3)==npk) then
                           impact_vel(aux3,bp_arr(npj)%body) =                 &
                              body_arr(bp_arr(npi)%body)%umax +                &
                              body_arr(bp_arr(npj)%body)%umax
                           cycle
                        endif
                     enddo
                  enddo
               endif  
               if (impact_vel(aux2,bp_arr(npj)%body)>0.d0) then 
                  k_masses = (bp_arr(npi)%mass * bp_arr(npj)%mass) /           &
                             (bp_arr(npi)%mass + bp_arr(npj)%mass)
                  f_coll_bp_bp(:) = - (2.d0 *                                  &
                                    (impact_vel(aux2,bp_arr(npj)%body) ** 2) / &
                                    r_per) * k_masses *                        &
                                    Gamma_boun(r_per,Domain%h) * (1 - r_par /  &
                                    (Domain%dx / dx_dxbodies)) *               &
                                    bp_arr(npj)%normal(:)
                  if (r_per<((Domain%dx/dx_dxbodies)/10.d0)) then   
                     write(file_name_test,"(a,a,i8.8,a,i8.8,a,i8.8,a)")        &
                        nomecaso(1:len_trim(nomecaso)),'_close_impact_',       &
                        on_going_time_step,'_',npi,'_',npj,".txt"
                     open (60,file=file_name_test,status="unknown",            &
                        form="formatted")
                     write (60,                                                &
              '((10x,a),4(11x,a),(2x,a),(9x,a),(6x,a),(4x,a),(9x,a),3(11x,a))')&
                       " time"," npi"," npj"," nbi"," nbj"," imp_vel(m/s)",    &
                       " r_per"," k_masses"," Gamma_boun"," r_par",            &
                       " n_x"," n_y"," n_z"
                     write (60,'((g14.7,1x),4(i14,1x),8(g14.7,1x))')           &
                        simulation_time,npi,npj,bp_arr(npi)%body,              &
                        bp_arr(npj)%body,impact_vel(aux2,bp_arr(npj)%body),    &
                        r_per,k_masses,Gamma_boun(r_per,Domain%h),r_par,       &
                        bp_arr(npj)%normal(1),bp_arr(npj)%normal(2),           &
                        bp_arr(npj)%normal(3)
                     close(60)
                  endif
                  Force(bp_arr(npi)%body,bp_arr(npj)%body,:) =                 &
                    Force(bp_arr(npi)%body,bp_arr(npj)%body,:) + f_coll_bp_bp(:)
! Body-body interactions: in the absence of neighburing fluid particles (dry 
! stage), sliding friction force depends on the interaction interface angle 
! (instead of the friction angle) and the normal reaction force under 
! sliding is correct. These forces exactly balance the gravity components.
                  if (body_arr(bp_arr(npi)%body)%pmax<1.d-5) then
                     aux_gravity(bp_arr(npi)%body,:) = 0.d0
                  endif
                  call Vector_Product(bp_arr(npi)%rel_pos,f_coll_bp_bp,temp,3)
                  Moment(bp_arr(npi)%body,bp_arr(npj)%body,:) =                &
                     Moment(bp_arr(npi)%body,bp_arr(npj)%body,:) + temp(:)
                  Force_mag_sum(bp_arr(npi)%body,bp_arr(npj)%body) =           &
                     Force_mag_sum(bp_arr(npi)%body,bp_arr(npj)%body) +        &
                     (1.d0 / r_per) * k_masses * Gamma_boun(r_per,Domain%h) *  &
                     (1.d0 - r_par / (Domain%dx / dx_dxbodies))
                  r_per_min(bp_arr(npi)%body,bp_arr(npj)%body) =               &
                     min(r_per,r_per_min(bp_arr(npi)%body,bp_arr(npj)%body)) 
               endif
            endif                      
         enddo
! Loop over boundaries (body-boundary impacts, partial contributions) (3D case) 
         if (simulation_time>time_max_no_body_frontier_impingements) then
         if (ncord==3) then
            do j=1,NumFacce
            if ((Tratto(BoundaryFace(j)%stretch)%tipo=="fixe") .or.            &
               (Tratto(BoundaryFace(j)%stretch)%tipo == "tapi")) then
               aux_vec2(:) = bp_arr(npi)%pos(:) - BoundaryFace(j)%Node(1)%GX(:)
               aux4 = dot_product(BoundaryFace(j)%T(:,3),aux_vec2)
               if (aux4>=0.d0) then
                  call dis_point_plane(bp_arr(npi)%pos,                        &
                     BoundaryFace(j)%Node(1)%GX,BoundaryFace(j)%Node(2)%GX,    &
                     BoundaryFace(j)%Node(3)%GX,r_per,normal_plane)
                  aux_mat(:,1) = BoundaryFace(j)%T(:,1)
                  aux_mat(:,2) = BoundaryFace(j)%T(:,2)
                  aux_mat(:,3) = BoundaryFace(j)%T(:,3)
                  if (BoundaryFace(j)%nodes==4) then
                     call reference_system_change(bp_arr(npi)%pos,             &
                        BoundaryFace(j)%Node(4)%GX,aux_mat,loc_pos)
                     elseif (BoundaryFace(j)%nodes==3) then
                        call reference_system_change(bp_arr(npi)%pos,          &
                           BoundaryFace(j)%Node(3)%GX,aux_mat,loc_pos)
                  endif
                  call point_inout_convex_non_degenerate_polygon(loc_pos(1:2), &
                     BoundaryFace(j)%nodes,BoundaryFace(j)%Node(1)%LX(1:2),    &
                     BoundaryFace(j)%Node(2)%LX(1:2),                          &
                     BoundaryFace(j)%Node(3)%LX(1:2),                          &
                     BoundaryFace(j)%Node(4)%LX(1:2),                          &
                     BoundaryFace(j)%Node(4)%LX(1:2),                          &
                     BoundaryFace(j)%Node(4)%LX(1:2),test)
                  if ((r_per>0.d0).and.(r_per<=(2.d0*Domain%h)).and.           &
                     (test==1) ) then     
                     aux_impact_vel =                                          &
                         - (dot_product(normal_plane,bp_arr(npi)%vel))
                     if ((impact_vel(aux2,n_bodies+j)==0.d0) .or.              &
                         (impact_vel(aux2,n_bodies+j)<aux_impact_vel) ) then
                        impact_vel(aux2,n_bodies+j) = aux_impact_vel
                     endif
                     if (impact_vel(aux2,n_bodies+j)>0.d0) then 
                        k_masses = bp_arr(npi)%mass
                        f_coll_bp_boun(:) = (2.d0 *                            &
                           (impact_vel(aux2,n_bodies+j) ** 2) / r_per) *       &
                           k_masses * Gamma_boun(r_per,Domain%h)               &
                           * normal_plane(:)
                        Force(bp_arr(npi)%body,n_bodies+j,:) =                 &
                           Force(bp_arr(npi)%body,n_bodies+j,:) +              &
                           f_coll_bp_boun(:)
! Body-frontier interactions. In the absence of neighburing fluid particles
! (dry stage), sliding friction force and explicit normal reaction force under 
! sliding.
                        if (body_arr(bp_arr(npi)%body)%pmax<1.d-5) then
                           if ((friction_angle>(0.d0-1.d-9)).and.              &
                              (simulation_time>time_max_no_body_gravity_force))&
                              then
! Sum of the normal vectors of the neighbouring frontiers
                              sum_normal_plane(bp_arr(npi)%body,:) =           &
                                 sum_normal_plane(bp_arr(npi)%body,:) +        & 
                                 normal_plane(:)
! To update the maximum tangential velocity
                              aux_vec(:) = dot_product(bp_arr(npi)%vel,        &
                                           normal_plane) * normal_plane(:)
                              aux_vec(:) = bp_arr(npi)%vel(:) - aux_vec(:)
                              aux_scalar = dsqrt(dot_product(aux_vec,aux_vec))
                              interface_sliding_vel_max(bp_arr(npi)%body) = max&
                                 (interface_sliding_vel_max(bp_arr(npi)%body), &
                                 aux_scalar)
                              else
! Contribution of the normal reaction force under sliding and of the sliding 
! friction force depending on the local slope angle (instead of the friction 
! angle)
                                 aux_gravity(bp_arr(npi)%body,:) = 0.d0
                           endif
                        endif
                        call Vector_Product(bp_arr(npi)%rel_pos,               &
                           f_coll_bp_boun,temp,3)
                        Moment(bp_arr(npi)%body,n_bodies+j,:) =                &
                           Moment(bp_arr(npi)%body,n_bodies+j,:) + temp(:)
                        Force_mag_sum(bp_arr(npi)%body,n_bodies+j) =           &
                           Force_mag_sum(bp_arr(npi)%body,n_bodies+j) +        &
                           (1.d0 / r_per) * k_masses *                         &
                           Gamma_boun(r_per,Domain%h) 
                        r_per_min(bp_arr(npi)%body,n_bodies+j) =               &
                           min(r_per,r_per_min(bp_arr(npi)%body,n_bodies+j))
                        if (Gamma_boun(r_per,Domain%h)<=0.d0)                  &
                           impact_vel(aux2,n_bodies+j) = 0.d0
! AA!!! test
                        aux_scal_test = aux_scal_test + 1
                        endif
                     endif
                  endif  
               endif
            enddo 
         endif
! AA!!! test start
         inter_front(bp_arr(npi)%body) = max(inter_front(bp_arr(npi)%body),    &
                                         aux_scal_test)
         aux_scal_test = 0
! AA!!! test end
! 2D case
         if (ncord==2) then
            do j=1,NumBSides
               if ((BoundarySide(j)%tipo=="fixe") .or.                         &
                  (BoundarySide(j)%tipo == "tapi")) then  
                  aux_vec2(:) = bp_arr(npi)%pos(:) -                           &
                                Vertice(1,BoundarySide(j)%Vertex(1))
                  aux4 = dot_product(BoundarySide(j)%T(:,3),aux_vec2)
                  if (aux4>=0.d0) then
                     pos_aux(1) =  Vertice(1,BoundarySide(j)%Vertex(1))
                     pos_aux(2) = 0.d0
                     pos_aux(3) = Vertice(3,BoundarySide(j)%Vertex(1))
                     call dis_point_plane(bp_arr(npi)%pos,                     &
                        Vertice(:,BoundarySide(j)%Vertex(1)),                  &
                        Vertice(:,BoundarySide(j)%Vertex(2)),pos_aux,r_per,    &
                        normal_plane)
                     aux_mat(:,1) = BoundarySide(j)%T(:,1)
                     aux_mat(2,1) = 0.d0
                     aux_mat(2,2) = 1.d0
                     aux_mat(2,3) = 0.d0
                     aux_mat(:,3) = BoundarySide(j)%T(:,3) 
                     call reference_system_change(bp_arr(npi)%pos,             &
                        Vertice(:,BoundarySide(j)%Vertex(1)),aux_mat,loc_pos)
                     call reference_system_change(                             &
                        Vertice(:,BoundarySide(j)%Vertex(2)),                  &
                        Vertice(:,BoundarySide(j)%Vertex(1)),aux_mat,          &
                        aux_locx_vert)
                     aux_locx_min = min(zero,aux_locx_vert(1))
                     aux_locx_max = max(zero,aux_locx_vert(1))
                     if ((loc_pos(1)>=aux_locx_min) .and.                      &
                        (loc_pos(1)<=aux_locx_max)) then
                        test = 1
                        else
                           test = 0
                     endif
                     if ((r_per>0.d0).and.(r_per<=(2.d0*Domain%h)).and.        &
                        (test==1) ) then     
                        aux_impact_vel =                                       &
                           - (dot_product(normal_plane,bp_arr(npi)%vel))
                     if ((impact_vel(aux2,n_bodies+j)==0.d0) .or.              &
                        (impact_vel(aux2,n_bodies+j)<aux_impact_vel)) then
                        impact_vel(aux2,n_bodies+j) = aux_impact_vel
                     endif
                     if (impact_vel(aux2,n_bodies+j)>0.d0) then 
                        k_masses = bp_arr(npi)%mass 
                        f_coll_bp_boun(:) = (2.d0*                             &
                           (impact_vel(aux2,n_bodies+j) ** 2) / r_per) *       &
                           k_masses * Gamma_boun(r_per,Domain%h)               &
                           * normal_plane(:)
                        Force(bp_arr(npi)%body,n_bodies+j,:) =                 &
                           Force(bp_arr(npi)%body,n_bodies+j,:) +              &
                           f_coll_bp_boun(:)
! Body-frontier interactions. In the absence of neighburing fluid particles
! (dry stage), sliding friction force and explicit normal reaction force under 
! sliding.
                        if (body_arr(bp_arr(npi)%body)%pmax<1.d-5) then
                           if ((friction_angle>(0.d0-1.d-9)).and.              &
                              (simulation_time>time_max_no_body_gravity_force))&
                              then
! Sum of the normal vectors of the neighbouring frontiers
                              sum_normal_plane(bp_arr(npi)%body,:) =           &
                                 sum_normal_plane(bp_arr(npi)%body,:) +        & 
                                 normal_plane(:)
! To update the maximum tangential velocity
                              aux_vec(:) = dot_product(bp_arr(npi)%vel,        &
                                           normal_plane) * normal_plane(:)
                              aux_vec(:) = bp_arr(npi)%vel(:) - aux_vec(:)
                              aux_scalar = dsqrt(dot_product(aux_vec,aux_vec))
                              interface_sliding_vel_max(bp_arr(npi)%body) = max&
                                 (interface_sliding_vel_max(bp_arr(npi)%body), &
                                 aux_scalar)
                              else
! Contribution of the normal reaction force under sliding and of the sliding 
! friction force depending on the local slope angle (instead of the friction 
! angle)
                                 aux_gravity(bp_arr(npi)%body,:) = 0.d0
                           endif
                        endif
                        call Vector_Product(bp_arr(npi)%rel_pos,               &
                           f_coll_bp_boun,temp,3)
                        Moment(bp_arr(npi)%body,n_bodies+j,:) =                &
                           Moment(bp_arr(npi)%body,n_bodies+j,:) + temp(:)
                        Force_mag_sum(bp_arr(npi)%body,n_bodies+j) =           &
                           Force_mag_sum(bp_arr(npi)%body,n_bodies+j) +        &
                           (1.d0 / r_per) * k_masses *                         & 
                           Gamma_boun(r_per,Domain%h) 
                        r_per_min(bp_arr(npi)%body,n_bodies+j) =               &
                           min(r_per,r_per_min(bp_arr(npi)%body,n_bodies+j))
                        if (Gamma_boun(r_per,Domain%h)<=0.d0)                  &
                           impact_vel(aux2,n_bodies+j) = 0.d0
                        endif
                     endif
                  endif
               endif
            enddo 
         endif
         endif
      endif
   endif
enddo
! Contributions to the moment of inertia
do npi=1,n_body_part
! No computations for imposed kinematics
   if (body_arr(bp_arr(npi)%body)%imposed_kinematics==0) then
      if (body_arr(bp_arr(npi)%body)%Ic_imposed==0) then
! Contributions to the moment of inertia coefficients (diagonal elements)
         if (ncord==3) then
            body_arr(bp_arr(npi)%body)%Ic(1,1) =                               &
               body_arr(bp_arr(npi)%body)%Ic(1,1) + bp_arr(npi)%mass *         &
               (bp_arr(npi)%rel_pos(2) ** 2 + bp_arr(npi)%rel_pos(3) ** 2)
            body_arr(bp_arr(npi)%body)%Ic(2,2) =                               &
               body_arr(bp_arr(npi)%body)%Ic(2,2) + bp_arr(npi)%mass *         &
               (bp_arr(npi)%rel_pos(1) ** 2 + bp_arr(npi)%rel_pos(3) ** 2)
            body_arr(bp_arr(npi)%body)%Ic(3,3) =                               &
               body_arr(bp_arr(npi)%body)%Ic(3,3) + bp_arr(npi)%mass *         &
               (bp_arr(npi)%rel_pos(1) ** 2 + bp_arr(npi)%rel_pos(2) ** 2)
! Contributions to the products of inertia (off-diagonal elements)
            body_arr(bp_arr(npi)%body)%Ic(1,2) =                               &
               body_arr(bp_arr(npi)%body)%Ic(1,2) - bp_arr(npi)%mass *         &
               (bp_arr(npi)%rel_pos(1) * bp_arr(npi)%rel_pos(2))
            body_arr(bp_arr(npi)%body)%Ic(1,3) =                               &
               body_arr(bp_arr(npi)%body)%Ic(1,3) - bp_arr(npi)%mass *         &
               (bp_arr(npi)%rel_pos(1) * bp_arr(npi)%rel_pos(3))
            body_arr(bp_arr(npi)%body)%Ic(2,3) =                               &
               body_arr(bp_arr(npi)%body)%Ic(2,3) - bp_arr(npi)%mass *         &
               (bp_arr(npi)%rel_pos(2) * bp_arr(npi)%rel_pos(3)) 
         endif
         if ((ncord==2).and.((it_start+1)==on_going_time_step)) then
            body_arr(bp_arr(npi)%body)%Ic(2,2) =                               &
               body_arr(bp_arr(npi)%body)%Ic(2,2) + bp_arr(npi)%mass *         &
               (bp_arr(npi)%rel_pos(1) ** 2 + bp_arr(npi)%rel_pos(3) ** 2)
         endif
      endif 
   endif
enddo     
! Loop over the transported bodies (global contributions from body-body and 
! boundary-body impacts; computation of Ic and its inverse)
! AA!!!test sub for inter_front
!$omp parallel do default(none) private(i,j,k_masses,alfa_boun,aux_int)        &
!$omp shared(n_bodies,body_arr,ncord,it_start,on_going_time_step,aux,r_per_min)&
!$omp shared(Domain,Force_mag_sum,Force,Moment,inter_front)
do i=1,n_bodies
   if (body_arr(i)%imposed_kinematics==0) then
! Forces and torques/moments
      do j=1,aux
         if (r_per_min(i,j)<1000000.d0) then
            if (j<=n_bodies) then
               k_masses = (body_arr(i)%mass * body_arr(j)%mass) /              & 
                          (body_arr(i)%mass + body_arr(j)%mass)
               else
                  k_masses = body_arr(i)%mass 
            endif
            if (Force_mag_sum(i,j)>0.d0) then
               alfa_boun = (Gamma_boun(r_per_min(i,j),Domain%h) /              &
                           r_per_min(i,j) * k_masses) / Force_mag_sum(i,j)
               if (j>n_bodies) alfa_boun = alfa_boun / inter_front(i)
               Force(i,j,:) = Force(i,j,:) * alfa_boun
               body_arr(i)%Force(:) = body_arr(i)%Force(:) + Force(i,j,:)
               Moment(i,j,:) = Moment(i,j,:) * alfa_boun
               body_arr(i)%Moment(:) = body_arr(i)%Moment(:) + Moment(i,j,:)
            endif
         endif
      enddo
! Ic and its inverse 
      if (ncord==3) then
         body_arr(i)%Ic(2,1) = body_arr(i)%Ic(1,2)
         body_arr(i)%Ic(3,1) = body_arr(i)%Ic(1,3)
         body_arr(i)%Ic(3,2) = body_arr(i)%Ic(2,3)
         call Matrix_Inversion_3x3(body_arr(i)%Ic,body_arr(i)%Ic_inv,aux_int)
      endif
      if ((ncord==2).and.((it_start+1)==on_going_time_step)) then
         body_arr(i)%Ic_inv = 0.d0 
         body_arr(i)%Ic_inv(2,2) = 1.d0/body_arr(i)%Ic(2,2)
      endif
   endif
enddo
!$omp end parallel do     
! Check body-body interactions and eventually zeroing impact velocities; 
! adding gravity 
!$omp parallel do default(none)                                                &
!$omp shared(n_bodies,Moment,Force,body_arr,impact_vel,aux,Domain,aux_gravity) &
!$omp shared(n_surf_body_part,surf_body_part,sum_normal_plane,friction_angle)  &
!$omp shared(inter_front,simulation_time,time_max_no_body_gravity_force)       &
!$omp shared(interface_sliding_vel_max,dtvel)                                  &
!$omp private(j,k,nbi,npk,nbk,n_interactions,aux2,aux_scalar,aux_vec)          &
!$omp private(sliding_friction_dir,aux_scalar_2,sliding_friction)              &
!$omp private(friction_limiter)  
do nbi=1,n_bodies
   if (body_arr(nbi)%imposed_kinematics==0) then
      if ((friction_angle>(0.d0-1.d-9)).and.(inter_front(nbi)>0).and.          &
         (simulation_time>time_max_no_body_gravity_force)) then
! Overall normal representing the neighbouring frontiers
         aux_scalar = dsqrt(dot_product(sum_normal_plane(nbi,:),               &
                      sum_normal_plane(nbi,:)))
         if (aux_scalar>1.d-9) then
            sum_normal_plane(nbi,:) = sum_normal_plane(nbi,:) / aux_scalar
            else
               sum_normal_plane(nbi,:) = 0.d0
         endif
! Gravity component normal to the overall normal representing the neighbouring 
! frontiers
         aux_vec(:) = dot_product(aux_gravity(nbi,:),sum_normal_plane(nbi,:))  &
                      * sum_normal_plane(nbi,:)
! Absolute value of the gravity component normal to the overall surface
         aux_scalar = dsqrt(dot_product(aux_vec,aux_vec))
! Gravity + normal reaction force under sliding (aligned with the overall normal
! of the neighbouring frontiers)
         aux_gravity(nbi,:) = aux_gravity(nbi,:) - aux_vec(:)
! Direction of the sliding friction force
         aux_vec(:) = dot_product(body_arr(nbi)%u_CM,sum_normal_plane(nbi,:)) *&
                      sum_normal_plane(nbi,:)
         sliding_friction_dir(:) = body_arr(nbi)%u_CM(:) - aux_vec(:)
         aux_scalar_2 = dsqrt(dot_product(sliding_friction_dir,                &
                        sliding_friction_dir))
         if (aux_scalar_2>1.d-9) then
            sliding_friction_dir(:) = - sliding_friction_dir(:) / aux_scalar_2
            else
               sliding_friction_dir(:) = 0.d0
         endif
! Sliding friction force: first estimation
         sliding_friction(:) = aux_scalar * dtan(friction_angle) *             &
                               sliding_friction_dir(:)
! Limiter to the sliding friction force
         aux_scalar = dsqrt(dot_product(sliding_friction(:),                   &
                      sliding_friction(:)))
         if (dtvel>1.d-9) then
            friction_limiter = interface_sliding_vel_max(nbi) / dtvel *        &
                               body_arr(nbi)%mass
            else
               friction_limiter = 0.d0
         endif
         aux_scalar_2 = min(aux_scalar,friction_limiter)
         if (aux_scalar>1.d-9) then
            sliding_friction(:) = sliding_friction(:) * aux_scalar_2 /         &
                                  aux_scalar
         endif
! Contribution of the sliding friction force to the global force
         body_arr(nbi)%Force(:) = body_arr(nbi)%Force(:) + sliding_friction(:)
      endif
! To update the resultant force with gravity + normal reaction force under 
! sliding
      body_arr(nbi)%Force(:) = body_arr(nbi)%Force(:) + aux_gravity(nbi,:)
!impact velocity
      n_interactions = 0
      do j=1,aux
         if (nbi==j) cycle
         if ((Moment(nbi,j,1)==0.d0).and.(Force(nbi,j,1)==0.d0).and.           &
             (Moment(nbi,j,2)==0.d0).and.(Force(nbi,j,2)==0.d0).and.           &
             (Moment(nbi,j,3)==0.d0).and.(Force(nbi,j,3)==0.d0)) then
            if (j<=n_bodies) then
               k = 0
               nbk = 1
               do while (nbk.ne.nbi)
                  k = k + body_arr(nbk)%npart
                  nbk = nbk + 1
               enddo
               do npk=(k+1),(k+body_arr(nbi)%npart)
                  do aux2=1,n_surf_body_part
                     if (surf_body_part(aux2)==npk) then
                        impact_vel(aux2,j) = zero
                        cycle
                     endif
                  enddo
               enddo
            endif
            else
               n_interactions = n_interactions + 1 
         endif
      enddo
   endif
enddo
!$omp end parallel do     
!------------------------
! Deallocations
!------------------------
deallocate(Force)
deallocate(Moment)
deallocate(Force_mag_sum)
deallocate(r_per_min)
deallocate(aux_gravity)
deallocate(sum_normal_plane)
! AA!!! test
deallocate(inter_front)
deallocate(interface_sliding_vel_max)
return
end subroutine RHS_body_dynamics

