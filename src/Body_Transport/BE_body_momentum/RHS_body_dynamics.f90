!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: RHS_body_dynamics  
! Description:  To estimate the RHS of the body dynamics equations (Amicarelli 
!               et al., 2015, CAF; Amicarelli et al., 2020, CPC; Amicarelli et 
!               al., 2022, IJCFD)
!-------------------------------------------------------------------------------
#ifdef SOLID_BODIES
subroutine RHS_body_dynamics(dtvel)
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
double precision,intent(in) :: dtvel
integer(4) :: npartint,i,j,npi,npj,aux,nbi,npk,k,nbk,alloc_stat
integer(4) :: n_interactions,aux2,aux3,test,aux_scal_test
#ifdef SPACE_3D
integer(4) :: aux_int
#endif
double precision :: k_masses,r_per,r_par,alfa_boun,dxbp,abs_det_thresh
double precision :: aux_impact_vel,aux4,aux_scalar,aux_scalar_2,Sum_W_vol,W_vol
double precision :: friction_limiter,pres_mir
#ifdef SPACE_2D
double precision :: aux_locx_max,aux_locx_min
#endif
double precision :: f_pres(3),r_par_vec(3),f_coll_bp_bp(3)          
double precision :: f_coll_bp_boun(3),normal_plane(3)
double precision :: u_rel(3),x_rel(3),aux_vec(3),aux_vec_2(3),rel_pos(3)
double precision :: loc_pos(3)
#ifdef SPACE_2D
double precision :: aux_locx_vert(3),pos_aux(3)
#endif
double precision :: sliding_friction(3),fluid_vel_npj(3)
double precision :: aux_mat(3,3)
! "inter_front": Number of neighbouring frontiers for a given body. It is 
! approximated by the maximum number of neighbouring frontiers of a single body 
! particle belonging to the body.
integer(4),dimension(:),allocatable :: inter_front
! Array of the number of "body particle - frontier" interactions for a given 
! body
integer(4),dimension(:),allocatable :: bp_bound_interactions
double precision,dimension(:),allocatable :: interface_sliding_vel_max
! Array of the hydrodynamic forces
double precision,dimension(:,:),allocatable :: Force_bod_flu
! "alfa_denom": denominator of the normalizing parameter alfa for body-body and 
! body-bundary interactions
double precision,dimension(:,:),allocatable :: alfa_denom
double precision,dimension(:,:),allocatable :: r_per_min
double precision,dimension(:,:),allocatable :: aux_gravity
! Average of the normal vectors of the neighbouring boundaries of a given body
double precision,dimension(:,:),allocatable :: mean_bound_normal
! Sliding force application point for any body. It is computed as the average 
! position of the body particles in the "body particle - frontier" interactions.
! Each interaction has the same weight in the averaging process (every particle
! might be considered more than once).
double precision,dimension(:,:),allocatable :: sliding_app_point
! Sliding direction. It is computed as the direction of the vector sum of the 
! sliding velocity of the body particles in the sliding region (every particle 
! might be considered more than once).
double precision,dimension(:,:),allocatable :: sliding_dir
! Array of the body-body and body-boundary forces
double precision,dimension(:,:,:),allocatable :: Force_bod_sol
! Array of the body-body and body-boundary torques
double precision,dimension(:,:,:),allocatable :: Moment_bod_sol
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
interface
   subroutine body_boundary_for_sliding_friction_normal_reaction(i_bp,         &
      bp_bound_interactions,normal_plane,bp_pos,interface_sliding_vel_max,     &
      mean_bound_normal,sliding_app_point,sliding_dir,aux_gravity)
      implicit none
      integer(4),intent(in) :: i_bp
      integer(4),intent(inout) :: bp_bound_interactions
      double precision,intent(in) :: normal_plane(3)
      double precision,intent(in) :: bp_pos(3)
      double precision,intent(inout) :: interface_sliding_vel_max
      double precision,intent(inout) :: mean_bound_normal(3)
      double precision,intent(inout) :: sliding_app_point(3)
      double precision,intent(inout) :: sliding_dir(3)
      double precision,intent(out) :: aux_gravity(3)
      double precision :: aux_scalar
      double precision :: aux_vec(3)
   end subroutine
end interface
interface
   subroutine Vector_Product(uu,VV,ww,SPACEDIM)
      implicit none
      integer(4),intent(in) :: SPACEDIM
      double precision,intent(in),dimension(SPACEDIM) :: uu,VV
      double precision,intent(inout),dimension(SPACEDIM) :: ww
      integer(4) :: i,j,k
      integer(4),dimension(3) :: iseg=(/2,3,1/)
   end subroutine
   subroutine MatrixProduct(AA,BB,CC,nr,nrc,nc)
      implicit none
      integer(4),intent(in) :: nr,nrc,nc
      double precision,intent(in),dimension(nr,nrc) :: AA
      double precision,intent(in),dimension(nrc,nc) :: BB
      double precision,intent(inout),dimension(nr,nc) :: CC
   end subroutine MatrixProduct
end interface
!------------------------
! Allocations
!------------------------
#ifdef SPACE_3D
   aux = n_bodies + NumFacce
#elif defined SPACE_2D
      aux = n_bodies + NumBSides
#endif
if (.not.allocated(Force_bod_flu)) then
   allocate(Force_bod_flu(n_bodies,3),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Allocation of "Force_bod_flu" failed in the program ',    &
         'unit "RHS_body_dynamics"; the execution terminates here.'
      stop
   endif
endif
if (.not.allocated(Force_bod_sol)) then
   allocate(Force_bod_sol(n_bodies,aux,3),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Allocation of "Force_bod_sol" failed in the program ',    &
         'unit "RHS_body_dynamics"; the execution terminates here.'
      stop
   endif
endif
if (.not.allocated(Moment_bod_sol)) then
   allocate(Moment_bod_sol(n_bodies,aux,3),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Allocation of "Moment_bod_sol" failed in the program ',   &
         'unit "RHS_body_dynamics"; the execution terminates here.'
      stop
   endif
endif
if (.not.allocated(alfa_denom)) then
   allocate(alfa_denom(n_bodies,aux),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Allocation of "alfa_denom" failed in the program unit ',  &
         '"RHS_body_dynamics"; the execution terminates here.'
      stop
   endif
endif
if (.not.allocated(r_per_min)) then
   allocate(r_per_min(n_bodies,aux),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Allocation of "r_per_min" failed in the program unit ',   &
         '"RHS_body_dynamics"; the execution terminates here.'
      stop
   endif
endif
if (.not.allocated(aux_gravity)) then
   allocate(aux_gravity(n_bodies,3),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Allocation of "aux_gravity" failed in the program unit ', &
         '"RHS_body_dynamics"; the execution terminates here.'
      stop
   endif
endif
if (.not.allocated(mean_bound_normal)) then
   allocate(mean_bound_normal(n_bodies,3),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Allocation of "mean_bound_normal" failed in the program ',&
         'unit "RHS_body_dynamics"; the execution terminates here.'
      stop
   endif
endif
if (.not.allocated(sliding_app_point)) then
   allocate(sliding_app_point(n_bodies,3),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Allocation of "sliding_app_point" failed in the program ',&
         'unit "RHS_body_dynamics"; the execution terminates here.'
      stop
   endif
endif
if (.not.allocated(sliding_dir)) then
   allocate(sliding_dir(n_bodies,3),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Allocation of "sliding_dir" failed in the program ',      &
         'unit "RHS_body_dynamics"; the execution terminates here.'
      stop
   endif
endif
if (.not.allocated(inter_front)) then
   allocate(inter_front(n_bodies),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Allocation of "inter_front" failed in the program ',      &
         'unit "RHS_body_dynamics"; the execution terminates here.'
      stop
   endif
endif
if (.not.allocated(bp_bound_interactions)) then
   allocate(bp_bound_interactions(n_bodies),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Allocation of "bp_bound_interactions" failed in the ',    &
         'program nit "RHS_body_dynamics"; the execution terminates here.'
      stop
   endif
endif
if (.not.allocated(interface_sliding_vel_max)) then
   allocate(interface_sliding_vel_max(n_bodies),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Allocation of "interface_sliding_vel_max" failed in the', &
         'program unit "RHS_body_dynamics"; the execution terminates here.'
      stop
   endif
endif
!------------------------
! Initializations
!------------------------
Force_bod_flu = 0.d0
Force_bod_sol = 0.d0
Moment_bod_sol = 0.d0
alfa_denom = 0.d0
r_per_min = 1.d6
aux2 = 0
inter_front(:) = 0
bp_bound_interactions(:) = 0
aux_scal_test = 0
mean_bound_normal(:,:) = 0.d0
sliding_app_point(:,:) = 0.d0
sliding_dir(:,:) = 0.d0
interface_sliding_vel_max(:) = 0.d0
abs_det_thresh = 1.d-9
!------------------------
! Statements
!------------------------
call body_p_max_limiter
! Contributions to fluid dynamics momentum equations
!$omp parallel do default(none)                                                &
!$omp shared(n_body_part,bp_arr,nPartIntorno_bp_f,NMAXPARTJ,PartIntorno_bp_f)  &
!$omp shared(pg,thin_walls,input_any_t,rag_bp_f,KerDer_bp_f_Gal)               &
!$omp shared(FSI_slip_conditions,body_minimum_pressure_limiter,body_arr)       &
!$omp shared(body_maximum_pressure_limiter,KerDer_bp_f_cub_spl)                &
!$omp private(npi,Sum_W_vol,j,npartint,npj,pres_mir,W_vol,aux_scalar)          &
!$omp private(aux_scalar_2,aux_vec,fluid_vel_npj)
do npi=1,n_body_part
   bp_arr(npi)%pres = 0.d0
   Sum_W_vol = 0.d0
   do j=1,nPartIntorno_bp_f(npi)
      npartint = (npi - 1) * NMAXPARTJ + j
      npj = PartIntorno_bp_f(npartint)
! Contribution to the "grad_p + ALE term": start
! Term for grad_p
      call body_pressure_mirror_interaction(npi,npj,npartint,pres_mir,W_vol)
      bp_arr(npi)%pres = bp_arr(npi)%pres + pres_mir * W_vol
      Sum_W_vol = Sum_W_vol + W_vol
      aux_scalar = (pres_mir - pg(npj)%pres) / (pg(npj)%dens ** 2)
! Term for ALE
      aux_scalar_2 = 2.d0 * pg(npj)%pres / (pg(npj)%dens * pg(npj)%dens)
      if (thin_walls) then
! Treatment for thin walls (to the coupling term on the pressure gradient)
         aux_scalar = aux_scalar * (1.d0 + (1.d0 - pg(npj)%sigma_fp -          &
                      pg(npj)%sigma_bp) / pg(npj)%sigma_bp)
         aux_scalar_2 = aux_scalar_2 * (1.d0 + (1.d0 - pg(npj)%sigma_fp -      &
                      pg(npj)%sigma_bp) / pg(npj)%sigma_bp)
      endif
      if (input_any_t%C1_BE) then
! grad_p (renormalization)
         call MatrixProduct(pg(npj)%B_ren_gradp,BB=rag_bp_f(1:3,npartint),     &
            CC=aux_vec,nr=3,nrc=3,nc=1)
!$omp critical (omp_RHS_bd_acc)
         pg(npj)%acc(1:3) = pg(npj)%acc(1:3) + pg(npj)%dens *                  &
                            bp_arr(npi)%volume * aux_scalar * (-aux_vec(1:3))  &
                            * KerDer_bp_f_Gal(npartint)
!$omp end critical (omp_RHS_bd_acc)
         else
! grad_p (no renormalization)
!$omp critical (omp_RHS_bd_acc)
            pg(npj)%acc(1:3) = pg(npj)%acc(1:3) + (-pg(npj)%dens *             &
                               bp_arr(npi)%volume * aux_scalar *               &
                               ( - rag_bp_f(1:3,npartint)) *                   &
                               KerDer_bp_f_Gal(npartint))
!$omp end critical (omp_RHS_bd_acc)
      endif
! ALE contribution to the acceleration
      if (.not.((pg(npj)%p0_neg_ALE).and.(pg(npj)%B_ren_gradp_stat==1))) then
         aux_vec(1:3) = -pg(npj)%dens * bp_arr(npi)%volume * aux_scalar_2 *    &
                        (-rag_bp_f(1:3,npartint)) * KerDer_bp_f_Gal(npartint)
!$omp critical (omp_RHS_bd_acc)
         pg(npj)%acc(1:3) = pg(npj)%acc(1:3) + aux_vec(1:3)
! Contribution to the ALE velocity increment (here it is still an acceleration)
         pg(npj)%dvel_ALE1(1:3) = pg(npj)%dvel_ALE1(1:3) + aux_vec(1:3)
!$omp end critical (omp_RHS_bd_acc)
      endif
! Contribution to the "grad_p + ALE term": end
      if ((FSI_slip_conditions==1).or.(FSI_slip_conditions==3)) then
! Body particle volume
         aux_scalar = bp_arr(npi)%volume
         if (thin_walls) then
! Treatment for thin walls (to the coupling term on the shear-stress gradient)
            aux_scalar = aux_scalar * (1.d0 + (1.d0 - pg(npj)%sigma_fp -       &
                         pg(npj)%sigma_bp) / pg(npj)%sigma_bp)
         endif
! Contribution to the shear stress gradient term
         if ((input_any_t%ALE3).and.(.not.(pg(npj)%p0_neg_ALE))) then
            fluid_vel_npj(1:3) = pg(npj)%vel_fluid(1:3)
            else
               fluid_vel_npj(1:3) = pg(npj)%vel(1:3)            
         endif
!$omp critical (omp_RHS_bd_acc)
! Only slip condition available: mirror fluid velocity as solid velocity
         pg(npj)%acc(:) = pg(npj)%acc(:) - 2.d0 * pg(npj)%kin_visc *           &
                          (bp_arr(npi)%vel(1:3) - fluid_vel_npj(1:3)) *        &
                          KerDer_bp_f_cub_spl(npartint) * aux_scalar
!$omp end critical (omp_RHS_bd_acc)
      endif
   enddo
! Unique value representative of the body-particle pressure (for body dynamics)  
   if (Sum_W_vol>1.d-3) bp_arr(npi)%pres = bp_arr(npi)%pres / Sum_W_vol
! Body-particle pressure limiters
   if (body_minimum_pressure_limiter) then
      if (bp_arr(npi)%pres<0.d0) bp_arr(npi)%pres = 0.d0 
   endif
   if (body_maximum_pressure_limiter) then
      if (bp_arr(npi)%pres>body_arr(bp_arr(npi)%body)%p_max_limiter) then
         bp_arr(npi)%pres = body_arr(bp_arr(npi)%body)%p_max_limiter
      endif
   endif
enddo
!$omp end parallel do
! Loop over the transported bodies (gravity forces and Ic initialization)
!$omp parallel do default(none) private(i)                                     &
!$omp shared(n_bodies,body_arr,Domain,it_start,on_going_time_step)             &
!$omp shared(aux_gravity,simulation_time,time_max_no_body_gravity_force)
do i=1,n_bodies
   if ((body_arr(i)%imposed_kinematics==0).or.                                 &
      (body_arr(i)%imposed_kinematics==2))  then
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
   if ((body_arr(bp_arr(npi)%body)%imposed_kinematics==0).or.                  &
      (body_arr(bp_arr(npi)%body)%imposed_kinematics==2)) then
! Only body particles at the body surface     
      if (bp_arr(npi)%surface) then  
! Computation of the ID of the surface body particles
         aux2 = aux2 + 1
! Fluid pressure forces
         f_pres(:) = bp_arr(npi)%pres * bp_arr(npi)%area * bp_arr(npi)%normal(:)
         body_arr(bp_arr(npi)%body)%Force(:) =                                 &
            body_arr(bp_arr(npi)%body)%Force(:) + f_pres(:)
         call Vector_Product(bp_arr(npi)%rel_pos(:),f_pres(:),aux_vec(:),3)
         body_arr(bp_arr(npi)%body)%Moment(:) =                                &
            body_arr(bp_arr(npi)%body)%Moment(:) + aux_vec(:)
! Loop over the neighbouring body particles, belonging to other bodies 
! (inter-body impacts, partial contributions) 
         do j=1,nPartIntorno_bp_bp(aux2)
            npartint = (aux2 - 1) * NMAXPARTJ + j
            npj = PartIntorno_bp_bp(npartint)
            r_per = abs(rag_bp_bp(1,npartint) * bp_arr(npj)%normal(1) +        &
                    rag_bp_bp(2,npartint) * bp_arr(npj)%normal(2) +            &
                    rag_bp_bp(3,npartint) * bp_arr(npj)%normal(3))
! rag=-r
            r_par_vec(:) = - rag_bp_bp(:,npartint) - r_per *                   &
                           bp_arr(npj)%normal(:)  
            r_par = dsqrt(r_par_vec(1) * r_par_vec(1) + r_par_vec(2) *         &
                    r_par_vec(2) + r_par_vec(3) * r_par_vec(3))
#ifdef SPACE_3D
            dxbp = bp_arr(npi)%volume ** (1.d0/3.d0)
#elif defined SPACE_2D
            dxbp = bp_arr(npi)%volume ** (1.d0/2.d0)
#endif
            if ((r_per>0.d0).and.((r_par/dxbp)<=1.d0)) then
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
                                    dxbp) * bp_arr(npj)%normal(:)
                  if (r_per<(dxbp/10.d0)) then   
                     write(file_name_test,"(a,a,i8.8,a,i8.8,a,i8.8,a)")        &
                        nomecaso(1:len_trim(nomecaso)),'_close_impact_',       &
                        on_going_time_step,'_',npi,'_',npj,".txt"
                     open(60,file=file_name_test,status="unknown",             &
                        form="formatted")
                     write(60,                                                &
              '((10x,a),4(11x,a),(2x,a),(9x,a),(6x,a),(4x,a),(9x,a),3(11x,a))')&
                       " time"," npi"," npj"," nbi"," nbj"," imp_vel(m/s)",    &
                       " r_per"," k_masses"," Gamma_boun"," r_par",            &
                       " n_x"," n_y"," n_z"
                     write(60,'((g14.7,1x),4(i14,1x),8(g14.7,1x))')           &
                        simulation_time,npi,npj,bp_arr(npi)%body,              &
                        bp_arr(npj)%body,impact_vel(aux2,bp_arr(npj)%body),    &
                        r_per,k_masses,Gamma_boun(r_per,Domain%h),r_par,       &
                        bp_arr(npj)%normal(1),bp_arr(npj)%normal(2),           &
                        bp_arr(npj)%normal(3)
                     close(60)
                  endif
                  Force_bod_sol(bp_arr(npi)%body,bp_arr(npj)%body,:) =         &
                    Force_bod_sol(bp_arr(npi)%body,bp_arr(npj)%body,:) +       &
                    f_coll_bp_bp(:)
! Body-body interactions: in the absence of neighburing fluid particles (dry 
! stage), sliding friction force depends on the interaction interface angle 
! (instead of the friction angle) and the normal reaction force under 
! sliding is correct. These forces exactly balance the gravity components.
                  if (body_arr(bp_arr(npi)%body)%pmax<1.d-5) then
                     aux_gravity(bp_arr(npi)%body,:) = 0.d0
                  endif
                  call Vector_Product(bp_arr(npi)%rel_pos,f_coll_bp_bp,        &
                     aux_vec,3)
                  Moment_bod_sol(bp_arr(npi)%body,bp_arr(npj)%body,:) =        &
                     Moment_bod_sol(bp_arr(npi)%body,bp_arr(npj)%body,:) +     &
                     aux_vec(:)
                  alfa_denom(bp_arr(npi)%body,bp_arr(npj)%body) =              &
                     alfa_denom(bp_arr(npi)%body,bp_arr(npj)%body) +           &
                     (1.d0 / r_per) * k_masses * Gamma_boun(r_per,Domain%h) *  &
                     (1.d0 - r_par / dxbp)
                  r_per_min(bp_arr(npi)%body,bp_arr(npj)%body) =               &
                     min(r_per,r_per_min(bp_arr(npi)%body,bp_arr(npj)%body)) 
               endif
            endif
         enddo
! Loop over boundaries (body-boundary impacts, partial contributions) (3D case) 
         if (simulation_time>time_max_no_body_frontier_impingements) then
#ifdef SPACE_3D
            do j=1,NumFacce
            if ((Tratto(BoundaryFace(j)%stretch)%tipo=="fixe") .or.            &
               (Tratto(BoundaryFace(j)%stretch)%tipo == "tapi")) then
               aux_vec(:) = bp_arr(npi)%pos(:) - BoundaryFace(j)%Node(1)%GX(:)
               aux4 = dot_product(BoundaryFace(j)%T(:,3),aux_vec)
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
                        Force_bod_sol(bp_arr(npi)%body,n_bodies+j,:) =         &
                           Force_bod_sol(bp_arr(npi)%body,n_bodies+j,:) +      &
                           f_coll_bp_boun(:)
                        call body_boundary_for_sliding_friction_normal_reaction&
                           (npi,bp_bound_interactions(bp_arr(npi)%body),       &
                           normal_plane(:),bp_arr(npi)%pos(:),                 &
                           interface_sliding_vel_max(bp_arr(npi)%body),        &
                           mean_bound_normal(bp_arr(npi)%body,:),              &
                           sliding_app_point(bp_arr(npi)%body,:),              &
                           sliding_dir(bp_arr(npi)%body,:),                    &
                           aux_gravity(bp_arr(npi)%body,:))
                        call Vector_Product(bp_arr(npi)%rel_pos,               &
                           f_coll_bp_boun,aux_vec,3)
                        Moment_bod_sol(bp_arr(npi)%body,n_bodies+j,:) =        &
                           Moment_bod_sol(bp_arr(npi)%body,n_bodies+j,:) +     &
                           aux_vec(:)
                        alfa_denom(bp_arr(npi)%body,n_bodies+j) =              &
                           alfa_denom(bp_arr(npi)%body,n_bodies+j) +           &
                           (1.d0 / r_per) * k_masses *                         &
                           Gamma_boun(r_per,Domain%h) 
                        r_per_min(bp_arr(npi)%body,n_bodies+j) =               &
                           min(r_per,r_per_min(bp_arr(npi)%body,n_bodies+j))
                        if (Gamma_boun(r_per,Domain%h)<=0.d0)                  &
                           impact_vel(aux2,n_bodies+j) = 0.d0
                        aux_scal_test = aux_scal_test + 1
                        endif
                     endif
                  endif  
               endif
            enddo
#endif
         inter_front(bp_arr(npi)%body) = max(inter_front(bp_arr(npi)%body),    &
                                         aux_scal_test)
         aux_scal_test = 0
! 2D case
#ifdef SPACE_2D
            do j=1,NumBSides
               if ((BoundarySide(j)%tipo=="fixe") .or.                         &
                  (BoundarySide(j)%tipo == "tapi")) then  
                  aux_vec(:) = bp_arr(npi)%pos(:) -                           &
                                Vertice(1,BoundarySide(j)%Vertex(1))
                  aux4 = dot_product(BoundarySide(j)%T(:,3),aux_vec)
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
                        Force_bod_sol(bp_arr(npi)%body,n_bodies+j,:) =         &
                           Force_bod_sol(bp_arr(npi)%body,n_bodies+j,:) +      &
                           f_coll_bp_boun(:)
                        call body_boundary_for_sliding_friction_normal_reaction&
                           (npi,bp_bound_interactions(bp_arr(npi)%body),       &
                           normal_plane(:),bp_arr(npi)%pos(:),                 &
                           interface_sliding_vel_max(bp_arr(npi)%body),        &
                           mean_bound_normal(bp_arr(npi)%body,:),              &
                           sliding_app_point(bp_arr(npi)%body,:),              &
                           sliding_dir(bp_arr(npi)%body,:),                    &
                           aux_gravity(bp_arr(npi)%body,:))
                        call Vector_Product(bp_arr(npi)%rel_pos,               &
                           f_coll_bp_boun,aux_vec,3)
                        Moment_bod_sol(bp_arr(npi)%body,n_bodies+j,:) =        &
                           Moment_bod_sol(bp_arr(npi)%body,n_bodies+j,:) +     &
                           aux_vec(:)
                        alfa_denom(bp_arr(npi)%body,n_bodies+j) =              &
                           alfa_denom(bp_arr(npi)%body,n_bodies+j) +           &
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
#endif
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
#ifdef SPACE_3D
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
#elif SPACE_2D
               if ((it_start+1)==on_going_time_step) then
                  body_arr(bp_arr(npi)%body)%Ic(2,2) =                         &
                     body_arr(bp_arr(npi)%body)%Ic(2,2) + bp_arr(npi)%mass *   &
                     (bp_arr(npi)%rel_pos(1) ** 2 + bp_arr(npi)%rel_pos(3) ** 2)
               endif
#endif
      endif 
   endif
enddo     
! Loop over the transported bodies (global contributions from body-body and 
! boundary-body interactions; computation of Ic and its inverse)
!$omp parallel do default(none) private(i,j,k_masses,alfa_boun)                &
#ifdef SPACE_3D
!$omp private(aux_int)                                                         &
#endif
!$omp shared(n_bodies,body_arr,it_start,on_going_time_step,aux,r_per_min)      &
!$omp shared(Domain,alfa_denom,Force_bod_sol,Moment_bod_sol,inter_front)       &
!$omp shared(Force_bod_flu,abs_det_thresh)
do i=1,n_bodies
   if ((body_arr(i)%imposed_kinematics==0).or.                                 &
      (body_arr(i)%imposed_kinematics==2)) then
! Saving the hydrodynamic forces 
      Force_bod_flu(i,:) = body_arr(i)%Force(:)
! Finalizing the assessment of the body-solid forces 
      do j=1,aux
         if (r_per_min(i,j)<1.d6) then
            if (j<=n_bodies) then
               k_masses = (body_arr(i)%mass * body_arr(j)%mass) /              & 
                          (body_arr(i)%mass + body_arr(j)%mass)
               else
                  k_masses = body_arr(i)%mass 
            endif
            if (alfa_denom(i,j)>0.d0) then
               alfa_boun = (Gamma_boun(r_per_min(i,j),Domain%h) /              &
                           r_per_min(i,j) * k_masses) / alfa_denom(i,j)
! The body-boundary forces are normalized by the number of neighbouring 
! frontiers
               if (j>n_bodies) alfa_boun = alfa_boun / inter_front(i)
! Body-solid forces and torques/moments
               Force_bod_sol(i,j,:) = Force_bod_sol(i,j,:) * alfa_boun
               body_arr(i)%Force(:) = body_arr(i)%Force(:) +                   &
                                      Force_bod_sol(i,j,:)
               Moment_bod_sol(i,j,:) = Moment_bod_sol(i,j,:) * alfa_boun
               body_arr(i)%Moment(:) = body_arr(i)%Moment(:) +                 &
                                       Moment_bod_sol(i,j,:)
            endif
         endif
      enddo
   endif
   if (body_arr(i)%imposed_kinematics==0) then
! Ic and its inverse
#ifdef SPACE_3D
         body_arr(i)%Ic(2,1) = body_arr(i)%Ic(1,2)
         body_arr(i)%Ic(3,1) = body_arr(i)%Ic(1,3)
         body_arr(i)%Ic(3,2) = body_arr(i)%Ic(2,3)
         call Matrix_Inversion_3x3(body_arr(i)%Ic,body_arr(i)%Ic_inv,          &
            abs_det_thresh,aux_int)
#elif defined SPACE_2D
            if ((it_start+1)==on_going_time_step) then
               body_arr(i)%Ic_inv(:,:) = 0.d0 
               body_arr(i)%Ic_inv(2,2) = 1.d0 / body_arr(i)%Ic(2,2)
            endif
#endif
   endif
enddo
!$omp end parallel do     
! Finalizing body-body interactions and body-boundary interactions; possible
! zeroing of the impact velocities; adding gravity. 
!$omp parallel do default(none)                                                &
!$omp shared(n_bodies,Moment_bod_sol,Force_bod_sol,body_arr,impact_vel,aux)    &
!$omp shared(n_surf_body_part,surf_body_part,mean_bound_normal,friction_angle) &
!$omp shared(inter_front,simulation_time,time_max_no_body_gravity_force)       &
!$omp shared(interface_sliding_vel_max,dtvel,aux_gravity,Domain,Force_bod_flu) &
!$omp shared(sliding_app_point,sliding_dir,bp_bound_interactions)              &
!$omp private(j,k,nbi,npk,nbk,n_interactions,aux2,aux_scalar,aux_vec,aux_vec_2)&
!$omp private(aux_scalar_2,sliding_friction,friction_limiter,rel_pos)
do nbi=1,n_bodies
   if ((body_arr(nbi)%imposed_kinematics==0).or.                               &
      (body_arr(nbi)%imposed_kinematics==2)) then
! Sliding friction, gravity and normal reaction: start
      if ((friction_angle>-1.d-9).and.(inter_front(nbi)>0).and.                &
         (simulation_time>time_max_no_body_gravity_force)) then
! Overall normal representing the neighbouring frontiers
         aux_scalar = dsqrt(dot_product(mean_bound_normal(nbi,:),              &
                      mean_bound_normal(nbi,:)))
         if (aux_scalar>1.d-9) then
            mean_bound_normal(nbi,:) = mean_bound_normal(nbi,:) / aux_scalar
            else
               mean_bound_normal(nbi,:) = 0.d0
         endif
! Normal force for sliding friction: start
! Gravity component normal to the overall normal representing the neighbouring 
! frontiers
         aux_vec(:) = dot_product(aux_gravity(nbi,:),mean_bound_normal(nbi,:)) &
                      * mean_bound_normal(nbi,:)
         if (body_arr(nbi)%pmax<1.d-5) then
! In case of dry body, normal boundary reaction under sliding is here 
! explicitly considered. Gravity force is replaced by the vector sum of the 
! gravity force and the normal reaction.
            aux_gravity(nbi,:) = aux_gravity(nbi,:) - aux_vec(:)
         endif
! Hydrodynamic force component normal to the overall normal representing the 
! neighbouring frontiers
         aux_vec_2(:) = dot_product(Force_bod_flu(nbi,:),                      &
                        mean_bound_normal(nbi,:)) * mean_bound_normal(nbi,:)
! Vector sum of the gravity (and normal reaction force under sliding) and the 
! hydrodynamic force, all aligned with the overall boundary normal.
         aux_vec(:) = aux_vec(:) + aux_vec_2(:)
! Absolute value of the normal force above
         aux_scalar = dsqrt(dot_product(aux_vec,aux_vec))
! Normal force for sliding friction: end
! Sliding direction
         aux_scalar_2 = dsqrt(dot_product(sliding_dir(nbi,:),sliding_dir(nbi,:)))
         if (aux_scalar_2>1.d-9) then
            sliding_dir(nbi,:) = - sliding_dir(nbi,:) / aux_scalar_2
            else
               sliding_dir(nbi,:) = 0.d0
         endif
! Sliding friction force: first estimation
         sliding_friction(:) = aux_scalar * dtan(friction_angle) *             &
                               sliding_dir(nbi,:)
! Limiter to the sliding friction force
         aux_scalar = dsqrt(dot_product(sliding_friction(:),                   &
                      sliding_friction(:)))
         if (dtvel>1.d-9) then
            friction_limiter = 2.d0 * interface_sliding_vel_max(nbi) / dtvel * &
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
! Contribution of the sliding friction to the global momentum
         sliding_app_point(nbi,:) = sliding_app_point(nbi,:) /                 &
                                    bp_bound_interactions(nbi)
         rel_pos(:) = sliding_app_point(nbi,:) - body_arr(nbi)%x_CM(:)
         call Vector_Product(rel_pos(:),sliding_friction(:),aux_vec(:),3)
         body_arr(nbi)%Moment(:) = body_arr(nbi)%Moment(:) + aux_vec(:)
      endif
! To update the resultant force with gravity + normal reaction force under 
! sliding
      body_arr(nbi)%Force(:) = body_arr(nbi)%Force(:) + aux_gravity(nbi,:)
! Sliding friction, gravity and normal reaction: end
! Impact velocity
      n_interactions = 0
      do j=1,aux
         if (nbi==j) cycle
         if ((Moment_bod_sol(nbi,j,1)==0.d0).and.                              &
             (Force_bod_sol(nbi,j,1)==0.d0).and.                               &
             (Moment_bod_sol(nbi,j,2)==0.d0).and.                              &
             (Force_bod_sol(nbi,j,2)==0.d0).and.                               &
             (Moment_bod_sol(nbi,j,3)==0.d0).and.                              &
             (Force_bod_sol(nbi,j,3)==0.d0)) then
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
if(allocated(Force_bod_flu)) then
   deallocate(Force_bod_flu,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Deallocation of "Force_bod_flu" in the program unit  ',   &
         '"RHS_body_dynamics" failed; the execution terminates here. '
      stop
   endif
endif
if(allocated(Force_bod_sol)) then
   deallocate(Force_bod_sol,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Deallocation of "Force_bod_sol" in the program unit  ',   &
         '"RHS_body_dynamics" failed; the execution terminates here. '
      stop
   endif
endif
if(allocated(Moment_bod_sol)) then
   deallocate(Moment_bod_sol,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Deallocation of "Moment_bod_sol" in the program unit  ',  &
         '"RHS_body_dynamics" failed; the execution terminates here. '
      stop
   endif
endif
if(allocated(alfa_denom)) then
   deallocate(alfa_denom,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Deallocation of "alfa_denom" in the program unit  ',      &
         '"RHS_body_dynamics" failed; the execution terminates here. '
      stop
   endif
endif
if(allocated(r_per_min)) then
   deallocate(r_per_min,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Deallocation of "r_per_min" in the program unit  ',       &
         '"RHS_body_dynamics.f90" failed; the execution terminates here. '
      stop
   endif
endif
if(allocated(aux_gravity)) then
   deallocate(aux_gravity,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Deallocation of "aux_gravity" in the program unit  ',     &
         '"RHS_body_dynamics" failed; the execution terminates here. '
      stop
   endif
endif
if(allocated(mean_bound_normal)) then
   deallocate(mean_bound_normal,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Deallocation of "mean_bound_normal" in the program unit ',&
         '"RHS_body_dynamics" failed; the execution terminates here. '
      stop
   endif
endif
if(allocated(sliding_app_point)) then
   deallocate(sliding_app_point,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Deallocation of "sliding_app_point" in the program unit ',&
         '"RHS_body_dynamics" failed; the execution terminates here. '
      stop
   endif
endif
if(allocated(sliding_dir)) then
   deallocate(sliding_dir,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Deallocation of "sliding_dir" in the program unit ',      &
         '"RHS_body_dynamics" failed; the execution terminates here. '
      stop
   endif
endif
if(allocated(inter_front)) then
   deallocate(inter_front,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Deallocation of "inter_front" in the program unit  ',     &
         '"RHS_body_dynamics" failed; the execution terminates here. '
      stop
   endif
endif
if(allocated(bp_bound_interactions)) then
   deallocate(bp_bound_interactions,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Deallocation of "bp_bound_interactions" in the program ', &
         'unit "RHS_body_dynamics" failed; the execution terminates here. '
      stop
   endif
endif
if(allocated(interface_sliding_vel_max)) then
   deallocate(interface_sliding_vel_max,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Deallocation of "interface_sliding_vel_max" in the ',     &
         'program unit "RHS_body_dynamics" failed; the execution terminates ', &
         'here. '
      stop
   endif
endif
return
end subroutine RHS_body_dynamics
#endif
