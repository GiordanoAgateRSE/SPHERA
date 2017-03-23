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
! Program unit: body_pressure_mirror
! Description: Computation of the body particle pressure (Amicarelli et al.,
!              2015, CAF).   
!-------------------------------------------------------------------------------
subroutine body_pressure_mirror
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
integer(4) :: npi,j,npartint,npj,i
double precision :: Sum_W_vol,W_vol,dis,pres_mir,aux,rho_ref,z_max,c_ref
double precision :: abs_u_max,abs_u,z_s_min_body,abs_gravity_acc,aux_scalar
double precision :: aux_acc(3)
double precision,external :: w
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
! Rough approximation of the maximum admissible pressure value on each body
if (body_maximum_pressure_limiter.eqv..true.) then
   rho_ref = maxval(Med(1:size(Med))%den0)
   z_max = maxval(pg(1:size(pg))%coord(3),mask=pg(1:nag)%cella/=0)
   c_ref = maxval(Med(1:size(Med))%celerita)
   abs_gravity_acc = dsqrt(dot_product(Domain%grav(:),Domain%grav(:)))
   abs_u_max = 0.d0
   do npi=1,nag
      if (pg(npi)%cella==0) cycle
      abs_u = pg(npi)%vel(1) * pg(npi)%vel(1) + pg(npi)%vel(2) *               &
              pg(npi)%vel(2) + pg(npi)%vel(3) * pg(npi)%vel(3)
      if (abs_u>abs_u_max) then
         abs_u_max = abs_u
      endif
   enddo
   do i=1,n_bodies
      z_s_min_body = minval(bp_arr(1:size(bp_arr))%pos(3),                     &
                     mask=bp_arr(1:size(bp_arr))%body==i)
      body_arr(i)%p_max_limiter = rho_ref * abs_gravity_acc * (z_max -         &
                         z_s_min_body + 2.d0 * Domain%h) + (1.05d0 * rho_ref)  &
                         * c_ref * (body_arr(i)%umax + abs_u_max)
      if (body_arr(i)%p_max_limiter<0.d0) body_arr(i)%p_max_limiter = 0.d0             
   enddo
endif
! Loop over body particles
! Draft for omp parallelization with critical section 
!$omp parallel do default(none)                                                &
!$omp private(npi,Sum_W_vol,j,npartint,npj,aux_acc,pres_mir,dis,W_vol,aux)     &
!$omp private(aux_scalar)                                                      &
!$omp shared(n_body_part,bp_arr,nPartIntorno_bp_f,NMAXPARTJ,PartIntorno_bp_f)  &
!$omp shared(Domain,pg,rag_bp_f,body_minimum_pressure_limiter)                 &
!$omp shared(body_maximum_pressure_limiter,body_arr,FSI_free_slip_conditions)
do npi=1,n_body_part
   bp_arr(npi)%pres = 0.d0
   Sum_W_vol = 0.d0
   do j=1,nPartIntorno_bp_f(npi)
      npartint = (npi - 1) * NMAXPARTJ + j
      npj = PartIntorno_bp_f(npartint)
      aux = dsqrt(dot_product(bp_arr(npi)%acc(:),bp_arr(npi)%acc(:)))
! Wall acceleration should be less than 100m/2^2, otherwise an impulse is 
! assumed to occur and the formulation with acc_body is not valid
      aux_scalar = 10.d0 * dsqrt(dot_product(Domain%grav(:),Domain%grav(:)))
      if (aux<=aux_scalar) then
         aux_acc(:) = Domain%grav(:) - bp_arr(npi)%acc(:)
         else
            aux_acc(:) = Domain%grav(:) - aux_scalar / aux * bp_arr(npi)%acc(:)
      endif
      if (FSI_free_slip_conditions.eqv..true.) then
         pres_mir = pg(npj)%pres + pg(npj)%dens *                              &
                    dot_product(aux_acc(:),bp_arr(npi)%normal(:)) *            &
                    dot_product(rag_bp_f(:,npartint),bp_arr(npi)%normal(:))
         else
            pres_mir = pg(npj)%pres + pg(npj)%dens *                           &
                       dot_product(aux_acc(:),rag_bp_f(:,npartint))
      endif
      dis = dsqrt(dot_product(rag_bp_f(:,npartint),rag_bp_f(:,npartint)))
      W_vol = w(dis,Domain%h,Domain%coefke) * (pg(npj)%mass / pg(npj)%dens)
      bp_arr(npi)%pres = bp_arr(npi)%pres + pres_mir * W_vol
      Sum_W_vol = Sum_W_vol + W_vol
   enddo
   if (Sum_W_vol>1.d-3) bp_arr(npi)%pres = bp_arr(npi)%pres / Sum_W_vol
   if (body_minimum_pressure_limiter.eqv..true.) then
      if (bp_arr(npi)%pres<0.d0) bp_arr(npi)%pres = 0.d0 
   endif
   if (body_maximum_pressure_limiter.eqv..true.) then
      if (bp_arr(npi)%pres>body_arr(bp_arr(npi)%body)%p_max_limiter) then
         bp_arr(npi)%pres = body_arr(bp_arr(npi)%body)%p_max_limiter
      endif
   endif
enddo
! Draft for omp parallelization with critical section 
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine body_pressure_mirror

