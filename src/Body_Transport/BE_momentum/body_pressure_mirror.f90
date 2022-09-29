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
! Program unit: body_pressure_mirror
! Description: Computation of the body particle pressure (Amicarelli et al.,
!              2015, CAF).   
!-------------------------------------------------------------------------------
#ifdef SOLID_BODIES
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
double precision :: Sum_W_vol,W_vol,pres_mir,rho_ref,z_max,c_ref
double precision :: abs_u_max,abs_u,z_s_min_body,abs_gravity_acc
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
   z_max = maxval(pg(1:size(pg))%coord(3),mask=pg(1:size(pg))%cella/=0)
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
!$omp parallel do default(none)                                                &
!$omp private(npi,Sum_W_vol,j,npartint,npj,pres_mir,W_vol)                     &
!$omp shared(n_body_part,bp_arr,nPartIntorno_bp_f,NMAXPARTJ,PartIntorno_bp_f)  &
!$omp shared(pg,body_minimum_pressure_limiter,body_maximum_pressure_limiter)   &
!$omp shared(body_arr)
do npi=1,n_body_part
   bp_arr(npi)%pres = 0.d0
   Sum_W_vol = 0.d0
   do j=1,nPartIntorno_bp_f(npi)
      npartint = (npi - 1) * NMAXPARTJ + j
      npj = PartIntorno_bp_f(npartint)
      call body_pressure_mirror_interaction(npi,npj,npartint,pres_mir,W_vol)
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
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine body_pressure_mirror
#endif
