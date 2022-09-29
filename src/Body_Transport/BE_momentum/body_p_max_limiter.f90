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
! Program unit: body_p_max_limiter
! Description: Limiter for the body-particle maximum pressure (Amicarelli et 
!              al., 2015, CAF).   
!-------------------------------------------------------------------------------
#ifdef SOLID_BODIES
subroutine body_p_max_limiter
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
integer(4) :: npi,i
double precision :: rho_ref,z_max,c_ref
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
if (body_maximum_pressure_limiter) then
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
!------------------------
! Deallocations
!------------------------
return
end subroutine body_p_max_limiter
#endif
