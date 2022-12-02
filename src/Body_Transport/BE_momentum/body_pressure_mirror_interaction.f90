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
! Program unit: body_pressure_mirror_interaction
! Description: Computation of the mirror pressure for a generic fluid-body 
!              inter-particle interaction (Amicarelli et al., 2015, CAF; 
!              Amicarelli et al., 2020, CPC; Amicarelli et al., 2022, IJCFD).   
!-------------------------------------------------------------------------------
#ifdef SOLID_BODIES
subroutine body_pressure_mirror_interaction(npi,npj,npartint,pres_mir,W_vol)
!------------------------
! Modules
!------------------------
use Hybrid_allocation_module
use Dynamic_allocation_module
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: npi,npj,npartint
double precision,intent(out) :: pres_mir,W_vol
double precision :: aux,aux_scalar,dis
double precision,dimension(3) :: aux_acc
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
aux = dsqrt(dot_product(bp_arr(proxy_normal_bp_f(npartint))%acc(:),            &
      bp_arr(proxy_normal_bp_f(npartint))%acc(:)))
! Wall acceleration should be less than 100m/2^2, otherwise an impulse is 
! assumed to occur and the formulation with acc_body is not valid
aux_scalar = 10.d0 * dsqrt(dot_product(Domain%grav(:),Domain%grav(:)))
if ((aux_scalar<1.d-9).or.(aux<=aux_scalar)) then
! Approximated second-order scheme
   aux_acc(:) = Domain%grav(:) - bp_arr(proxy_normal_bp_f(npartint))%acc(:) -  &
                2.d0 * (bp_arr(npi)%acc(:) -                                   &
                bp_arr(proxy_normal_bp_f(npartint))%acc(:))
   else
! In the presence of the limiter the scheme remains first-order
      aux_acc(:) = Domain%grav(:) - aux_scalar / aux *                         &
                   bp_arr(proxy_normal_bp_f(npartint))%acc(:)
endif
if ((FSI_slip_conditions==0).or.(FSI_slip_conditions==2)) then
   pres_mir = pg(npj)%pres + pg(npj)%dens * dot_product(aux_acc(:),            &
              bp_arr(proxy_normal_bp_f(npartint))%normal(:)) *                 &
              dot_product(rag_bp_f(:,npartint),                                &
              bp_arr(proxy_normal_bp_f(npartint))%normal(:))
   else
      pres_mir = pg(npj)%pres + pg(npj)%dens * dot_product(aux_acc(:),         &
                 rag_bp_f(:,npartint))
endif
dis = dsqrt(dot_product(rag_bp_f(:,npartint),rag_bp_f(:,npartint)))
W_vol = w(dis,Domain%h,Domain%coefke) * (pg(npj)%mass / pg(npj)%dens)
!------------------------
! Deallocations
!------------------------
return
end subroutine body_pressure_mirror_interaction
#endif
