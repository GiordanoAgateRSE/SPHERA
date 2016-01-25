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
! Program unit: body_pressure_mirror
! Description: Computation of the body particle pressure (Amicarelli et al., 2015, CAF).   
!----------------------------------------------------------------------------------------------------------------------------------

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
integer(4) :: npi,j,npartint,npj
double precision :: Sum_W_vol,W_vol,dis,pres_mir,aux
double precision :: aux_acc(3)
double precision, external :: w
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
! Loop over body particles 
! Draft for omp parallelization with critical section 
!$omp parallel do default(none)                                                &
!$omp private(npi,Sum_W_vol,j,npartint,npj,aux_acc,pres_mir,dis,W_vol,aux)     &
!$omp shared(n_body_part,bp_arr,nPartIntorno_bp_f,NMAXPARTJ,PartIntorno_bp_f)  &
!$omp shared(Domain,pg,rag_bp_f)
do npi=1,n_body_part
   bp_arr(npi)%pres = 0.
   Sum_W_vol = 0.  
   do j=1,nPartIntorno_bp_f(npi)
      npartint = (npi - 1) * NMAXPARTJ + j
      npj = PartIntorno_bp_f(npartint)
      aux = dsqrt(dot_product(bp_arr(npi)%acc,bp_arr(npi)%acc))
! Wall acceleration should be less than 100m/2^2, otherwise an impulse is 
! assumed to occur and the formulation with acc_body is not valid
      if (aux<=100.) then
         aux_acc(:) = Domain%grav(:) - bp_arr(npi)%acc(:)
         else
            aux_acc(:) = Domain%grav(:) - (100./aux) * bp_arr(npi)%acc(:)
      endif
      pres_mir = pg(npj)%pres - pg(npj)%dens *                                 &
         dot_product(aux_acc,-rag_bp_f(:,npartint))
      dis = dsqrt(dot_product(rag_bp_f(:,npartint),rag_bp_f(:,npartint)))
      W_vol = w(dis,Domain%h,Domain%coefke) * (pg(npj)%mass / pg(npj)%dens)
      bp_arr(npi)%pres = bp_arr(npi)%pres + pres_mir * W_vol
      Sum_W_vol = Sum_W_vol + W_vol  
   end do
   if (Sum_W_vol>0.) bp_arr(npi)%pres = bp_arr(npi)%pres / Sum_W_vol      
end do
! Draft for omp parallelization with critical section 
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine body_pressure_mirror

