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
! Program unit: body_particles_to_continuity
! Description: Contributions of the body particles to the continuity equation, 
!              included explicit BODY BC ALE terms and ALE implicit terms in CE. 
!-------------------------------------------------------------------------------
#ifdef SOLID_BODIES
subroutine body_particles_to_continuity
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
double precision :: temp_dden,dis,W_vol,sum_W_vol,ALE1_CE_BODY,aux_scalar
double precision :: ALE2_CE_BODY
double precision :: dvar(3),aux_vec(3),rag_bp_f_aux(3),aux_vec_ALE1(3),tau_s(3)
double precision :: delta_dvel_ALE1(3),aux_vec_ALE2(3)
double precision,external :: w  
!------------------------
! Explicit interfaces
!------------------------
interface
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
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
! Loop over the body particles
!$omp parallel do default(none)                                                &
!$omp shared(n_body_part,bp_arr,nPartIntorno_bp_f,NMAXPARTJ,PartIntorno_bp_f)  &
!$omp shared(KerDer_bp_f_cub_spl,rag_bp_f,pg,Domain,FSI_slip_conditions)       &
!$omp shared(surf_body_part,thin_walls,proxy_normal_bp_f,input_any_t)          &
!$omp private(npi,sum_W_vol,W_vol,j,npartint,npj,temp_dden,dis,dvar,aux_vec)   &
!$omp private(rag_bp_f_aux,aux_vec_ALE1,delta_dvel_ALE1,aux_vec_ALE2)          &
!$omp private(ALE1_CE_BODY,ALE2_CE_BODY,tau_s,aux_scalar)
do npi=1,n_body_part
   bp_arr(npi)%vel_mir(:) = 0.d0
   sum_W_vol = 0.d0
! Loop over fluid particles (contributions to fluid particles) 
   do j=1,nPartIntorno_bp_f(npi)
      npartint = (npi - 1) * NMAXPARTJ + j
      npj = PartIntorno_bp_f(npartint)
      temp_dden = 0.d0
! Continuity equation
! Relative velocity for the continuity equation
      select case (FSI_slip_conditions)
         case(0)
! free-slip conditions
            aux_vec(1:3) = bp_arr(proxy_normal_bp_f(npartint))%vel(1:3) -      &
                           pg(npj)%vel(1:3)
            dvar(1:3) = bp_arr(proxy_normal_bp_f(npartint))%normal(1:3) *      &
                        2.d0 * dot_product(aux_vec,                            &
                        bp_arr(proxy_normal_bp_f(npartint))%normal)
         case(1)
! no-slip conditions
            aux_vec(1:3) = bp_arr(proxy_normal_bp_f(npartint))%vel(1:3) -      &
                           pg(npj)%vel(1:3)
            dvar(1:3) = 2.d0 * aux_vec(1:3)
         case(2,3)
! mirror velocity as solid velocity
            dvar(1:3) = bp_arr(npi)%vel(1:3) - pg(npj)%vel(:)
      endselect
      if ((input_any_t%ALE3).and.(.not.(pg(npj)%p0_neg_ALE))) then
! For the ALE2-CE term (valid for any slip condition)
         delta_dvel_ALE1(1:3) = -2.d0 * pg(npj)%dvel_ALE1(1:3)
         if (FSI_slip_conditions==0) then
! Correction for the velocity divergence (free-slip conditions)
            aux_vec(1:3) = pg(npj)%dvel_ALE1(1:3) + pg(npj)%dvel_ALE3(1:3)
            tau_s(1:3) = pg(npj)%vel_fluid(1:3) -                              &
                         bp_arr(proxy_normal_bp_f(npartint))%normal(1:3) *     &
                         dot_product(pg(npj)%vel_fluid,                        &
                         bp_arr(proxy_normal_bp_f(npartint))%normal)
            aux_scalar = dsqrt(dot_product(tau_s,tau_s))
            if (aux_scalar>1.d-9) then
               tau_s(1:3) = tau_s(1:3) / aux_scalar
               else
                  tau_s(1:3) = 0.d0
            endif
            dvar(1:3) = dvar(1:3) - 2.d0 * dot_product(aux_vec,tau_s) *        &
                        tau_s(1:3)
            elseif (FSI_slip_conditions>1) then
! Correction for the velocity divergence (mirror velocity as solid velocity); 
! no-slip conditions require no correction.
               dvar(1:3) = dvar(1:3) - (pg(npj)%dvel_ALE1(1:3) +               &
                           pg(npj)%dvel_ALE3(1:3))
         endif
      endif
      dis = dsqrt(dot_product(rag_bp_f(:,npartint),rag_bp_f(:,npartint)))
      W_vol = w(dis,Domain%h,Domain%coefke) * pg(npj)%mass / pg(npj)%dens
      bp_arr(npi)%vel_mir(:) = bp_arr(npi)%vel_mir(:) + (dvar(:) +             &
                               pg(npj)%vel(:)) * W_vol
      sum_W_vol = sum_W_vol + W_vol
! Contributions to the continuity equation
      rag_bp_f_aux(1:3) = rag_bp_f(1:3,npartint)
      if ((input_any_t%ALE3).and.(.not.(pg(npj)%p0_neg_ALE))) then
! BODY BC CE ALE2 term (it has not be renormalized)
         aux_vec_ALE2(1:3) = bp_arr(npi)%volume * KerDer_bp_f_cub_spl(npartint)&
                             * rag_bp_f_aux(1:3)
         ALE2_CE_BODY = pg(npj)%dens * dot_product(delta_dvel_ALE1,aux_vec_ALE2)
! Contribution to CE
         temp_dden = temp_dden + ALE2_CE_BODY
      endif
      if (input_any_t%C1_BE) then
! 1st-order consistency
         call MatrixProduct(pg(npj)%B_ren_divu,BB=rag_bp_f(1:3,npartint),      &
            CC=aux_vec,nr=3,nrc=3,nc=1)
         rag_bp_f_aux(1:3) = -aux_vec(1:3)
      endif
! Body BC contribution to CE (velocity-divergence term)
      temp_dden = temp_dden + pg(npj)%dens * bp_arr(npi)%volume *              &
                  KerDer_bp_f_cub_spl(npartint) *                              &
                  dot_product(dvar,rag_bp_f_aux)
      if ((input_any_t%ALE3).and.(.not.(pg(npj)%p0_neg_ALE))) then
! BODY BC CE ALE1 term
!!!test: start
!         aux_vec_ALE1(1:3) = bp_arr(npi)%volume * KerDer_bp_f_cub_spl(npartint)&
!                             * rag_bp_f_aux(1:3)
         aux_vec_ALE1(1:3) = bp_arr(npi)%volume * KerDer_bp_f_cub_spl(npartint)&
                             * rag_bp_f(1:3,npartint)
!!!test: end
         ALE1_CE_BODY = -pg(npj)%dens * dot_product(delta_dvel_ALE1,           &
                        aux_vec_ALE1)
! Contribution to CE
         temp_dden = temp_dden + ALE1_CE_BODY
         if (thin_walls) then
            ALE1_CE_BODY = ALE1_CE_BODY * (1.d0 + (1.d0 - pg(npj)%sigma_fp -   &
                           pg(npj)%sigma_bp) / pg(npj)%sigma_bp)
            ALE2_CE_BODY = ALE2_CE_BODY * (1.d0 + (1.d0 - pg(npj)%sigma_fp -   &
                           pg(npj)%sigma_bp) / pg(npj)%sigma_bp)
         endif
! ALE terms are collected (for the volume equation)
!$omp critical (omp_dden_ALE12)
         pg(npj)%dden_ALE12 = pg(npj)%dden_ALE12 + ALE1_CE_BODY + ALE2_CE_BODY
!$omp end critical (omp_dden_ALE12)
      endif
      if (thin_walls) then
! Treatment for thin walls (to all the coupling terms for the continuity 
! equation)
         temp_dden = temp_dden * (1.d0 + (1.d0 - pg(npj)%sigma_fp -            &
                     pg(npj)%sigma_bp) / pg(npj)%sigma_bp)
      endif
!$omp critical (omp_body_particles_to_continuity)
      pg(npj)%dden = pg(npj)%dden + temp_dden
!$omp end critical (omp_body_particles_to_continuity)
   enddo
   if (sum_W_vol>1.d-18) bp_arr(npi)%vel_mir(:) = bp_arr(npi)%vel_mir(:) /     &
                                                  sum_W_vol
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine body_particles_to_continuity
#endif
