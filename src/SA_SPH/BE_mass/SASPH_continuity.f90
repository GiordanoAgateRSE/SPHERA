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
! Program unit: SASPH_continuity
! Description: SASPH boundary term for the continuity equation.
!              Inversion of the renormalization matrices for the 3D and 2D 
!              velocity-divergence term and pressure-gradient term (also in the 
!              absence of SASPH neighbouring frontiers).
!              Renormalization of the SASPH velocity-divergence term of the 
!              continuity equation.
!              ALE3-LC: explicit ALE1 and ALE2 SASPH terms of CE (both C0 and 
!              C1 consistency).
!-------------------------------------------------------------------------------
subroutine SASPH_continuity(npi                                                &
#ifdef SPACE_3D
                               ,Ncbf_Max                                       &
#endif
                            )
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Dynamic_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: npi
#ifdef SPACE_3D
integer(4),intent(inout) :: Ncbf_Max
#endif
#ifdef SPACE_3D
integer(4) :: Ncbf
#elif defined SPACE_2D
integer(4) :: Ncbs,IntNcbs
#endif
! Explicit SASPH ALE1 and ALE2 terms in CE
double precision :: ALE1_CE_SA,ALE2_CE_SA
! SASPH terms for the gradients of the velocity components
double precision,dimension(3) :: grad_u_SA,grad_v_SA,grad_w_SA
double precision,dimension(3) :: grad_rhod1u_SA,grad_rhod1v_SA,grad_rhod1w_SA
double precision,dimension(3) :: aux_vec
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
grad_u_SA(1:3) = 0.d0
grad_v_SA(1:3) = 0.d0
grad_w_SA(1:3) = 0.d0
grad_rhod1u_SA(1:3) = 0.d0
grad_rhod1v_SA(1:3) = 0.d0
grad_rhod1w_SA(1:3) = 0.d0
!------------------------
! Statements
!------------------------
! Seaching for the neighbouring faces/sides of the particle "npi"
#ifdef SPACE_3D
   Ncbf = BoundaryDataPointer(1,npi)
! Detecting the faces with actual contributions
   if (Ncbf>0) then
!$omp critical (omp_Ncbf_Max_2)
      Ncbf_Max = max(Ncbf_Max,Ncbf)
!$omp end critical (omp_Ncbf_Max_2)
      call AddBoundaryContribution_to_CE3D(npi,Ncbf,grad_u_SA,grad_v_SA,       &
         grad_w_SA,grad_rhod1u_SA,grad_rhod1v_SA,grad_rhod1w_SA)
   endif
#elif defined SPACE_2D
      Ncbs = BoundaryDataPointer(1,npi)
      IntNcbs = BoundaryDataPointer(2,npi)
! Detecting the sides with actual contributions
      if ((Ncbs>0).and.(IntNcbs>0)) then
         call AddBoundaryContribution_to_CE2D(npi,IntNcbs,grad_u_SA,grad_w_SA, &
            grad_rhod1u_SA,grad_rhod1w_SA)
      endif
#endif
if ((input_any_t%ALE3).and.(.not.(pg(npi)%p0_neg_ALE))) then
   ALE2_CE_SA = -(grad_rhod1u_SA(1) + grad_rhod1v_SA(2) + grad_rhod1w_SA(3))
endif
if (input_any_t%C1_BE) then
! The renormalization matrices are inverted just after all its components are 
! collected and before the 1st-order consistency scheme applies to the 
! summation of all the particle-boundary contributions (even in the absence of 
! SASPH neighbouring frontiers).
   call B_ren_divu_inversion(npi)
   call B_ren_gradp_inversion(npi)
endif
! SASPH term of the momentum divergence: start
#ifdef SPACE_3D
if (Ncbf>0) then
#elif defined SPACE_2D
if ((Ncbs>0).and.(IntNcbs>0)) then
#endif
   if (input_any_t%C1_BE) then
! 1st-order consistency for the SASPH terms
! For the velocity-divergence SASPH term: start
      call MatrixProduct(pg(npi)%B_ren_divu,BB=grad_u_SA,CC=aux_vec,nr=3,      &
         nrc=3,nc=1)
      grad_u_SA(1:3) = -aux_vec(1:3)
#ifdef SPACE_3D
      call MatrixProduct(pg(npi)%B_ren_divu,BB=grad_v_SA,CC=aux_vec,nr=3,      &
         nrc=3,nc=1)
      grad_v_SA(1:3) = -aux_vec(1:3)
#endif
      call MatrixProduct(pg(npi)%B_ren_divu,BB=grad_w_SA,CC=aux_vec,nr=3,      &
         nrc=3,nc=1)
      grad_w_SA(1:3) = -aux_vec(1:3)
   endif
! For the velocity-divergence SASPH term: end
! For the explicit ALE1 SASPH term: start
!!!test: start
!   if ((input_any_t%C1_BE).and.(.not.(pg(npi)%p0_neg_ALE))) then
!      call MatrixProduct(pg(npi)%B_ren_divu,BB=grad_rhod1u_SA,CC=aux_vec,      &
!         nr=3,nrc=3,nc=1)
!      grad_rhod1u_SA(1:3) = -aux_vec(1:3)
!#ifdef SPACE_3D
!      call MatrixProduct(pg(npi)%B_ren_divu,BB=grad_rhod1v_SA,CC=aux_vec,      &
!         nr=3,nrc=3,nc=1)
!      grad_rhod1v_SA(1:3) = -aux_vec(1:3)
!#endif
!      call MatrixProduct(pg(npi)%B_ren_divu,BB=grad_rhod1w_SA,CC=aux_vec,      &
!         nr=3,nrc=3,nc=1)
!      grad_rhod1w_SA(1:3) = -aux_vec(1:3)
!   endif
!!!test: end
! For the explicit ALE1 SASPH term: end
! Adding the SASPH boundary term of the momentum divergence to the continuity 
! equation
   if (pg(npi)%koddens==0) then
      pg(npi)%dden = pg(npi)%dden - pg(npi)%dens * (grad_u_SA(1) + grad_v_SA(2)&
                     + grad_w_SA(3))
      if ((input_any_t%ALE3).and.(.not.(pg(npi)%p0_neg_ALE))) then
         ALE1_CE_SA = grad_rhod1u_SA(1) + grad_rhod1v_SA(2) + grad_rhod1w_SA(3)
         pg(npi)%dden = pg(npi)%dden + ALE1_CE_SA + ALE2_CE_SA
         pg(npi)%dden_ALE12 = pg(npi)%dden_ALE12 + ALE1_CE_SA + ALE2_CE_SA
      endif
   endif
endif
! SASPH term of the momentum divergence: end
!------------------------
! Deallocations
!------------------------
return
end subroutine SASPH_continuity
