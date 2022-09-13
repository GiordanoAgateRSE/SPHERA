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
! Program unit: AddBoundaryContribution_to_CE3D                                 
! Description: Computation of the SASPH boundary terms for the 3D continuity 
!              equation (Di Monaco et al., 2011, EACFM).
!              SASPH contributions to the renormalization matrix for the 3D 
!              velocity-divergence term (only for the first step).
!              ALE3-LC: 3D auxiliary vectors for the explicit ALE1 SASPH term 
!              of CE.
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine AddBoundaryContribution_to_CE3D(npi,Ncbf,grad_u_SA,grad_v_SA,       &
   grad_w_SA,grad_rhod1u_SA,grad_rhod1v_SA,grad_rhod1w_SA)
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
integer(4),intent(in)    :: npi,Ncbf
double precision,dimension(3),intent(inout) :: grad_u_SA,grad_v_SA,grad_w_SA
double precision,dimension(3),intent(inout) :: grad_rhod1u_SA,grad_rhod1v_SA
double precision,dimension(3),intent(inout) :: grad_rhod1w_SA
integer(4) :: sd,sdj,icbf,iface,ibdt,ibdp,stretch,ii,test
double precision :: aux_scalar
double precision,dimension(1:SPACEDIM) :: LocPi,one_Loc,dvel
double precision,dimension(1:SPACEDIM) :: one_vec_dir,B_ren_aux_Loc,aux_vec_3
double precision,dimension(1:SPACEDIM) :: B_ren_aux_Glo,aux_vec,aux_vec_2,tau_s
character(4) :: boundtype
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
if (Ncbf<=0) return
ibdt = BoundaryDataPointer(3,npi)
one_vec_dir(1:3) = 1.d0 / dsqrt(3.d0)
!------------------------
! Statements
!------------------------
do icbf=1,Ncbf
   ibdp = ibdt + icbf - 1
   LocPi(1:SPACEDIM) = BoundaryDataTab(ibdp)%LocXYZ(1:SPACEDIM)
   iface = BoundaryDataTab(ibdp)%CloBoNum
   stretch = BoundaryFace(iface)%stretch
   boundtype = Tratto(stretch)%tipo
   if (boundtype=="fixe".or.boundtype=="tapi") then
! The face "iface" interacts with particle Pi
      if (LocPi(3)>zero) then
! Contributions of the neighbouring SASPH frontiers to the inverse of the 
! renormalization matrix for div_u_ and grad_p: start
         if (input_any_t%C1_BE) then
! Local components explicitly depending on the unit vector of the unity vector
            do sd=1,SPACEDIM
               one_Loc(sd) = 0.d0
               do sdj=1,SPACEDIM
                  one_Loc(sd) = one_Loc(sd) + BoundaryFace(iface)%T(sdj,sd) *  &
                     one_vec_dir(sdj)
               enddo
            enddo
! Local components (second assessment)
! "IntGiWrRdV", or equivalently "J_3,w" is always computed using the 
! beta-spline cubic kernel, no matter about the renormalization
            call MatrixProduct(BoundaryDataTab(ibdp)%IntGiWrRdV,BB=one_Loc,    &
               CC=B_ren_aux_Loc,nr=3,nrc=3,nc=1)
! Global components
            call MatrixProduct(BoundaryFace(iface)%T,BB=B_ren_aux_Loc,         &
               CC=B_ren_aux_Glo,nr=3,nrc=3,nc=1)
! Contribution to the renormalization matrix
            do ii=1,3
               pg(npi)%B_ren_divu(ii,1:3) = pg(npi)%B_ren_divu(ii,1:3) +       &
                                            B_ren_aux_Glo(ii)
               pg(npi)%B_ren_gradp(ii,1:3) = pg(npi)%B_ren_gradp(ii,1:3) +     &
                                             B_ren_aux_Glo(ii)
            enddo
         endif
! Contributions of the neighbouring SASPH frontiers to the inverse of the 
! renormalization matrix for div_u_ and grad_p: end
! Summations of the SASPH terms for grad_u_SA, grad_v_SA and grad_w_SA: start
         aux_vec_2(1:3) = two * (BoundaryFace(iface)%velocity(1:3) -           &
                          pg(npi)%var(1:3))
         aux_vec(1:3) = BoundaryFace(iface)%T(1:3,3)
         dvel(1:3) = dot_product(aux_vec_2,aux_vec) * aux_vec(1:3)
         if (input_any_t%ALE3) then
! Correction for the velocity divergence
! Free-slip conditions always apply to the 3D SASPH CE term
            aux_vec_3(1:3) = pg(npi)%dvel_ALE1(1:3) + pg(npi)%dvel_ALE3(1:3)
            tau_s(1:3) = pg(npi)%vel(1:3) - BoundaryFace(iface)%T(1:3,3) *     &
                         dot_product(pg(npi)%vel,BoundaryFace(iface)%T(1:3,3))
            aux_scalar = dsqrt(dot_product(tau_s,tau_s))
            if (aux_scalar>1.d-9) then
               tau_s(1:3) = tau_s(1:3) / aux_scalar
               else
                  tau_s(1:3) = aux_vec_3(1:3) -                                &
                               BoundaryFace(iface)%T(1:3,3)                    &
                               * dot_product(aux_vec_3,                        &
                               BoundaryFace(iface)%T(1:3,3))
                  aux_scalar = dsqrt(dot_product(tau_s,tau_s))
                  if (aux_scalar>1.d-9) then
                     tau_s(1:3) = tau_s(1:3) / aux_scalar
                     else
                        tau_s(1:3) = 0.d0
                  endif
            endif
            dvel(1:3) = dvel(1:3) - 2.d0 * dot_product(aux_vec_3,tau_s) *      &
                        tau_s(1:3)
         endif
         call MatrixProduct(BoundaryFace(iface)%T,                             &
            BB=BoundaryDataTab(ibdp)%BoundaryIntegral(4:6),CC=aux_vec,nr=3,    &
            nrc=3,nc=1)
         grad_u_SA(1:3) = grad_u_SA(1:3) - dvel(1) * aux_vec(1:3)
         grad_v_SA(1:3) = grad_v_SA(1:3) - dvel(2) * aux_vec(1:3)
         grad_w_SA(1:3) = grad_w_SA(1:3) - dvel(3) * aux_vec(1:3)
! Summations of the SASPH terms for grad_u_SA, grad_v_SA and grad_w_SA: end
         if (input_any_t%ALE3) then
! Auxiliary vectors for the SASPH ALE1 explicit term in CE
            grad_rhod1u_SA(1:3) = grad_rhod1u_SA(1:3) + 2.d0 * pg(npi)%dens *  &
                                  pg(npi)%dvel_ALE1(1) *                       &
                                  BoundaryDataTab(ibdp)%BoundaryIntegral(4:6)
#ifdef SPACE_3D
            grad_rhod1v_SA(1:3) = grad_rhod1v_SA(1:3) + 2.d0 * pg(npi)%dens *  &
                                  pg(npi)%dvel_ALE1(2) *                       &
                                  BoundaryDataTab(ibdp)%BoundaryIntegral(4:6)
#endif
            grad_rhod1w_SA(1:3) = grad_rhod1w_SA(1:3) + 2.d0 * pg(npi)%dens *  &
                                  pg(npi)%dvel_ALE1(3) *                       &
                                  BoundaryDataTab(ibdp)%BoundaryIntegral(4:6)
         endif
      endif
      elseif (boundtype=="velo".or.boundtype=="flow".or.boundtype=="sour") then
         if ((Domain%time_stage==1).or.(Domain%time_split==1)) then 
            pg(npi)%koddens = 2
            grad_u_SA(1:3) = 0.d0
            grad_v_SA(1:3) = 0.d0
            grad_w_SA(1:3) = 0.d0
         endif
         return
   endif
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine AddBoundaryContribution_to_CE3D
#endif
