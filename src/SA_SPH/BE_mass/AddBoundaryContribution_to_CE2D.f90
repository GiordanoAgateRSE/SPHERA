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
! Program unit: AddBoundaryContribution_to_CE2D                                
! Description: Computation of the SASPH terms for the 2D continuity equation 
!              of a generic "npi" SPH particle (Di Monaco et al., 2011, EACFM).
!              SASPH contributions to the renormalization matrix for the 2D 
!              velocity-divergence term (only for the first step).
!              ALE3-LC: 2D auxiliary vectors for the explicit ALE1 SASPH term 
!              of CE; ALE implicit terms in CE.
!-------------------------------------------------------------------------------
#ifdef SPACE_2D
subroutine AddBoundaryContribution_to_CE2D(npi,IntNcbs,grad_u_SA,grad_w_SA,    &
   grad_rhod1u_SA,grad_rhod1w_SA)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use Memory_I_O_interface_module
use SA_SPH_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: npi,IntNcbs
double precision,dimension(3),intent(inout) :: grad_u_SA,grad_w_SA
double precision,dimension(3),intent(inout) :: grad_rhod1u_SA,grad_rhod1w_SA
integer(4) :: pd,icbs,iside,sidestr,ibdt,ibdp,ii
double precision :: IntWds,roi,IntWdV
integer(4),dimension(1:PLANEDIM) :: acix
double precision,dimension(1:PLANEDIM) :: IntLocXY,nnlocal,dvel
type (TyBoundarySide) :: RifBoundarySide
character(4) :: strtype
double precision,dimension(1:SPACEDIM,1:SPACEDIM) :: B_ren_aux
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
if (IntNcbs<=0) return
! Active coordinate indices
acix(1) = 1        
acix(2) = 3
roi = pg(npi)%dens
ibdt = BoundaryDataPointer(3,npi)
!------------------------
! Statements
!------------------------
do icbs=1,IntNcbs
   ibdp = ibdt + icbs - 1
   IntLocXY(1:PLANEDIM) = BoundaryDataTab(ibdp)%LocXYZ(1:PLANEDIM)
   iside = BoundaryDataTab(ibdp)%CloBoNum
   RifBoundarySide = BoundarySide(iside)
   sidestr = RifBoundarySide%stretch
   strtype = Tratto(sidestr)%tipo
! SASPH contributions to the renormalization matrix for div_u_ and grad_p: start
   if (input_any_t%C1_BE) then
! "IntWdV", or equivalently "J_3,w" (2D version) is always computed using the 
! beta-spline cubic kernel, no matter about the renormalization
      IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(3)
      B_ren_aux(1:3,1:3) = -RifBoundarySide%RN(1:3,1:3) * IntWdV
! Renormalization at SASPH frontiers
      if (input_any_t%C1_divu) then
         pg(npi)%B_ren_divu(1:3,1:3) = pg(npi)%B_ren_divu(1:3,1:3) +           &
                                       B_ren_aux(1:3,1:3)
      endif
      pg(npi)%B_ren_gradp(1:3,1:3) = pg(npi)%B_ren_gradp(1:3,1:3) +            &
                                     B_ren_aux(1:3,1:3)
   endif
! SASPH contributions to the renormalization matrix for div_u_ and grad_p: end
   if (strtype=="fixe".or.strtype=="tapi".or.strtype=="velo".or.               &
      strtype=="flow".or.strtype=="sour") then 
      IntWdS = BoundaryDataTab(ibdp)%BoundaryIntegral(1)
      do pd=1,PLANEDIM
         nnlocal(pd) = RifBoundarySide%T(acix(pd),acix(2))
      enddo
      select case (strtype)
! No-slip conditions always apply to the 2D SASPH CE term. Thus, no correction 
! is requested if (ALE3).
         case ("fixe")
            do pd=1,PLANEDIM
               dvel(pd) = 2.d0 * (-pg(npi)%var(acix(pd)))
            enddo
         case ("tapi")
            do pd=1,PLANEDIM
               dvel(pd) = 2.d0 * (RifBoundarySide%velocity(acix(pd)) -         &
                          pg(npi)%var(acix(pd)))
            enddo
         case ("velo","flow","sour")
            if ((Domain%time_stage==1).or.(Domain%time_split==1)) then 
               pg(npi)%koddens = 2
               grad_u_SA(1:3) = 0.d0
               grad_w_SA(1:3) = 0.d0
            endif
            return
      endselect
 ! Summations of the SASPH terms for grad_u_SA and grad_w_SA: start
      do ii=1,PLANEDIM
         grad_u_SA(acix(ii)) = grad_u_SA(acix(ii)) - dvel(1) *                 &
                               RifBoundarySide%T(acix(ii),acix(2)) * IntWdS
         grad_w_SA(acix(ii)) = grad_w_SA(acix(ii)) - dvel(2) *                 &
                               RifBoundarySide%T(acix(ii),acix(2)) * IntWdS
      enddo
! Summations of the SASPH terms for grad_u_SA and grad_w_SA: end
      if ((input_any_t%ALE3).and.(.not.(pg(npi)%p0_neg_ALE))) then
! Auxiliary vectors for the SASPH ALE1 explicit term in CE
         do ii=1,PLANEDIM
            grad_rhod1u_SA(ii) = grad_rhod1u_SA(ii) + 2.d0 * pg(npi)%dens *    &
                                 pg(npi)%dvel_ALE1(ii) *                       &
                                 RifBoundarySide%T(acix(ii),acix(2)) * IntWdS
            grad_rhod1w_SA(ii) = grad_rhod1w_SA(ii) + 2.d0 * pg(npi)%dens *    &
                                 pg(npi)%dvel_ALE1(ii) *                       &
                                 RifBoundarySide%T(acix(ii),acix(2)) * IntWdS
         enddo
      endif
   endif
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine AddBoundaryContribution_to_CE2D
#endif
