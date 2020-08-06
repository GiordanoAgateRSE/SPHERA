!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2020 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: velocity_smoothing_SA_SPH_2D 
! Description: To calculate a corrective term for velocity.    
!-------------------------------------------------------------------------------
#ifdef SPACE_2D
subroutine velocity_smoothing_SA_SPH_2D(npi)
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
integer(4),intent(in) :: npi
integer(4) :: i,icbs,Ncbs,IntNcbs,ibdt,ibdp,iside,sidestr
double precision :: IntWdV,u_t_0,slip_coefficient,aux_scalar
integer(4),dimension(1:PLANEDIM) :: acix
double precision,dimension(1:PLANEDIM) :: sss,nnn,DVLoc,DVGlo,BCLoc,BCGlo
double precision,dimension(1:SPACEDIM) :: u_t_0_vector
character(4) :: strtype
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
acix(1) = 1
acix(2) = 3
!------------------------
! Statements
!------------------------
Ncbs = BoundaryDataPointer(1,npi)
IntNcbs = BoundaryDataPointer(2,npi)
ibdt = BoundaryDataPointer(3,npi)
if (IntNcbs>0) then
   do icbs=1,IntNcbs
      ibdp = ibdt + icbs - 1
      iside = BoundaryDataTab(ibdp)%CloBoNum
      sidestr = BoundarySide(iside)%stretch
      strtype = Tratto(sidestr)%tipo
      if ((strtype=='sour').or.(strtype=='velo').or.(strtype=='flow'))&
         then
         pg(npi)%var(:) = zero   
         exit  
      endif
      do i=1,PLANEDIM
         sss(i) = BoundarySide(iside)%T(acix(i),1)
         nnn(i) = BoundarySide(iside)%T(acix(i),3)
         DVGlo(i) = two * (Tratto(sidestr)%velocity(acix(i)) -                 &
                    pg(npi)%vel(acix(i)))
      enddo
      DVLoc(1) = sss(1) * DVGlo(1) + sss(2) * DVGlo(2)
      DVLoc(2) = nnn(1) * DVGlo(1) + nnn(2) * DVGlo(2)
      IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(3)
      if ((strtype=='fixe').or.(strtype=='tapi')) then
         if (Partz(Tratto(sidestr)%zone)%slip_coefficient_mode==1) then
! Slip coefficient from input
            slip_coefficient = Partz(Tratto(sidestr)%zone)%BC_shear_stress_input
            elseif (Partz(Tratto(sidestr)%zone)%slip_coefficient_mode==2) then
! Slip coefficient computed
! Particle tangential (relative) velocity (vector)
! Both "DVLoc" and "T" are defined with an opposite direction
               u_t_0_vector(:) = (Tratto(sidestr)%velocity(:) - pg(npi)%vel(:))&
                                 - (0.5d0 * DVLoc(2) *                         &
                                 BoundarySide(iside)%T(:,3))
! Particle tangential (relative) velocity (absolute value)
               u_t_0 = dsqrt(dot_product(u_t_0_vector(:),u_t_0_vector(:)))
! To assess the slip coefficient
               call wall_function_for_SASPH(u_t_0,                             &
                  Partz(Tratto(sidestr)%zone)%BC_shear_stress_input,           &
                  BoundaryDataTab(ibdp)%LocXYZ(2),slip_coefficient,aux_scalar)
         endif
         BCLoc(1) = DVLoc(1) * IntWdV * slip_coefficient
         BCLoc(2) = DVLoc(2) * IntWdV
         BCGlo(1) = sss(1) * BCLoc(1) + nnn(1) * BCLoc(2)
         BCGlo(2) = sss(2) * BCLoc(1) + nnn(2) * BCLoc(2)
         pg(npi)%var(1) = pg(npi)%var(1) + BCGlo(1)   
         pg(npi)%var(3) = pg(npi)%var(3) + BCGlo(2)   
      endif
   enddo
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine velocity_smoothing_SA_SPH_2D
#endif
