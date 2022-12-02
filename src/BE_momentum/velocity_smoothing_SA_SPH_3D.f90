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
! Program unit: velocity_smoothing_SA_SPH_3D 
! Description: SASPH 3D contribution to the partial velocity smoothing    
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine velocity_smoothing_SA_SPH_3D(npi)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use Neighbouring_Search_interface_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: npi
integer(4) :: i,j,ibdt,ibdp,Ncbf,icbf,iface,facestr,ix,iy,aux_int,aux_int_2
double precision :: IntWdV,u_t_0,aux_scalar,slip_coefficient,aux_scalar_2
double precision,dimension(1:SPACEDIM) :: DVLoc,DVGlo,BCLoc,BCGlo,LocX
double precision,dimension(1:SPACEDIM) :: u_t_0_vector
character(4) :: strtype
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine wall_function_for_SASPH(u_t_0,d_50,r_0w,slip_coefficient_0w,     &
      ni_T_0w)
      implicit none
         double precision,intent(in) :: u_t_0,d_50,r_0w
         double precision,intent(out) :: slip_coefficient_0w,ni_T_0w
   end subroutine wall_function_for_SASPH
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
Ncbf = BoundaryDataPointer(1,npi)
ibdt = BoundaryDataPointer(3,npi)
if (Ncbf>0) then
   do icbf=1,Ncbf
      ibdp = ibdt + icbf - 1
      LocX(1:SPACEDIM) = BoundaryDataTab(ibdp)%LocXYZ(1:SPACEDIM)
      iface = BoundaryDataTab(ibdp)%CloBoNum
      facestr = BoundaryFace(iface)%stretch
      strtype = Tratto(facestr)%tipo
      if ((strtype=='sour').or.(strtype=='velo').or.(strtype=='flow')) then
         pg(npi)%var(:) = zero
         exit
      endif
      DVGlo(:) = two * (Tratto(facestr)%velocity(:) - pg(npi)%vel(:))
      do i=1,SPACEDIM
         DVLoc(i) = zero
         do j=1,SPACEDIM
            DVLoc(i) = DVLoc(i) + BoundaryFace(iface)%T(j,i) * DVGlo(j)
         enddo
      enddo
      IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(2)
      if ((strtype=='fixe').or.(strtype=='tapi')) then
         if (Partz(Tratto(facestr)%zone)%slip_coefficient_mode==1) then
! Slip coefficient from input
            slip_coefficient = Partz(Tratto(facestr)%zone)%BC_shear_stress_input
            elseif (Partz(Tratto(facestr)%zone)%slip_coefficient_mode==2) then
! Slip coefficient computed
! Particle tangential (relative) velocity (vector)
! Both "DVLoc" and "T" are defined with an opposite direction
               u_t_0_vector(:) = 0.5d0 * (DVGlo(:) - DVLoc(3) *                &
                                 BoundaryFace(iface)%T(:,3))
! Particle tangential (relative) velocity (absolute value)
               u_t_0 = dsqrt(dot_product(u_t_0_vector(:),u_t_0_vector(:)))
! To assess the slip coefficient and the turbulent viscosity
               if (CLC_flag.eqv..true.) then
                  aux_int_2 = CellIndices(pg(npi)%cella,ix,iy,aux_int)
                  aux_scalar_2 = CLC%z0(ix,iy) * 10.d0
                  call wall_function_for_SASPH(u_t_0,aux_scalar_2,LocX(3),     &
                     slip_coefficient,aux_scalar)
                  else
                     aux_scalar_2 =                                            &
                        Partz(Tratto(facestr)%zone)%BC_shear_stress_input *    &
                        10.d0
                     call wall_function_for_SASPH(u_t_0,aux_scalar_2,LocX(3),  &
                        slip_coefficient,aux_scalar)
               endif
         endif
         BCLoc(1:2) = DVLoc(1:2) * IntWdV * slip_coefficient
         BCLoc(3) = DVLoc(3) * IntWdV
         do i=1,SPACEDIM
            BCGlo(i) = zero
            do j=1,SPACEDIM
               BCGlo(i) = BCGlo(i) + BoundaryFace(iface)%T(i,j) * BCLoc(j)
            enddo
         enddo
         pg(npi)%var(:) = pg(npi)%var(:) + BCGlo(:)   
      endif
   enddo
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine velocity_smoothing_SA_SPH_3D
#endif
