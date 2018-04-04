!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2018 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: velocity_smoothing_SA_SPH_3D 
! Description: To calculate a corrective term for velocity.    
!-------------------------------------------------------------------------------
subroutine velocity_smoothing_SA_SPH_3D(npi)
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
integer(4) :: i,j,ibdt,ibdp,Ncbf,icbf,iface,facestr
double precision :: IntWdV
double precision,dimension(1:SPACEDIM) :: DVLoc,DVGlo,BCLoc,BCGlo,LocX
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
!------------------------
! Statements
!------------------------
ncbf = BoundaryDataPointer(1,npi)
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
         BCLoc(1) = DVLoc(1) * IntWdV * Tratto(facestr)%ShearCoeff
         BCLoc(2) = DVLoc(2) * IntWdV * Tratto(facestr)%ShearCoeff
         BCLoc(3) = DVLoc(3) * IntWdV
         do i=1,SPACEDIM
            BCGlo(i) = zero
            do j=1,SPACEDIM
               BCGlo(i) = BCGlo(i) + BoundaryFace(iface)%T(i,j) *     &
                          BCLoc(j)
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

