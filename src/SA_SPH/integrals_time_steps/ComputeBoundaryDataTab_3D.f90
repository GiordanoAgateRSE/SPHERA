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
! Program unit: ComputeBoundaryDataTab_3D                                 
! Description: To calculate the array to store close boundaries and integrals
!              in 3D (Di Monaco et al., 2011, EACFM)
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine ComputeBoundaryDataTab_3D
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use Memory_I_O_interface_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: npi,Ncbf,icbf,ibdt,Nfzn
double precision :: IntWdV,IntdWrm1dV,IntGWZrm1dV
character(len=lencard) :: nomsub = "ComputeBoundaryDataTab"
integer(4),dimension(1:input_any_t%MAXCLOSEBOUNDFACES) :: Cloboface
double precision,dimension(1:SPACEDIM) :: IntGWdV
double precision,dimension(1:SPACEDIM,1:input_any_t%MAXCLOSEBOUNDFACES) :: LocX
double precision,dimension(1:SPACEDIM,1:SPACEDIM) :: IntGWrRdV
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
BoundaryDataPointer = 0
!------------------------
! Statements
!------------------------
! Zeroing the counter of particles close to boundaries 
BoundaryFace(:)%CloseParticles = 0
BoundaryFace(:)%CloseParticles_maxQuota = const_m_9999
!$omp parallel do default(none)                                                &
!$omp private(npi,Ncbf,Cloboface,LocX,Nfzn,icbf,ibdt)                          &
!$omp private(IntWdV,IntdWrm1dV,IntGWZrm1dV,IntGWdV,IntGWrRdV)                 &
!$omp shared(nag,pg,BoundaryDataTab,BoundaryDataPointer,EpCount,MaxNcbf)       &
!$omp shared(nomsub,Domain,input_any_t)
loop_particle: do npi=1,nag
   if (pg(npi)%cella<=0.or.pg(npi)%vel_type/="std") cycle loop_particle
! Searching for the boundary faces, which are the nearest the current          
! particle "npi"
   call FindCloseBoundaryFaces3D(npi,Ncbf,Cloboface,LocX,Nfzn)
   if (Ncbf==0) then
      BoundaryDataPointer(:,npi) = 0
      else
! Some nearest boundaries have been detected
! EpCount counts particles exited from the frontier, but still within the grid
         if (Nfzn==Ncbf) then
!$omp critical (omp_ComputeBoundaryDataTab)
            EpCount(pg(npi)%imed) = EpCount(pg(npi)%imed) + 1
!$omp end critical (omp_ComputeBoundaryDataTab)
         endif
         BoundaryDataPointer(1,npi) = Ncbf
         BoundaryDataPointer(2,npi) = 0
         ibdt = input_any_t%MAXCLOSEBOUNDFACES * (npi - 1)
         BoundaryDataPointer(3,npi) = ibdt + 1
! Check array sizes
         if (ibdt>MaxNcbf) then
            call diagnostic(arg1=8,arg2=2,arg3=nomsub)
         endif
         do icbf=1,Ncbf
            call ComputeBoundaryVolumeIntegrals_P0(icbf,Cloboface,LocX,IntWdV, &
               IntdWrm1dV,IntGWZrm1dV,IntGWdV,IntGWrRdV)
            ibdt = ibdt + 1
            BoundaryDataTab(ibdt)%CloBoNum  = Cloboface(icbf)
            BoundaryDataTab(ibdt)%LocXYZ(1:SPACEDIM) = LocX(1:SPACEDIM,icbf)
            BoundaryDataTab(ibdt)%BoundaryIntegral(1) = zero
            BoundaryDataTab(ibdt)%BoundaryIntegral(2) = IntWdV
            BoundaryDataTab(ibdt)%BoundaryIntegral(3) = IntdWrm1dV
            BoundaryDataTab(ibdt)%BoundaryIntegral(4:6) = IntGWdV(1:3)
            BoundaryDataTab(ibdt)%BoundaryIntegral(7) = IntGWZrm1dV
            BoundaryDataTab(ibdt)%BoundaryIntegral(8) = zero
            BoundaryDataTab(ibdt)%IntGiWrRdV(:,:) = IntGWrRdV(:,:)
         enddo
   endif
enddo loop_particle
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine ComputeBoundaryDataTab_3D
#endif
