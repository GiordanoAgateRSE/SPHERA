!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
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
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: ComputeBoundaryDataTab                                   
! Description: To calculate the array to store close boundaries and integrals.
!              (Di Monaco et al., 2011, EACFM)                        
!-------------------------------------------------------------------------------
subroutine ComputeBoundaryDataTab 
!------------------------
! Modules
!------------------------ 
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: npi,Ncbs,IntNcbs,icbs,Ncbf,icbf,ibdt,Nfzn,Ncols
double precision :: IntWdS,IntWdV,IntWdV1,IntdWrm1dV,IntGWZrm1dV,IntWd1s0
double precision :: IntWd3s0,IntWd1s2,deltai,ypi,xpmin,xpmax,interlen  
character(len=lencard) :: nomsub = "ComputeBoundaryDataTab"
integer(4),dimension(1:NUMCOLS_BIT) :: Colmn
integer(4),dimension(1:MAXCLOSEBOUNDSIDES) :: Cloboside,Intboside
integer(4),dimension(1:Domain%MAXCLOSEBOUNDFACES) :: Cloboface
double precision,dimension(1:PLANEDIM) :: IntDpWdV
double precision,dimension(1:SPACEDIM) :: IntGWdV
double precision,dimension(1:NUMCOLS_BIT) :: Func
double precision,dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES) :: LocXY,IntLocXY
double precision,dimension(1:SPACEDIM,1:Domain%MAXCLOSEBOUNDFACES) :: LocX
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
! 2D computation: start
if (ncord==2) then
! Zeroing the counter of particles close to boundaries 
   BoundarySide(:)%CloseParticles = 0
   BoundarySide(:)%CloseParticles_maxQuota = const_m_9999
!$omp parallel do default(none)                                                &
!$omp private(npi,Ncbs,Cloboside,LocXY,IntNcbs,Intboside,IntLocXY,ibdt,icbs)   &
!$omp private(xpmin,xpmax,interlen,Ncols,Colmn,deltai,Func,ypi)                &
!$omp private(IntWdS,IntWdV,IntDpWdV,IntWdV1,IntWd1s0,IntWd3s0,IntWd1s2)       &
!$omp shared(nag,pg,Domain,BoundaryDataTab,BoundaryDataPointer,MaxNcbs,nomsub) &
!$omp shared(nout,nscr,BoundarySide,squareh)
   do npi=1,nag
      if (pg(npi)%cella==0.or.pg(npi)%vel_type/="std") cycle
! Searching for the boundary sides, which are the nearest the current          
! particle "npi"
      call FindCloseBoundarySides2D (npi,Ncbs, Cloboside, LocXY)
! Some nearest boundaries have been detected
      if (Ncbs>0) then
! To select the boundary sides that effectively contribute to the 
! governing equations
         call SelectCloseBoundarySides2D(npi,Ncbs,Cloboside,LocXY,IntNcbs,     &
            Intboside,IntLocXY)
         if (IntNcbs>0) then
            BoundaryDataPointer(1,npi) = Ncbs
            BoundaryDataPointer(2,npi) = IntNcbs
            ibdt = MAXCLOSEBOUNDSIDES * (npi-1)
            BoundaryDataPointer(3,npi) = ibdt+1
            do icbs=1,IntNcbs
               ibdt = ibdt + 1
! Check array sizes
               if (ibdt>MaxNcbs) then
                  call diagnostic(arg1=8,arg2=1,arg3=nomsub)
               endif
               BoundaryDataTab(ibdt)%CloBoNum  = Intboside(icbs)
               BoundaryDataTab(ibdt)%LocXYZ(1:PLANEDIM) =                      &
                  IntLocXY(1:PLANEDIM,icbs)
               BoundaryDataTab(ibdt)%LocXYZ(3) = zero
               Func = zero
               if (IntNcbs==2) then  
! Numerical computation of the integrals "IntWds" and "IntWdV"
! To find the intersection between the boundary and the kernel support
                  call FindBoundaryIntersection2D(icbs,Intboside,IntLocXY,     &
                     BoundarySide,xpmin,xpmax,interlen)
! Computation of the 1D integrals
                  call ComputeSurfaceIntegral_WdS2D(icbs,IntLocXY,xpmin,       &
                     interlen,IntWdS)   
! Computation of the 2D integrals
                  call BoundaryVolumeIntegrals2D(icbs,IntLocXY,xpmin,xpmax,    &
                     interlen,IntWdV,IntDpWdV)
                  call ComputeVolumeIntegral_WdV2D(icbs,IntNcbs,Intboside,     &
                     IntLocXY,BoundarySide,xpmin,xpmax,interlen,IntWdV1)
! Interpolation of the integrals "IntWd1s0", "IntWd3s0" and "IntWd1s2"  
! (from tables)
                  Ncols = 3
                  Colmn(1) = 3
                  Colmn(2) = 4
                  Colmn(3) = 5
                  ypi = IntLocXY(2,icbs)
                  deltai = ypi / Domain%h
                  call InterpolateBoundaryIntegrals2D(deltai,Ncols,Colmn,Func)
                  IntWd1s0 = Func(1) / squareh
                  IntWd3s0 = Func(2) / squareh
                  IntWd1s2 = Func(3) / squareh
                  BoundaryDataTab(ibdt)%BoundaryIntegral(1) = IntWdS
                  BoundaryDataTab(ibdt)%BoundaryIntegral(2) = IntWdV
                  BoundaryDataTab(ibdt)%BoundaryIntegral(3) = IntWdV1
                  BoundaryDataTab(ibdt)%BoundaryIntegral(4:5) = IntDpWdV(1:2)
                  BoundaryDataTab(ibdt)%BoundaryIntegral(6) = IntWd1s0
                  BoundaryDataTab(ibdt)%BoundaryIntegral(7) = IntWd3s0
                  BoundaryDataTab(ibdt)%BoundaryIntegral(8) = IntWd1s2
                  elseif (IntNcbs==1) then
! Interpolation of the integrals "IntWdS", "IntWdV", "IntWd1s0", "IntWd3s0" and 
! "IntWd1s2" (from tables)
                     Ncols=5
                     Colmn(1) = 1
                     Colmn(2) = 2
                     Colmn(3) = 3
                     Colmn(4) = 4
                     Colmn(5) = 5
                     ypi = IntLocXY(2, icbs)
                     deltai = ypi / Domain%h
                     call InterpolateBoundaryIntegrals2D(deltai,Ncols,Colmn,   &
                        Func)
                     IntWdS = Func(1) / Domain%h
                     IntWdV = Func(2)
                     IntWd1s0 = Func(3) / squareh
                     IntWd3s0 = Func(4) / squareh
                     IntWd1s2 = Func(5) / squareh
                     BoundaryDataTab(ibdt)%BoundaryIntegral(1) = IntWds
                     BoundaryDataTab(ibdt)%BoundaryIntegral(2) = IntWdV
                     BoundaryDataTab(ibdt)%BoundaryIntegral(3) = IntWdV
                     BoundaryDataTab(ibdt)%BoundaryIntegral(4:5) = zero
                     BoundaryDataTab(ibdt)%BoundaryIntegral(6) = IntWd1s0
                     BoundaryDataTab(ibdt)%BoundaryIntegral(7) = IntWd3s0
                     BoundaryDataTab(ibdt)%BoundaryIntegral(8) = IntWd1s2
               endif
            enddo
         endif
      endif
   enddo
!$omp end parallel do
! 3D computation: start
   else
! Zeroing the counter of particles close to boundaries 
      BoundaryFace(:)%CloseParticles = 0
      BoundaryFace(:)%CloseParticles_maxQuota = const_m_9999
!$omp parallel do default(none)                                                &
!$omp private(npi,Ncbf,Cloboface,LocX,Nfzn,icbf,ibdt)                          &
!$omp private(IntWdV,IntdWrm1dV,IntGWZrm1dV,IntGWdV,IntGWrRdV)                 &
!$omp shared(nag,pg,BoundaryDataTab,BoundaryDataPointer,EpCount,MaxNcbf,nout)  &
!$omp shared(nscr,nomsub,Domain)
      loop_particle: do npi=1,nag
         if (pg(npi)%cella==0.or.pg(npi)%vel_type/="std") cycle loop_particle
! Searching for the boundary faces, which are the nearest the current          
! particle "npi"
         call FindCloseBoundaryFaces3D(npi,Ncbf,Cloboface,LocX,Nfzn)
         if (Ncbf==0) then
            BoundaryDataPointer(:,npi) = 0
            else
! Some nearest boundaries have been detected
! EpCount counts particles exited from the frontier, but still within the grid
               if (Nfzn==Ncbf) then
                  EpCount(pg(npi)%imed) = EpCount(pg(npi)%imed) + 1
               endif
               BoundaryDataPointer(1,npi) = Ncbf
               BoundaryDataPointer(2,npi) = 0
               ibdt = Domain%MAXCLOSEBOUNDFACES * (npi - 1)
               BoundaryDataPointer(3,npi) = ibdt + 1
! Check array sizes
               if (ibdt>MaxNcbf) then
                  call diagnostic(arg1=8,arg2=2,arg3=nomsub)
               endif
               do icbf=1,Ncbf
                  call ComputeBoundaryVolumeIntegrals_P0(icbf,Cloboface,LocX,  &
                     IntWdV,IntdWrm1dV,IntGWZrm1dV,IntGWdV,IntGWrRdV)
                  ibdt = ibdt + 1
                  BoundaryDataTab(ibdt)%CloBoNum  = CloboFace(icbf)
                  BoundaryDataTab(ibdt)%LocXYZ(1:SPACEDIM) =                   &
                     LocX(1:SPACEDIM,icbf)
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
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ComputeBoundaryDataTab

