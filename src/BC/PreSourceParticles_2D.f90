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
! Program unit: PreSourceParticles_2D
! Description: To generate new source particles at the inlet section (only in 2D
!              and with one inlet section). 
!-------------------------------------------------------------------------------
subroutine PreSourceParticles_2D
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use I_O_file_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nt,nA,isi,sd,ip,i_source
double precision :: deltapart,linedist,sidelen,eps
double precision,dimension(1:SPACEDIM) :: A,ss
integer(4), external :: ParticleCellNumber
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
! Searching the ID of the source side
SourceSide = 0
SpCount = 0
i_source=0
!------------------------
! Statements
!------------------------
do isi=1,NumBSides
   if (BoundarySide(isi)%tipo=="sour") then
      SourceSide = isi
      i_source=i_source+1
      nt = BoundarySide(SourceSide)%stretch
      irz = Tratto(nt)%zone
      mat = partz(irz)%Medium 
      nA = BoundarySide(SourceSide)%Vertex(1)
      do sd=1,SPACEDIM
         A(sd) = Vertice(sd,nA)
         ss(sd) = BoundarySide(SourceSide)%T(sd,1)
         nn(sd) = BoundarySide(SourceSide)%T(sd,3)
      end do
      deltapart = Domain%dx
      sidelen = BoundarySide(SourceSide)%length
      NumPartperLine(i_source) = Int(sidelen / deltapart + 0.01d0)
      eps = -half
      yfila = eps * deltapart 
      linedist = -half * deltapart
      do ip=1,NumPartperLine(i_source)
         linedist = linedist + deltapart
         do sd=1,SPACEDIM
            PartLine(i_source, ip, sd) = A(sd) + linedist * ss(sd)
         end do
      end do
      ParticleVolume = Domain%PVolume
      RowPeriod = ParticleVolume * NumPartperLine(i_source) / Tratto(nt)%FlowRate 
      RowVelocity(i_source) = Domain%dx / RowPeriod
      Tratto(nt)%NormVelocity = RowVelocity(i_source)
      partz(irz)%vel(1) = RowVelocity(i_source) * nn(1)
      partz(irz)%vel(2) = RowVelocity(i_source) * nn(2)
      partz(irz)%vel(3) = RowVelocity(i_source) * nn(3)
      pinttimeratio = -1
   end if 
end do
!------------------------
! Deallocations
!------------------------
return
end subroutine PreSourceParticles_2D

