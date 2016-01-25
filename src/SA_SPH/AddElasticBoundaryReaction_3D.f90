!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-Ricerca di Sistema)



! SPHERA authors and email contact are provided on SPHERA documentation.

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
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: AddElasticBoundaryReaction_3D                                
! Description: To add supplementariìy normal boundary reaction to support eventual insufficient pressure gradient boundary term. 
!              in case of few neighbouring particles and presence of normal component of mass force (gravity).
!              The normal reaction is computed with the formula R=(c0^2/d) ln(zi/d) [for zi<d], stemming from the compressible 
!              reaction of the fluid, where:
!                 c0^2 = E/ro0 is the square of the sound speed within the fluid;
!                 zi is the distance of the particle Pi from the boundary face;
!                 d is a reference distance from which the reaction is added.
!              Check that the elastic boundary reaction never works.
!              (Di Monaco et al., 2011, EACFM).                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine AddElasticBoundaryReaction_3D(npi,Ncbf,BoundReaction)
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
integer(4),intent(IN) :: npi, Ncbf
double precision,intent(INOUT),dimension(1:SPACEDIM) :: BoundReaction
double precision,parameter :: zmincoeff = 0.25d0
double precision,parameter :: reafactor = 1.0d0
integer(4) :: sd,icbf,iface,nt,ibdt,ibdp,mate,fkod,ne,NCloseEdgeF
double precision :: zi,zimin,celer02,vin,normreact 
double precision :: scaprod,edgelen2,tau,edgedist2
logical,dimension(1:Domain%MAXCLOSEBOUNDFACES) :: ReaFace
double precision,dimension(1:SPACEDIM) :: PXLoc,csi,XQ,QP,QPcosdir
logical,external :: IsPointInternal
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
BoundReaction = zero
mate = pg(npi)%imed
zimin = zmincoeff * Domain%dd
celer02 = Med(mate)%eps / Med(mate)%den0
ibdt = BoundaryDataPointer(3,npi)
!------------------------
! Statements
!------------------------
do icbf=1,Ncbf
   ibdp = ibdt + icbf - 1
   iface = BoundaryDataTab(ibdp)%CloBoNum
   nt = BoundaryFace(iface)%stretch
   ReaFace(icbf) = .false.
   if (Tratto(nt)%tipo=="fixe".or.Tratto(nt)%tipo=="tapi") then
      PXLoc(:) = BoundaryDataTab(ibdp)%LocXYZ(:)
      zi = PXLoc(3)
      if (zi<zimin) then
         call LocalNormalCoordinates (PXLoc, csi, iface)
         fkod = BoundaryFace(iface)%nodes - 2
         if (IsPointInternal(fkod, csi)) then 
! The projection of the particle position, which is normal to the plane of the 
! face "iface", is internal to "iface".
            vin = zero
            do sd=1,SPACEDIM
               vin = vin + pg(npi)%var(sd) * BoundaryFace(iface)%T(sd,3)
            enddo
            if (vin<zero) then
               normreact = -reafactor * celer02 * DLog((Domain%h + zi - zimin) &
                  / Domain%h) / Domain%h
               BoundReaction(:) = BoundReaction(:) + normreact *               &
                  BoundaryFace(iface)%T(:, 3)
            endif
            ReaFace(icbf) = .true.
            else
               ReaFace(icbf) = .false.
         endif
      endif
   endif
enddo
! Reaction from possible close edges
do ne=1,NumBEdges
! To check if the particle is close to at least one of the faces to which the 
! side belongs 
   NCloseEdgeF = 0
   ibdt = BoundaryDataPointer(3,npi)
   do icbf=1,Ncbf
      ibdp = ibdt + icbf - 1
      iface = BoundaryDataTab(ibdp)%CloBoNum
      if ((iface==BoundaryConvexEdge(ne)%face(1)).and.(.Not. ReaFace(icbf)))   &
         then
         NCloseEdgeF = NCloseEdgeF + 1
         elseif ((iface==BoundaryConvexEdge(ne)%face(2)).and.                  &
            (.not.ReaFace(icbf))) then
            NCloseEdgeF = NCloseEdgeF + 1
      endif
   enddo
   if (NCloseEdgeF/=2) return
! Distance "zi" between the particle "npi" and the side "ne"
   scaprod = zero
   do sd=1,SPACEDIM
      scaprod = scaprod + (pg(npi)%Coord(sd) -                                 &
         BoundaryConvexEdge(ne)%node(1)%GX(sd)) *                              &
         BoundaryConvexEdge(ne)%component(sd)
   enddo
   edgelen2 = BoundaryConvexEdge(ne)%length * BoundaryConvexEdge(ne)%length
   tau = scaprod / edgelen2
   if ((tau>=zero).and.(tau<=one)) then
      edgedist2 = zero
      XQ(:) = BoundaryConvexEdge(ne)%node(1)%GX(:) +                           &
         BoundaryConvexEdge(ne)%component(:) * tau
      QP(:) = pg(npi)%Coord(:) - XQ(:)
      do sd=1,SPACEDIM
         edgedist2 = edgedist2 + QP(sd) * QP(sd)
      enddo
      zi = Dsqrt(edgedist2)
      if (zi<zimin) then
         vin = zero
         QPcosdir(:) = QP(:) / zi
         do sd=1,SPACEDIM
            vin = vin + pg(npi)%var(sd) * QPcosdir(sd)
         enddo
         if (vin<zero) then
            normreact = -reafactor * celer02 * DLog(zi / zimin) / zimin
            BoundReaction(:) = BoundReaction(:) +  normreact * QPcosdir(:)
         endif
      endif
   endif
enddo
!------------------------
! Deallocations
!------------------------
return 
end subroutine AddElasticBoundaryReaction_3D

