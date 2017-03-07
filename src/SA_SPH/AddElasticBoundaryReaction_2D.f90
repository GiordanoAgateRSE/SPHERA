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
! Program unit: AddElasticBoundaryReaction_2D                                
! Description: To add supplementary normal boundary reaction to support eventual
!              insufficient pressure gradient boundary term. 
!              In case of few neighbouring particles and presence of normal 
!              component of mass force (gravity).
!              The normal reaction is computed with the formula 
!              R=(c0^2/d) ln(zi/d) [for zi<d], stemming from the compressible 
!              reaction of the fluid, where:
!                c0^2 = E/ro0 is the square of the sound speed within the fluid;
!                zi is the distance of the particle Pi from the boundary face;
!                d is a reference distance from which the reaction is added.
!              Check that the elastic boundary reaction never works.
!              To compute the boundary integral IntWdS 
!              (Di Monaco et al., 2011, EACFM).                        
!-------------------------------------------------------------------------------
subroutine AddElasticBoundaryReaction_2D(npi,Ncbs,BoundReaction)
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
integer(4),intent(IN) :: npi,Ncbs
double precision,intent(INOUT),dimension(1:SPACEDIM) :: BoundReaction
double precision,parameter :: ymincoeff = 0.25d0
double precision,parameter :: reafactor = 1.0d0
integer(4) :: sd,icbs,iside,nt,ibdt,ibdp,mate
double precision :: xpi,ypi,ypimin,celer02,vin,normreact
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
mate = pg(npi)%imed
ypimin = ymincoeff * Domain%dx
celer02 = Med(mate)%eps / Med(mate)%den0
ibdt = BoundaryDataPointer(3,npi)
!------------------------
! Statements
!------------------------
do icbs=1,Ncbs
   ibdp = ibdt + icbs - 1
   iside = BoundaryDataTab(ibdp)%CloBoNum
   nt = BoundarySide(iside)%stretch
   if (Tratto(nt)%tipo=="fixe".or.Tratto(nt)%tipo=="tapi") then
      xpi = BoundaryDataTab(ibdp)%LocXYZ(1)
      ypi = BoundaryDataTab(ibdp)%LocXYZ(2)
      if (ypi<ypimin) then
         if (xpi>zero.and.xpi<BoundarySide(iside)%Length) then  
! The projection of the particle position, which is normal to the plane of the 
! side "iside", is internal to "iside".
            vin = zero
            do sd=1,SPACEDIM
               vin = vin + pg(npi)%var(sd) * BoundarySide(iside)%T(sd,3)
            enddo
            if (vin<zero) then
               normreact = -reafactor * celer02 * DLog((Domain%h + ypi -       &
                  ypimin) / Domain%h) / Domain%h
               do sd=1,SPACEDIM
                  BoundReaction(sd) = BoundReaction(sd) + normreact *          &
                     BoundarySide(iside)%T(sd,3)
               enddo
            endif
         endif
      endif
   endif
enddo
!------------------------
! Deallocations
!------------------------
return 
end subroutine AddElasticBoundaryReaction_2D

