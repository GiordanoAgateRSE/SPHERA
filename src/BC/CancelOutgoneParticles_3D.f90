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
! Program unit: CancelOutgoneParticles_3D
! Description: To count and delete the outgoing particles on boundaries of type 
!              "leve", "flow", "velo", "crit", "open". Deletion occurs in 2 
!              different ways:
!              a) If the particle belongs to a particle zone (maxzone) with the 
!                 highest index (the only zone where both particle number 
!                 reduction and increase are allowed), then the outgoing 
!                 particle (npi) is replaced by the last particle (nag) in the 
!                 particle array pg, and the total number of particle becomes 
!                 nag=nag-1; simultaneously, the index of the last particle of 
!                 the zone is changed (Partz(maxzone)%limit(2)).   
!              b) Otherwise, simply pg(npi)%cella = 0 (particle out of the 
!                 domain boundaries).
!-------------------------------------------------------------------------------
subroutine CancelOutgoneParticles_3D  
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
logical :: esci
integer(4) :: nfc,iof,npi,sdi,sdj,fkod,nodes,i 
double precision :: deltax21,deltay21,deltaz21,deltax31,deltay31,deltaz31
double precision :: deltax,deltay,deltaz,DetPnew,DetPold
double precision,dimension(1:SPACEDIM) :: LocXY,LocXYZnew,LocXYZold,csi
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
! Loop on the boundary opened faces
do iof=1,NumOpenFaces
   nfc = OpenFace(iof)
   nodes = BoundaryFace(nfc)%nodes
   deltax21 = BoundaryFace(nfc)%Node(2)%GX(1) - BoundaryFace(nfc)%Node(1)%GX(1)
   deltax31 = BoundaryFace(nfc)%Node(3)%GX(1) - BoundaryFace(nfc)%Node(1)%GX(1)
   deltay21 = BoundaryFace(nfc)%Node(2)%GX(2) - BoundaryFace(nfc)%Node(1)%GX(2)
   deltay31 = BoundaryFace(nfc)%Node(3)%GX(2) - BoundaryFace(nfc)%Node(1)%GX(2)
   deltaz21 = BoundaryFace(nfc)%Node(2)%GX(3) - BoundaryFace(nfc)%Node(1)%GX(3)
   deltaz31 = BoundaryFace(nfc)%Node(3)%GX(3) - BoundaryFace(nfc)%Node(1)%GX(3)
!$omp parallel do default(none) &
!$omp private(npi,i,deltax,deltay,deltaz,DetPnew,DetPold,LocXY,LocXYZnew)      &
!$omp private(locXYZold,sdi,sdj,esci,fkod,csi) &
!$omp shared(nag,pg,BoundaryFace,nfc,nodes,deltax21,deltax31,deltay21,deltay31)&
!$omp shared(deltaz21,deltaz31,OpCount)
! Loop on the particles in the domain
   do npi = 1,nag
      if (pg(npi)%cella==0) cycle
! Consider the three nodes that identify the plane of the current face and 
! evaluates the current and old volume of the tetrahedrons having the three 
! points as base and the current and old particle positions as vertices
      deltax = pg(npi)%coord(1) - BoundaryFace(nfc)%Node(1)%GX(1)
      deltay = pg(npi)%coord(2) - BoundaryFace(nfc)%Node(1)%GX(2)
      deltaz = pg(npi)%coord(3) - BoundaryFace(nfc)%Node(1)%GX(3)
      DetPnew = deltax21 * deltay31 * deltaz + deltax31 * deltay * deltaz21 +  &
                deltax * deltay21 * deltaz31 - deltax21 * deltay * deltaz31 -  &
                deltax31 * deltay21 * deltaz - deltax * deltay31 * deltaz21
! If the current determinant is smaller than zero (volume oriented), the particle
! might be out of the boundary face, since the reference normal is oriented 
! inside the domain
      if (DetPnew<=zero) then
! It verifies the signes of the determinants: if different, the particle passed 
! through the plane of the boundary face
         deltax = pg(npi)%CoordOld(1) - BoundaryFace(nfc)%Node(1)%GX(1)
         deltay = pg(npi)%CoordOld(2) - BoundaryFace(nfc)%Node(1)%GX(2)
         deltaz = pg(npi)%CoordOld(3) - BoundaryFace(nfc)%Node(1)%GX(3)
         DetPold = deltax21 * deltay31 * deltaz + deltax31 * deltay * deltaz21 &
                   + deltax * deltay21 * deltaz31 - deltax21 * deltay *        &
                   deltaz31 - deltax31 * deltay21 * deltaz - deltax * deltay31 &
                   * deltaz21
         if (sign(one,DetPnew)==sign(one,DetPold)) cycle
! It verifies if the particle path crosses the boundary face area. As first 
! step, it evaluates all the coordinates in the local system of the face plane
! [x'.y',z'] = Tgl*[x-x3,y-y3,z-z3]       
         do sdi = 1, spacedim
            LocXYZnew(sdi) = zero
            LocXYZold(sdi) = zero
            do sdj = 1, spacedim
               LocXYZnew(sdi) = LocXYZnew(sdi) + BoundaryFace(nfc)%T(sdj,sdi)  &
                                * (pg(npi)%coord(sdj) -                        &
                                BoundaryFace(nfc)%Node(nodes)%GX(sdj))
               LocXYZold(sdi) = LocXYZold(sdi) + BoundaryFace(nfc)%T(sdj,sdi)  &
                                * (pg(npi)%CoordOld(sdj) -                     &
                                BoundaryFace(nfc)%Node(nodes)%GX(sdj))
            enddo
         enddo
! It evaluates the coordinates of the intersection point in the local coordinate
! system; by definition, the z' value is equal to the constant of the plane 
! equation z' = c
         deltaz =  - LocXYZold(3) / (LocXYZnew(3) - LocXYZold(3))
         locXY(1) =  deltaz * (LocXYZnew(1) - LocXYZold(1))
         locXY(2) =  deltaz * (LocXYZnew(2) - LocXYZold(2))
! It transforms the local coordinates on the face in the XYZ system into the 
! local coordinates (csi1,csi2,csi3) in the local system
         Call LocalNormalCoordinates (LocXY,csi,nfc)
! It checks if the projected point falls inside the face (triangular=3 nodes, 
! rectangular=4 nodes)
         esci = .true.
         fkod = 6 - BoundaryFace(nfc)%nodes 
         do i=1,fkod
            if ((csi(i)<zero).or.(csi(i)>one)) then
               esci = .false.
               exit
            endif
         enddo
! The particle projection falls outside the face and therefore must be deleted 
         if (esci) then  
            OpCount(pg(npi)%imed) = OpCount(pg(npi)%imed) + 1    
            pg(npi)%cella = -1
         endif
      endif
   enddo
!$omp end parallel do
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine CancelOutgoneParticles_3D

