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
! Program unit: FindCloseBoundaryFaces3D                                       
! Description: To finds the "close" boundary faces, i.e. those faces located at a distance from the particle npi smaller 
!              than or equal to 2h. It returns:
!                 Ncbf: number of close boundary faces
!                 Clobface(1 to Ncbf): list of close boundary faces
!                 LocX(1:SPACEDIM,Ncbf): local coordinates of particle npi with respect each boundary side
!              The algorithm looks for the boundary faces intersected by the cell boxes of the reference frame located 
!              all around particle npi, and cancels the repeated ones.
!              (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine FindCloseBoundaryFaces3D(npi,Ncbf,Clobface,LocX,Nfzn)
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
integer(4),intent(IN) :: npi
integer(4),intent(INOUT) :: Ncbf, Nfzn
integer(4),intent(INOUT),dimension(1:Domain%MAXCLOSEBOUNDFACES) :: Clobface
double precision,intent(INOUT),dimension(1:SPACEDIM,1:Domain%MAXCLOSEBOUNDFACES) :: LocX
logical :: Thereis
integer(4) :: nc,ic,jc,kc,i,j,k,sdi,sdj,nodes,irestocell,fkod  
integer(4) :: flpini,flp,flpfin,nfpercell,intbf,icbf,nbface,stretch
double precision :: pin,pinmin,pinmax
double precision,dimension(1:SPACEDIM) :: PXLoc,csi
character(len=lencard) :: nomsub = "FindCloseBoundaryFaces3D"
logical,external :: IsPointInternal
integer(4), external :: CellNumber,ParticleCellNumber,CellIndices
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
Clobface = 0
Ncbf = 0
Nfzn = 0
LocX = zero
pg(npi)%CloseBcOut = 0
!------------------------
! Statements
!------------------------
! To find the cell ID of the current particle
nc = ParticleCellNumber(pg(npi)%coord)
if (nc<=0) return
! To find the cell indices in the background positioning grid
irestocell = CellIndices(nc,ic,jc,kc)
! Loop over adjacent cells
do i=(ic-1),(ic+1)
   do j=(jc-1),(jc+1)
      do k=(kc-1),(kc+1)
! To find the cell ID of the adjacent cell
         nc = CellNumber(i, j, k)
         if (nc==0) cycle
! To load the number of boundary faces cutting the cell and the initial and 
! final pointers
         nfpercell = GCBFPointers(nc,1)
         if (nfpercell>0) then
            flpini = GCBFPointers(nc, 2)
            flpfin = flpini + nfpercell - 1
! Loop over the cutting faces
            do flp=flpini,flpfin
! To load the face index
               intbf = GCBFVector(flp)
               thereis = .false.
               pinmin = zero
               pinmax = doubleh
               stretch = BoundaryFace(intbf)%stretch
               if (Tratto(stretch)%tipo=="sour") pinmin = -doubleh   
! To check if the face "intbf" is already included in the array "Clobface"
               if (Ncbf>0) thereis = any(Clobface(1:Ncbf)==intbf)
! The face is not yet considered
               if (.Not. Thereis) then         
! Normal to the face "intbf"
                  nodes = BoundaryFace(intbf)%nodes
                  do sdi=1,SPACEDIM
                     PXLoc(sdi) = zero
                     do sdj=1,SPACEDIM
                        PXLoc(sdi) = PXLoc(sdi) +                              &
                           BoundaryFace(intbf)%T(sdj,sdi) *                    &
                           (pg(npi)%coord(sdj) -                               &
                           BoundaryFace(intbf)%Node(nodes)%GX(sdj))
                     enddo
                  enddo
                  pin = PXLoc(3)
                  call LocalNormalCoordinates(PXLoc,csi,intbf)
                  fkod = BoundaryFace(intbf)%nodes - 2
! Distance between the particle and the face is bigger than zero and
! smaller then 2h 
! The face is considered and it is added to the array "Clobface"
                  if (pin>=pinmin.and.pin<pinmax) then 
                     if ((Tratto(stretch)%tipo=="sour").or.                    &
                        (Tratto(stretch)%tipo=="velo").or.                     &
                        (Tratto(stretch)%tipo=="flow")) then
                        if (IsPointInternal(fkod,csi)) then  
! The projection of particle "npi" on the plane containing the face "iface"  
! is internal to the face itself 
                           Ncbf = ncbf + 1 
                           if (ncbf<=Domain%MAXCLOSEBOUNDFACES) then
                              Clobface(Ncbf) = intbf
                              LocX(3, Ncbf) = pin
                              pg(npi)%CloseBcOut = 1
                              else
                                 call diagnostic (arg1=8,arg2=6,arg3=nomsub)
                           endif
                        endif
                        else
                           Ncbf = ncbf + 1 
                           if (ncbf<=Domain%MAXCLOSEBOUNDFACES) then
                              Clobface(Ncbf) = intbf
                              LocX(3, Ncbf) = pin
                              else
                                 call diagnostic (arg1=8,arg2=7,arg3=nomsub)
                           endif
                     endif
                     elseif (pin<pinmin) then   
! One may add a test in case the projection belongs to the face 
                        if (IsPointInternal(fkod,csi)) then    
! The projection of particle "npi" on the plane containing the face "iface"  
! is internal to the face itself 
                           Nfzn = Nfzn + 1
                        endif
                  endif
               endif
            enddo
         endif
      enddo
   enddo
enddo
! To compute the particle coordinates (r,s) in the local reference system 
! (along the two axis, which belong to the face plane)
if (Ncbf>0) then
! Loop over the found faces
   do icbf=1,Ncbf
      nbface = Clobface(icbf)
      nodes = BoundaryFace(nbface)%nodes
! To increase the number of particles close to boundaries 
! and to estimate the maximum height 
!$omp critical (numpa)
      BoundaryFace(nbface)%CloseParticles =                                    &
         BoundaryFace(nbface)%CloseParticles + 1
      if (BoundaryFace(nbface)%CloseParticles_maxQuota<pg(npi)%coord(3))       &
         BoundaryFace(nbface)%CloseParticles_maxQuota = pg(npi)%coord(3)
!$omp end critical (numpa)
      do sdi=1,PLANEDIM
         LocX(sdi,icbf) = zero
         do sdj=1,SPACEDIM
            LocX(sdi,icbf) = LocX(sdi,icbf) + BoundaryFace(nbface)%T(sdj,sdi)  &
               * (pg(npi)%coord(sdj) - BoundaryFace(nbface)%Node(nodes)%GX(sdj))
         enddo
      enddo
   enddo
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine FindCloseBoundaryFaces3D

