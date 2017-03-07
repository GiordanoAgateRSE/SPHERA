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
! Program unit: GridCellBoundaryFacesIntersections3D                                       
! Description: To find the boundary faces intercepted by each frame cell of the 
!              grid nc[1,NumCells]. In the generic row nc of the vector
!              GCBFPointers(1:NumCells,1:2), it sets:
!                 in the first column: the number of the intercepted faces
!                 in the second column: the pointer to CFBFVector, where the 
!                                       list of intercepted faces begins 
!              Searching is based on a principle of exclusion and is carried out
!              in two phases: 
!                 First phase: for every cell, it excludes (as possibly 
!                              intercepted) the faces, whose vertices all lie in
!                              one of the semispaces (defined by the planes 
!                              containing the cell faces), which do not include
!                              the cell itself.  
!                 Second phase: for every remaining face, it verifies if all the
!                               8 cell vertices belong to one of the semispaces 
!                               defined by the plane containing the face. In the
!                               positive case, the face is excluded.  
!              (Di Monaco et al., 2011, EACFM)                        
!-------------------------------------------------------------------------------
subroutine GridCellBoundaryFacesIntersections3D(NumCellmax)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(IN) :: NumCellmax
logical :: Found
integer(4) :: nc,nf,kf,i,j,k,i0,j0,k0,flpointer,nodes,no,sd,ier,flpointer_cell 
integer(4) :: i_flpointer,irestocell,ii,jj,kk,nv
double precision :: XYZnod,CellNodeZita,deltaXYZ
integer(4),dimension(-2:2) :: Signcount
integer(4),dimension(:),allocatable :: GCBFVector_aux,GCBFVector_cell
double precision,dimension(1:SPACEDIM,1:2) :: CellXYZ
double precision,dimension(1:8,1:SPACEDIM) :: CellNodeXYZ
character(len=lencard) :: nomsub = "GridCellBoundaryFacesIntersections3D"
integer(4),external :: CellIndices
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
allocate(GCBFVector_aux(GCBFVecDim))    
if (NumCellmax<Grid%nmax) call diagnostic(arg1=7,arg3=nomsub)            
!------------------------
! Initializations
!------------------------
flpointer = 0
!------------------------
! Statements
!------------------------
! Loop over all the cells of the grid
!$omp parallel do default(none)                                                &
!$omp shared(Grid,GCBFPointers,NumFacce,BFaceList,Tratto,BoundaryFace)         &
!$omp shared(flpointer,GCBFVector_aux,nomsub,Domain,nout)                      &
!$omp private(nc,irestocell,i,i0,j,j0,k,k0,CellXYZ,nv,ii,jj,kk,CellNodeXYZ,kf) &
!$omp private(nf,nodes,Found,sd,Signcount,no,XYZnod,CellNodeZita,deltaXYZ)     &
!$omp private(i_flpointer,flpointer_cell,GCBFVector_cell)
do nc=1,Grid%nmax
   GCBFPointers(nc,1) = zero
! To detect the cell IDs 
   irestocell = CellIndices(nc,i,j,k)        
! To assess the minimum and maximum coordinates of the cell 
   i0 = i - 1
   j0 = j - 1
   k0 = k - 1
   CellXYZ(1,1) = Grid%extr(1,1) + Grid%dcd(1)*i0
   CellXYZ(1,2) = Grid%extr(1,1) + Grid%dcd(1)*i
   CellXYZ(2,1) = Grid%extr(2,1) + Grid%dcd(2)*j0
   CellXYZ(2,2) = Grid%extr(2,1) + Grid%dcd(2)*j
   CellXYZ(3,1) = Grid%extr(3,1) + Grid%dcd(3)*k0
   CellXYZ(3,2) = Grid%extr(3,1) + Grid%dcd(3)*k
! To assess the coordinates of the cell vertices
   nv = 0
   do ii=1,2
      do jj=1,2
         do kk=1,2
            nv = nv + 1
            CellNodeXYZ(nv,1) = CellXYZ(1,ii)
            CellNodeXYZ(nv,2) = CellXYZ(2,jj)
            CellNodeXYZ(nv,3) = CellXYZ(3,kk)
         enddo
      enddo
   enddo
! Allocating "GCBFVector_cell" and initializing the private variable 
! "flpointer_cell"
   allocate(GCBFVector_cell(Domain%MAXCLOSEBOUNDFACES))
   flpointer_cell = 0
! Loop over all the boundary faces 
   do kf=1,NumFacce                                    
      nf = BFaceList(kf)
! To skip in case of boundary type equal to "perimeter" or "pool"
      if ((Tratto(BoundaryFace(nf)%stretch)%tipo/="peri").and.                 &
         (Tratto(BoundaryFace(nf)%stretch)%tipo/="pool")) then
         nodes = BoundaryFace(nf)%nodes
! First phase of exclusion 
         Found = .True.
         do sd=1,SPACEDIM
            Signcount(-2:2) = 0
! Loop over the face vertices (3 or 4)
            do no=1,nodes
               XYZnod = BoundaryFace(nf)%Node(no)%GX(sd)
               if (XYZnod<CellXYZ(sd,1)) then
                  Signcount(-2) = Signcount(-2) + 1
                  elseif (XYZnod==CellXYZ(sd,1)) then
                     Signcount(-1) = Signcount(-1) + 1
                     elseif (XYZnod<CellXYZ(sd,2)) then
                        Signcount(0) = Signcount(0) + 1
                        elseif (XYZnod==CellXYZ(sd,2)) then
                           Signcount(1) = Signcount(1) + 1
                           elseif (XYZnod>CellXYZ(sd,2)) then
                              Signcount(2) = Signcount(2) + 1
               endif
            enddo
            if (((Signcount(-2)+Signcount(-1))==nodes).and.                    &
               (.not.(Signcount(-1)==nodes))) then
               Found = .False.
               exit
               elseif (((Signcount(2)+Signcount(1))==nodes).and.               &
                  (.not.(Signcount(1)==nodes))) then
                  Found = .False.
                  exit
            endif
         enddo
! Second phase of exclusion
         if (Found) then
            SignCount(-2:2) = 0
            do nv=1,8
! Normal vectors 
               CellNodeZita = zero
               do sd=1,SPACEDIM
                  deltaXYZ = CellNodeXYZ(nv,sd) -                              &
                     BoundaryFace(nf)%node(nodes)%GX(sd)
                  CellNodeZita = CellNodeZita + BoundaryFace(nf)%T(sd,3) *     &
                     deltaXYZ
               enddo
               if (CellNodeZita<zero) then
                  SignCount(-1) = SignCount(-1) + 1
                  elseif (CellNodeZita==zero) then
                     SignCount(0) = SignCount(0) + 1
                     else
                        SignCount(1) = SignCount(1) + 1
               endif
            enddo
            if (((SignCount(-1)+SignCount(0))==8).and.                         &
               (.not.(SignCount(0)==4))) then
               Found = .False.
               elseif (((SignCount(1)+SignCount(0))==8).and.                   &
                  (.not.(SignCount(0)==4))) then
                  Found = .False.
            endif
! Final check. If "found" is equal to ".true.", then the face "nf" cuts the 
! cell "nc" and is counted.
            if (Found) then        
               flpointer_cell = flpointer_cell + 1
               GCBFVector_cell(flpointer_cell) = nf
               GCBFPointers(nc,1) = GCBFPointers(nc,1) + 1
               if (flpointer_cell>Domain%MAXCLOSEBOUNDFACES) then
                  write(nout,'(1x,a)')                                         &
" Too many faces crossing a given cell. Please increase the parameter MAXCLOSEBOUNDFACES. "
                  call diagnostic(arg1=4,arg3=nomsub)
               endif
            endif
         endif
      endif
   enddo
!$omp critical (cs_GridCellBoundaryFacesIntersections3D)
   do i_flpointer=1,flpointer_cell    
      flpointer = flpointer + 1
      GCBFVector_aux(flpointer) = GCBFVector_cell(i_flpointer)
      if (i_flpointer==1) GCBFPointers(nc,2) = flpointer
   enddo
!$omp end critical (cs_GridCellBoundaryFacesIntersections3D)
! Deallocate the auxiliary array
   deallocate(GCBFVector_cell)
enddo
!$omp end parallel do
! To save the size of the GCBFVector array
GCBFVecDim = flpointer
! Allocating "GCBFVector" 
allocate(GCBFVector(GCBFVecDim),stat=ier)    
if (ier/=0) then
   write(nout,'(1x,a,i2)') "   Array GCBFVector not allocated. Error code: ",ier
   call diagnostic(arg1=4,arg3=nomsub)
   else
      write(nout,'(1x,a)') "   Array GCBFVector successfully allocated "
endif
! Loop over the reference array "GCBFVector"
!$omp parallel do default(none)                                                &
!$omp shared(GCBFVector,GCBFVector_aux,GCBFVecDim)                             &
!$omp private(i)
do i=1,GCBFVecDim
! Copy the auxiliary array in the reference array
   GCBFVector(i) = GCBFVector_aux(i)
enddo
!$omp end parallel do  
!------------------------
! Deallocations
!------------------------
deallocate(GCBFVector_aux)
return
end subroutine GridCellBoundaryFacesIntersections3D

