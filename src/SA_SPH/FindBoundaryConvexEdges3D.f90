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
! Program unit: FindBoundaryConvexEdges3D                                       
! Description: To look for possible edges with an associated convex geometry. 
!              Their geometrical data are saved in BoundaryConvexEdge as 
!              TyBoundaryConvexEdge.
!              (Di Monaco et al., 2011, EACFM)                        
!-------------------------------------------------------------------------------
subroutine FindBoundaryConvexEdges3D
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
logical :: EdgeFound
integer(4) :: kf,nf,sd,nodes,nodrif,nt,fk,i,j,nodi,nodj
integer(4) :: kf1,nf1,fk1,nt1,ii,jj,kk,nodii,nodjj,nodkk,NBFfin,NBFini,nodes1
double precision :: delta,length2,zitakk
integer(4),dimension(2,4) :: iseg
character(len=lencard) :: nomsub = "FindBoundaryConvexEdges3D"
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
iseg(1,1) = 2
iseg(1,2) = 3
iseg(1,3) = 1
iseg(1,4) = 0
iseg(2,1) = 2
iseg(2,2) = 3
iseg(2,3) = 4
iseg(2,4) = 1
NumBEdges = 0
NBFfin = NumFacce - 1
!------------------------
! Statements
!------------------------
!$omp parallel do default(none)                                                &
!$omp shared(NBFfin,BFaceList,BoundaryFace,NumFacce,Tratto,Vertice,NumBEdges)  &
!$omp shared(Domain,nomsub,BoundaryConvexEdge,iseg)                            &
!$omp private(kf,nf,nt,nodes,fk,NBFini,kf1,nf1,nt1,EdgeFound,nodes1,fk1,i,j)   &
!$omp private(nodi,nodj,ii,jj,nodii,nodjj,kk,nodkk,zitakk,length2,nodrif,sd)   &
!$omp private(delta)
do kf=1,NBFfin
   nf = BFaceList(kf)
   nt = BoundaryFace(nf)%stretch
   if ((Tratto(nt)%tipo=="fixe").or.(Tratto(nt)%tipo=="tapi")) then
      nodes = BoundaryFace(nf)%nodes
      fk = nodes - 2
      NBFini = kf + 1
      do kf1 = NBFini,NumFacce
         nf1 = BFaceList(kf1)
         nt1 = BoundaryFace(nf1)%stretch
         EdgeFound = .False.
         if ((Tratto(nt1)%tipo=="fixe").or.(Tratto(nt1)%tipo=="tapi")) then
            nodes1 = BoundaryFace(nf1)%nodes
            fk1 = nodes1 - 2
            do i=1,nodes
               j = iseg(fk, i)
               nodi = BoundaryFace(nf)%node(i)%name 
               nodj = BoundaryFace(nf)%node(j)%name 
               do ii=1,nodes1
                  jj = iseg(fk1, ii)
                  nodii = BoundaryFace(nf1)%node(ii)%name 
                  nodjj = BoundaryFace(nf1)%node(jj)%name 
                  if (((nodi==nodii).and.(nodj==nodjj)).or.((nodi==nodjj).and. &
                     (nodj==nodii))) then
! Found a common side for both faces "nf" and "nf1"
! To check if the side is "convex" 
! Node of face "nf1" different from "nodii" and "nodjj"
                     kk = iseg(fk1,jj)             
                     nodkk = BoundaryFace(nf1)%node(kk)%name
! Normal with respect to face "nf" 
                     zitakk = zero              
                     nodrif = BoundaryFace(nf)%nodes
                     do sd=1,SPACEDIM
                        zitakk = zitakk + BoundaryFace(nf)%T(sd,3) *           &
                           (Vertice(sd,nodkk) -                                &
                           BoundaryFace(nf)%node(nodrif)%GX(sd))   
                     enddo
                     if (zitakk<zero) then          
! The side of vertices (nodi,nodj) is "convex"
!$omp critical (omp_FBCE3D)
                        NumBEdges = NumBEdges + 1
                        if (NumBEdges>Domain%MAXNUMCONVEXEDGES) call           &
                           diagnostic(arg1=8,arg2=10,arg3=nomsub)
                        BoundaryConvexEdge(NumBEdges)%face(1) = nf
                        BoundaryConvexEdge(NumBEdges)%face(2) = nf1
                        BoundaryConvexEdge(NumBEdges)%node(1)%name = nodi
                        BoundaryConvexEdge(NumBEdges)%node(2)%name = nodj
                        length2 = zero
                        do sd=1,SPACEDIM
                           BoundaryConvexEdge(NumBEdges)%node(1)%GX(sd) =      &
                              Vertice(sd,nodi)
                           BoundaryConvexEdge(NumBEdges)%node(2)%GX(sd) =      &
                              Vertice(sd,nodj)
                           delta = Vertice(sd,nodj) - Vertice(sd,nodi)
                           BoundaryConvexEdge(NumBEdges)%component(sd) = delta
                           length2 = length2 + delta * delta
                        enddo
                        BoundaryConvexEdge(NumBEdges)%length = Dsqrt(length2)
!$omp end critical (omp_FBCE3D)
                        EdgeFound = .True.
                        exit 
                     endif
                  endif
               enddo
               if (EdgeFound) exit
            enddo
            if (EdgeFound) exit
         endif
      enddo
   endif
enddo
!$omp end parallel do  
!------------------------
! Deallocations
!------------------------
return
end subroutine FindBoundaryConvexEdges3D

