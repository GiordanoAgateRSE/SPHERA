!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2021 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: DefineLocalSystemVersors                                      
! Description:  To define the directional cosines of the local reference system 
!               (Di Monaco et al., 2011, EACFM). Further modifications take into
!               account pentagon and hexagon faces (only for complex "perimeter"
!               zones / fluid reseroirs, not for SASPH frontiers).                       
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine DefineLocalSystemVersors(Nf)
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
integer(4) :: Nf
integer(4) :: i,j,n,nnodes,sidek,refnode,nod
double precision :: U1len,W12len,LocX
double precision,dimension(SPACEDIM) :: RR,ss,nnlocal,U1,U2,W12
integer(4),dimension(2,3) :: nindex 
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
nindex(1,1) = 1
nindex(1,2) = 2
nindex(1,3) = 3
nindex(2,1) = 1
nindex(2,2) = 3
nindex(2,3) = 4
nnodes = 6
if (BoundaryFace(nf)%Node(6)%name<=0) nnodes = 5
if (BoundaryFace(nf)%Node(5)%name<=0) nnodes = 4
if (BoundaryFace(Nf)%Node(4)%name<=0) nnodes = 3
BoundaryFace(Nf)%nodes = nnodes
!------------------------
! Statements
!------------------------
! sidek=1 (for triangles), sidek=2 (for quadrilaterals), sidek=3 (for pentagons)
! sidek=4 (for hexagons)
sidek = nnodes - 2                   
do n=1,nnodes                     
    nod = BoundaryFace(nf)%Node(n)%name
    do i=1,SPACEDIM
       BoundaryFace(nf)%Node(n)%GX(i) = Vertice(i,nod)
    enddo
enddo
! To define the components of the main sides
U1(1:SPACEDIM) = BoundaryFace(Nf)%Node(nindex(sidek,1))%GX(1:SPACEDIM) -       &
   BoundaryFace(Nf)%Node(nindex(sidek,3))%GX(1:SPACEDIM)
U2(1:SPACEDIM) = BoundaryFace(Nf)%Node(nindex(sidek,2))%GX(1:SPACEDIM) -       &
   BoundaryFace(Nf)%Node(nindex(sidek,3))%GX(1:SPACEDIM)
! Length of side U1
U1len = dsqrt(U1(1) * U1(1) + U1(2) * U1(2) + U1(3) * U1(3)) 
! To compute directional cosines of side W12
RR(1:SPACEDIM) = U1(1:SPACEDIM) / U1len     
! Vector product W12=W12xU2 
call Vector_Product(U1,U2,W12,SPACEDIM) 
! Length of vector W12
W12len = dsqrt(W12(1) * W12(1) + W12(2) * W12(2) + W12(3) * W12(3))
! Area of the face "nf" (denominator=2 for triangles, =1 for parallelograms) 
! The estimation of pentagon and hexagon face areas is not mandatory as they are
! only considered for "perimeter" zones, not for SASPH frontiers.
if (nnodes<=4) BoundaryFace(Nf)%Area = W12len / float(3 - sidek)
! Directional cosines of the normal to the face "nf"     
nnlocal(1:SPACEDIM) = W12(1:SPACEDIM) / W12len
! Vector product ss=rrxnn; ss is the unity vector of the third local axis      
call Vector_Product(nnlocal,RR,ss,SPACEDIM)   
! Directional cosine matrix of face "nf"
BoundaryFace(Nf)%T(1:SPACEDIM,1) = RR(1:SPACEDIM)
BoundaryFace(Nf)%T(1:SPACEDIM,2) = ss(1:SPACEDIM)
BoundaryFace(Nf)%T(1:SPACEDIM,3) = nnlocal(1:SPACEDIM)
! Local coordinates of face nodes
refnode = nnodes        
do n=1,(nnodes-1)
   do i=1,PLANEDIM
      LocX = zero
      do j=1,SPACEDIM
         LocX = LocX + BoundaryFace(Nf)%T(j,i) *                               &
               (BoundaryFace(Nf)%Node(n)%GX(j) -                               &
               BoundaryFace(Nf)%Node(refnode)%GX(j))
      enddo
      BoundaryFace(nf)%Node(n)%LX(i) = LocX
   enddo
enddo
do i=1,PLANEDIM
   BoundaryFace(nf)%Node(refnode)%LX(i) = zero
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine DefineLocalSystemVersors
#endif
