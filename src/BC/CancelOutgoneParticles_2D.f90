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
! Program unit: CancelOutgoneParticles_2D
! Description: To count and delete the outgoing particles on boundaries of type "leve", "flow", "velo", "crit", "open".
!----------------------------------------------------------------------------------------------------------------------------------

subroutine CancelOutgoneParticles_2D
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
integer(4) :: ios,v1,pd,npi,isi,v2
double precision :: detV1PnewV2,detV1PoldV2,detPoldV1Pnew,detPoldV2Pnew
integer(4),dimension(1:PLANEDIM) :: acix
double precision,dimension(1:PLANEDIM) :: OP1,op2
double precision,dimension(1:PLANEDIM) :: Plocalnew,Plocalold
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
! 2D scheme uses the X and Z axes
acix(1)=1  
acix(2)=3
!------------------------
! Statements
!------------------------
! All the opened boundary sides are considered and the relative position of 
! all the particles is checked
do ios=1,NumOpenSides
   isi = OpenSide(ios)
   v1 = BoundarySide(isi)%Vertex(1)
   v2 = BoundarySide(isi)%Vertex(2)
! Load the reference vertex coordinates and local side versor (tangent and 
! normal)
   do pd=1,PLANEDIM
      OP1(pd) = Vertice(acix(pd),v1)
      OP2(pd) = Vertice(acix(pd),v2)
   end do
!$omp parallel do default(none) &
!$omp private(npi,ios,Plocalnew,Plocalold,detV1PnewV2,detV1PoldV2)             &
!$omp private(detPoldV1Pnew,detPoldV2Pnew)                                     &
!$omp shared(nag,pg,acix,OP1,OP2,OpCount)
! Loop over all the particles
   do npi = 1,nag
      if (pg(npi)%cella == 0) cycle
! Load the particle coordinates
      Plocalnew(1:planedim) = pg(npi)%coord(acix(1:planedim))
      Plocalold(1:planedim) = pg(npi)%CoordOld(acix(1:planedim))
! Evaluates the determinant of the outside triplet V1-Pnew-V2
      detV1PnewV2 = (Plocalnew(2) - op1(2)) * (op2(1) - op1(1)) -              &
                    (plocalnew(1) - op1(1)) * (op2(2) - op1(2))
! If the normal component is smaller than zero, then the particle might be out 
! of the boundary side, since the reference normal is oriented inside the domain
      if (detV1PnewV2 <= zero) then
! So the crossing point between the last path of the particle and the boundary
! segment is tested
         detV1PoldV2 = (Plocalold(2) - op1(2)) * (op2(1) - op1(1)) -           &
                       (plocalold(1) - op1(1)) * (op2(2) - op1(2))
        if (sign(one,detV1PnewV2)==sign(one,detV1PoldV2)) cycle
        detPoldV1Pnew = (op1(2) - plocalold(2)) * (plocalnew(1) - plocalold(1))&
                        - (op1(1) - plocalold(1)) * (plocalnew(2) -            &
                        plocalold(2))
        detPoldV2Pnew = (op2(2) - plocalold(2)) * (plocalnew(1) - plocalold(1))&
                        - (op2(1) - plocalold(1)) * (plocalnew(2) -            &
                        plocalold(2))
        if (sign(one,detPoldV1Pnew)==sign(one,detPoldV2Pnew)) cycle
        OpCount(pg(npi)%imed) = OpCount(pg(npi)%imed) + 1 
        pg(npi)%cella = -1
      end if
   end do
!$omp end parallel do
end do
!------------------------
! Deallocations
!------------------------
return
end subroutine CancelOutgoneParticles_2D

