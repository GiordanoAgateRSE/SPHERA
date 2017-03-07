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
! Program unit: FindFrame
! Description: It finds extremes of the rectangular frame which contains the 
!              boundary mib. 
!-------------------------------------------------------------------------------
subroutine FindFrame(Xmin,Xmax,Nt)
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
integer(4),intent(IN) :: Nt
double precision,intent(INOUT),dimension(SPACEDIM,NumFacce) :: Xmin,Xmax
integer(4) :: i,n,iv,nf,nod
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
do iv=Tratto(Nt)%iniface,(Tratto(Nt)%iniface+Tratto(Nt)%numvertices-1)
   nf = BFaceList(iv)
   do n=1,6
      nod = BoundaryFace(nf)%Node(n)%name
      if (nod<=0) cycle
      do i=1,Ncord
         if (Vertice(i,nod)<Xmin(i,Nt)) Xmin(i,Nt) = Vertice(i,nod)
         if (Vertice(i,nod)>Xmax(i,Nt)) Xmax(i,Nt) = Vertice(i,nod)
      end do
   end do
end do
!------------------------
! Deallocations
!------------------------
return
end subroutine FindFrame

