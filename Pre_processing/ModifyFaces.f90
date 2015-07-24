!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-; SPHERA has been authored for RSE SpA by 
!    Andrea Amicarelli, Antonio Di Monaco, Sauro Manenti, Elia Bon, Daria Gatti, Giordano Agate, Stefano Falappi, 
!    Barbara Flamini, Roberto Guandalini, David Zuccalà).
! Main numerical developments of SPHERA: 
!    Amicarelli et al. (2015,CAF), Amicarelli et al. (2013,IJNME), Manenti et al. (2012,JHE), Di Monaco et al. (2011,EACFM). 
! Email contact: andrea.amicarelli@rse-web.it

! This file is part of SPHERA.
! SPHERA is free software: you can redistribute it and/or modify
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
! Program unit: ModifyFaces                   
! Description: To generate triangles from quadrilaterals (partitioning along the shortest diagonal)                 
!----------------------------------------------------------------------------------------------------------------------------------

subroutine ModifyFaces (NumberEntities)
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
integer(4),intent(IN),dimension(20) :: NumberEntities
integer(4) :: n,i,new
double precision :: d13, d24
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
new = NumberEntities(11)
!------------------------
! Statements
!------------------------
do n=1,NumberEntities(11)
   if (BoundaryFace(n)%Node(4)%name==0) cycle
   new = new + 1
   d13 = zero
   d24 = zero
   do i=1,SPACEDIM
      d13 = d13 + (Vertice(i,BoundaryFace(n)%Node(1)%name) -                   &
            Vertice(i,BoundaryFace(n)%Node(3)%name)) *                         &
            (Vertice(i,BoundaryFace(n)%Node(1)%name) -                         &
            Vertice(i,BoundaryFace(n)%Node(3)%name))
      d24 = d24 + (Vertice(i,BoundaryFace(n)%Node(2)%name) -                   &
            Vertice(i,BoundaryFace(n)%Node(4)%name)) *                         &
            (Vertice(i,BoundaryFace(n)%Node(2)%name) -                         &
            Vertice(i,BoundaryFace(n)%Node(4)%name))
   enddo
   if (d13<d24) then
      BoundaryFace(new) = BoundaryFace(n)
      BoundaryFace(n)%Node(4)%name = -BoundaryFace(n)%Node(4)%name
      BoundaryFace(new)%Node(2)%name = BoundaryFace(new)%Node(3)%name
      BoundaryFace(new)%Node(3)%name = BoundaryFace(new)%Node(4)%name
      BoundaryFace(new)%Node(4)%name =-99999999
      else
         BoundaryFace(new) = BoundaryFace(n)
         i = BoundaryFace(n)%Node(1)%name
         BoundaryFace(n)%Node(1)%name = BoundaryFace(n)%Node(2)%name
         BoundaryFace(n)%Node(2)%name = BoundaryFace(n)%Node(3)%name
         BoundaryFace(n)%Node(3)%name = BoundaryFace(n)%Node(4)%name
         BoundaryFace(n)%Node(4)%name =-i
         BoundaryFace(new)%Node(3)%name = BoundaryFace(new)%Node(4)%name
         BoundaryFace(new)%Node(4)%name =-99999999
   endif
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ModifyFaces

