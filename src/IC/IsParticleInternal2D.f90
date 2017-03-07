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
! Program unit: IsParticleInternal2D
! Description: To check whether a particle is internal to the 2D domain. 
!-------------------------------------------------------------------------------
logical function IsParticleInternal2D (Nt,PX)
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
integer(4),intent(IN) :: Nt
double precision,intent(IN),dimension(SPACEDIM) :: PX
integer(4) :: inizio,fine,iv,ni,n,n2
double precision :: xa,za,xba,zba,xi,zs
character(len=lencard) :: nomsub = "IsParticleInternal2D"
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
IsParticleInternal2D = .FALSE.
ni = 0
inizio = Tratto(Nt)%inivertex
fine   = Tratto(Nt)%inivertex + Tratto(Nt)%numvertices - 2
!------------------------
! Statements
!------------------------
! Loop on the line vertices
do iv=inizio,fine
! To Check for the storage limits
   if (iv>NumBVertices) then
      call diagnostic(arg1=10,arg2=2,arg3=nomsub)
   end if
! To set the pointer to the current vertex
   n = BoundaryVertex(iv)
   if (n>NumVertici) then
     call diagnostic(arg1=10,arg2=3,arg3=nomsub)
   end if
! To set the coordinates of the current vertex
   xa = Vertice(1,n)
   za = Vertice(3,n)
! To set the next vertex and its coordinates
   n2 = BoundaryVertex(iv+1)    
   xba = Vertice(1,n2)
   zba = Vertice(3,n2)
! The current segment is not vertical
   if (abs(xa-xba)>=xyz_tolerance) then
! The current segment is not horizontal
      if (abs(za-zba)>=xyz_tolerance) then
! It evaluates the x coordinate of the particle projection on the segment 
! along X
         xi = xa + (xba - xa) * (px(3) - za) / (zba - za)
         else
! The segment is horizontal: the X value of the mean point is assumed
            xi = half * (xa+xba)
      end if
      else
! The segment is vertical: the X value of the vertices is assumed
         xi = xa
   end if
! To order the vertices to have the first one with the lower Z value
   if (za>zba) then
      zs=za
      za=zba
      zba=zs
   end if
! The Z value of the particle is inside the segment Z values
   if (((PX(3)-za)>xyz_tolerance).AND.((PX(3)-zba)<xyz_tolerance)) then
      if (xi>PX(1)) then
         ni = ni + 1
      end if
   end if
end do
if (MOD(ni,2)==1) then
   IsParticleInternal2D = .TRUE.
end if
!------------------------
! Deallocations
!------------------------
return
end function IsParticleInternal2D

