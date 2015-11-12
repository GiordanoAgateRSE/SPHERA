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
! Program unit: SearchforParticleZone_3D              
! Description: It returns in "partizone" the highest index of wet cells. In case
!              no cell is wet, "partzone = sourzone" ("sourzone"is the inlet section cell).
!----------------------------------------------------------------------------------------------------------------------------------

subroutine SearchforParticleZone_3D(partizone)
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
integer(4),intent(INOUT) :: partizone
integer(4) :: iz,sourzone,mate
character(4) :: tipo
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
partizone = 0
sourzone = 0
!------------------------
! Statements
!------------------------
do iz=NPartZone,1,-1
   tipo = Partz(iz)%tipo
   if (tipo/="sour") then
      mate = Partz(iz)%Medium
      if (mate>0) then
         partizone = iz
         exit
      endif
      else
         sourzone = iz
   endif
enddo
if (partizone==0) partizone = sourzone
!------------------------
! Deallocations
!------------------------
return
end subroutine SearchforParticleZone_3D

