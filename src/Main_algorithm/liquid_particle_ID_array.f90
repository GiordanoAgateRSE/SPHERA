!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2020 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: liquid_particle_ID_array
! Description: To count the liquid particles and update the associated ID array    
!-------------------------------------------------------------------------------
subroutine liquid_particle_ID_array
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Dynamic_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: npi
character(len=lencard) :: nomsub = "liquid_particle_ID_array"
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
indarrayFlu = 0
do npi=1,nag
   if ((pg(npi)%cella==0).or.(pg(npi)%vel_type/="std")) cycle
   if (pg(npi)%state=="flu") then
      indarrayFlu = indarrayFlu + 1
! To check the maximum dimension of the array and possible resizing
      if (indarrayFlu>PARTICLEBUFFER) then
         call diagnostic(arg1=9,arg2=1,arg3=nomsub)
      endif
      Array_Flu(indarrayFlu) = npi
   endif
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine liquid_particle_ID_array
