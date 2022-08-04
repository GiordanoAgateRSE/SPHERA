!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: ReadCheck
! Description: Check on a generic reading (oldest checks)             
!-------------------------------------------------------------------------------
logical function ReadCheck(IoErr,Ier,Nrighe,ainp,listadati,ninp,ulog)
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: IoErr,Ier,Nrighe,ninp,ulog
character(*) :: ainp,listadati
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
if (IoErr==0) then
   Ier = 0
   ReadCheck = .true.
   else
      Ier = 4
      ReadCheck = .false.
      write(ulog,"(1x,a)")    ">>>>>>>>>>>>>> Warning:"
      write(ulog,"(1x,a)")  
      write(ulog,"(1x,a,i5)") "Error reading unit:  ",ninp
      write(ulog,"(1x,a,a)")  "Expected data:       ",listadati
      write(ulog,"(1x,a,i8,a)")                                                &
         "Last input line read:"//trim(ainp)//"(line number:",Nrighe,")"
endif
!------------------------
! Deallocations
!------------------------
return
end function ReadCheck
