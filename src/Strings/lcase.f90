!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)

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
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Program unit: lcase                                           
! Description: 
!-------------------------------------------------------------------------------
character(100) function lcase(ainp)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: n,ia
character(*) :: ainp
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
do n=1,100
   lcase(n:n) = " "
enddo
!------------------------
! Statements
!------------------------
do n=1,len_trim(ainp)
   ia = iachar(ainp(n:n))
   if ((ia>=65).and.(ia<=90)) then
      lcase(n:n) = char(ia+32)
      else
         lcase(n:n) = ainp(n:n)
   endif
enddo
!------------------------
! Deallocations
!------------------------
return
end function lcase

