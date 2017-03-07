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
! Program unit: GetToken                                           
! Description: 
!-------------------------------------------------------------------------------

character(100) function GetToken(ainp,itok,ioerr)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: itok,ioerr
character(*) :: ainp
logical :: blank
integer(4) :: n,number_token
integer(4),dimension(2,20) :: index_token
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
number_token = 0
index_token  = 0
blank = .TRUE.
!------------------------
! Statements
!------------------------
do n=1,len_trim(ainp)
   if ((blank).and.(ainp(n:n)/=" ")) then
      number_token = number_token + 1
      index_token(1,number_token) = n
      index_token(2,number_token) = n
      blank = .FALSE.
      elseif ((.not.blank).and.(ainp(n:n)/=" ")) then
         index_token(2,number_token) = n
         elseif (ainp(n:n)==" ") then 
         blank = .TRUE.
   endif
end do
if (itok<=number_token) then
   ioerr = 0
   GetToken = ainp(index_token(1,itok):index_token(2,itok))
   else
   ioerr = itok
   GetToken = ""
endif
!------------------------
! Deallocations
!------------------------
return
end function GetToken

