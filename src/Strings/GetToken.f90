!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2021 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)

! SPHERA authors and email contact are provided on SPHERA documentation.

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
! Program unit: GetToken                                           
! Description: To extract the "itok"-th token from the input string "ainp"
!-------------------------------------------------------------------------------
character(100) function GetToken(itok,ainp,io_err)
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: itok
character(*),intent(in) :: ainp
integer(4),intent(out) :: io_err
logical :: blank
integer(4) :: n,number_token
integer(4),dimension(2,820) :: index_token
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
blank = .true.
!------------------------
! Statements
!------------------------
do n=1,len_trim(ainp)
   if ((blank).and.(ainp(n:n)/=" ")) then
      number_token = number_token + 1
      index_token(1,number_token) = n
      index_token(2,number_token) = n
      blank = .false.
      elseif ((.not.blank).and.(ainp(n:n)/=" ")) then
         index_token(2,number_token) = n
         elseif (ainp(n:n)==" ") then 
            blank = .true.
   endif
enddo
if (itok<=number_token) then
   io_err = 0
   GetToken = ainp(index_token(1,itok):index_token(2,itok))
   else
      io_err = itok
      GetToken = ""
endif
!------------------------
! Deallocations
!------------------------
return
end function GetToken
