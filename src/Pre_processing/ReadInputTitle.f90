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
! Program unit: ReadInputTitle                             
! Description:                        
!-------------------------------------------------------------------------------
subroutine ReadInputTitle(ainp,comment,nrighe,ier,ninp,ulog)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nrighe,ier,ninp,ulog
character(1) :: comment
character(len=lencard) :: ainp
integer(4) :: n,ioerr
character(100),external :: lcase
logical,external :: ReadCheck
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
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"TITLE DATA",ninp,ulog)) return 
n = 0
do while (trim(lcase(ainp))/="##### end title #####" )
   n = n + 1
   if (n<=maxtit) then
      title(n) = ainp
      if ((input_second_read.eqv..true.).and.(ulog>0)) write(ulog,"(1x,a)") title(n)
   endif
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"TITLE DATA",ninp,ulog)) return 
enddo
if ((input_second_read.eqv..true.).and.(ulog>0)) write(ulog,"(1x,a)") " "
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputTitle
