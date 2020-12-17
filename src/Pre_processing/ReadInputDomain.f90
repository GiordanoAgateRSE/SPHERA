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
! Program unit: ReadInputDomain                       
! Description:                        
!-------------------------------------------------------------------------------
subroutine ReadInputDomain(ainp,comment,nrighe,ier,ninp,ulog)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nrighe,ier, ninp,ulog
character(1) :: comment
character(len=lencard) :: ainp
integer(4) :: ioerr
double precision :: dx, trunc
character(100) :: token
logical,external :: ReadCheck
character(100),external :: lcase,GetToken
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
! In case of restart, this input file section is read just once (not twice as 
! for regular runs). restart=.false. during the first reading of the main input
! file, even for restarted simulations.
if (restart) then
   do while (trim(lcase(ainp))/="##### end domain #####")
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DOMAIN DATA",ninp,ulog))       &
         return
   enddo
   return
endif
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DOMAIN DATA",ninp,ulog)) return
do while (trim(lcase(ainp))/="##### end domain #####")
   token = lcase(GetToken(ainp,1,ioerr))
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DOMAIN TYPE",ninp,ulog)) return 
   select case (token(1:4))
      case ("bsph","semi") 
         Domain%tipo = token(1:4)
         if ((input_second_read.eqv..true.).and.(ulog>0)) then
            write(ulog,"(1x,a,1x,a)"  ) "Domain Type            : ",           &
               trim(token)
         endif      
      case default
         ier = 3
         return
   endselect
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (ioerr==0) read(ainp,*,iostat=ioerr) dx,trunc
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DD & TRUNC",ninp,ulog)) return
   Domain%dx = dx
   input_any_t%trunc = trunc
   token = lcase(GetToken(ainp,3,ioerr))
   if (token(1:1)=='r') then
      Domain%RandomPos = 'r'
      else
         Domain%RandomPos = 'n'
   endif
   if ((input_second_read.eqv..true.).and.(ulog>0)) then
      write(ulog,"(1x,a,1pe12.4)") "dx                     : ",dx
      write(ulog,"(1x,a,1pe12.4)") "Trunc                  : ",trunc
      write(ulog,"(1x,a,1x,a)") "Random Initial Position: ",Domain%RandomPos
      write(ulog,"(1x,a)")  " "
   endif
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DOMAIN DATA",ninp,ulog)) return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputDomain
