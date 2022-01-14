!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2021 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: ReadInputRestart                           
! Description: To read the restart parameters from the main input file.                        
!-------------------------------------------------------------------------------
subroutine ReadInputRestart(ainp,comment,nrighe,ier)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use I_O_file_module
use Memory_I_O_interface_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nrighe,ier
character(1) :: comment
character(len=lencard) :: ainp
logical :: restartOK
integer(4) :: ioerr,i_tok
character(100) :: token,file_name
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
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
if (.not.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART DATA",ninp,ulog)) return
do while (trim(lcase(ainp))/="##### end restart #####")
   select case (trim(lcase(GetToken(ainp,1,ioerr))))
      case ("step")
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"RESTART STEP ",ninp,ulog))  &
            return
         token = lcase(GetToken(ainp,1,ioerr))
         read(token,*,iostat=ioerr) Domain%istart
         if (.not.ReadCheck (ioerr,ier,nrighe,ainp,"Restart step",             &
            ninp,ulog)) return
         if ((input_second_read.eqv..true.).and.(ulog>0)) then
            write(ulog,"(1x,a,i12)") "Restart from step: ",Domain%istart
            if (Domain%istart<0) write(ulog,"(1x,a)") "Negative restart step!"
! The time mode is deactivated
            if (Domain%start>zero) then
               write(ulog,"(1x,a,f20.12,a)") "Restart from time: ",            &
                  Domain%start," option ignored!"
               Domain%start = zero
            endif
         endif
      case ("time")
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"RESTART TIME ",ninp,ulog))  &
            return
         token = lcase(GetToken(ainp,1,ioerr))
         read(token,*,iostat=ioerr) Domain%start
         if (.not.ReadCheck (ioerr,ier,nrighe,ainp,"Restart time ",ninp,ulog)) &
            return
         if ((input_second_read.eqv..true.).and.(ulog>0)) then
            write(ulog,"(1x,a,f20.12)") "Restart from time: ",Domain%start
            if (Domain%start<zero) write(ulog,"(1x,a)")                        &
               "Negative restart time!"
! The step mode is deactivated
            if (Domain%istart>0) then
               write(ulog,"(1x,a,i12,a)") "Restart from step: ",Domain%istart, &
                  " option ignored!"
               Domain%istart = 0 
            endif
         endif
      case ("last_time")
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"RESTART PATH ",ninp,ulog))  &
            return
         token = GetToken(ainp,1,ioerr)
         read(token,*,iostat=ioerr) Domain%restart_path
         if (.not.ReadCheck (ioerr,ier,nrighe,ainp,"Restart path ",ninp,ulog)) &
            return
! Replace "*" with "/"
         do i_tok=1,100
            if (iachar(Domain%restart_path(i_tok:i_tok))==iachar("*")) then
               Domain%restart_path(i_tok:i_tok) = "/"
            endif
         enddo
! Read the last restart time
         file_name = trim(adjustl(Domain%restart_path)) //                     &
                     "/restart_last_time.rist"
         call open_close_file(.true.,max_file_unit_booked+1,file_name)
         read(max_file_unit_booked+1,*,iostat=ioerr) Domain%start
! Restart input time has to be smaller than the actual restart time 
         Domain%start = Domain%start * (1.d0 - 1.d-9)
         if (.not.ReadCheck(ioerr,ier,1,file_name,                             &
            "Domain%start from Domain%restart_path",max_file_unit_booked+1,    &
            ulog)) then
            write(uerr,*) "Error in reading the file ",file_name,              &
               ". The execution stops here."
            stop
         endif
         call open_close_file(.false.,max_file_unit_booked+1,file_name)
         if ((input_second_read.eqv..true.).and.(ulog>0)) then
            write(ulog,"(1x,a,a)") "File for restart time: ",file_name
            write(ulog,"(1x,a,f20.12)") "Restart from time: ",Domain%start
! The step mode is deactivated
            if (Domain%istart>0) then
               write(ulog,"(1x,a,i12,a)") "Restart from step: ",Domain%istart, &
                  " option ignored!"
               Domain%istart = 0 
            endif
         endif
      case default
! Input line for the restart file name
         Domain%file = ainp
         file_name = trim(adjustl(Domain%restart_path)) // "/" //              &
            trim(adjustl(Domain%file))
         if ((input_second_read.eqv..true.).and.(ulog>0)) then
            inquire(file=file_name,exist=restartOK)
            if (restartOK) then
               write(ulog,"(1x,2a)") "Restart file: ",trim(adjustl(file_name))
               else
                  write(ulog,"(1x,3a)") "Restart file: ",                      &
                     trim(adjustl(file_name)), " not found!"
            endif
         endif
   endselect
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"RESTART DATA",ninp,ulog)) return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputRestart
