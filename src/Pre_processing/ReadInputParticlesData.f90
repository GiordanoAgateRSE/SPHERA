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
! Program unit: ReadInputParticlesData                          
! Description:                        
!-------------------------------------------------------------------------------
subroutine ReadInputParticlesData(NumberEntities,Medium,icolor,bends,move,slip,&
                                  npointv,valuev,values3,pressu,valp,ainp,     &
                                  comment,nrighe,ier,ninp,ulog)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(inout) :: Medium,icolor,npointv,nrighe,ier,ninp,ulog
double precision,intent(inout) :: valp
integer(4),dimension(20),intent(inout) :: NumberEntities
double precision,dimension(3),intent(inout) :: values3
double precision,dimension(0:3,maxpointsvlaw),intent(inout) :: valuev
character(1),intent(inout) :: bends,slip,comment
character(2),intent(inout) :: pressu
character(3),intent(inout) :: move
character(len=lencard),intent(inout) :: ainp
integer(4) :: ioerr,i,n,icord
character(100) :: token
character(6) :: token_color
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
read(ainp,*,iostat=ioerr) Medium
if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"MEDIUM INDEX",ninp,ulog)) return
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
token = GetToken(ainp,1,ioerr)
if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"COLOURING TYPE",ninp,ulog)) return
bends = trim(lcase(token(1:1)))
token = GetToken(ainp,2,ioerr)
! Saving colours 
token_color(1:2) = token(5:6)
token_color(3:4) = token(3:4)
token_color(5:6) = token(1:2) 
read(token_color,'(Z6)',iostat=ioerr) icolor
if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"COLOUR",ninp,ulog)) return
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
move = trim(GetToken(ainp,1,ioerr))
if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL MOVE STATUS TYPE",ninp,ulog))&
   return
values3(:) = zero
slip = " "                             
npointv = 0
select case (lcase(move))
   case ("std")
      npointv = 0
      do n=1,3
         token = GetToken(ainp,(n+1),ioerr)
         read(token,*,iostat=ioerr) values3(n)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL VELOCITY",ninp,ulog)&
            ) return
      enddo
   case ("fix")
      npointv = 1
      valuev = zero
      do n=1,ncord
         token = GetToken(ainp,(n+1),ioerr)
         read(token,*,iostat=ioerr) values3(n)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL VELOCITY",ninp,ulog &
            )) return
      enddo
! No-slip / free-slip
      token = GetToken(ainp,(ncord+2),ioerr)
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"NO_SLIP/FREE_SLIP",ninp,ulog)) &
         return
      if (lcase(token(1:7))=="no_slip") then
         slip = "n"
         elseif (lcase(token(1:9))=="free_slip") then
            slip = "f"
            elseif (lcase(token(1:9))=="cont_slip") then  
               slip = "c"
               NumberEntities(19) = 1
               else
                  slip = "?"
                  ioerr= 1
                  if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                    &
                     "NO_SLIP/FREE_SLIP/CONT_SLIP",ninp,ulog)) return
      endif
   case ("law")
      ioerr   = 0
      npointv = 0
      valuev  = zero  
      token = GetToken(ainp,2,ioerr)
      read(token,*,iostat=ioerr) npointv
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"POINTS OF VELOCITY LAW",ninp,  &
         ulog)) return
      slip = "n"
      do i=1,npointv
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"VALUES OF VELOCITY LAW",    &
            ninp,ulog)) return
         if ((ioerr/=0).or.(input_second_read.eqv..false.)) cycle
         do n=0,ncord
            icord = icoordp(n,ncord-1)
            token = GetToken(ainp,(n+1),ioerr)
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"VALUES OF VELOCITY LAW", &
               ninp,ulog)) return
            read(token,*,iostat=ioerr) valuev(n,i)
         enddo
         if (i==1) values3(1:3) = valuev(1:3,1)
      enddo
   case default
      if (ulog>0) write(ulog,*) "Unknown option: ",trim(ainp)
      stop
endselect
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
pressu = trim(GetToken(ainp,1,ioerr))
if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL PRESSURE TYPE",ninp,ulog))   &
   return
! IC uniform pressure ("pa") is assigned 
if (pressu=="pa") then  
   token = lcase(GetToken(ainp,2,ioerr))
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"PRESSURE VALUES",ninp,ulog)) return
   read(token,*,iostat=ioerr) valp
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"PRESSURE VALUES",ninp,ulog)) return
! IC hydrostatic pressure based on a reference free surface height ("qp")
   elseif (pressu=="qp") then 
      token = lcase(GetToken(ainp,2,ioerr))
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL PIEZO LINE",ninp,ulog))&
         return
      read(token,*,iostat=ioerr) valp
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL PIEZO LINE",ninp,ulog))&
         return
! IC hydrostatic pressure based on the maximum level ("pl") of
! an assigned fluid 
         elseif (pressu=="pl") then
            token = lcase(GetToken(ainp,2,ioerr))
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"FREE LEVEL LINE",ninp,   &
               ulog)) return
            read(token,*,iostat=ioerr) valp
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"FREE LEVEL LINE",ninp,   &
               ulog)) return
            else
               if (ulog>0) write(ulog,*) "Unknown option: ",trim(ainp)
               stop
endif
if (input_second_read.eqv..true.) then
   if ((bends=="b").and.(icolor>5)) then
      icolor = 5
      if (ulog>0) write(ulog,*) "Maximum number of bends is 5!"
   endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputParticlesData
