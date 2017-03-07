!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
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
! Program unit: ReadInputParticlesData                          
! Description:                        
!-------------------------------------------------------------------------------
subroutine ReadInputParticlesData(NumberEntities,Medium,icolor,bends,move,slip,&
                                  npointv,valuev,values3,pressu,valp,ainp,     &
                                  comment,nrighe,ier,ninp,nout)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4):: nrighe,ier,ninp,nout,Medium,npointv,icolor
integer(4),dimension(20) :: NumberEntities
double precision :: valp
character(1) :: comment
character(100) :: ainp
character(1) :: bends,slip
character(3) :: move
character(2) :: pressu
double precision,dimension(3) :: values3
double precision,dimension(0:3,maxpointsvlaw) :: valuev
integer(4) :: ioerr,i,n,icord
character(100) :: token
character(6) :: token_color
logical,external :: ReadCheck
character(100),external :: lcase, GetToken
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
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"MEDIUM INDEX",ninp,nout)) return
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
token = GetToken(ainp,1,ioerr)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"COLOURING TYPE",ninp,nout)) return
bends = lcase(token(1:1))
token = GetToken(ainp,2,ioerr)
! Saving colours 
token_color(1:2) = token(5:6)
token_color(3:4) = token(3:4)
token_color(5:6) = token(1:2) 
read(token_color,'(Z6)',iostat=ioerr) icolor
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"COLOUR",ninp,nout)) return
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
move = GetToken(ainp,1,ioerr)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL MOVE STATUS TYPE",ninp,nout))&
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
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL VELOCITY",ninp,nout)&
            ) return
      enddo
   case ("fix")
      npointv = 1
      valuev  = zero  
      do n=1,NumberEntities(1)
         token = GetToken(ainp,(n+1),ioerr)
         read(token,*,iostat=ioerr) values3(n)
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL VELOCITY",ninp,nout &
            )) return
      enddo
! No-slip / free-slip 
      token = GetToken(ainp,(NumberEntities(1)+2),ioerr)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"NO_SLIP/FREE_SLIP",ninp,nout)) &
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
                  if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                    &
                     "NO_SLIP/FREE_SLIP/CONT_SLIP",ninp,nout)) return
      endif
   case ("law")
      ioerr   = 0
      npointv = 0
      valuev  = zero  
      token = GetToken(ainp,2,ioerr)
      read(token,*,iostat=ioerr) npointv
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"POINTS OF VELOCITY LAW",ninp,  &
         nout)) return
      slip = "n"
      do i=1,npointv
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VALUES OF VELOCITY LAW",    &
            ninp,nout)) return
         if ((ioerr/=0).or.(ncord<=0)) cycle
         do n=0,NumberEntities(1)
            icord = icoordp(n,ncord-1)
            token = GetToken(ainp,(n+1),ioerr)
            if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VALUES OF VELOCITY LAW", &
               ninp,nout)) return
            read(token,*,iostat=ioerr) valuev(n,i)
         enddo
         if (i==1) values3(1:3) = valuev(1:3,1)
      enddo
   case default
      if (nout>0) write(nout,*) "Unknown option: ",trim(ainp)
      stop
endselect
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
pressu = GetToken(ainp,1,ioerr)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL PRESSURE TYPE",ninp,nout))   &
   return
! IC uniform pressure ("pa") is assigned 
if (pressu=="pa") then  
   token = lcase(GetToken(ainp,2,ioerr))
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"PRESSURE VALUES",ninp,nout)) return
   read(token,*,iostat=ioerr) valp
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"PRESSURE VALUES",ninp,nout)) return
! IC hydrostatic pressure based on a reference free surface height ("qp")
   elseif (pressu=="qp") then 
      token = lcase(GetToken(ainp,2,ioerr))
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL PIEZO LINE",ninp,nout))&
         return
      read(token,*,iostat=ioerr) valp
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL PIEZO LINE",ninp,nout))&
         return
! IC hydrostatic pressure based on the maximum level ("pl") of
! an assigned fluid 
         elseif (pressu=="pl") then      
            token = lcase(GetToken(ainp,2,ioerr))
            if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FREE LEVEL LINE",ninp,   &
               nout)) return
            read(token,*,iostat=ioerr) valp
            if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FREE LEVEL LINE",ninp,   &
               nout)) return
            else
               if (nout>0) write(nout,*) "Unknown option: ",trim(ainp)
               stop
endif
if (ncord>0) then
   if ((bends=="b").AND.(icolor>5)) then
      icolor = 5
      if (nout>0) write(nout,*) "Maximum number of bends is 5!"
   endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputParticlesData

