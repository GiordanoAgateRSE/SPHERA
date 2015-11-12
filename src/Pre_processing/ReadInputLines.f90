!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-Ricerca di Sistema)



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
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: ReadInputLines                          
! Description:                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine ReadInputLines(NumberEntities,BoundaryVertex,Tratto,ainp,comment,   &
                          nrighe,ier,ninp,nout)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nrighe,ier, ninp,nout
integer(4),dimension(20) :: NumberEntities
integer(4),dimension(NumBVertices) :: BoundaryVertex
type (TyBoundaryStretch),dimension(NumTratti) :: Tratto
character(1) :: comment
character(80) :: ainp
integer(4),parameter :: MAXLINENODES = 20
integer(4) :: n,i,ioerr,i1,index,numv,numv_line,ipointer
character(5) :: txt
character(80) :: token
logical,external :: ReadCheck
character(80),external :: lcase, GetToken
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
! In case of restart, input data are not read
if (restart) then
   do while (TRIM(lcase(ainp))/="##### end lines #####")
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"LINES DATA",ninp,nout)) return
   enddo
   return
endif
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"LINES DATA",ninp,nout)) return
do while (TRIM(lcase(ainp))/="##### end lines #####")
   select case (TRIM(Domain%tipo))
      case ("semi","bsph") 
! Reading the boundary vertices 
         numv = 0
! Reading the line index  
         numv_line = 1
         token = GetToken(ainp,numv_line,ioerr)
         if (ioerr==0) read (token,*,iostat=ioerr) index
! "NumTratti"
         NumberEntities(8) = max(NumberEntities(8),index)  
         VERTEX_LOOP: do while (ioerr==0)
            numv_line = numv_line + 1
            token = GetToken(ainp,numv_line,ioerr)
! It exits when either it finds a comment beginning, or the EOR 
! (no more data on the input line)
            if ((ioerr/=0).OR.(trim(token)=="").OR.(ichar(trim(token(1:1)))==9)&
               .OR.(token(1:1)=="!")) exit VERTEX_LOOP
            if ((token(1:1)=="&").OR.(token(1:1)==">")) then
               call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
               if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "VERTICES LIST (continue...)",ninp,nout)) return
               numv_line = 0
               cycle VERTEX_LOOP
            endif
            numv  = numv + 1
            if (numv>MAXLINENODES) then
               stop 'ERRORE in ReadInputLines numv>MAXLINENODES'
            endif
! Counter of the boundary vertices
            NumberEntities(9) = NumberEntities(9) + 1  
            read (token,*,iostat=ioerr) i
            write(txt,"(i5)") i
            if (numv==1) i1 = i
            if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VERTEX n."//txt,ninp,nout&
               )) return
            if (ncord>0) then
               if (numv==1) ipointer = NumberEntities(9)
               BoundaryVertex(NumberEntities(9)) = i
            endif
         enddo VERTEX_LOOP
! To count the number of "BoundarySide"
         NumberEntities(10) = NumberEntities(10) + numv - 1
         if (ncord>0) then
            Tratto(index)%numvertices = numv
            Tratto(index)%inivertex   = ipointer
         endif
      case default
         if (nout>0) then
            write (nout,*) "Unknown Domain Type: ",Domain%tipo
         endif
         ier = 2
         return
   end select
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"LINES DATA",ninp,nout)) return
enddo
if ((ncord>0).AND.(nout>0)) then
   write (nout,"(1x,a)") "List of lines"
   write (nout,*)
   do n=1,NumberEntities(8)
      write (nout,"(1x,a,i3,1x,a)") "Line: ",n
      write (nout,"(1x,a,i3,1x,a)") "Number of Vertices:  ",                   &
         Tratto(n)%numvertices
      write (nout,"(1x,a,i3,1x,a)") "Vertices Pointer:    ",Tratto(n)%inivertex
      write (nout,"(1x,a,i3,1x,a)") "Vertices List"
      write (nout,"(1x,10i5)")                                                 &
BoundaryVertex(Tratto(n)%inivertex:Tratto(n)%inivertex+Tratto(n)%numvertices-1)
      write (nout,"(1x,a)") " "
   enddo
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputLines

