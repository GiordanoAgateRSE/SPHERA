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
! Program unit: ReadInputLines                          
! Description:                        
!-------------------------------------------------------------------------------
#ifdef SPACE_2D
subroutine ReadInputLines(NumberEntities,BoundaryVertex,Tratto,ainp,comment,   &
                          nrighe,ier,ninp,ulog)
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
integer(4),dimension(20) :: NumberEntities
integer(4),dimension(NumBVertices) :: BoundaryVertex
type (TyBoundaryStretch),dimension(NumTratti) :: Tratto
character(1) :: comment
character(len=lencard) :: ainp
integer(4),parameter :: MAXLINENODES = 820
integer(4) :: n,i,ioerr,i1,index,numv,numv_line,ipointer
character(5) :: txt
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
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"LINES DATA",ninp,ulog)) return
do while (trim(lcase(ainp))/="##### end lines #####")
   select case (trim(Domain%tipo))
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
            if ((ioerr/=0).or.(trim(token)=="").or.(ichar(trim(token(1:1)))==9)&
               .or.(token(1:1)=="!")) exit VERTEX_LOOP
            if ((token(1:1)=="&").or.(token(1:1)==">")) then
               call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "VERTICES LIST (continue...)",ninp,ulog)) return
               numv_line = 0
               cycle VERTEX_LOOP
            endif
            numv = numv + 1
            if (numv>MAXLINENODES) then
               stop 'ERRORE in ReadInputLines numv>MAXLINENODES'
            endif
! Counter of the boundary vertices
            NumberEntities(9) = NumberEntities(9) + 1  
            read (token,*,iostat=ioerr) i
            write(txt,"(i5)") i
            if (numv==1) i1 = i
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"VERTEX n."//txt,ninp,ulog&
               )) return
            if (input_second_read.eqv..true.) then
               if (numv==1) ipointer = NumberEntities(9)
               BoundaryVertex(NumberEntities(9)) = i
            endif
         enddo VERTEX_LOOP
         NumberEntities(10) = NumberEntities(10) + numv - 1
         if (input_second_read.eqv..true.) then
            Tratto(index)%numvertices = numv
            Tratto(index)%inivertex = ipointer
         endif
      case default
         if (ulog>0) then
            write(ulog,*) "Unknown Domain Type: ",Domain%tipo
         endif
         ier = 2
         return
   endselect
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"LINES DATA",ninp,ulog)) return
enddo
if ((input_second_read.eqv..true.).and.(ulog>0)) then
   write(ulog,"(1x,a)") "List of lines"
   write(ulog,*)
   do n=1,NumberEntities(8)
      write(ulog,"(1x,a,i3,1x,a)") "Line: ",n
      write(ulog,"(1x,a,i3,1x,a)") "Number of Vertices:  ",                    &
         Tratto(n)%numvertices
      write(ulog,"(1x,a,i3,1x,a)") "Vertices Pointer:    ",Tratto(n)%inivertex
      write(ulog,"(1x,a,i3,1x,a)") "Vertices List"
      write(ulog,"(1x,10i5)")                                                  &
BoundaryVertex(Tratto(n)%inivertex:Tratto(n)%inivertex+Tratto(n)%numvertices-1)
      write(ulog,"(1x,a)") " "
   enddo
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputLines
#endif
