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
! Program unit: ReadInputControlLines                    
! Description: Reading monitoring lines.                      
!-------------------------------------------------------------------------------
subroutine ReadInputControlLines(NumberEntities,Control_Points,Control_Lines,  &
                                 ainp,comment,nrighe,ier,ninp,nout)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nrighe,ier,ninp,nout
integer(4),dimension(20) :: NumberEntities
type (TyCtlPoint),dimension(NPointst) :: Control_Points
type (TyCtlLine),dimension(NLines) :: Control_Lines
character(1) :: comment
character(100) :: ainp
integer(4) :: n,i,ndiv,icord,ioerr,npts
double precision :: vp
double precision,dimension(3) :: values1,values2,values3
character(5) :: txt
character(8) :: label
logical,external :: ReadCheck 
character(100), external :: lcase
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
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL LINES DATA",ninp,nout))      &
   return
npts = npoints
do while (TRIM(lcase(ainp))/="##### end control lines #####")
   values1 = zero
   values2 = zero
   values3 = zero
   NumberEntities(5) = NumberEntities(5) + 1
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL LINE LABEL",ninp,nout))   &
      return
   label(1:8) = ainp(1:8)
   write(txt,"(i5)") NumberEntities(5)
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   read(ainp,*,iostat=ioerr) values1(1:ncord)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                                   &
      "CONTROL LINE"//txt//" - FIRST POINT",ninp,nout)) return
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   read(ainp,*,iostat=ioerr) values2(1:ncord)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                                   &
      "CONTROL LINE"//txt//" - SECOND POINT",ninp,nout)) return
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   read(ainp,*,iostat=ioerr) ndiv
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                                   &
      "CONTROL LINE"//txt//" - POINTS NUMBER",ninp,nout)) return
   NumberEntities(6) = NumberEntities(6) + ndiv
   if (ncord>0) then
      control_lines(NumberEntities(5))%label = label
      control_lines(NumberEntities(5))%icont(1) = npts + 1
      control_lines(NumberEntities(5))%icont(2) = npts + ndiv
      npts = npts + ndiv
      values3(:) = (values2(:) - values1(:)) / (ndiv - 1)
      vp = Dsqrt(values3(1) * values3(1) + values3(2) * values3(2) +           &
           values3(3) * values3(3))
      if (nout>0) then
         write (nout,"(1x,a,i3,1x,a)") "Control line      ",NumberEntities(5), &
            "("//control_lines(NumberEntities(5))%label//")"
         write (nout,"(1x,a,i12)") "First Point:      ",                       &
            control_lines(NumberEntities(5))%icont(1)
         write (nout,"(1x,a,i12)") "Last  Point:      ",                       &
            control_lines(NumberEntities(5))%icont(2)
      endif
      do i=control_lines(NumberEntities(5))%icont(1),                          &
         control_lines(NumberEntities(5))%icont(2)
         do n=1,ncord
            icord = icoordp(n,ncord-1)
            control_points(i)%coord(icord) = values1(n)
         enddo
         if (i==control_lines(NumberEntities(5))%icont(1)) then
            control_points(i)%dist = zero
            else
               control_points(i)%dist = control_points(i-1)%dist + vp
         endif
         values1 = values1 + values3
         if (nout>0) then
            write (nout,"(1x,a,i5,1pe12.4,3(3x,a,e12.4))") "Point ",i,         &
               control_points(i)%dist,(xyzlabel(icoordp(n,ncord-1))//" = ",    &
               control_points(i)%coord(icoordp(n,ncord-1)),n=1,ncord)
         endif
      enddo
      write (nout,"(1x,a)") " "
   endif
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL LINES DATA",ninp,nout))   &
      return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputControlLines

