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
! Program unit: ReadInputControlLines                    
! Description: Reading monitoring lines.                      
!-------------------------------------------------------------------------------
subroutine ReadInputControlLines(NumberEntities,Control_Points,Control_Lines,  &
                                 ainp,comment,nrighe,ier,ninp,ulog)
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
integer(4),dimension(20) :: NumberEntities
type (TyCtlPoint),dimension(npointst) :: Control_Points
type (TyCtlLine),dimension(NLines) :: Control_Lines
character(1) :: comment
character(len=lencard) :: ainp
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
interface
   subroutine ReadRiga(ninp,ainp,io_err,comment_sym,lines_treated)
      implicit none
      integer(4),intent(in) :: ninp
      character(*),intent(inout) :: ainp
      integer(4),intent(out) :: io_err
      character(1),intent(in),optional :: comment_sym
      integer(4),intent(inout),optional :: lines_treated
   end subroutine ReadRiga
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL LINES DATA",ninp,ulog))      &
   return
npts = npoints
do while (trim(lcase(ainp))/="##### end control lines #####")
   values1 = zero
   values2 = zero
   values3 = zero
   NumberEntities(5) = NumberEntities(5) + 1
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL LINE LABEL",ninp,ulog))   &
      return
   label(1:8) = ainp(1:8)
   write(txt,"(i5)") NumberEntities(5)
   call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
   read(ainp,*,iostat=ioerr) values1(1:ncord)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                                   &
      "CONTROL LINE"//txt//" - FIRST POINT",ninp,ulog)) return
   call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
   read(ainp,*,iostat=ioerr) values2(1:ncord)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                                   &
      "CONTROL LINE"//txt//" - SECOND POINT",ninp,ulog)) return
   call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
   read(ainp,*,iostat=ioerr) ndiv
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                                   &
      "CONTROL LINE"//txt//" - POINTS NUMBER",ninp,ulog)) return
   NumberEntities(6) = NumberEntities(6) + ndiv
   if (input_second_read.eqv..true.) then
      control_lines(NumberEntities(5))%label = label
      control_lines(NumberEntities(5))%icont(1) = npts + 1
      control_lines(NumberEntities(5))%icont(2) = npts + ndiv
      npts = npts + ndiv
      values3(:) = (values2(:) - values1(:)) / (ndiv - 1)
      vp = dsqrt(values3(1) * values3(1) + values3(2) * values3(2) +           &
           values3(3) * values3(3))
      if (ulog>0) then
         write(ulog,"(1x,a,i3,1x,a)") "Control line      ",NumberEntities(5),  &
            "("//control_lines(NumberEntities(5))%label//")"
         write(ulog,"(1x,a,i12)") "First Point:      ",                        &
            control_lines(NumberEntities(5))%icont(1)
         write(ulog,"(1x,a,i12)") "Last  Point:      ",                        &
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
         if (ulog>0) then
            write(ulog,"(1x,a,i5,1pe12.4,3(3x,a,e12.4))") "Point ",i,          &
               control_points(i)%dist,(xyzlabel(icoordp(n,ncord-1))//" = ",    &
               control_points(i)%coord(icoordp(n,ncord-1)),n=1,ncord)
         endif
      enddo
      write(ulog,"(1x,a)") " "
   endif
   call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL LINES DATA",ninp,ulog))   &
      return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputControlLines
