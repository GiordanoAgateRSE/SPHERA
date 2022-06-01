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
! Program unit: ReadInputGeneralPhysical                          
! Description:                        
!-------------------------------------------------------------------------------
subroutine ReadInputGeneralPhysical(ainp,comment,nrighe,ier,ninp,ulog)
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
integer(4) :: n,icord,ioerr
double precision :: prif
double precision,dimension(3) :: values1
character(100),external :: lcase
logical,external :: ReadCheck
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
if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"GENERAL PHYSICAL PROPERTIES DATA",   &
   ninp,ulog)) return
do while (trim(lcase(ainp))/="##### end general physical properties #####")
   read(ainp,*,iostat=ioerr) values1(1:ncord)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"GRAVITAL ACCELERATION VECTOR",ninp&
      ,ulog)) return
   call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
   read(ainp,*,iostat=ioerr) prif
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"REFERENCE PRESSURE",ninp,ulog))   &
      return
   call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"GENERAL PHYSICAL PROPERTIES DATA",&
      ninp,ulog)) return
enddo
if (input_second_read.eqv..true.) then
   Domain%grav(:) = zero            
   do n=1,ncord
      icord = icoordp(n,ncord-1)
      Domain%grav(icord) = values1(n)
      if (ulog>0) write(ulog,"(1x,a,a,1p,e12.4)") xyzlabel(icord),             &
         "gravity acceler. :",Domain%grav(icord)
   enddo
   Domain%prif = prif
   if (ulog>0) write(ulog,"(1x,a,1p,e12.4)") "P rif:            :",Domain%prif
   if (ulog>0) write(ulog,"(1x,a)") " "
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputGeneralPhysical
