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
! Program unit: ReadInputGeneralPhysical                          
! Description:                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine ReadInputGeneralPhysical(NumberEntities,ainp,comment,nrighe,ier,ninp,nout)
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
character(1) :: comment
character(80) :: ainp
integer(4) :: n,icord,ioerr
double precision :: prif
double precision,dimension(3) :: values1
character(80),external :: lcase
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
! In case of restart, input data are not read
if (restart) then
   do while (TRIM(lcase(ainp))/="##### end general physical properties #####")
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                                &
         "GENERAL PHYSICAL PROPERTIES DATA",ninp,nout)) return
   enddo
   return
endif
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"GENERAL PHYSICAL PROPERTIES DATA",   &
   ninp,nout)) return
do while (TRIM(lcase(ainp))/="##### end general physical properties #####")
   read (ainp,*,iostat=ioerr) values1(1:NumberEntities(1))
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"GRAVITAL ACCELERATION VECTOR",ninp&
      ,nout)) return
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   read (ainp,*,iostat=ioerr) prif
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"REFERENCE PRESSURE",ninp,nout))   &
      return
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"GENERAL PHYSICAL PROPERTIES DATA",&
      ninp,nout)) return
enddo
if (ncord>0) then
   Domain%grav(:) = zero            
   do n=1,NumberEntities(1)
      icord = icoordp(n,ncord-1)
      Domain%grav(icord) = values1(n)
      if (nout>0) write (nout,"(1x,a,a,1p,e12.4)") xyzlabel(icord),            &
         "gravity acceler. :",Domain%grav(icord)
   enddo
   Domain%prif = prif
   if (nout>0) write (nout,"(1x,a,1p,e12.4)") "P rif:            :",Domain%prif
   if (nout>0) write (nout,"(1x,a)") " "
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputGeneralPhysical

