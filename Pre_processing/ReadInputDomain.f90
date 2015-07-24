!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-; SPHERA has been authored for RSE SpA by 
!    Andrea Amicarelli, Antonio Di Monaco, Sauro Manenti, Elia Bon, Daria Gatti, Giordano Agate, Stefano Falappi, 
!    Barbara Flamini, Roberto Guandalini, David Zuccalà).
! Main numerical developments of SPHERA: 
!    Amicarelli et al. (2015,CAF), Amicarelli et al. (2013,IJNME), Manenti et al. (2012,JHE), Di Monaco et al. (2011,EACFM). 
! Email contact: andrea.amicarelli@rse-web.it

! This file is part of SPHERA.
! SPHERA is free software: you can redistribute it and/or modify
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
! Program unit: ReadInputDomain                       
! Description:                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine ReadInputDomain(NumberEntities,ainp,comment,nrighe,ier,ninp,nout,   &
                           nscr)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nrighe,ier, ninp,nout,nscr
integer(4),dimension(20) :: NumberEntities
character(1) :: comment
character(80) :: ainp
integer(4) :: ioerr
double precision :: dd, trunc
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
   do while (TRIM(lcase(ainp))/="##### end domain #####")
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DOMAIN DATA",ninp,nout))       &
         return
   enddo
   return
endif
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DOMAIN DATA",ninp,nout)) return
do while (TRIM(lcase(ainp))/="##### end domain #####")
   token = lcase(GetToken(ainp,1,ioerr))
   read(token,*,iostat=ioerr) NumberEntities(1)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DOMAIN COORDINATES NUMBER",ninp,  &
      nout)) return 
   if ((ncord>0).and.(nout>0)) then
      write(nout,"(1x,a,i3,1x,a)") "Domain Dimension       : ",ncord,          &
         ncordlabel(ncord)
   endif
   token = lcase(GetToken(ainp,2,ioerr))
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DOMAIN TYPE",ninp,nout)) return 
   select case (token(1:4))
      case ("bsph","semi") 
         if ((.not.DBSPH_Module).and.(token(1:4)=="bsph")) then
            write(nscr,"(1x,a)") " "
            write(nout,"(1x,a)") " "
            write(nscr,"(1x,a)")                                               &
               " >>WARNING! - The DBSPH module is not available."
            write(nout,"(1x,a)")                                               &
               " >>WARNING! - The DBSPH module is not available."
            ier = 40
            return
         endif
         Domain%tipo = token(1:4)
         if (ncord>0.and.nout>0) then
            write(nout,"(1x,a,1x,a)"  ) "Domain Type            : ",           &
               trim(token)
         endif      
      case default
         ier = 3
         return
   end select
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (ioerr==0) read(ainp,*,iostat=ioerr) dd,trunc
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DD & TRUNC",ninp,nout)) return
   Domain%dd = dd
   Domain%trunc = trunc
   token = lcase(GetToken(ainp,3,ioerr))
   if (token(1:1)=='r') then
      Domain%RandomPos = 'r'
      else
         Domain%RandomPos = 'n'
   endif
   if ((ncord>0).and.(nout>0)) then
      write(nout,"(1x,a,1pe12.4)") "Dd                     : ",dd
      write(nout,"(1x,a,1pe12.4)") "Trunc                  : ",trunc
      write(nout,"(1x,a,1x,a)") "Random Initial Position: ",Domain%RandomPos
      write(nout,"(1x,a)")  " "
   endif
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DOMAIN DATA",ninp,nout)) return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputDomain

