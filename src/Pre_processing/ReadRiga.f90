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
! Program unit: ReadRiga                               
! Description:                         
!----------------------------------------------------------------------------------------------------------------------------------

subroutine ReadRiga(ainp,comment,nrighe,ier,ninp)
!------------------------
! Modules
!------------------------ 
implicit none
integer(4) :: ier,ninp,nrighe
character(1) :: comment
character(*) :: ainp
integer(4) :: ioerr,n,l
!------------------------
! Declarations
!------------------------
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
ioerr = 0
READ_LOOP: do while (ioerr==0)
   read(ninp,"(a)",iostat=ioerr) ainp
   nrighe = nrighe + 1
   if ((ioerr==0).AND.(trim(ainp)/="")) then
! Replacing tabs with blank spaces 
      if (ainp(1:1)/=comment ) then
         do n=1,len(trim(ainp))
            if (iachar(ainp(n:n))==9) ainp(n:n) = " "
         enddo
         l = index(ainp,comment)
         if (l>0) then
            do n = l,len(trim(ainp))
               ainp(n:n)=" "
            enddo
         endif
         exit READ_LOOP 
      endif
   endif
enddo  READ_LOOP
ier = ioerr
!------------------------
! Statements
!------------------------
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadRiga

