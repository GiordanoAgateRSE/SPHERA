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
! Program unit: ReadInputExternalFile                        
! Description:                        
!-------------------------------------------------------------------------------
subroutine ReadInputExternalFile(NumberEntities,ainp,comment,nrighe,ier,       &
                                 OnlyTriangle,ninp,nout,ninp2)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nrighe,ier,ninp,nout,ninp2
integer(4),dimension(20) :: NumberEntities
character(1) :: comment
character(100) :: ainp
logical :: OnlyTriangle
integer(4) :: ioerr
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
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"GEOMETRY FILE",ninp,nout)) return
OnlyTriangle = .TRUE.
do while (TRIM(lcase(ainp))/="##### end geometry file #####")
   open(ninp2,file=trim(ainp),form="formatted",status="old",iostat=ioerr)
   if (nout>0) then
      if (ioerr==0) then
         write (nout,"(1x,3a)") "Geometry File: ",trim(ainp)
         else
            write (nout,"(1x,3a)") "Geometry File: ",trim(ainp)," not found!"
            return
      endif
   endif
! To read the first line of the file 
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp2)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"GEOMETRY FILE",ninp2,nout)) return
   SECTION_LOOP: do while (ioerr==0)
      select case (TRIM(lcase(trim(ainp))))
         case ("##### vertices #####")
            call ReadInputVertices (NumberEntities,Vertice,ainp,comment,       &
                                    nrighe,ier,.FALSE.,ninp2,nout)
         case ("##### lines #####")
            call ReadInputLines(NumberEntities,BoundaryVertex,Tratto,ainp,     &
                                comment,nrighe,ier,ninp2,nout)
         case ("##### faces #####")
            call ReadInputFaces(NumberEntities,ainp,comment,nrighe,ier,.FALSE.,&
                                ninp2,nout)
         case default
      endselect
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp2)
! In case of EOF, then it exits, otherwise it checks the error 
      if (ioerr==-1) cycle SECTION_LOOP
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"GEOMETRY FILE",ninp,nout))     &
         return
   enddo SECTION_LOOP
   close (ninp2)
   if (nout>0) then
      write (nout,"(1x,3a)") "End Reading Geometry File"
   endif
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"GEOMETRY FILE",ninp,nout)) return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputExternalFile

