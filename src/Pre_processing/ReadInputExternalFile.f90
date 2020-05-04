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
! Program unit: ReadInputExternalFile                       
! Description:                        
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine ReadInputExternalFile(NumberEntities,ainp,comment,nrighe,ier,       &
                                 OnlyTriangle,ninp,ulog,ninp2)
#elif defined SPACE_2D
subroutine ReadInputExternalFile(NumberEntities,ainp,comment,nrighe,ier,       &
                                 ninp,ulog,ninp2)
#endif
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
integer(4) :: nrighe,ier,ninp,ulog,ninp2
integer(4),dimension(20) :: NumberEntities
character(1) :: comment
character(len=lencard) :: ainp
#ifdef SPACE_3D
logical :: OnlyTriangle
#endif
integer(4) :: ioerr
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
if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"GEOMETRY FILE",ninp,ulog)) return
#ifdef SPACE_3D
OnlyTriangle = .true.
#endif
do while (trim(lcase(ainp))/="##### end geometry file #####")
   open(ninp2,file=trim(ainp),form="formatted",status="old",iostat=ioerr)
   if (ulog>0) then
      if (ioerr==0) then
         write(ulog,"(1x,3a)") "Geometry File: ",trim(ainp)
         else
            write(ulog,"(1x,3a)") "Geometry File: ",trim(ainp)," not found!"
            return
      endif
   endif
! To read the first line of the file
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp2)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"GEOMETRY FILE",ninp2,ulog)) return
   SECTION_LOOP: do while (ioerr==0)
      select case (trim(lcase(trim(ainp))))
         case ("##### vertices #####")
            call ReadInputVertices (NumberEntities,Vertice,ainp,comment,       &
                                    nrighe,ier,.false.,ninp2,ulog)
#ifdef SPACE_2D
         case ("##### lines #####")
            call ReadInputLines(NumberEntities,BoundaryVertex,Tratto,ainp,     &
                                comment,nrighe,ier,ninp2,ulog)
#elif defined SPACE_3D
         case ("##### faces #####")
            call ReadInputFaces(NumberEntities,ainp,comment,nrighe,ier,.false.,&
                                ninp2,ulog)
#endif
         case default
      endselect
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp2)
! In case of EOF, then it exits, otherwise it checks the error 
      if (ioerr==-1) cycle SECTION_LOOP
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"GEOMETRY FILE",ninp,ulog))     &
         return
   enddo SECTION_LOOP
   close (ninp2)
   if (ulog>0) then
      write(ulog,"(1x,3a)") "End Reading Geometry File"
   endif
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"GEOMETRY FILE",ninp,ulog)) return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputExternalFile
