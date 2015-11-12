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
! Program unit: ReadInputFaces                         
! Description:                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine ReadInputFaces(NumberEntities,ainp,comment,nrighe,ier,prtopt,ninp,nout)
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
integer(4) :: nrighe,ier, ninp,nout
logical(4) :: prtopt
integer(4),dimension(20) :: NumberEntities
character(1) :: comment
character(80) :: ainp
integer(4) :: n,i,ioerr,stretch
integer(4),dimension(4) :: ivalues
character(8) :: label
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
   do while (TRIM(lcase(ainp))/="##### end faces #####")
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FACES DATA",ninp,nout)) return
   enddo
   return
endif
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FACES DATA",ninp,nout)) return
do while (TRIM(lcase(ainp))/="##### end faces #####")
   select case (TRIM(Domain%tipo))
      case ("semi","bsph") 
         ivalues = 0
! ivalues(4) is the 4-th vertex (null in case of trinagles)
         read(ainp,*,iostat=ioerr) i,ivalues,stretch  
         write(label,"(i8)") i
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FACE n."//label,ninp,nout)) &
            return
         NumberEntities(11) = max(i,NumberEntities(11))
! To count quadrilaterals to be splitted in triangles (if required)
         if (ivalues(4)>0) NumberEntities(18) = NumberEntities(18) + 1
         if (ncord>0) then
            if (BoundaryFace(i)%Node(1)%name==0) then 
               BoundaryFace(i)%Node(1:MAXFACENODES)%name =                     &
                  ivalues(1:MAXFACENODES)
               BoundaryFace(i)%stretch = stretch
               else
                  if (nout>0) then
                     write(nout,*) "Face definition: ",trim(ainp)
                     write(nout,*) "Face already defined: ",i,                 &
                        BoundaryFace(i)%Node(1:MAXFACENODES)%name
                  endif
                  ier = 104
                  return
            endif
         endif
      case default
         if (nout>0) then
            write(nout,*) "Unknown Domain Type: ",Domain%tipo
         endif
         ier = 2
         return
   end select
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FACES DATA",ninp,nout)) return
enddo
if ((ncord>0).AND.(nout>0).AND.(prtopt)) then
   write(nout,*)
   write(nout,"(1x,a)") "List of faces:"
   do n=1,NumberEntities(11)
      write(nout,"(i10,' - ',4i10,' - ',i8)") n,BoundaryFace(n)%Node(1)%name,  &
         BoundaryFace(n)%Node(2)%name,BoundaryFace(n)%Node(3)%name,            &
         BoundaryFace(n)%Node(4)%name,BoundaryFace(n)%stretch
   enddo
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputFaces

