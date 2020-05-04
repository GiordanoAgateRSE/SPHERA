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
! Program unit: ReadInputVertices                              
! Description:                        
!-------------------------------------------------------------------------------
subroutine ReadInputVertices(NumberEntities,Vertice,ainp,comment,nrighe,ier,   &
                             prtopt,ninp,ulog)
!------------------------
! Modules
!------------------------
use Static_allocation_module                              
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
logical(4) :: prtopt
integer(4) :: nrighe,ier,ninp,ulog
integer(4),dimension(20) :: NumberEntities
double precision,dimension(1:SPACEDIM,NumVertici) :: Vertice
character(1) :: comment
character(len=lencard) :: ainp
integer(4) :: n,i,icord,ioerr
double precision,dimension(3) :: values1
character(8) :: label
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
! In case of restart, input data are not read
if (restart) then
   do while (trim(lcase(ainp))/="##### end vertices #####")
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
      if (.not.ReadCheck (ioerr,ier,nrighe,ainp,"VERTICES DATA",ninp,ulog))    &
         return
   enddo
  return
endif
call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
if (.not.ReadCheck (ioerr,ier,nrighe,ainp,"VERTICES DATA",ninp,ulog)) return
if ((input_second_read.eqv..true.).and.(ulog>0).and.(prtopt)) then
   write(ulog,"(1x,a)") "List of vertices:"
endif
do while (trim(lcase(ainp))/="##### end vertices #####")
   select case (TRIM(Domain%tipo))
      case ("semi","bsph") 
         read (ainp,*,iostat=ioerr) i, values1(1:ncord)
         write(label,"(i8)") i
         if (.not.ReadCheck (ioerr,ier,nrighe,ainp,"VERTEX n."//label,ninp,    &
            ulog)) return
         NumberEntities(7) = max(i,NumberEntities(7))
         if (input_second_read.eqv..true.) then
            do n=1,ncord
               icord = icoordp(n,ncord-1)
               if (NumberEntities(7)==1) then
                  Domain%coord(icord,1) = values1(n)
                  Domain%coord(icord,2) = values1(n)
               endif
               Vertice(icord,i) = values1(n)
               Domain%coord(icord,1) = min(values1(n),Domain%coord(icord,1))
               Domain%coord(icord,2) = max(values1(n),Domain%coord(icord,2))
            enddo
         endif
      case default
         if (ulog>0) then
            write(ulog,*) "Unknown Domain Type: ",Domain%tipo
         endif
         ier = 2
         return
   endselect
   if ((input_second_read.eqv..true.).and.(ulog>0).and.(prtopt)) then
      write(ulog,"(i6,1p,3(2x,a,e12.4))") i,(xyzlabel(icoordp(n,ncord-1)),     &
         Vertice(icoordp(n,ncord-1),i),n=1,ncord)
   endif
   call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
   if (.not.ReadCheck (ioerr,ier,nrighe,ainp,"VERTICES DATA",ninp,ulog)) return
enddo
if ((input_second_read.eqv..true.).and.(ulog>0)) then
   do n=1,ncord
      icord = icoordp(n,ncord-1)
      write(ulog,"(1x,a,a,1p,e12.4)") xyzlabel(icord)," coordinate min. ",     &
         Domain%coord(icord,1)
      write(ulog,"(1x,a,a,1p,e12.4)") xyzlabel(icord)," coordinate max. ",     &
         Domain%coord(icord,2)
   enddo
   write(ulog,"(1x,a)") " "
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputVertices
