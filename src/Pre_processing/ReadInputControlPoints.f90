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
! Program unit: ReadInputControlPoints                     
! Description: Reading monitoring points.                      
!-------------------------------------------------------------------------------
subroutine ReadInputControlPoints(NumberEntities,Control_Points,ainp,comment,  &
                                  nrighe,ier,ninp,nout)
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
character(1) :: comment
character(100) :: ainp
integer(4) :: n,i,icord,ioerr
double precision,dimension(3) :: values1
logical,external :: ReadCheck
character(100),external :: lcase
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
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL POINTS DATA",ninp,nout))     &
   return
do while (TRIM(lcase(ainp))/="##### end control points #####" )
   NumberEntities(4) = NumberEntities(4) + 1
   read(ainp,*,iostat=ioerr) values1(1:ncord)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL POINT COORDINATES",ninp,  &
      nout)) return
   if (ncord>0) then
      Control_Points(NumberEntities(4))%coord(1:3) = zero
      Control_Points(NumberEntities(4))%dist = zero
      do n=1,ncord
         icord = icoordp(n,ncord-1)
         Control_Points(NumberEntities(4))%coord(icord) = values1(n)
      enddo
      if (nout>0) then
         i = NumberEntities(4)
         write (nout,"(1x,a,i5,3(3x,a,1p,e12.4))") "Control point ",i,         &
            (xyzlabel(icoordp(n,ncord-1))//" = ",                              &
            Control_Points(i)%coord(icoordp(n,ncord-1)),n=1,ncord)
       endif
    endif
    call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
    if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL POINTS DATA",ninp,nout)) &
       return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputControlPoints

