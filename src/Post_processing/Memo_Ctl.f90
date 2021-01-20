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
! Program unit: Memo_Ctl               
! Description: Post-processing for monitoring lines and points.       
!-------------------------------------------------------------------------------
subroutine Memo_Ctl
!------------------------
! Modules
!------------------------
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: i,j
character(255) :: nomefilectl
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
if (npoints>0) then
   write(nomefilectl,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_',       &
      on_going_time_step,".cpt"
   open(ncpt,file=nomefilectl,status="unknown",form="formatted")
   write(ncpt,*) "Control points "
   write(ncpt,'(1x,2(a,10x),3(a,8x),3(a,5x),a,7x,a)') " Time","Iter",         &
      "X Coord","Y Coord","Z Coord","X Velocity","Y Velocity","Z Velocity",    &
      " Pressure","Density "
  flush(ncpt)
  do i=1,npoints
     if (control_points(i)%cella==0) then
        write(ncpt,'(a,i10,a,3g14.7)') "control point ",i,                    &
           " is outside. Coord=",Control_Points(i)%coord(:)
        else
           write(ncpt,'(g14.7,i14,8(1x,g14.7))') simulation_time,             &
              on_going_time_step,Control_Points(i)%coord(:),                   &
              Control_Points(i)%vel(:),Control_Points(i)%pres,                 &
              Control_Points(i)%dens
     endif
  enddo
  close(ncpt)
endif
! Printing monitoring line data
if (nlines>0) then
   write(nomefilectl,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_',       &
      on_going_time_step,".cln"
   open(ncpt,file=nomefilectl,status="unknown",form="formatted")
   write(ncpt,*) "Control lines "
   write(ncpt,'(5(11x,a),3(9x,a),(10x,a),(5x,a))') "t(s)","step","x(m)",      &
      "y(m)","z(m)","u(m/s)","v(m/s)","w(m/s)","p(Pa)","rho(kg/m3)"
  flush(ncpt)
  do i=1,nlines
     write(ncpt,*) "line #", i,"    Label ",Control_Lines(i)%label
     do j=Control_Lines(i)%icont(1),Control_Lines(i)%icont(2)
        if (control_points(j)%cella==0) then
           write(ncpt,'(a,i10,a,3(g14.7))') "control point ",j,               &
              " is outside. Coord=",Control_Points(j)%coord(:)
           else
              write(ncpt,'(g14.7,i14,8(1x,g14.7))') simulation_time,          &
                 on_going_time_step,Control_Points(j)%coord(:),                &
                 Control_Points(j)%vel(:),Control_Points(j)%pres,              &
                 Control_Points(j)%dens
        endif
     enddo
  enddo
  close(ncpt)
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine Memo_Ctl
