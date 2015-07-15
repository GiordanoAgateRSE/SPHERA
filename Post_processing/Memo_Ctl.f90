!cfile Memo_Ctl.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : Memo_Ctl
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!
!************************************************************************************
! Module purpose : Module to write results in the control lines and control points
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
!

subroutine Memo_Ctl
!
!.. assign modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Scalars ..
!
integer(4)     :: i,j
character(255) :: nomefilectl
!
!.. Executable Statements ..
!
!.. stampa control points
!
if (Npoints > 0) then
  write(nomefilectl,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_',it_corrente,".cpt"
  open ( ncpt, file=nomefilectl, status="unknown", form="formatted" )
  write (ncpt,*) "Control points "
  write (ncpt,'(1x,2(a,10x),3(a,8x),3(a,5x),a,7x,a)') &
  " Time","Iter","X Coord","Y Coord","Z Coord","X Velocity","Y Velocity","Z Velocity"," Pressure","Density "
  flush(ncpt)
  do i = 1,Npoints
    if (control_points(i)%cella == 0) then
      write (ncpt,'(a,i10,a,3g14.7)') "control point ",i," is outside. Coord=",Control_Points(i)%coord(:)
    else
      write (ncpt,'(g14.7,i14,8(1x,g14.7))') tempo, it_corrente &
                             ,Control_Points(i)%coord(:) &
                             ,Control_Points(i)%vel(:)   &
                             ,Control_Points(i)%pres     &
                             ,Control_Points(i)%dens
    end if
  end do
  close (ncpt)
end if
!
!.. stampa control lines
!
if (Nlines > 0) then
  write(nomefilectl,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_',it_corrente,".cln"
  open ( ncpt, file=nomefilectl, status="unknown", form="formatted" )
  write (ncpt,*) "Control lines "
  write (ncpt,'(1x,2(a,10x),3(a,8x),3(a,5x),a,7x,a)') &
  " Time","Iter","X Coord","Y Coord","Z Coord","X Velocity","Y Velocity","Z Velocity"," Pressure","Density "
  flush(ncpt)
  do i = 1,Nlines
    write (ncpt,*) "line #", i,"    Label ",Control_Lines(i)%label
    do j = Control_Lines(i)%icont(1),Control_Lines(i)%icont(2)
      if (control_points(j)%cella == 0) then
        write (ncpt,'(a,i10,a,g14.7)') "control point ",j," is outside. Coord=",Control_Points(j)%coord(:)
      else
        write (ncpt,'(g14.7,i14,8(1x,g14.7))') tempo, it_corrente &
                               ,Control_Points(j)%coord(:) &
                               ,Control_Points(j)%vel(:)   &
                               ,Control_Points(j)%pres     &
                               ,Control_Points(j)%dens
      end if
    end do
  end do
  close (ncpt)
end if
!AA601 rm part
!.. stampa control sections
!
if (Nsections > 0) then
  write(nomefilectl,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_',it_corrente,".csc"
  open ( ncpt, file=nomefilectl, status="unknown", form="formatted" )
  write (ncpt,*) "Control sections "
  write (ncpt,'(1x,2(a,10x),3(a,8x),3(a,5x),a,7x,a)') &
  " Time","Iter","X Coord","Y Coord","Z Coord","X Velocity","Y Velocity","Z Velocity"," Pressure","Density "
  flush(ncpt)
  do i = 1,Nsections
    write (ncpt,*) "section #", i,"    Label ",Control_sections(i)%label,"    Type ",Control_sections(i)%Tipo
    do j = Control_sections(i)%icont(1),Control_sections(i)%icont(2)
      if (control_points(j)%cella == 0) then
        write (ncpt,'(a,i10,a,g14.7)') "control point ",j," is outside. Coord=",Control_Points(j)%coord(:)
      else
        write (ncpt,'(g14.7,i14,8(1x,g14.7))') tempo, it_corrente &
                               ,Control_Points(j)%coord(:) &
                               ,Control_Points(j)%vel(:)   &
                               ,Control_Points(j)%pres     &
                               ,Control_Points(j)%dens
      end if
    end do
  end do
  close (ncpt)
end if

return

end subroutine Memo_Ctl
!---split

