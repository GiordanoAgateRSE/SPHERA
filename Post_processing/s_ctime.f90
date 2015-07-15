!cfile s_ctime.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine s_ctime ( nout )

implicit none

!include "use_dflib.f90"
!use DFPORT
!PANE_WINNT      integer time,time_

integer(4)    :: nout
character(24) :: str
integer(4),    external :: time
character(24), external :: ctime

!PANE_WINNT      time = time_()
!PANE_WINNT      call ctime_(str,time)
 str = ctime(time())
 write (nout,'(a,a)') 's_ctime routine --> ',str

return
end subroutine s_ctime
!---split

