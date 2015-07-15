!cfile writime2.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine writime2 (ti,tf,nout)
!
use language_writime2
!
implicit none
!
double precision, dimension(2) :: ti,tf
integer(4)                     :: nout
!
! local arrays and scalars
logical, save          :: first=.true.
double precision       :: telat, tcput, tcpup, telap
double precision, save :: tcpus, telas
!
if ( first ) then
   tcpus = ti(1)
   telas = ti(2) !- ti(1)  ! toglie tcpus per aggiustare tempo elapsed (granularita' in secondi)
   tcpup = ti(1)
   telap = ti(2)
end if
first = .false.
!
tcput = tf(1)
telat = tf(2) - telas
!
write(nout,1001) " "
write(nout,1001) cpulbl,tcput,totlbl,tf(1)-ti(1),przlbl

write(nout,1001) elalbl,max(1.0d0,telat),totlbl,max(1.0d0,tf(2)-ti(2)),przlbl
!
tcpup = tf(1)
telap = tf(2)
!
return
 1001 format(1x,a,f20.4,1x,a,f20.4,1x,a)
end subroutine writime2
!---split

