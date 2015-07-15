!cfile inidt2.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : inidt2
!
! Last updating : November 23, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Agate             2011           varie
! 04  Amicarelli        23/11/2011     multiple inlet
!
!************************************************************************************
! Module purpose : Module to act on pl(0) where are pre-fixed x,y,z and define the
!                  particle variables
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
  subroutine inidt2 

!Versione AdM del 19/11/05 - Da usarsi INSIEME a RUNDT2

!Calcola il passo dt iniziale secondo la procedura euristica proposta da AdM
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Local Parameters ..
!double precision,parameter :: k1 = 2.33333333333333d0  ! 7/3
!double precision,parameter :: Skmaxq = 0.55d0 * 0.55d0
!double precision,parameter :: Ckmax = 1.33333333333333d0  ! 4/3
!
!.. Local Scalars ..
!AA401
!  integer(4)       :: j
!  double precision :: celmedq, viscmed, viscequi, dtmed, dtmin, dtdiff1, diffmax, Denom, TermB, TermC  !, hh, viscterm
  integer(4)       :: j,npi,mate,ii
  double precision :: dtmin,dt_CFL,dt_dif,dt_vis,diffmax,U
!
!.. Executable Statements ..
!
!.. initializations
! 
  dtmin   = 1.0d+30                                                        
  diffmax = zero
!
!.. evaluates the time step minimum value depending on the different celerities in the different materials
!
!  do j = 1,nmedium
!!
!    celmedq  = Med(j)%celerita * Med(j)%celerita
!    viscmed  = max ( Med(j)%visc, Med(j)%numx )
!!    viscequi = 0.75d0 * Med(j)%alfaMon * Med(j)%celerita * Domain%h + 2.3333333d0 * viscmed
!!    viscterm = 3.287311d0 * viscequi / celmedq
!!    dtmed    = Dsqrt( viscterm * viscterm + 3.698224d0 * Domain%h * Domain%h / celmedq ) - viscterm
!    viscequi = Med(j)%alfaMon * Med(j)%celerita * Domain%h + k1 * viscmed
!    Denom    = Skmaxq * celmedq
!    termB    = Ckmax * viscequi / (Denom + Denom)
!    termC    = squareh / Denom
!    dtmed    = Dsqrt(termB * termB + termC) - termB
!!
!    diffmax  = max ( diffmax,Med(j)%codif)
!    dtdiff1  = half * squareh / (diffmax+0.000000001d0)
!    dtmin    = min ( dtmed, dtmin , dtdiff1)
!!
!  end do
!!
!!.. evaluates the time step current value applying the time coefficient to the found value
!!
!!  dt = Domain%cote * dtmin
!  dt = Domain%CFL * dtmin
!
!AA401 loop to compute the time step, according to 3 conditions:
!1) the CFL condition: dt_CFL=min(2h/(c+U))
!2) viscous stability condition dt_vis=min(rho*h**2/(0.5*mu)) 
!3) interface diffusion condition dt_diff=(h**2/2*teta)
  do ii = 1,indarrayFlu
    npi = Array_Flu(ii)

    mate = pg(npi)%imed
    U = sqrt(pg(npi)%vel(1)**2+pg(npi)%vel(2)**2+pg(npi)%vel(3)**2)
    dt_CFL = 2.*Domain%h/(Med(mate)%celerita+U)
!    dt_vis = pg(npi)%dens*Domain%h**2/(0.5*pg(npi)%mu)
    dt_vis = pg(npi)%dens*squareh/(half*pg(npi)%mu)
    dtmin  = min(dtmin,dt_CFL,dt_vis)
  end do
  do j = 1,nmedium
    diffmax  = max(diffmax,Med(j)%codif)
  end do
  dt_dif  = half * squareh / (diffmax+0.000000001d0)
  dtmin   = min (dtmin,dt_dif)
!AA401 end
!.. evaluates the time step current value applying the time coefficient to the found value
!
!AA405 
!initial dt for a jet
  if (indarrayFlu == 0) dtmin = 2.*Domain%h/(Med(1)%celerita)
!
!AA401
! CFL is used as a constant for every condition (the CFL, the viscous and the diffusive one)
  dt = Domain%CFL * dtmin
!
!.. initializes the medium time value for step to the initial value
!                             
  dt_average = dt
!
  return
  end subroutine inidt2
!---split

