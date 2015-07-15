!cfile BoundaryVolumeIntegrals2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : BoundaryVolumeIntegrals2D
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
! Module purpose : Module to compute the boundary volume integrals IntWdV only
!
! Calling routine: ComputeBoundaryDataTab
!
! Called routines: 
!
!************************************************************************************
!
subroutine BoundaryVolumeIntegrals2D (icbs, LocXY, xpmin, xpmax, interlen, IntWdV, IntDpWdV)
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
integer(4),      parameter :: Nalfadiv = 10
double precision,parameter :: eps = 0.05d0
!
!.. Formal Arguments ..
integer(4),      intent(IN)    :: icbs
double precision,intent(IN),   dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES) :: LocXY
double precision,intent(IN)    :: xpmin
double precision,intent(IN)    :: xpmax
double precision,intent(IN)    :: interlen
double precision,intent(INOUT) :: IntWdV
double precision,intent(INOUT),dimension(1:2) :: IntDpWdV
!
!.. Local Scalars ..
!integer(4)       :: NumBSides, iside
integer(4)       :: k, ndiv
double precision :: xpi, ypi, yplimite, &             !ypiq, 
                    dalfarif, Intalfa, dalfa, &
                    alfaA, alfaB, alfa_k, csiPA, csiPB, sinalfa, cosalfa, rb, rob, mult
!
!.. Local Arrays ..
!
! External functions and subrotuines
double precision,external :: WIntegr, J2Wro2
!
!.. Executable Statements ..
!
  dalfarif = PIGRECO / Nalfadiv
!
  IntWdV = zero
  IntDpWdV = zero
!
  if (interlen <= zero) return
!
  yplimite = eps * Domain%h
  xpi = LocXY(1, icbs)
  ypi = LocXY(2, icbs)
!
  if (ypi < yplimite) ypi = yplimite
!
!  ypiq = ypi * ypi
  csiPA = xpmin - xpi
  csiPB = xpmax - xpi
  alfaA = Atan2(csiPA, ypi)
  alfaB = Atan2(csiPB, ypi)
  Intalfa = alfaB - alfaA
  ndiv = Int(intalfa / dalfarif + half)
  if ( ndiv < 1 ) ndiv = 1
  dalfa = intalfa / ndiv
  alfa_k = alfaA - half * dalfa
!
  do k = 1, ndiv
    alfa_k = alfa_k + dalfa
    sinalfa = sin(alfa_k)
    cosalfa = cos(alfa_k)
!    if (cosalfa == zero) then 
!      rb = zero
!    else
      rb = ypi / cosalfa
!    end if
    rob = rb / Domain%h
    mult = Domain%h * J2Wro2(rob) * dalfa
    IntDpWdV(1) = IntDpWdV(1) + sinalfa * mult
    IntDpWdV(2) = IntDpWdV(2) - cosalfa * mult
    IntWdV = IntWdV + WIntegr(rb, Domain%h) * dalfa
  end do
!
! controllo se la particella è molto vicina allo spigolo interno 
!
  if (ypi == yplimite) then
!    if ( ) then ! la proiezione della particella è interna alla faccia
      IntWdV = half
!    else
!      IntWdV = zero
  end if
!
return
end subroutine BoundaryVolumeIntegrals2D
!---split

