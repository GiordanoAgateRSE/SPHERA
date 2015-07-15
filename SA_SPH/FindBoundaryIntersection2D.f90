!cfile FindBoundaryIntersection2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : FindBoundaryIntersection2D
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
! Module purpose : Module for Find the intersection segment between the circle of
!                  influence of particle i, whose local coordinates are
!                  xpi=LocXY(1, icbs) and ypi=LocXY(2, icbs), and the straight
!                  boundary side iside=Cloboside(icbs), which lies on the local
!                  x-axis and extends from x=0 to bsidelen = BoundarySide(iside)%Length 
!
! Calling routine: ComputeBoundaryDataTab
!
! Called routines: 
!
!************************************************************************************
!
subroutine FindBoundaryIntersection2D (icbs, Cloboside, LocXY, BoundarySide, &
                                       xpmin, xpmax, interlen)

!Finds the intersection segment between the circle of influence of particle i,
!whose local coordinates are xpi=LocXY(1, icbs) and ypi=LocXY(2, icbs), and the
!straight boundary side iside=Cloboside(icbs), which lies on the local x-axis and extends
!from x=0 to bsidelen = BoundarySide(iside)%Length 
!Returns:
!xpmin =        Minimum abscissa of intersected segment
!xpmax =        Maximum abscissa of intersected segment
!interlen =        Length of intersected segment
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
integer(4)                                                       :: icbs
integer(4),           dimension(1:MAXCLOSEBOUNDSIDES)            :: Cloboside
double precision,     dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES) :: LocXY
type (TyBoundarySide),dimension(1:NumBSides)                     :: BoundarySide
double precision                                                 :: xpmin
double precision                                                 :: xpmax
double precision                                                 :: interlen
!
!.. Local Scalars ..
integer(4)       :: iside
double precision :: xpi, ypi, ypiq, bsidelen, halfchord, minlen, eps, XA, XB
!
!.. Local Arrays ..
!integer(4),dimension(1:PLANEDIM) :: acix
!
!.. Executable Statements ..
!
!  acix(1)=1        !active coordinate indexes
!  acix(2)=3
!
  eps = 0.01d0
  minlen = eps*Domain%h
  interlen = zero
!
  xpi = LocXY(1, icbs)
  ypi = LocXY(2, icbs)
!
  if (ypi < zero) then    !if the particle is out of the boundary operates as if the particle
    ypi = zero            !where on the boundary
  end if    
  ypiq = ypi * ypi
  halfchord = Dsqrt(doublesquareh - ypiq)
  XA = xpi - halfchord
  XB = xpi + halfchord
  iside = Cloboside(icbs)
  bsidelen = BoundarySide(iside)%Length
  xpmin = zero
  if (XA > xpmin) then
    xpmin = XA
  end if
  xpmax = bsidelen
  if (XB < xpmax) then
    xpmax = XB
  end if
  interlen = xpmax - xpmin
  if (interlen <= minlen) then
    interlen = zero               !Intersection is too small
  end if
!
return
end subroutine FindBoundaryIntersection2D
!---split

