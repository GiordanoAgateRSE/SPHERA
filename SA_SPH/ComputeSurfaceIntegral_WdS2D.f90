!cfile ComputeSurfaceIntegral_WdS2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ComputeSurfaceIntegral_WdS2D
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
! Module purpose : Module for Computing the surface integral of kernel W along the
!                  segments intercepted by the influence circle (radius=2h) of the
!                  particle i, whose local coordinates are xpi=LocXY(1, icbs) and
!                  ypi=LocXY(2, icbs), on the adjacent boundary side icbs.
 
!
! Calling routine: ComputeBoundaryDataTab
!
! Called routines: 
!
!************************************************************************************
!
subroutine ComputeSurfaceIntegral_WdS2D (icbs, LocXY, xpmin, interlen, SIntWds) !Cloboside, xpmax, 

!Computes the surface integral of kernel W along the segments intercepted
!by the influence circle (radius=2h) of the particle i, whose local coordinates
!are xpi=LocXY(1, icbs) and ypi=LocXY(2, icbs), on the adjacent boundary side icbs.
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
integer(4), parameter :: ndivrif = 2
!
!.. Formal Arguments ..
integer(4),      intent(IN)    :: icbs
!integer(4),      intent(IN),dimension(1:MAXCLOSEBOUNDSIDES)            :: Cloboside
double precision,intent(IN),dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES) :: LocXY
double precision,intent(IN)    :: xpmin
!double precision,intent(IN)    :: xpmax
double precision,intent(IN)    :: interlen
double precision,intent(INOUT) :: SIntWds
!
!.. Local Scalars ..
integer(4)       :: ndiv, ipt    !iside, 
double precision :: xpi, ypi, ypiq, dsrif, deltas, xpip, dxpip, ris, risq, Wds
!
!.. Local Arrays ..
!integer(4),dimension(1:PLANEDIM)  :: acix
!
! External functions and subrotuines
double precision,    external    :: w
!
!.. Executable Statements ..
!
!  acix(1) = 1        !active coordinate indexes
!  acix(2) = 3

  SIntWds = zero
!  iside = Cloboside(icbs)
  dsrif = Domain%h / ndivrif

  if (interlen > zero) then
    xpi = LocXY(1, icbs)
    ypi = LocXY(2, icbs)
    if (ypi < zero) then      !if the particle is out of the boundary operates as if the particle
        ypi = zero            !where on the boundary
    end if    
    ypiq = ypi * ypi
    ndiv = Int(interlen / dsrif + half)
    if (ndiv < ndivrif) then
        ndiv = ndivrif
    end if
    deltas = interlen / ndiv
    xpip = xpmin - half * deltas
    do ipt = 1, ndiv
      xpip = xpip + deltas
      dxpip = xpip - xpi
      risq = dxpip * dxpip + ypiq
      ris = Dsqrt(risq)
      Wds = w(ris, Domain%h, Domain%coefke) * deltas
      SIntWds = SIntWds + Wds
    end do
  end if
!
return
end subroutine ComputeSurfaceIntegral_WdS2D
!---split

