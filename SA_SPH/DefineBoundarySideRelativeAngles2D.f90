!cfile DefineBoundarySideRelativeAngles2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : DefineBoundarySideRelativeAngles2D
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
! Module purpose : Module for calculation of previous adjacent side and 
!                  relative angle for each boundary side
!
! Calling routine: DefineBoundarySideGeometry2D
!
! Called routines: 
!
!************************************************************************************
!
  subroutine DefineBoundarySideRelativeAngles2D
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Local Scalars ..
  integer(4)       :: nti, ntj, isi, jsi, ksi
  double precision :: sinangle, cosangle
!
!.. local arrays ..
  double precision, dimension(1:SPACEDIM, 1:SPACEDIM) ::    Tmx,  PTmx
!
!.. Executable Statements ..
!
!.. loops on all the domain boundary sides
!
  do isi = 1, NumBSides
!
    nti = BoundarySide(isi)%stretch
!
!.. skips the perimeter and pool types
!
    if (Tratto(nti)%tipo /= "peri" .AND. Tratto(nti)%tipo /= "pool" ) then
!
!.. loops on all the other sides
!
      ksi=0
      do jsi = 1, NumBSides
!
       ntj = BoundarySide(jsi)%stretch
       if (Tratto(ntj)%tipo == "peri" .or. Tratto(ntj)%tipo == "pool" ) cycle
!
!.. checks if the sides are adjacents; !ksi is the previous adjacent side
!
        if (BoundarySide(jsi)%Vertex(2) == BoundarySide(isi)%Vertex(1)) then
          ksi = jsi                    
          exit
        end if
!
      end do
!
      BoundarySide(isi)%angle = zero
!      BoundarySide(isi)%previous_side = 0
!
!.. an adjacent side has been found
!
      if (ksi > 0) then
!
!.. evaluates the angle between the two sides in radiants
!
        ntj = BoundarySide(ksi)%stretch
!
! errore corretto 25set08       if (Tratto(nti)%tipo /= "peri" .AND. Tratto(nti)%tipo /= "pool" ) then
        if (Tratto(ntj)%tipo /= "peri" .AND. Tratto(ntj)%tipo /= "pool" ) then
!            
          Tmx  = BoundarySide(isi)%T
! errore corretto 25set08         PTmx = BoundarySide(jsi)%T
          PTmx = BoundarySide(ksi)%T
!
          sinangle = PTmx(1, 1) * Tmx(3, 1) - PTmx(3, 1) * Tmx(1, 1)
          cosangle = PTmx(1, 1) * Tmx(1, 1) + PTmx(3, 1) * Tmx(3, 1)
          BoundarySide(isi)%angle = Atan2(sinangle, cosangle)
!
        end if
! errore corretto 06ott08
!        BoundarySide(isi)%previous_side = ksi
!
      end if
!
      BoundarySide(isi)%previous_side = ksi
!
    end if
!
  end do
!
  return
  end subroutine DefineBoundarySideRelativeAngles2D
!---split

