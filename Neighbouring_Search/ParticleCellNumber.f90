!cfile ParticleCellNumber.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ParticleCellNumber
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
! Module purpose : Module Returns the number of the grid cell where particle np is 
!                  at the present instant 
!                  if particle is outside of the reference grids returns -1
!
! Calling routine: 
!
! Called routines: CellNumber
!
!************************************************************************************
!
!AA406 sub
  Integer(4) function ParticleCellNumber (coord)

!.. Returns the number of the grid cell where particle np is at the present instant
!.. if particle is outside of the reference grids returns -1
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
!
!AA406 sub
  double precision,intent(in)  :: coord(3)
!
!.. Local Scalars ..
  integer(4)       :: i, j, k
  double precision :: xp, yp, zp
!
!.. External Routines ..
  integer(4), external :: CellNumber
!
!.. Executable Statements ..
!
!.. evaluates the progressive cell pointer in all the directions with respect the lower grid coordinates
!
!AA406 sub
  xp = coord(1) - Grid%extr(1,1)
  yp = coord(2) - Grid%extr(2,1)
  zp = coord(3) - Grid%extr(3,1)
!
  i = ceiling(xp / Grid%dcd(1))
  k = ceiling(zp / Grid%dcd(3)) 
!
  if (ncord == 3) then
    j = ceiling(yp / Grid%dcd(2))
  else
    j = 1
  end if
!
!.. return the cell number (-1 if the particle is out of the grid)
!
  ParticleCellNumber = CellNumber(i, j, k)
!
  return
  End Function ParticleCellNumber
!---split

