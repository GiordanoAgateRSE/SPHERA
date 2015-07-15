!cfile CalcVarp.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : CalcVarp
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
! Module purpose : Module to calculate p on point (xpu,ypu)
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: GetVarPart
!
!************************************************************************************
!
subroutine CalcVarp 
!* calcola grandezze nei punti di monitor
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
integer(4) :: i,ii,jj,kk,pointcellnumber
double precision xp,yp,zp
type (TyCtlPoint) :: pglocal
!
!.. External routines ..
  integer(4), external   :: CellNumber
!
!.. Executable Statements ..
!
  if (Npointst < 1) return

  do i = 1, Npointst

    pglocal%coord(:) = control_points(i)%coord(:)
    pglocal%pres     = zero
    pglocal%dens     = zero
    pglocal%vel(:)   = zero
    pglocal%uni      = zero
    xp = pglocal%coord(1) - Grid%extr(1,1)
    yp = pglocal%coord(2) - Grid%extr(2,1)
    zp = pglocal%coord(3) - Grid%extr(3,1)
    ii = ceiling(xp / Grid%dcd(1))
    jj = ceiling(yp / Grid%dcd(2))
    kk = ceiling(zp / Grid%dcd(3)) 
    pointcellnumber = CellNumber(ii, jj, kk)
    pglocal%cella   = pointcellnumber

    call GetVarPart (pglocal)

    control_points(i)%pres   = pglocal%pres
    control_points(i)%dens   = pglocal%dens
    control_points(i)%vel(:) = pglocal%vel(:)
    control_points(i)%uni    = pglocal%uni
    control_points(i)%cella  = pglocal%cella

  end do

return
end subroutine CalcVarp
!---split

