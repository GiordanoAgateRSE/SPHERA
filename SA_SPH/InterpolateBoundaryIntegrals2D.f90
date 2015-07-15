!cfile InterpolateBoundaryIntegrals2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : InterpolateBoundaryIntegrals2D
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
! Module purpose : Module for Interpolation in table "BoundIntegralTab(:,:)", 
!                  defined in module "BoundIntegralTab_Module", the values in 
!                  columns "Colmn(nc), nc=1, Ncols" corresponding to the input
!                  value "x" to be interpolated, in turn, in column 0.
!
! Calling routine: ComputeBoundaryDataTab
!
! Called routines: 
!
!************************************************************************************
!
subroutine InterpolateBoundaryIntegrals2D (x, Ncols, Colmn, Func)

!Interpolates in table "BoundIntegralTab(:,:)", defined in module "BoundIntegralTab_Module",  
!the values in columns "Colmn(nc), nc=1, Ncols" corresponding to the input value "x" to be
!interpolated, in turn, in column 0.
!Returns:
!Func(nc), nc=1, Ncols    =    Values interpolated in columns Col(nc), nc=1, Ncols
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use BoundIntegralTab_Module
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
double precision                          :: x 
integer(4)                                :: Ncols
integer(4),      dimension(1:NUMCOLS_BIT) :: Colmn
double precision,dimension(1:NUMCOLS_BIT) :: Func
!
!.. Local Scalars ..
integer(4)       :: nc, i, j
double precision :: xi, fi, fip1   !xip1, 
!
!.. Executable Statements ..
!
  if (x < BoundIntegralTab2D(1,0)) x = BoundIntegralTab2D(1,0)
!
  i = Int((x-BoundIntegralTab2D(1,0))/DELTAX_BIT)+1
  xi = BoundIntegralTab2D(i,0)
!  xip1 = BoundIntegralTab2D(i+1,0)

  do nc = 1, Ncols
    j = Colmn(nc)
    fi = BoundIntegralTab2D(i,j)
    fip1 = BoundIntegralTab2D(i+1,j)
    Func(nc) = fi+(fip1-fi)*(x-xi)/DELTAX_BIT
  end do
!
return
end subroutine InterpolateBoundaryIntegrals2D
!---split

