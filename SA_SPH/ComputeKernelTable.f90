!cfile ComputeKernelTable.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ComputeKernelTable
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
! Module purpose : Module to compute boundary contributions to rodivV
!
! Calling routine: Gest_Trans
!
! Called routines: diagnostic
!
!************************************************************************************
!
subroutine ComputeKernelTable
!
!Precomputes and stores in kerneltab(0:ktrows, 0:ktcols) the following values
!kerneltab(0:ktrows, 0) = rob = rb/h
!kerneltab(0:ktrows, 1) = Int W* ro2 dro         (from rob to 2)
!kerneltab(0:ktrows, 2) = Int dW*/dro ro dro     (from rob to 2)
!kerneltab(0:ktrows, 3) = Int dW*/dro ro^2 dro   (from rob to 2)
!kerneltab(0:ktrows, 4) = Int dW*/dro ro^3 dro   (from rob to 2)
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Scalars ..
integer(4)        :: nr
double precision  :: rob
character(len=lencard) :: nomsub = "ComputeKernelTable"
!
!.. External routines ..
double precision, external :: IWro2dro
double precision, external :: JdWsRn
!
!.. Executable Statements ..
!
ktdelta = two / INT_KERNELTABLE
!
if (ncord == 3) then
  do nr = 0, ktrows
    rob = ktdelta * nr
    kerneltab(nr,0) = rob
    kerneltab(nr,1) = IWro2dro(rob)
    kerneltab(nr,2) = JdWsRn(rob,3,1,1) * Unosusquareh
    kerneltab(nr,3) = JdWsRn(rob,3,2,1) * Unosuh
    kerneltab(nr,4) = JdWsRn(rob,3,3,1)
  end do
!
else if (ncord == 2) then
  call diagnostic (arg1=8,arg2=3,arg3=nomsub)
!  do nr = 0,ktrows
!    rob = ktdelta * nr
!    kerneltab(nr,0) = rob
!    kerneltab(nr,1) = IWro2dro(rob)
!    kerneltab(nr,2) = JdWsRn(rob,2,1,1) * Unosusquareh
!    kerneltab(nr,3) = JdWsRn(rob,2,2,1) * Unosuh
!    kerneltab(nr,4) = JdWsRn(rob,2,3,1)
!  end do
end if
!
return
end subroutine ComputeKernelTable
!---split

