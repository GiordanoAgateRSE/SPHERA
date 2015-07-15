!cfile EvaluateBER_TimeStep.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : EvaluateBER_TimeStep
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
! Called routines: 
!
!************************************************************************************
!
  subroutine EvaluateBER_TimeStep
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
  integer(4) :: imed
  double precision :: DTlocal
!
!.. Executable Statements ..
!
  DTminBER = 1.0d30
  do imed = 1,NMedium
    DTlocal = Domain%h / Med(imed)%celerita
    if (DTlocal < DTminBER) DTminBER = DTlocal
  end do
  DTminBER = sqrttwo * DTminBER
!
  return
  end subroutine EvaluateBER_TimeStep
!---split

