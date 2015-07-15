!cfile BoundaryMassForceMatrix2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : BoundaryMassForceMatrix2D
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
! Module purpose : Module for generation of generalized boundary mass force matrix
!
! Calling routine: DefineBoundarySideGeometry2D
!
! Called routines: 
!
!************************************************************************************
!
subroutine BoundaryMassForceMatrix2D (T, RN, FiS, FiN) 

!Generates the generalised boundary mass force matrix RN(1 To SPACEDIM, 1 To SPACEDIM)
!on the base of the cosine matrix T(1 To SPACEDIM, 1 To SPACEDIM) and the parameter Fi
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
double precision                                   :: FiS, FiN 
double precision,dimension(1:SPACEDIM, 1:SPACEDIM) :: T
double precision,dimension(1:SPACEDIM, 1:SPACEDIM) :: RN
!
!.. Executable Statements ..
!
RN(1, 1) = FiS * T(1, 1) * T(1, 1) + FiN * T(3, 1) * T(3, 1)
RN(1, 3) = (FiS - FiN) * T(1, 1) * T(3, 1) 
RN(3, 1) = RN(1, 3)         
RN(3, 3) = FiS * T(3, 1) * T(3, 1) + FiN * T(1, 1) * T(1, 1)
!
return
!
end subroutine BoundaryMassForceMatrix2D
!---split

