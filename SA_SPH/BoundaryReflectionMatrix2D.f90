!cfile BoundaryReflectionMatrix2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : BoundaryReflectionMatrix2D
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
! Module purpose : Module for generation of generalized reflection matrix
!
! Calling routine: DefineBoundarySideGeometry2D
!
! Called routines: 
!
!************************************************************************************
!
subroutine BoundaryReflectionMatrix2D (T, R, PsiS, PsiN)

!Generates the generalised reflection matrix R(1 To SPACEDIM, 1 To SPACEDIM) on the base of
!the cosine matrix T(1 To SPACEDIM, 1 To SPACEDIM) and the parameters PsiS, PsiN
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
double precision   :: PsiS, PsiN 
double precision,dimension(1:SPACEDIM, 1:SPACEDIM)  :: T
double precision,dimension(1:SPACEDIM, 1:SPACEDIM)  :: R
!
!.. Executable Statements ..
!


R(1, 1) = PsiS * T(1, 1) * T(1, 1) + PsiN * T(3, 1) * T(3, 1) 
R(1, 3) = (PsiS - PsiN) * T(1, 1) * T(3, 1)
R(3, 1) = R(1, 3)
R(3, 3) = PsiS * T(3, 1) * T(3, 1) + PsiN * T(1, 1) * T(1, 1)
!
return
!
end subroutine BoundaryReflectionMatrix2D
!---split

