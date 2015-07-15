!cfile BoundaryPressureGradientMatrix3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : BoundaryPressureGradientMatrix3D
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
! Module purpose : Module to generate the pressure gradient matrix
!                  RRP(1 To SPACEDIM, 1 To SPACEDIM) on the base of the cosine
!                  matrix T(1 To SPACEDIM, 1 To SPACEDIM) and the parameters Psi()
!
! Calling routine: DefineBoundaryFaceGeometry3D
!
! Called routines: MatrixTransposition
!                  MatrixProduct
!
!************************************************************************************
!
subroutine BoundaryPressureGradientMatrix3D ( T, RGP, Psi )
!Generates the pressure gradient matrix RRP(1 To SPACEDIM, 1 To SPACEDIM) on the base of
!the cosine matrix T(1 To SPACEDIM, 1 To SPACEDIM) and the parameters Psi()
!
!.. assign modules
use GLOBAL_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
double precision,intent(INOUT),dimension(SPACEDIM,SPACEDIM) :: T
double precision,intent(INOUT),dimension(SPACEDIM,SPACEDIM) :: RGP
double precision,intent(INOUT),dimension(SPACEDIM)          :: Psi
!
!.. Local Scalars ..
Integer(4) :: i
!
!.. Local Arrays ..
double precision,dimension(SPACEDIM,SPACEDIM) :: Diag, PsiR, TTR
!
!.. Executable Statements ..
!
 Diag = zero
 do i = 1, SPACEDIM
   Diag(i, i) = Psi(i)
 end do

 call MatrixTransposition ( T, TTR, SPACEDIM, SPACEDIM )

 call MatrixProduct ( Diag, TTR, PsiR, SPACEDIM, SPACEDIM, SPACEDIM )

 call MatrixProduct ( T, PsiR, RGP, SPACEDIM, SPACEDIM, SPACEDIM )
!
return
!
end subroutine BoundaryPressureGradientMatrix3D
!---split

