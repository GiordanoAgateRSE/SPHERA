!cfile BoundaryMassForceMatrix3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : BoundaryMassForceMatrix3D
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
! Module purpose : Module to generate the boundary mass force matrix
!                  RN(1 To SPACEDIM, 1 To SPACEDIM) on the base of the cosine
!                  matrix T(1 To SPACEDIM, 1 To SPACEDIM) and the parameters Fi()
!
! Calling routine: DefineBoundaryFaceGeometry3D
!
! Called routines: MatrixTransposition
!                  MatrixProduct
!
!************************************************************************************
!
subroutine BoundaryMassForceMatrix3D ( T, RMF, Fi )
!Generates the boundary mass force matrix RN(1 To SPACEDIM, 1 To SPACEDIM)
!on the base of the cosine matrix T(1 To SPACEDIM, 1 To SPACEDIM) and the parameters Fi()
!
!.. assign modules
use GLOBAL_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
double precision,intent(INOUT),dimension(SPACEDIM,SPACEDIM) :: T
double precision,intent(INOUT),dimension(SPACEDIM,SPACEDIM) :: RMF
double precision,intent(INOUT),dimension(SPACEDIM)          :: Fi
!
!.. Local Scalars ..
integer(4) :: i
!
!.. Local Arrays ..
double precision, dimension(SPACEDIM,SPACEDIM) :: Diag, FiR, TTR
!
!.. Executable Statements ..
!
 Diag = zero
 do i = 1, SPACEDIM
   Diag(i, i) = Fi(i)
 end do
!
 call MatrixTransposition ( T, TTR, SPACEDIM, SPACEDIM )
!
 call MatrixProduct ( Diag, TTR, FiR, SPACEDIM, SPACEDIM, SPACEDIM )
!
 call MatrixProduct ( T, FiR, RMF, SPACEDIM, SPACEDIM, SPACEDIM )
!
return
!
end subroutine BoundaryMassForceMatrix3D
!---split

