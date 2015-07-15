!cfile DefineBoundaryFaceGeometry3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : DefineBoundaryFaceGeometry3D
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
! Module purpose : Module to define boundary face from geometry 3D
!
! Calling routine: Gest_Input
!
! Called routines: BoundaryMassForceMatrix3D
!                  BoundaryPressureGradientMatrix3D
!                  DefineLocalSystemVersors
!
!************************************************************************************
!
subroutine DefineBoundaryFaceGeometry3D
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
integer(4) :: Kf, Nf, Nt
!
!.. Local Arrays ..
!double precision, dimension(SPACEDIM)          :: Psi
double precision, dimension(SPACEDIM)          :: Fi
double precision, dimension(SPACEDIM,SPACEDIM) :: TT, RGP, RMF

! 20081028  Data Psi /1.0d0, 1.0d0, 0.0d0/
! 20081028  Data Fi /0.0d0, 0.0d0, 1.0d0/
!Data Psi /0.0d0, 0.0d0, 0.0d0/
Data Fi /1.0d0, 1.0d0, 1.0d0/
!
!.. Executable Statements ..
!
  do Kf = 1, NumFacce, 1
!
    Nf = BFaceList(Kf)
!
    if ( Nf == 0 ) cycle
!
    Nt = BoundaryFace(Nf)%stretch

!AA504 sub
    call DefineLocalSystemVersors (Nf)
!
    !Psi(1:SPACEDIM) = Tratto(nt)%PsiCoeff(1:SPACEDIM)
    !Fi (1:SPACEDIM) = Tratto(nt)%FiCoeff(1:SPACEDIM)
    TT (1:SPACEDIM, 1:SPACEDIM) = BoundaryFace(nf)%T(1:SPACEDIM, 1:SPACEDIM)
!
    RGP = zero
    RMF = zero
! 20081028     call BoundaryPressureGradientMatrix3D ( TT, RGP, Psi )
    call BoundaryMassForceMatrix3D ( TT, RMF, Fi )
!
! 20081028    BoundaryFace(nf)%RPsi(1:SPACEDIM, 1:SPACEDIM) = RGP(1:SPACEDIM, 1:SPACEDIM)
    BoundaryFace(nf)%RFi (1:SPACEDIM, 1:SPACEDIM) = RMF(1:SPACEDIM, 1:SPACEDIM)
!
    if ( Tratto(nt)%tipo == "tapi" ) then
       BoundaryFace(nf)%velocity(1:SPACEDIM) = Tratto(nt)%velocity(1:SPACEDIM)
    else
       BoundaryFace(nf)%velocity(1:SPACEDIM) = zero
    end if
   !nv = nv + 1
!
  end do
!
return
end subroutine DefineBoundaryFaceGeometry3D
!---split

