!cfile LocalNormalCoordinates.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : LocalNormalCoordinates
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
! Module purpose : Module that Given the local coordinates PX(1 to 2) of a point P
!                  laing on the plane of the boundary face nf, the procedure
!                  computes in csi(1 to 3) the normal coordinates of the point Q
!                  corresponding to P in the inverse linear tranformation
!
! Calling routine: AddElasticBoundaryReaction_3D, ComputeBoundaryVolumeIntegrals_P0,
!                  FindCloseBoundaryFaces3D, CancelOutgoneParticles_3D, IsParticleInternal3D
! 
! Called routines: 
!
!************************************************************************************
!
subroutine LocalNormalCoordinates ( PX, csi, nf )
!Given the local coordinates PX(1 to 2) of a point P laing on the plane
!of the boundary face nf, the procedure
!computes in csi(1 to 3) the normal coordinates of the point Q
!corresponding to P in the inverse linear tranformation
!
!.. assign modules
use GLOBAL_MODULE
use AdM_User_Type
use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
double precision,intent(IN),   dimension(1:SPACEDIM) :: PX
double precision,intent(INOUT),dimension(1:SPACEDIM) :: csi
integer(4),      intent(IN)                          :: nf
!
!.. Local Scalars ..
integer(4)       :: i, j, k, nodes, fkod
double precision :: AA, BB, CC, DueArea, UsuDueArea, xj, yj, xk, yk
!
!.. Local Arrays ..
integer(4), dimension(3)   :: iseg = (/ 2,3,1 /)
integer(4), dimension(2,3) :: mainod ! = (/ 1,1, 2,3, 3,4 /) ! modifica per compatibilita xlf90
!
!.. Executable Statements ..
!
! modifica per compatibilita xlf90
  mainod(1, 1) = 1
  mainod(1, 2) = 2
  mainod(1, 3) = 3
  mainod(2, 1) = 1
  mainod(2, 2) = 3
  mainod(2, 3) = 4
! fine modifica
!
  nodes = 4
  if ( BoundaryFace(nf)%Node(4)%name <= 0 ) nodes = 3
  fkod = nodes - 2                    !=1 (triangolo), =2 (parallelogramma)
  DueArea = (3 - fkod) * BoundaryFace(nf)%Area
  UsuDueArea = one / DueArea
  do i = 1, 2
    j = iseg(i)
    k = iseg(j)
    xj = BoundaryFace(nf)%Node(mainod(fkod, j))%LX(1)
    yj = BoundaryFace(nf)%Node(mainod(fkod, j))%LX(2)
    xk = BoundaryFace(nf)%Node(mainod(fkod, k))%LX(1)
    yk = BoundaryFace(nf)%Node(mainod(fkod, k))%LX(2)
!
    AA = xj * yk - xk * yj
    BB = yj - yk
    CC = xk - xj
!
    csi(i) = (AA + BB * PX(1) + CC * PX(2)) * UsuDueArea
  end do
!
  csi(3) = one - (csi(1) + csi(2))
!
return
end subroutine LocalNormalCoordinates

!---split

