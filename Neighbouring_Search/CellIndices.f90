!cfile CellIndices.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
  integer(4) function CellIndices (nc, i, j, k)

!.Returns indices i, j, k, of the cell nc
!..in a 3D box frame with ni, nj, nk, cells on each side
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
  integer(4),intent(IN)  :: nc
  integer(4),intent(OUT) :: i,j,k
!
!.. Local Scalars ..
  integer(4) :: ncij, nucellsij,ni,nj   !,nk
!
!.. Executable Statements ..
!
  ni = Grid%ncd(1)
  nj = Grid%ncd(2)
!  nk = Grid%ncd(3)
!
  nucellsij = ni * nj
!
!.. grid index in the Z direction
!
  k = int((nc - 1) / nucellsij) + 1
  ncij = nc - nucellsij * (k - 1)
!
!.. grid index in the Y direction
!
  j = int((ncij - 1) / ni) + 1
!
!.. grid index in the X direction
!
  i = ncij - ni * (j - 1)
!
  CellIndices = ncij
!
  End function CellIndices
!---split

