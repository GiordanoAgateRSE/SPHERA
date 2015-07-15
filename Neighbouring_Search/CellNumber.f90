!cfile CellNumber.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
  integer(4) function CellNumber (i, j, k)

!Returns the number of the cell sited at integer coordinates i, j, k
!in a 3D box frame with ni, nj, nk cells on each side
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
  integer(4),intent(IN) :: i,j,k
!
!.. scalars
  integer(4) :: ni,nj,nk
!
!.. Executable Statements ..
!
  ni = Grid%ncd(1)
  nj = Grid%ncd(2)
  nk = Grid%ncd(3)
!
  if ( i<1 .or. i>ni .or. j<1 .or. j>nj .or. k<1 .or. k>nk) then
!
!.. the cell is outside the grid limits
!
    CellNumber = 0  !-1
!
  else
!
!.. return the cell number
!
    CellNumber = ((k - 1) * nj + (j - 1)) * ni + i
!
  end if
!
  End Function CellNumber
!---split

