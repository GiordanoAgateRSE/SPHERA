!cfile MatrixProduct.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : MatrixProduct
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
! Module purpose : Module that Returns in CC() the product between matrices
!                  AA() and BB()
!
! Calling routine: BoundaryMassForceMatrix3D, BoundaryPressureGradientMatrix3D,
!AA501b
!                  vector_rotation
!
! Called routines: 
!
!************************************************************************************
!
subroutine MatrixProduct ( AA, BB, CC, nr, nrc, nc )
!Returns in CC() the product between matrices AA() and BB()
!nr =   number of rows of AA() and CC()
!nc =   number of columns of BB() and CC()
!nrc =  number of columns of AA() = number of rows of BB()
!
!.. assign modules
use GLOBAL_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
integer(4),      intent(IN) :: nr
integer(4),      intent(IN) :: nrc
integer(4),      intent(IN) :: nc
double precision,intent(IN),   dimension(nr,nrc) :: AA
double precision,intent(IN),   dimension(nrc,nc) :: BB
double precision,intent(INOUT),dimension(nr, nc) :: CC
!
!.. Local Scalars ..
integer(4)       :: i,j,k
double precision :: sum
!
!.. Executable Statements ..
!
 do i = 1, nr
    do j = 1, nc
       sum = zero
       do k = 1, nrc
          sum = sum + AA(i, k) * BB(k, j)
       end do
       CC(i, j) = sum
    end do
 end do
!
return
end subroutine MatrixProduct
!---split

