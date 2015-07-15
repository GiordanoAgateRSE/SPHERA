!cfile MatrixTransposition.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : MatrixTransposition
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
! Module purpose : Module that Returns in AAT(n, m) the transponse of AA(m, n)
!
! Calling routine: BoundaryMassForceMatrix3D, BoundaryPressureGradientMatrix3D
!
! Called routines: 
!
!************************************************************************************
!
subroutine MatrixTransposition (AA, AAT, m, n )
!Returns in AAT(n, m) the transponse of AA(m, n)
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
integer(4),      intent(IN)   :: m
integer(4),      intent(IN)   :: n
double precision,intent(IN),   dimension(m,n) :: AA
double precision,intent(INOUT),dimension(n,m) :: AAT
!
!.. Local Scalars ..
integer(4) :: i,j
!
!.. Executable Statements ..
!
 do i = 1, n
    do j = 1, m
        AAT(i, j) = AA(j, i)
    end do
 end do

return
end subroutine MatrixTransposition
!---split

