!AA501b the whole subroutine
!AA504 sub
!cfile Body_dynamics.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : Matrix_Inversion_2x2
!
! Creation : Amicarelli-Agate, 19 July 2012
!
!************************************************************************************
! Module purpose : computation of the inverse (inv) of a provided 3x3 matrix (mat)
!
! Calling routine: RHS_body_dynamics
!
! Called routines: /
!
!************************************************************************************

subroutine Matrix_Inversion_2x2(mat,inv)

! Declarations
 implicit none
 double precision,intent(IN) :: mat(2,2) 
 double precision,intent(INOUT) :: inv(2,2)
 double precision :: det

! Statements
 det = mat(1,1)*mat(2,2) - mat(2,1)*mat(1,2)
 inv(1,1) = mat(2,2)
 inv(1,2) = -mat(2,1)
 inv(2,1) = -mat(1,2)
 inv(2,2) = mat(1,1)
 inv(:,:) = inv(:,:) / det

return
end subroutine Matrix_Inversion_2x2
!---split

