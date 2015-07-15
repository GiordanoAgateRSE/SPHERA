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
! Module name     : Matrix_Inversion_3x3
!
!AA504 sub start
! Versions: 
! 01  Amicarelli      19Jul12       (creation) 
! 02  Amicarelli      08apr14       (v5.04) treatment of bad conditioned matrices      
!AA504 sub end
!
!************************************************************************************
! Module purpose : computation of the inverse (inv) of a provided 3x3 matrix (mat)
!
! Calling routine: RHS_body_dynamics
!
! Called routines: /
!
!************************************************************************************

!AA504 sub
subroutine Matrix_Inversion_3x3(mat,inv,test)

! Declarations
 implicit none
 double precision,intent(IN) :: mat(3,3) 
 double precision,intent(INOUT) :: inv(3,3)
!AA504
! integer(4),intent(inout),optional :: test  
 integer(4),intent(inout) :: test  
 double precision :: det

! Statements
 det = mat(1,1)*mat(2,2)*mat(3,3) + mat(2,1)*mat(3,2)*mat(1,3) + mat(3,1)*mat(1,2)*mat(2,3) &
      -mat(1,1)*mat(3,2)*mat(2,3) - mat(3,1)*mat(2,2)*mat(1,3) - mat(2,1)*mat(1,2)*mat(3,3)
!AA504 start
 if (det.ne.0) then
 test = 1
!AA504 end
 inv(1,1) = mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)
 inv(1,2) = mat(1,3)*mat(3,2) - mat(1,2)*mat(3,3)
 inv(1,3) = mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2)
 inv(2,1) = mat(2,3)*mat(3,1) - mat(2,1)*mat(3,3)
 inv(2,2) = mat(1,1)*mat(3,3) - mat(1,3)*mat(3,1)
 inv(2,3) = mat(1,3)*mat(2,1) - mat(1,1)*mat(2,3)
 inv(3,1) = mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1)
 inv(3,2) = mat(1,2)*mat(3,1) - mat(1,1)*mat(3,2)
 inv(3,3) = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)
 inv(:,:) = inv(:,:) / det
!AA504
 else
 test = 0
 endif

return
end subroutine Matrix_Inversion_3x3
!---split

