!AA501b the whole subroutine
!AA504 sub
!cfile vector_rotation.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : vector_rotation
!
! Creation : Amicarelli-Agate, 2Aug12
!
!************************************************************************************
! Module purpose : Rotation of a given vector provided the vector of the rotation 
!                  angles (3D)
!
! Calling routine: time_integration_body_dynamics,Input_Body_Dynamics
!
! Called routines: MatrixProduct
!
!************************************************************************************

subroutine vector_rotation(vector,angle)

! Declarations
 implicit none
 double precision,intent(IN) :: angle(3) 
 double precision,intent(INOUT) :: vector(3)
 double precision               :: vec_temp(3)
 double precision               :: mat1_temp(3,3),mat2_temp(3,3),mat3_temp(3,3),mat4_temp(3,3),cos_dir(3,3)
 
! Statements
  mat1_temp(1,1) = 1.
  mat1_temp(1,2) = 0.
  mat1_temp(1,3) = 0.
  mat1_temp(2,1) = 0.
  mat1_temp(2,2) = cos(angle(1))
  mat1_temp(2,3) = sin(-angle(1))
  mat1_temp(3,1) = 0.
  mat1_temp(3,2) = -sin(-angle(1))
  mat1_temp(3,3) = cos(-angle(1))
  mat2_temp(1,1) = cos(-angle(2))
  mat2_temp(1,2) = 0.
  mat2_temp(1,3) = -sin(-angle(2))
  mat2_temp(2,1) = 0. 
  mat2_temp(2,2) = 1.
  mat2_temp(2,3) = 0.
  mat2_temp(3,1) = sin(-angle(2)) 
  mat2_temp(3,2) = 0.
  mat2_temp(3,3) = cos(-angle(2))
  mat3_temp(1,1) = cos(-angle(3))
  mat3_temp(1,2) = sin(-angle(3))
  mat3_temp(1,3) = 0.
  mat3_temp(2,1) = -sin(-angle(3))
  mat3_temp(2,2) = cos(-angle(3))
  mat3_temp(2,3) = 0.
  mat3_temp(3,1) = 0. 
  mat3_temp(3,2) = 0.
  mat3_temp(3,3) = 1.     
  call MatrixProduct(mat1_temp,mat2_temp,mat4_temp,3,3,3)
  call MatrixProduct(mat4_temp,mat3_temp,cos_dir,3,3,3)
  vec_temp(:) = vector(:) 
  vector(1) = dot_product(cos_dir(1,:),vec_temp)
  vector(2) = dot_product(cos_dir(2,:),vec_temp)
  vector(3) = dot_product(cos_dir(3,:),vec_temp)

return
end subroutine vector_rotation
!---split

