!cfile body_to_smoothing_vel.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: body_to_smoothing_vel
!
!************************************************************************************
! Module purpose : contribution of body particles to pressure smoothing
!
! Creation: Amicarelli-Agate 13nov12
!
! Calling routines: inter_SmoothVelo_2D,inter_SmoothVelo_3D
!
! Called subroutines: /
!
!************************************************************************************
!
  subroutine body_to_smoothing_vel(dervel_mat,unity_vec)

! Assigning modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE

! Declarations
  implicit none
  integer(4) :: npi,j,npartint,npj
  double precision :: W_vol,dis
  double precision,dimension(nag,3),intent(inout) :: dervel_mat
  double precision,dimension(nag),intent(inout) :: unity_vec

! External functions
  double precision, external :: w

! Statements
! Loop over body particles (neighbours: fluid particles)
  do npi=1,n_body_part
! Loop over the neighbouring fluid particles 
     do j=1,nPartIntorno_bp_f(npi)
        npartint = (npi-1)* NMAXPARTJ + j
        npj = PartIntorno_bp_f(npartint)
        dis = dsqrt(dot_product(rag_bp_f(:,npartint),rag_bp_f(:,npartint)))
        W_vol = w(dis,Domain%h,Domain%coefke) * ((Domain%dd/dx_dxbodies)**ncord)
        unity_vec(npj) = unity_vec(npj) + W_vol
        dervel_mat(npj,:) = dervel_mat(npj,:) + (bp_arr(npi)%vel_mir(:)-pg(npj)%vel(:)) * W_vol
     enddo
  enddo
  
  return
  end subroutine body_to_smoothing_vel
!---split

