!cfile body_to_smoothing_pres.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: body_to_smoothing_pres
!
! Creation: Amicarelli-Agate 18October12
!
!************************************************************************************
! Module purpose : contribution of body particles to pressure smoothing
!
! Calling routines: PressureSmoothing_2D,PressureSmoothing_3D
!
! Called subroutines: /
!
! Creation: Amicarelli-Agate 13nov12
!
!************************************************************************************
!
  subroutine body_to_smoothing_pres(sompW_vec,AppUnity_vec)

! Assigning modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE

! Declarations
  implicit none
  integer(4) :: npi,j,npartint,npj
  double precision :: W_vol,dis
  double precision,dimension(nag),intent(inout) :: sompW_vec,AppUnity_vec

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
        AppUnity_vec(npj) = AppUnity_vec(npj) + W_vol
        sompW_vec(npj) = sompW_vec(npj) + (bp_arr(npi)%pres-pg(npj)%pres) * W_vol
     enddo
  enddo
  
  return
  end subroutine body_to_smoothing_pres
!---split

