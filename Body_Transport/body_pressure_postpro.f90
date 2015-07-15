!cfile body_pressure_postpro.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: body_pressure_postpro
!
! Creation: Amicarelli-Agate 13nov12
!
!************************************************************************************
! Module purpose : Computation of pressure for body partuicles
!
! Calling routines: Loop_Irre_2D,Loop_Irre_3D
!
! Called subroutines: /
!
!************************************************************************************
!
  subroutine body_pressure_postpro

! Assigning modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE

! Declarations
  implicit none
  integer(4) :: npi,j,npartint,npj
  double precision :: Sum_W_vol,W_vol,dis,mod_normal
  double precision :: dis_vec(3)
  double precision,dimension(:),allocatable :: aux_pres
  integer,dimension(:),allocatable :: wet
  double precision :: aux

! External functions
  double precision, external :: w

! Allocations
  allocate(aux_pres(n_body_part))
  allocate(wet(n_body_part))

!Initializations
  wet = 0.

! Statements

!AA501btest start
! Loop over the body particles (neighbours: body particles at the surface of the same body)
!$omp parallel do default(none) private(Sum_W_vol,npi,npj,dis_vec,dis,mod_normal,W_vol) &
!$omp shared(aux_pres,n_body_part,bp_arr,Domain,wet)
  do npi=1,n_body_part
     aux_pres(npi) = 0.  
     mod_normal = dsqrt(dot_product(bp_arr(npi)%normal,bp_arr(npi)%normal))

!AA501btest     
!     if ((mod_normal>0.).and.(wet(npi)==1)) then
      if (mod_normal>0.) then

        Sum_W_vol = 0.
! Loop over the body particles at the surface of the same body
        do npj=1,n_body_part
           mod_normal = dsqrt(dot_product(bp_arr(npj)%normal,bp_arr(npj)%normal))
           if ( (bp_arr(npi)%body == bp_arr(npj)%body) .and. (mod_normal>0.) ) then
              dis_vec(:) = bp_arr(npj)%pos(:) - bp_arr(npi)%pos(:) 
              dis = dsqrt(dot_product(dis_vec,dis_vec))
              W_vol = w(dis,Domain%h,Domain%coefke) * bp_arr(npj)%mass
              aux_pres(npi) = aux_pres(npi) + bp_arr(npj)%pres * W_vol
              Sum_W_vol = Sum_W_vol + W_vol 
           endif
        end do
        if (Sum_W_vol > 0.000001) aux_pres(npi) = aux_pres(npi) / Sum_W_vol
     endif
  enddo
!$omp end parallel do
! Loop over the body particles: update of pressure values 
!$omp parallel do default(none) private(npi) shared(bp_arr,aux_pres,n_body_part)
  do npi=1,n_body_part

!AA501btest
     if (aux_pres(npi).ne.0.) bp_arr(npi)%pres = aux_pres(npi)

  enddo
!$omp end parallel do
!AA501btest end

! Deallocations
  deallocate(aux_pres)
  deallocate(wet)

  return
  end subroutine body_pressure_postpro
!---split

