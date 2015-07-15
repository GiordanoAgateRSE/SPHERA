!cfile body_pressure_mirror.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: body_pressure_mirror
!
! Creation: Amicarelli-Agate 13nov12
!
!************************************************************************************
! Module purpose : Computation of pressure for body particles (a mirror particle technique)
!
! Calling routines: Loop_Irre_2D,Loop_Irre_3D,RHS_body_dynamics
!
! Called subroutines: /
!
!************************************************************************************
!
  subroutine body_pressure_mirror

! Assigning modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE

! Declarations
  implicit none
  integer(4) :: npi,j,npartint,npj
  double precision :: Sum_W_vol,W_vol,dis,pres_mir,aux
  double precision :: aux_acc(3)

! External functions
  double precision, external :: w

! Statements
! Loop over body particles (discretized semi-analytic approach: a mirror particle technique)
!$omp parallel do default(none) &
!$omp private(npi,Sum_W_vol,j,npartint,npj,aux_acc,pres_mir,dis,W_vol,aux) &
!$omp shared(n_body_part,bp_arr,nPartIntorno_bp_f,NMAXPARTJ,PartIntorno_bp_f,Domain,pg,rag_bp_f)
  do npi=1,n_body_part
     bp_arr(npi)%pres = 0.
     Sum_W_vol = 0.  
     do j=1,nPartIntorno_bp_f(npi)
        npartint = (npi-1)* NMAXPARTJ + j
        npj = PartIntorno_bp_f(npartint)

!AA501btest start        
        aux = dsqrt(dot_product(bp_arr(npi)%acc,bp_arr(npi)%acc))
! wall acceleration should be less than 100m/2^2, otherwise an impulse is assumed to occur &
! and the formulation with acc_body is not valid

!AA501btest
        if (aux <= 100.) then
!if (aux <= 10.) then

           aux_acc(:) = Domain%grav(:) - bp_arr(npi)%acc(:)
           else

!AA501btest           
              aux_acc(:) = Domain%grav(:) - (100./aux) * bp_arr(npi)%acc(:)
!   aux_acc(:) = Domain%grav(:) - (10./aux) * bp_arr(npi)%acc(:)

        endif
!AA501btest end        

        pres_mir = pg(npj)%pres - pg(npj)%dens * dot_product(aux_acc,-rag_bp_f(:,npartint))
        dis = dsqrt(dot_product(rag_bp_f(:,npartint),rag_bp_f(:,npartint)))
        W_vol = w(dis,Domain%h,Domain%coefke) * (pg(npj)%mass/pg(npj)%dens)
        bp_arr(npi)%pres = bp_arr(npi)%pres + pres_mir * W_vol
        Sum_W_vol = Sum_W_vol + W_vol  
                  
     end do
     if (Sum_W_vol > 0.) bp_arr(npi)%pres = bp_arr(npi)%pres / Sum_W_vol      
  end do
 !$omp end parallel do
 
  return
  end subroutine body_pressure_mirror
!---split

