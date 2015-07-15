!cfile time_integration.f90
!AA401 all the subroutine
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : time_integration
!
! Last updating : April 18, 2013
!
! Improvement traceback:
! 00  Amicarelli       30Jun11   First development
! 01  Amicarelli/Agate 30Jun11   Adaptation for RK2 scheme
! 02  Amicarelli/Agate 18apr13   Adaptation for sloshing case
!
!************************************************************************************
! Module purpose : Explicit Runge-Kutta time integration
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: start_and_stop, Euler, Heun, inter_CoefDif, 
!                  aggdens, calcpre, diagnostic
!
!************************************************************************************
!
  subroutine time_integration  
!
!.. assign modules
  use FILES_ENTITIES
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
  use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Scalar declaration
  integer(4) :: ii,npi
  double precision :: appo1,appo2,appo3
  character(len=lencard)  :: nomsub = "time_integration"
!
!.. Executable Statements ..
!
! Explicit Runge-Kutta time integration scheme (Euler-RK1-, Heun-RK2, RK3c -classic RK3-, RK4c -classic RK4-)
  call start_and_stop(2,17)
  select case (Domain%RKscheme)
    case (1) 
      call Euler
!AA402 start
    case (2) 
      if (Domain%time_stage == 1) then
        call Euler
      else
        call Heun
      end if
!AA402 end
  end select
  call start_and_stop(3,17)
! 
! diffusive coefficient update
  if (diffusione) then
    call start_and_stop(2,15)
!.. update the diffusion coefficient (diffused solid-liquid intefrace)
!$omp parallel do default(none) private(npi,ii,appo1,appo2,appo3) shared(nag,Pg,Med,indarrayFlu,Array_Flu)
    do ii = 1,indarrayFlu
      npi = Array_Flu(ii)
      if (pg(npi)%VolFra == VFmx .and. pg(npi)%visc == Med(pg(npi)%imed)%mumx / pg(npi)%dens) then
        pg(npi)%coefdif = zero
      else
        call inter_CoefDif (npi)
        if (pg(npi)%uni > zero) pg(npi)%veldif = pg(npi)%veldif / pg(npi)%uni  !§
        appo1 = (pg(npi)%veldif(1)-pg(npi)%var(1)) * (pg(npi)%veldif(1)-pg(npi)%var(1))
        appo2 = (pg(npi)%veldif(2)-pg(npi)%var(2)) * (pg(npi)%veldif(2)-pg(npi)%var(2))
        appo3 = (pg(npi)%veldif(3)-pg(npi)%var(3)) * (pg(npi)%veldif(3)-pg(npi)%var(3))
        pg(npi)%coefdif = pg(npi)%coefdif * Dsqrt(appo1+appo2+appo3)
      end if
    end do
!$omp end parallel do
    call start_and_stop(3,15)
  end if
!solid-liquid interface diffusion model
  if (diffusione) then
    call start_and_stop(2,16)
    call aggdens
    call start_and_stop(3,16)
  end if
! State equation! 
  call start_and_stop(2,13)
!.. state equation according to the speed of sound definition and eventually to the liquid-solid diffusion model (Sphera_Tools.f90)
       call calcpre  
  call start_and_stop(3,13)
!
  return
  end subroutine time_integration
!---split

