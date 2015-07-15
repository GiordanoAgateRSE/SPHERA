!cfile calcpre.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : calcpre
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
! Module purpose : Module calculation of pressure in every particles of the general 
!                  field
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D, time_integration
!
! Called routines: 
!
!************************************************************************************
!
  subroutine CalcPre
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Local Scalars ..
  integer(4)       :: npi
  double precision :: rhorif,c2,crhorif,wrhorif,wc2,cc2 !,maxpExpl   !,presc,presw
!
!.. Executable Statements ..
!
!  maxpExpl = max_negative_number
!
!.. modello bifluido
  if (diffusione) then
!
!$omp parallel do default(none) private(npi,crhorif,wrhorif,wc2,cc2) shared(nag,Pg,Med)
!
    do npi = 1,nag
!
      if (pg(npi)%koddens /= 0) cycle
!
      if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!
      crhorif = Med(2)%den0
      wrhorif = Med(1)%den0
      wc2     = Med(1)%celerita * Med(1)%celerita
      cc2     = Med(2)%celerita * Med(2)%celerita
!§
      if (pg(npi)%imed == 1)then
!        presc = cc2 * (pg(npi)%rhoc * pg(npi)%VolFra - crhorif * VFmn)
!        presw = wc2 * (pg(npi)%rhow  * (one-pg(npi)%VolFra) - wrhorif * (one-VFmn))
        pg(npi)%pres = wc2 * (pg(npi)%dens - (crhorif * VFmn + wrhorif * (1-VFmn))) !ALTERNATIVO
      else if (pg(npi)%imed == 2)then
!        presc = cc2 * (pg(npi)%rhoc * pg(npi)%VolFra - crhorif * VFmx)
!        presw = wc2 * (pg(npi)%rhow  * (one-pg(npi)%VolFra) - wrhorif * (one-VFmx))
        pg(npi)%pres = cc2 * (pg(npi)%dens - (crhorif * VFmx + wrhorif * (1-VFmx))) !ALTERNATIVO
      end if
!      pg(npi)%pres = presc + presw
!§
!      presc = cc2 * (pg(npi)%rhoc - crhorif)
!      presw = wc2 * (pg(npi)%rhow - wrhorif)
!      pg(npi)%pres = pg(npi)%VolFra * presc + (one - pg(npi)%VolFra) * presw
!
    end do
!
!$omp end parallel do
!
  else
!
    if (it_corrente > 0) then  !inutile ??
!
!AA501 sub
!$omp parallel do default(none) private(npi,rhorif,c2) shared(nag,pg,Domain,Med,esplosione) !,maxpExpl)
!
!.. loops on all the particles
!
      do npi = 1,nag
!
!.. skips the outgone particles
!.. skips the particles with velocity type different from "standard"
!
        if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!
        if (esplosione) then
!...................................... 2011 mar 08
!.. modify for Specific Internal Energy
          pg(npi)%pres = (Med(pg(npi)%imed)%gamma - one) * pg(npi)%IntEn * pg(npi)%dens
!          maxpExpl = max(maxpExpl,pg(npi)%pres)
!......................................
        else
!.. evaluates the new pressure condition for the particle
!ç
!!!!        if (it_corrente > 0) then  !.and. pg(npi)%state /= 'sol') then
          rhorif     = Med(pg(npi)%imed)%den0
          c2         = Med(pg(npi)%imed)%eps / rhorif
!
!AA406
           if ((pg(npi)%dens - rhorif) /= 0.) pg(npi)%pres = c2 * (pg(npi)%dens - rhorif)
!!!!        end if
!
!AA406test
!
!AA601rm
        end if
!
      end do
!
!......................................
!!!!      if (maxpExpl < pg(??)%pres) esplosione = .false.
!      if (maxpExpl < 0.99e+9) esplosione = .false.
!......................................
!
!$omp end parallel do
!
    end if
!
  end if
!
  return
  end subroutine CalcPre
!---split



!cfile contrmach.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
!subroutine contrmach (k,celmax)
!
!use GLOBAL_MODULE
!use AdM_USER_TYPE
!use ALLOC_MODULE
!
!implicit none
!
!integer(4) :: k,
!
!double precision    :: celmax,amaxnmach,amaxnmach2,amachnumb2
!!double precision   :: vmod2,cel2,amachnumb
!
! amaxnmach  = 0.2d0
! cel2       = celmax*celmax
! amaxnmach2 = amaxnmach*amaxnmach
!
! controllo numero mach - scrive su nout
! vmod2 = pg(k)%vel(1)*pg(k)%vel(1) + pg(k)%vel(3)*pg(k)%vel(3)
! amachnumb2 = vmod2 / cel2
!
! if ( amachnumb2 > amaxnmach2 )then    ! -------------------------------
!    amachnumb = Dsqrt(amachnumb2)
!!   pg(k)%vel(:)=(amaxnmach*celmax/Dsqrt(vmod2))*pg(k)%vel(:) 
! end if        ! -------------------------------
!
!return
!end
!---split
