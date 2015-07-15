!cfile inter_CoefDif.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : inter_CoefDif
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
! Module purpose : Module to calculate the corrective term of the velocity
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D, time_integration
!
! Called routines: 
!
!************************************************************************************
!
subroutine inter_CoefDif (npi)
! ex inter6
!* implementa il meccanismo di ricerca delle particelle che agiscono
!* su quella i-esima.
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Local Parameters ..
!integer(4), parameter :: local_d = 500  ! num max part entro 2h
!
!.. Formal Arguments ..
integer(4),intent(IN) :: npi
!
!.. Local Scalars ..
integer(4)       :: npj,contj,npartint
double precision :: unity,rhoj,amassj,pesoj   !moddervel,  rhoi,
!
!.. Executable Statements ..
!
!* azzeramento quantita generali
 unity = zero
 pg(npi)%veldif(:) = zero
!
!*_______________________________________________________________
!*prima passata per trovare celle interagenti e memorizzazione
!
 do contj = 1, nPartIntorno(npi)
!
   npartint = (npi-1)* NMAXPARTJ + contj
   npj = PartIntorno(npartint)
!
!   rhoi   = pg(npi)%dens
   rhoj   = pg(npj)%dens
   amassj = pg(npj)%mass
!
!============= CALCOLO VELOCITA' PER COEFFICIENTE DIFFUSIVO ===================

   pesoj = amassj * PartKernel(4,npartint) / rhoj

   unity = unity + pesoj  

! calcolo velocita' per modello diffusivo
   pg(npi)%veldif(:) = pg(npi)%veldif(:) + pg(npj)%vel(:) * pesoj  

!============= FINE CALCOLO VELOCITA' PER COEFFICIENTE DIFFUSIVO ================
!
!prova!
!   pg(npi)%veldif = floor(pg(npi)%veldif * azzeramento) / azzeramento
!   where (dabs(pg(npi)%veldif) < arrotondamento) pg(npi)%veldif = zero 
!prova!
!
 end do
!
 pg(npi)%uni = unity
!
return
end subroutine inter_CoefDif
!---split


!cfile inter_SmoothVF.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : inter_SmoothVF
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!
!************************************************************************************
! Module purpose : Module to calculate the corrective term of the volume fraction
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
!
subroutine inter_SmoothVF (npi,appo1,unity)
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
integer(4),      intent(IN)    :: npi
double precision,intent(INOUT) :: appo1, unity
!
!.. Local Scalars ..
integer(4)       :: npj,contj,npartint
double precision :: voli,volj,rhoj,amassj,pesoj
!
!.. Executable Statements ..
!
!* azzeramento quantita generali
 unity = zero
 appo1 = zero
!
 do contj = 1, nPartIntorno(npi)
!
   npartint = (npi-1)* NMAXPARTJ + contj
   npj = PartIntorno(npartint)
!
   voli   = pg(npi)%VolFra
   rhoj   = pg(npj)%dens    
   volj   = pg(npj)%VolFra    
   amassj = pg(npj)%mass
!
   if ( pg(npj)%vel_type /= "std" ) cycle          !non part fix o altro
!
!.. calcola i pesi
!
   pesoj = amassj * PartKernel(4,npartint) / rhoj
!
   unity = unity + pesoj  
!
   appo1 = appo1 + (volj - voli ) * pesoj  
!
!============= FINE CALCOLO VELOCITA' PER COEFFICIENTE DIFFUSIVO ================
!
 end do
!
 pg(npi)%uni = unity
!
return
end subroutine inter_SmoothVF
!---split



!cfile AggDens.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name    : AggDens
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
!
!************************************************************************************
! Module purpose : 
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D, time_integration
!
! Called routines: 
!
!************************************************************************************
!
subroutine AggDens
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
integer(4) :: npi,ii
double precision :: tirhoc, tirhow, vdens      !,tirhoc1,tirhow1,tirhoc2,tirhow2
double precision :: a1, b1, c1     !,gam,qut,qu,appo, 
!double precision :: dop, denom, del, a, b, vol_temp, controllo,controllo2    !,gam,qut,qu,
!
!.. Executable Statements ..
!appo,
!$omp parallel do default(none) &
!$omp private(npi,vdens,tirhoc,tirhow,a1,b1,c1,ii) &  
!$omp shared(nag,Pg,Med,Domain,dt,it_corrente,ncord,indarrayFlu,Array_Flu)
!
!!!  do npi = 1,nag
!!!!
!!!    if (pg(npi)%koddens /= 0) cycle
!!!!
!!!!!    if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!!!!!gio12mar2009
!!!    if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" .or. pg(npi)%state == "sol") cycle
!$$$$$
  do ii = 1,indarrayFlu
    npi = Array_Flu(ii)
    if (pg(npi)%koddens /= 0) cycle
!$$$$$
!
    vdens  = pg(npi)%dden / pg(npi)%dens
!
!§
    if(it_corrente <= 1) then
      pg(npi)%tiroc = pg(npi)%rhoc * pg(npi)%VolFra
      tirhoc = pg(npi)%tiroc + dt * (pg(npi)%tiroc * vdens + pg(npi)%rhoc * pg(npi)%diffu - Med(2)%settlingcoef)
    else    
      tirhoc = pg(npi)%tiroc + dt * (pg(npi)%tiroc * vdens + pg(npi)%rhoc * pg(npi)%diffu - Med(2)%settlingcoef)
      pg(npi)%tiroc = tirhoc     
    end if
!§
!    tirhoc = pg(npi)%rhoc * pg(npi)%VolFra           ! rhograintilde     
!    tirhoc = tirhoc + dt * (tirhoc * vdens + pg(npi)%rhoc * pg(npi)%diffu - med(2)%settlingcoef)      ! new rhograintilde
!
    if (tirhoc >= (pg(npi)%dens / pg(npi)%VolFra)) then
      tirhoc = pg(npi)%dens
      tirhow = zero
      pg(npi)%VolFra = VFmx
      pg(npi)%rhoc = pg(npi)%dens
      pg(npi)%mass = pg(npi)%dens * (Domain%dd**ncord)
    else if (tirhoc <= zero) then
      tirhoc = zero
      tirhow = pg(npi)%dens
      pg(npi)%rhow = pg(npi)%dens
      pg(npi)%VolFra = VFmn
      pg(npi)%mass = pg(npi)%dens * (Domain%dd**ncord)
    else
!
!§
!      if(pg(npi)%imed==1)then
!        tirhow = (Med(2)%den0 * VFmn + Med(1)%den0 * (1-VFmn)) - tirhoc
!      else if(pg(npi)%imed==2)then
!        tirhow = (Med(2)%den0 * VFmx + Med(1)%den0 * (1-VFmx)) - tirhoc
!      end if
!§
      tirhow = pg(npi)%dens - tirhoc
!
      a1 = med(2)%den0 * med(2)%celerita*med(2)%celerita - med(1)%den0 * med(1)%celerita*med(1)%celerita
      b1 = - (med(2)%celerita*med(2)%celerita) * (tirhoc + med(2)%den0) - (med(1)%celerita*med(1)%celerita) * (tirhow - med(1)%den0)
      c1 = (med(2)%celerita*med(2)%celerita) * tirhoc 
!      appo = b1 * b1 - 4.d0 * a1 * c1
!
!      vol_temp = pg(npi)%VolFra
!if(( b1 * b1 - 4. * a1 * c1) <= zero) then
!continue
!end if
      pg(npi)%VolFra = (- b1 - Dsqrt( b1 * b1 - 4.d0 * a1 * c1)) / (two * a1)   ! nuovo calcolo della frazione di volume
!      pg(npi)%VolFra = vol_temp + (pg(npi)%VolFra - vol_temp) / 10
!
!§
!      if (pg(npi)%VolFra >= one) then
      if (pg(npi)%VolFra >= VFmx) then
        pg(npi)%VolFra  = VFmx
!§
        pg(npi)%rhoc = pg(npi)%dens 
        pg(npi)%rhow = Med(1)%den0 
        pg(npi)%mass = pg(npi)%dens * (Domain%dd**ncord)
!§
!      else if(pg(npi)%VolFra <= zero)then
      else if(pg(npi)%VolFra <= VFmn)then
        pg(npi)%VolFra = VFmn
!§
        pg(npi)%rhoc = Med(2)%den0 
        pg(npi)%rhow = pg(npi)%dens
        pg(npi)%mass = pg(npi)%dens * (Domain%dd**ncord)
      else
        pg(npi)%rhoc = tirhoc / pg(npi)%VolFra            ! new dens grain
        pg(npi)%rhow = tirhow / (one-pg(npi)%VolFra)      ! new dens wat
        pg(npi)%mass = (pg(npi)%VolFra * pg(npi)%rhoc + (one-pg(npi)%VolFra)*pg(npi)%rhow) * (Domain%dd**ncord)
      end if
!
    end if
!
  end do       
!
!$omp end parallel do
!
return
end subroutine AggDens
!---split

