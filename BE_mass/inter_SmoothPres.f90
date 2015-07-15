!cfile inter_SmoothPres.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : inter_SmoothPres
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
! Module purpose : Module to calculate the corrective term of the pressure
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
!
subroutine inter_SmoothPres
! ex inter7
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
integer(4)       :: npi,npj,contj,npartint,ii
double precision :: unity,presi,presj,rhoj,amassj,pesoj,appo1,appo2,TetaP1
!
!.. Executable Statements ..
!
!$omp parallel do default(none) &
!$omp private(npi,ii,unity,appo1,appo2,contj,npartint,npj,presi,rhoj,presj,amassj,pesoj) &
!$omp shared(nag,pg,Med,Domain,nPartIntorno,NMAXPARTJ,PartIntorno,PartKernel,indarrayFlu,Array_Flu)
!
!.. loops on all the particles
!
!!!!!!  do npi = 1,nag
!!!!!!!
!!!!!!    if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" .or. pg(npi)%state == "sol") cycle
!$$$$$
  do ii = 1,indarrayFlu
    npi = Array_Flu(ii)
!$$$$$
!
!* azzeramento quantita generali
    unity = zero
    appo1 = zero
    appo2 = zero
!
    do contj = 1, nPartIntorno(npi)
!
      npartint = (npi-1)* NMAXPARTJ + contj
      npj = PartIntorno(npartint)
!
      if ( pg(npj)%vel_type /= "std" ) cycle          !non part fix o altro
!
      presi  = pg(npi)%pres
      rhoj   = pg(npj)%dens    
      presj  = pg(npj)%pres    
      amassj = pg(npj)%mass
!
!.. calcola i pesi
      pesoj = amassj * PartKernel(4,npartint) / rhoj
!
      unity = unity + pesoj  
!
      appo1 = appo1 + (presj - presi) * pesoj  
      appo2 = appo2 - Domain%grav(3)*Med(pg(npi)%imed)%den0 * (pg(npj)%coord(3) - pg(npi)%coord(3)) * pesoj     !aggiunto SaMa
!
    end do
!
!--------------- SaMa -----------------
    if (Domain%Psurf /= 's') then
      pg(npi)%vpres = appo1
    else if (unity > 0.8d0) then
      pg(npi)%vpres = appo1
    else
      pg(npi)%vpres = appo1 + appo2
    end if
!
    pg(npi)%uni = unity
!--------------- SaMa -----------------
!
  end do
!
!$omp end parallel do
!
!  call cpu_time(cpu_loop1a)
!  call cpu_time(cpu_loop2a)
!
!$omp parallel do default(none) private(npi,ii,TetaP1) shared(nag,Pg,Med,Domain,dt,indarrayFlu,Array_Flu,esplosione)
!
!.. applies the correction to all the particles
!
!!!!!!         do npi = 1,nag
!!!!!!!
!!!!!!           if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" .or. pg(npi)%state == "sol") cycle
!$$$$$
  do ii = 1,indarrayFlu
    npi = Array_Flu(ii)
!$$$$$
!
    if (esplosione) then
!.. con Csound al posto di Celerita e' circa uguale
      TetaP1 = Domain%TetaP * pg(npi)%Csound * dt / Domain%h
    else
! calcolo TetaP adeguato al passo temporale
      TetaP1 = Domain%TetaP * Med(pg(npi)%imed)%Celerita * dt / Domain%h
    end if

!
!    if (pg(npi)%densass == 0) then
!.. updates the pressure 
    pg(npi)%pres = pg(npi)%pres + TetaP1 * pg(npi)%vpres / pg(npi)%uni
!
!.. updates the density depending on the local pressure, reference medium density and the comprimibility eps
    pg(npi)%dens = Med(pg(npi)%imed)%den0 * (one + pg(npi)%pres / Med(pg(npi)%imed)%eps)
!    end if
  end do

!$omp end parallel do
!
return
end subroutine inter_SmoothPres
!---split

