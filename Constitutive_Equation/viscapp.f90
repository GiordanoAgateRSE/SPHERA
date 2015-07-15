!cfile viscapp.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : viscapp
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
! Module purpose : Module to define for each particle the viscosity based on 
!                  reology model
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
!
subroutine viscapp
!* definisce per ogni particella la viscosita' in base 
!* al modello reologico di appartenenza
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
double precision,parameter :: alfa = 1.0d0
!
!.. Local Scalars ..
integer(4)       :: npi
double precision :: mu, mumax, secinv,cuin, smalen, smalenq, visc1, visc2 
!
!.. External routines ..
character(80),external :: lcase
!
!.. Executable Statements ..
!
if (.not. diffusione) then
!
!$omp parallel do default(none) &
!$omp private(npi,visc1,smalen,smalenq,visc2,secinv,cuin,mu,mumax) &
!$omp shared(nag,pg,Med,Domain,it_corrente)
!
  do npi = 1,nag
    if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!
    select case ( Med(pg(npi)%imed)%tipo )
!
    case ( "liquid  " )
       pg(npi)%visc = Med(pg(npi)%imed)%visc
!   
    case ( "gas     " )
       pg(npi)%visc = Med(pg(npi)%imed)%visc
!   
    case ( "smagorin" )  
       visc1   = Med(pg(npi)%imed)%visc 
       smalen  = Med(pg(npi)%imed)%Cs * Domain%h               
       smalenq = smalen * smalen  
       visc2   = smalenq * two * pg(npi)%secinv

       pg(npi)%visc = visc1 + visc2 
!
    case ( "general " )
       secinv = two * pg(npi)%secinv
       cuin   = Med(pg(npi)%imed)%cuin
       mu     = Med(pg(npi)%imed)%taucri/(secinv+0.0001d0) + Med(pg(npi)%imed)%cons*((secinv+0.0001d0)**(cuin-one))
       mumax  = Med(pg(npi)%imed)%mumx
       pg(npi)%visc = min(mumax,mu) / pg(npi)%dens
!
!!!!!VERIFICARE
       if (pg(npi)%state == "sol") then
!.. controllo prime NIterSol iterazioni mezzo granulare non si muove
         if ((pg(npi)%visc /= mumax/pg(npi)%dens .and. it_corrente > Med(pg(npi)%imed)%NIterSol) .or. pg(npi)%kodvel == 2) then
           pg(npi)%state = "flu"
         end if
       end if
!!!!!VERIFICARE
!
    case ( "granular" )
!       secinv = two * pg(npi)%secinv        
!!       cuin   = Med(pg(npi)%imed)%cuin
!       pre    = (max(zero,pg(npi)%pres))
!       coes   = Med(pg(npi)%imed)%coes
!       coeff  = sin (Med(pg(npi)%imed)%phi)
!       mu     = (coes * cos(Med(pg(npi)%imed)%phi) + pre * coeff) / (secinv+0.00001d0)    !11-4 coesione
!       mumax  = Med(pg(npi)%imed)%mumx
!       pg(npi)%visc = min(mumax,mu) / pg(npi)%dens
       pg(npi)%visc = pg(npi)%mu / pg(npi)%dens
!
!.. verifies the "state" of a granular particle
!
!ç       if (pg(npi)%state == "sol") then
!.. controllo prime NIterSol iterazioni mezzo granulare non si muove
!AGGR         if (pg(npi)%visc /= mumax/pg(npi)%dens .and. it_corrente > Med(pg(npi)%imed)%NIterSol) then
!AGGR         if ((pg(npi)%visc /= mumax/pg(npi)%dens .and. it_corrente > Med(pg(npi)%imed)%NIterSol) .or. pg(npi)%kodvel == 2) then
!ç         if ((pg(npi)%visc /= mumax/pg(npi)%dens .and. it_corrente > Med(pg(npi)%imed)%NIterSol) .or. pg(npi)%kodvel == 2) then
!!!!         if (it_corrente > Med(pg(npi)%imed)%NIterSol) then
!ç           pg(npi)%state = "flu"
!ç         end if
!ç       end if
!
    case default
!
    end select
!
  end do
!
!$omp end parallel do
!
!
else if (diffusione) then
!!! modello bifluido
!
!$omp parallel do default(none) private(npi) shared(nag,pg,Med)
!
  do npi = 1,nag
    if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle

    if (index(Med(1)%tipo,"liquid") > 0 .and. index(Med(2)%tipo,"liquid") > 0) then

      pg(npi)%visc = pg(npi)%VolFra*Med(2)%visc + (one-pg(npi)%VolFra)*Med(1)%visc

    else if (index(Med(1)%tipo,"liquid") > 0 .and. index(Med(2)%tipo,"granular") > 0) then

!     if (pg(npi)%VolFra > 0.1d0) then
!§
      if (pg(npi)%VolFra >= VFmn .and. pg(npi)%VolFra <= VFmx) then
!       secinv = two * pg(npi)%secinv        
!       cuin   = Med(2)%cuin
!       pre    = (max(zero,pg(npi)%pres))
!       coes   = Med(2)%coes
!       coeff  = sin (Med(2)%phi)
!       mu     = pg(npi)%VolFra * (coes * cos(Med(2)%phi) + pre*coeff) / (secinv+0.00001d0) + &
!                (one-pg(npi)%VolFra) * Med(1)%visc
!       mumax  = Med(2)%mumx
!       pg(npi)%visc = min(mumax,mu) / pg(npi)%dens
        pg(npi)%visc = pg(npi)%VolFra*Med(2)%visc + (one-pg(npi)%VolFra)*Med(1)%visc
     
      else if (pg(npi)%VolFra < VFmn) then

        pg(npi)%visc = Med(1)%visc

      else if (pg(npi)%VolFra < VFmx) then

        pg(npi)%visc = Med(2)%visc
       
      end if
!§
!
!§
!!!!!VERIFICARE
!     if (pg(npi)%state == "sol") then
!       if (pg(npi)%visc /= (mumax/pg(npi)%dens) .or. pg(npi)%VolFra <= VFmx) then
!.. controllo prime NIterSol iterazioni mezzo granulare non si muove
!         if ( it_corrente > Med(pg(npi)%imed)%NIterSol) then
!           pg(npi)%state = "flu"
!         end if
!       end if
!     end if
!!!!!VERIFICARE
!§
    end if
!
  end do
!
!$omp end parallel do
!
end if
!
return
end subroutine viscapp
!---split

