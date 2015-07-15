!cfile MohrC.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : MohrC
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  G. Agate, S. Manenti      Initial development of the code
!
!************************************************************************************
! Module purpose : MohrC Erosion Model.
! 
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: diagnostic
!
!************************************************************************************
!
subroutine MohrC 
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
!.. Local Parameters ..
!
!.. Formal Arguments ..
!
!.. Local Scalars ..
integer(4) :: npi, j, k, m
integer(4) :: intpl_id, imed, pl_imed, intliq_id, intsol_id, ncelcorr, igridi, jgridi, kgridi, iappo 
integer(4) :: numcell, nnlocal, flag, nsp, npartint, npj
double precision :: PartLiq_pres, peloloc, interf_liq, interf_sol, preidro, pretot, preeff, Velocity2, appo1
double precision :: secinv, coeff1, coeff2, mu, mumax
character(len=lencard)  :: nomsub = "Mohr-Coulomb"
!
!.. External Routines ..
integer(4),external :: ParticleCellNumber, CellIndices, CellNumber
!
!.. Executable Statements ..
!
!$omp parallel do default(none) &
!$omp private(npi,imed,ncelcorr,iappo,igridi,jgridi,kgridi,flag,nsp,appo1,j,npartint,npj,k,numcell,nnlocal,m) &
!$omp private(intpl_id,peloloc,pl_imed,interf_liq,intliq_id,interf_sol,intsol_id) &
!$omp private(secinv,PartLiq_pres,preidro,pretot,preeff,coeff1,coeff2,mu,mumax,Velocity2) &
!$omp shared(ind_interfaces,nomsub,nout,nscr) &
!$omp shared(nag,pg,Med,Domain,Grid,Icont,npartord,nPartIntorno,PartIntorno,NMAXPARTJ,diffusione,esplosione,it_corrente)
!
  do npi = 1,nag
!
    if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!
!.. controllo movimento granulare e calcolo viscosita'
!
    imed = pg(npi)%imed
    if (index(Med(imed)%tipo,"granular") <= 0) cycle
!
!.. controllo se la particella si muove in caso di esplosione non applico il modello di erosione  
    if (esplosione) then
      Velocity2 = pg(npi)%vel(1)*pg(npi)%vel(1) + pg(npi)%vel(2)*pg(npi)%vel(2) + pg(npi)%vel(3)*pg(npi)%vel(3)
      if (Velocity2 > 1.0e-3) cycle
    end if
!
    ncelcorr = ParticleCellNumber(pg(npi)%coord)
    iappo = CellIndices(ncelcorr,igridi,jgridi,kgridi)
!       
    intpl_id  = ind_interfaces(igridi,jgridi,1)
    intliq_id = ind_interfaces(igridi,jgridi,2)
    intsol_id = ind_interfaces(igridi,jgridi,3)
!
    if (intliq_id == 0) then
!.. localizzazione dell'interfaccia liquida su tutte le celle in direzione i e j del dominio
      interf_liq = max_positive_number
      intliq_id = 0
!.. ciclo sulle celle sopra e quella corrente compresa
      do k = 1,Grid%ncd(3)
        numcell = CellNumber(igridi,jgridi,k)
        if (numcell == 0) cycle
!.. ciclo su tutte le particelle della cella numcell
        do m = Icont(numcell),Icont(numcell+1)-1
          nnlocal = npartord(m)
          if (index(Med(pg(nnlocal)%imed)%tipo,"liquid") > 0) then
!.. ricerca interfaccia liquid/granular (particella liquid flu)
            if (pg(nnlocal)%Coord(3) < interf_liq) then
              interf_liq = pg(nnlocal)%Coord(3)
              intliq_id = nnlocal
            end if
          end if
        end do
      end do
    end if
!
    if (intsol_id == 0) then
!.. localizzazione dell'interfaccia solida su tutte le celle in direzione i e j del dominio
      interf_sol = max_negative_number
      intsol_id = 0
!.. ciclo sulle celle sopra e quella corrente compresa
      do k = 1,Grid%ncd(3)
        numcell = CellNumber(igridi,jgridi,k)
        if (numcell == 0) cycle
!.. ciclo su tutte le particelle della cella numcell
        do m = Icont(numcell),Icont(numcell+1)-1
          nnlocal = npartord(m)
          if (index(Med(pg(nnlocal)%imed)%tipo,"granular") > 0) then
!.. ricerca interfaccia liquid/granular (particella sol)
            if (pg(nnlocal)%Coord(3) > interf_sol) then
              interf_sol = pg(nnlocal)%Coord(3)
              intsol_id = nnlocal
            end if
          end if
        end do
      end do
    end if
!
    if (intliq_id == 0 .and. intsol_id == 0) then
      call diagnostic (arg1=11,arg2=2,arg3=nomsub)
    end if
!
    if (intpl_id == 0) then
!.. localizzazione del pelo libero su tutte le celle in direzione i e j del dominio
      peloloc = max_negative_number
      intpl_id = 0
!.. ciclo sulle celle sopra e quella corrente compresa
      do k = 1,Grid%ncd(3)
        numcell = CellNumber(igridi,jgridi,k)
        if (numcell == 0) cycle
!.. ciclo su tutte le particelle della cella numcell
        do m = Icont(numcell),Icont(numcell+1)-1
          nnlocal = npartord(m)
          if (index(Med(pg(nnlocal)%imed)%tipo,"liquid") > 0) then
!.. ricerca interfaccia liquid/granular (particella liquid flu)
            if (pg(nnlocal)%Coord(3) > peloloc) then
              peloloc = pg(nnlocal)%Coord(3)
              intpl_id = nnlocal
            end if
          end if
        end do
      end do
    end if
!
    if (intliq_id == 0) intliq_id  = intsol_id
    interf_liq = pg(intliq_id)%coord(3)
!
    if (intsol_id == 0) intsol_id  = intliq_id
    interf_sol = pg(intsol_id)%coord(3)
!
    if (intpl_id == 0) intpl_id = intliq_id
    peloloc = pg(intpl_id)%coord(3)
    pl_imed = pg(intpl_id)%imed
!
!.. calcolo viscosita' apparente
!
    secinv = two * pg(npi)%secinv
!.. modifica invariante secondo rapportandolo alla viscosita' del mezzo dove è situato il pelo libero locale
!!terzaprova ok
    secinv = secinv * Med(imed)%mumx / (Med(pl_imed)%visc*med(pl_imed)%den0)
!
!
!!se x<=xtubo peloloc = peloloc
!!primaprova
!!     if(pg(npi)%Coord(1) >= 3.0) peloloc = pg(pl_id)%pres / (-Domain%grav(3) * med(pl_imed)%den0) + pg(pl_id)%Coord(3)
!
    preidro = -Domain%grav(3) * med(pl_imed)%den0 * (peloloc-pg(npi)%Coord(3))
!
! ( caso0 )
!    pretot  = -Domain%grav(3) * ( &
!              (med(pg(intsol_id)%imed)%den0 * (interf_sol - pg(npi)%Coord(3))) + &
!              (med(pl_imed)%den0 * (peloloc - interf_sol)) )
!!    preeff  = pretot - preidro
!    preeff  = -Domain%grav(3) * (med(pg(intsol_id)%imed)%den0 - med(pg(intliq_id)%imed)%den0) * &
!              (interf_sol - pg(npi)%Coord(3))
!
! ( caso1 )
!    appo1 = (interf_liq + interf_sol) * half
!    pretot  = -Domain%grav(3) * ( &
!              (med(pg(intsol_id)%imed)%den0 * (appo1 - pg(npi)%Coord(3))) + &
!              (med(pl_imed)%den0 * (peloloc - appo1)) )
!!    preeff  = pretot - preidro
!    preeff  = -Domain%grav(3) * (med(pg(intsol_id)%imed)%den0 - med(pg(intliq_id)%imed)%den0) * &
!              (appo1 - pg(npi)%Coord(3))
!
! ( caso2 )
!.. calcolo pressione
!    PartLiq_pres = pg(intliq_id)%pres
!    pretot  = PartLiq_pres  + (-Domain%grav(3) * med(pg(intsol_id)%imed)%den0 * (interf_liq - pg(npi)%Coord(3)) )
!!    preeff  = pretot - preidro
!    preeff  = -Domain%grav(3) * (med(pg(intsol_id)%imed)%den0 - med(pg(intliq_id)%imed)%den0) * &
!              (interf_liq - pg(npi)%Coord(3))
!
! ( caso3 )
!.. calcolo pressione
!    PartLiq_pres = pg(intliq_id)%pres
!    pretot  = PartLiq_pres + (-Domain%grav(3)) * ( &
!              med(pg(intsol_id)%imed)%den0 * (interf_sol + interf_liq - two*pg(npi)%Coord(3)) * half + &
!              med(pg(intliq_id)%imed)%den0 * (interf_liq - interf_sol) * half )
!!    preeff  = pretot - preidro
!    preeff  = -Domain%grav(3) * (med(pg(intsol_id)%imed)%den0 - med(pg(intliq_id)%imed)%den0) * &
!              ((interf_sol + interf_liq - two*pg(npi)%Coord(3)) * half)
!
! ( caso4 )
!.. calcolo pressione
     Velocity2 = pg(intliq_id)%vel(1)*pg(intliq_id)%vel(1) + pg(intliq_id)%vel(2)*pg(intliq_id)%vel(2) + &
                 pg(intliq_id)%vel(3)*pg(intliq_id)%vel(3)
     PartLiq_pres = (interf_liq - interf_sol) * med(pg(intliq_id)%imed)%den0 * (-Domain%grav(3)) + &
                    pg(intliq_id)%pres + Velocity2 * half * med(pg(intliq_id)%imed)%den0
     pretot  = PartLiq_pres + &
               med(pg(intsol_id)%imed)%den0 * (-Domain%grav(3)) * (interf_sol - pg(npi)%Coord(3))
     preeff  = pretot - preidro
!---------- provare-------------
!
    coeff1  = cos (Med(imed)%phi)
    coeff2  = sin (Med(imed)%phi)
!
    mu      = (Med(imed)%coes * coeff1 + preeff * coeff2) / (secinv + 0.00001d0)
    mumax   = Med(imed)%mumx
!
    flag = 0  ! vale 1 se esiste una particella 'liquid' nell'intorno della particella 'granular' corrente
!
!.. controllo per eventuale cambio dello stato alle particelle solide
!
    if (mu < mumax) then ! .and. it_corrente > Med(pg(npi)%imed)%NIterSol) then
!
!---------- provare-------------
      Nsp = nPartIntorno(npi) 
! ricerca fluide vicino per individuare particelle solide di strati vicino all'interfaccia fluido-solido
      appo1 = peloloc - interf_sol
      do j = 1,Nsp
        npartint = (npi-1)* NMAXPARTJ + j
        npj = PartIntorno(npartint)
!
!.. modifica test per includere anche le particelle granulari ma con stato flu
!!! secondaprova ok
        if (index(Med(pg(npj)%imed)%tipo,"liquid") > 0 .or. (index(Med(pg(npj)%imed)%tipo,"granular") > 0 .and. &
            pg(npi)%state == "flu")) then
!!!         if (index(Med(pg(npj)%imed)%tipo,"liquid") > 0) then
          if (appo1 >= Domain%dd) then
            flag = 1
            exit
          end if
        end if
      end do
!---------- provare-------------
!      appo1 = pg(intsol_id)%Coord(1) - pg(npi)%Coord(1)
!      appo2 = pg(intsol_id)%Coord(2) - pg(npi)%Coord(2)
!      appo3 = pg(intsol_id)%Coord(3) - pg(npi)%Coord(3)
!      ragtemp = Dsqrt((appo1*appo1) + (appo2*appo2) + (appo3*appo3))
!      if (ragtemp <= two*Domain%h) then
!        if ((peloloc - interf_sol) >= Domain%dd) then !<<<superfluo????
!          flag = 1
!        end if
!      end if
!---------- provare-------------
!
! ricerca completata                         V crea l'instabilita' iniziale per inizio moto
      if (flag == 0 .and. it_corrente > Med(imed)%NIterSol) then
        if (pg(npi)%CloseBcOut == 0) then
          pg(npi)%state = "sol"
          pg(npi)%vel = zero
          pg(npi)%var = zero
!! inutile perche' assegnata in setparticle e poi mantenuta costante
!!!!           pg(npi)%mu = mumax
! viene assegnata una densita' coerente con la pressione idrostatica
          if (.not. diffusione) pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita*Med(imed)%celerita))
        end if
!
      else if (flag == 1) then
!! inutile perche' assegnata in setparticle e poi mantenuta costante
!!         if (pg(npi)%state /= "flu") pg(npi)%mu = med(imed)%visc * med(imed)%den0
        pg(npi)%state = "flu"
      end if
!
    else
      if (pg(npi)%CloseBcOut == 0) then
        pg(npi)%state = "sol"
        pg(npi)%vel = zero
        pg(npi)%var = zero
!! inutile perche' assegnata in setparticle e poi mantenuta costante
!!!!         pg(npi)%mu = mumax
! viene assegnata una densita' coerente con la pressione idrostatica
        if (.not. diffusione) pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita*Med(imed)%celerita))
      end if
!
    end if
!
  end do
!
!$omp end parallel do
!
return
end subroutine MohrC
!---split

