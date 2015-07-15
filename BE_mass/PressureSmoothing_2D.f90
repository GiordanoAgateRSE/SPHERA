!cfile PressureSmoothing_2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : PressureSmoothing_2D
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!AA501b
! 03  Amicarelli-Agate  13nov12        Body dynamics
!AA601
! 04  Amicarelli        26Jan15        DBSPH-input (AA601). Pressure smoothing considering DBSPH surface elements.
!
!************************************************************************************
! Module purpose : Module to compute pressure smoothig
!
! Calling routine: Loop_Irre_2D
!
!AA501b modified
! Called routines: body_pressure_postpro,body_to_smoothing_vel
!
!************************************************************************************
!
subroutine PressureSmoothing_2D
!Regolarizza la pressione, e quindi la densita' delle particelle
!Attenzione: non si adatta parfettamente al caso di mezzi con diversa densita' a contatto
!In questo caso sarebbe meglio mediare sulle particelle vicine solo dello stesso mezzo
!in modo da rispettare la discontinuita' di densita'
!
!.. assign modules
  use FILES_ENTITIES
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Local Scalars ..
  integer(4)       :: Ncbs, IntNcbs
  integer(4)       :: npi, npj, mati, nsp, icbs, iside, sidestr, izonelocal, ibdt, ibdp  !!!!!!!!matj, 
!20081028 da riattivare  integer(4)       :: nrows
  integer(4)       :: j, i, ii, npartint
  double precision :: ro0i, p0i, pi, sompW, pesoj, TetaP1 
  double precision :: VIntWdV_FT, VIntWdV_SO, VIntWdV_OSB, VIntWdV_OSP, press_so, press_osb, vin
  double precision :: Appunity,smoothpi,IntWdV,IntEnerg

  character(4)     :: strtype
!  Logical          :: PressLocked
!
!.. Local Arrays ..
  integer(4),      dimension(1:PLANEDIM) :: acix

!20081028 da riattivare  integer(4),      dimension(1:NUMCOLS_BIT)   :: Rows
  double precision,dimension(1:2)        :: IntDpWdV
  double precision,dimension(1:PLANEDIM) :: IntLocXY

!20081028 da riattivare  double precision,dimension(1:NUMCOLS_BIT)      :: Integral
  double precision,dimension(1:PLANEDIM) :: sss, nnn, massforce
!
  character(1), parameter :: SmoothVersion = "b"          ! = "a"  (SPHERA), "b" (FreeSurf), "c"
!

!AA501b start
  double precision,dimension(:),allocatable :: sompW_vec,AppUnity_vec    

  if (n_bodies > 0) then  
!Allocations
     allocate(sompW_vec(nag))
     allocate(AppUnity_vec(nag))
!Initializations
     sompW_vec = zero
     AppUnity_vec = zero
  endif
!AA501b end

!.. Executable Statements ..
!
!!!!  allocate (PartPress(1:nag), PartDens(1:nag))
!
!!!  call cpu_time(cpu_loop1a)
!
  acix(1) = 1
  acix(2) = 3
!

!AA501b
! Body particle contributions to pressure smoothing
  if (n_bodies > 0) then
!     call start_and_stop(3,14)
     call start_and_stop(2,19)
     call body_to_smoothing_pres(sompW_vec,AppUnity_vec)
     call start_and_stop(3,19)
!     call start_and_stop(2,14)
  endif
!AA501b end 

!matj,sidelen,
!AA501b modified
!$omp parallel do default(none) &
!$omp private(npi,ii,Appunity,TetaP1,Ncbs,IntNcbs,ibdt,icbs,ibdp,iside,sidestr,Nsp) &
!$omp private(mati,ro0i,p0i,pi,SompW,j,npartint,npj,pesoj,IntEnerg) &
!$omp private(VIntWdV_FT,VIntWdV_SO,VIntWdV_OSB,VIntWdV_OSP,press_so,press_osb) &
!$omp private(IntLocXY,strtype,sss,nnn,massforce,IntWdV,IntDpWdV,izonelocal,vin,smoothpi) &
!$omp shared(nag,Pg,Med,Tratto,Partz,Domain,nPartIntorno,PartIntorno,NMAXPARTJ,PartKernel,kernel_fw) &
!$omp shared(BoundarySide,BoundaryDataPointer,BoundaryDataTab,acix,dt,indarrayFlu,Array_Flu,esplosione,sompW_vec,n_bodies,AppUnity_vec) 
!
!!!!!  do npi = 1,Nag
!!!!!!
!!!!!    pg(npi)%vpres = pg(npi)%pres
!!!!!    pg(npi)%vden  = pg(npi)%dens
!!!!!    if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" .or. pg(npi)%state == "sol") cycle
!!!!!!    if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!$$$$$
  do ii = 1,indarrayFlu
    npi = Array_Flu(ii)
   
!$$$$$
!
   
!.. esclusione particelle vicino alla faccia con condizione 'flow', 'velo' e 'sour'
    if (pg(npi)%koddens == 0) then 
!
      Nsp = nPartIntorno(npi)
!
      if (Nsp > 0) then
!
        Appunity = zero
!
        mati = pg(npi)%imed
        ro0i = Med(mati)%den0
        p0i = Domain%prif
        pi = pg(npi)%pres
        sompW = zero
!
!.. valore mediato delle differenze di pressione rispetto alla pressione della particella
!
        do j = 1,Nsp
!
          npartint = (npi-1)* NMAXPARTJ + j
          npj = PartIntorno(npartint)
!          
!$$$
!          matj = pg(npj)%imed                   !!!!!!!!!!!!!!
!.. prova da fare   smoothing sullo stesso tipo di mezzo della particella corrente
!!          if (Med(mati)%tipo == Med(matj)%tipo) then                   !!!!!!!!!!!!!!
!
            pesoj = pg(npj)%mass * PartKernel(4,npartint) / pg(npj)%dens
            Appunity = Appunity + pesoj
            sompW = sompW + (pg(npj)%pres - pi) * pesoj
!!          end if                                                        !!!!!!!!!!!!!!
!$$$
!
        end do
!AA501b
        if (n_bodies > 0) then
           sompW = sompW + sompW_vec(npi)
           AppUnity = AppUnity + AppUnity_vec(npi)
        endif
!AA406
        if (Domain%tipo == "bsph") then
!AA601
           TetaP1 = Domain%TetaP * Med(pg(npi)%imed)%Celerita * dt / Domain%h
            pg(npi)%vpres = pi + TetaP1 * sompW / AppUnity
        else
!
        VIntWdV_FT = zero
        VIntWdV_SO = zero
        VIntWdV_OSB = zero
        VIntWdV_OSP = zero
        press_so = pi
        press_osb = pi
!
        Ncbs = BoundaryDataPointer(1,npi)
        IntNcbs = BoundaryDataPointer(2,npi)
        ibdt = BoundaryDataPointer(3,npi)
!
        do icbs = 1,IntNcbs
!
          ibdp = ibdt + icbs - 1
          IntLocXY(1:PLANEDIM) = BoundaryDataTab(ibdp)%LocXYZ(1:PLANEDIM)
          iside = BoundaryDataTab(ibdp)%CloBoNum
!
          sidestr = BoundarySide(iside)%stretch
          strtype = Tratto(sidestr)%tipo
!          sidelen = BoundarySide(iside)%Length
          do i = 1, PLANEDIM
            sss(i) = BoundarySide(iside)%T(acix(i),1)
            nnn(i) = BoundarySide(iside)%T(acix(i),3)
          end do
          massforce(1) = sss(1) * Domain%grav(acix(1)) + sss(2) * Domain%grav(acix(2))
          massforce(2) = nnn(1) * Domain%grav(acix(1)) + nnn(2) * Domain%grav(acix(2))

!20081028 da riattivare         if (IntNcbs == 2) then
!.. Calcolo numerico dell'integrale IntWdV
!
          IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(2)
!          IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(3)
          IntDpWdV(1:2) = BoundaryDataTab(ibdp)%BoundaryIntegral(4:5)
!
!20081028            call ComputeVolumeIntegral_WdV2D (icbs, IntNcbs, Intboside, IntLocXY, BoundarySide, xpmin, xpmax, interlen, IntWdV)
                                                          
!20081028 da riattivare          else if (IntNcbs == 1) then
!.. Interpolazione da tabella dell'integrale IntWdV,
                    
!20081028 da riattivare            deltai = ypi / Domain%h
!20081028 da riattivare            Nrows = 1
!20081028 da riattivare            Rows(1) = 2
                    
!20081028 da riattivare            call InterpolateBoundaryIntegrals2D (deltai, Nrows, Rows, Integral)
                    
!20081028 da riattivare            IntWdV = Integral(1)

!20081028 da riattivare          end if
!          end if
!
          if (strtype == "fixe" .OR. strtype == "tapi") then

            VIntWdV_FT = VIntWdV_FT + IntWdV
            sompW = sompW + ro0i * (massforce(1) * IntDpWdV(1) + massforce(2) * IntDpWdV(2))
!
          else if (strtype == "sour") then
!
            VIntWdV_SO = VIntWdV_SO + IntWdV
            izonelocal = pg(npi)%izona
            if (partz(izonelocal)%pressure == "pa") then
              press_so = partz(izonelocal)%valp
            else if (partz(izonelocal)%pressure == "qp" .or. partz(izonelocal)%pressure == "pl") then
              press_so = ro0i * Domain%grav(3) * (Pg(npi)%coord(3) - partz(izonelocal)%valp)
            end if
!
          else if (strtype == "crit") then
!
            VIntWdV_OSP = VIntWdV_OSP + IntWdV
!
!.. Implicitamente press_osp = pi
          else if (strtype == "leve") then
!
            VIntWdV_OSB = VIntWdV_OSB + IntWdV
            izonelocal = pg(npi)%izona
            vin = BoundarySide(iside)%T(acix(1), 3) * pg(npi)%Vel(acix(1)) + &
                  BoundarySide(iside)%T(acix(2), 3) * pg(npi)%Vel(acix(2))
            if (vin < zero) then
              press_osb = ro0i * Domain%grav(3) * (pg(npi)%coord(3) - partz(izonelocal)%valp)
            else
              press_osb = pi
            end if
!
          end if
!
        end do
!
        if (esplosione) then
!.. con Csound al posto di Celerita e' circa uguale
          TetaP1 = Domain%TetaP * pg(npi)%Csound * dt / Domain%h
        else
!.. calcolo TetaP adeguato al passo temporale
          TetaP1 = Domain%TetaP * Med(pg(npi)%imed)%Celerita * dt / Domain%h
        end if
!
        AppUnity = AppUnity + VIntWdV_FT + VIntWdV_SO + VIntWdV_OSP + VIntWdV_OSB
!
        select case (SmoothVersion)
          case ("a")
              smoothpi = pi + TetaP1 * (sompW + (press_so - pi) * VIntWdV_SO + (press_osb - pi) * VIntWdV_OSB) / AppUnity
          case ("b")
              smoothpi = pi + TetaP1 * (sompW + (press_so - pi) * VIntWdV_SO + (press_osb - pi) * VIntWdV_OSB + &
                         (p0i - pi) * (one - AppUnity))
          case ("c")
              smoothpi = pi + TetaP1 * (sompW + (press_so - pi) * VIntWdV_SO + (press_osb - pi) * VIntWdV_OSB)
          case default
              smoothpi = Pg(npi)%pres       !no smoothing
        end select
!
        pg(npi)%vpres = smoothpi
!
!AA406
      endif !Domain%tipo
!
      end if  !(nsp > 0)
!
    end if  !(koddens == 0)
!
  end do
!
!$omp end parallel do
!
!!!    call cpu_time(cpu_loop1b)
!!!    call cpu_time(cpu_loop2a)
!
!.. I nuovi valori di densita' e pressione sono stati "parcheggiati" nei vettori PartDens() e PartPress()
!.. Nel ciclo che segue tali valori vengono copiati nella riga k
!
!$omp parallel do default(none) &
!$omp private(npi,ii,IntEnerg) &
!$omp shared(nag,Pg,Domain,Med,indarrayFlu,Array_Flu,esplosione)
!
!!!!  do npi = 1,Nag
!$$$$$
  do ii = 1,indarrayFlu
    npi = Array_Flu(ii)
!$$$$$
!
!.. esclusione particelle vicino alla faccia con condizione 'flow', 'velo' e 'sour'
    if (pg(npi)%koddens == 0) then 
!
      pg(npi)%pres = pg(npi)%vpres
!
      if (esplosione) then
!...................................... 2011 mar 11
!..VERIFICARE??
!.. modify for Specific Internal Energy
        IntEnerg = pg(npi)%pres / ((Med(pg(npi)%imed)%gamma - one) * pg(npi)%Dens)
        pg(npi)%IntEn = half * (pg(npi)%IntEn + IntEnerg)
        pg(npi)%dens = pg(npi)%pres / ((Med(pg(npi)%imed)%gamma - one) * pg(npi)%IntEn)
!        IntEnerg = pg(npi)%pres / (Med(pg(npi)%imed)%gamma - one)
!        pg(npi)%IntEn = half * (pg(npi)%IntEn + IntEnerg / pg(npi)%Dens)
!        pg(npi)%dens = IntEnerg / pg(npi)%IntEn
!......................................
      else
        pg(npi)%dens = Med(pg(npi)%imed)%den0 * (one + (pg(npi)%vpres - Domain%Prif) / Med(pg(npi)%imed)%eps)
      end if
!
    end if
  end do
!
!$omp end parallel do

!AA501b start
!Deallocations
  if (n_bodies > 0) then
     deallocate(sompW_vec)
     deallocate(AppUnity_vec)
  endif
!AA501b end

!
!!!  call cpu_time(cpu_loop2b)
!
return
end subroutine PressureSmoothing_2D
!---split

