!cfile PressureSmoothing_3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : PressureSmoothing_3D
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
! Calling routine: Loop_Irre_3D
!
!AA501b modified
! Called routines: body_pressure_postpro,body_to_smoothing_vel
!
!************************************************************************************
!
subroutine PressureSmoothing_3D
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
!.. Local Parameters ..
  double precision, parameter :: MinTotUnit = 0.95d0
  character(10),    parameter :: SmoothingNormalisation = "complete  "          ! = "incomplete"
!
!.. Local Scalars ..
  integer(4)       :: npi, npj, nsp, icbf, Ncbf, ibdt, ibdp    !!!!!!!!!matj, mati, 
  integer(4)       :: ii, j, npartint
  double precision :: p0i, pi, sompW, pesoj, DiffP, TetaP1 
  double precision :: Appunity, smoothpi
  double precision :: IntWdV
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

!!!  call cpu_time(cpu_loop1a)
!
!matj,mati,
!
!AA601 sub
!$omp parallel do default(none) &
!$omp private(npi,ii,Nsp,DiffP,p0i,pi,sompW,Appunity) &
!$omp private(j,npartint,npj,pesoj,Ncbf,ibdt,icbf,ibdp,intWdV) &
!$omp shared(nag,Pg,Med,Domain,nPartIntorno,PartIntorno,NMAXPARTJ,PartKernel) &
!$omp shared(BoundaryDataPointer,BoundaryDataTab,indarrayFlu,Array_Flu,TetaP1,sompW_vec,AppUnity_vec,n_bodies,dt)
!
!!!!!!  do npi = 1,Nag
!!!!!!!
!!!!!!    pg(npi)%vpres = zero
!!!!!!!
!!!!!!    if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" .or. pg(npi)%state == "sol") cycle
!$$$$$
  do ii = 1,indarrayFlu
    npi = Array_Flu(ii)
!$$$$$
!
!.. esclusione particelle vicino alla faccia con condizione 'flow', 'velo' e 'sour'
    if (pg(npi)%koddens == 0) then 
!
    DiffP = zero
!
    Nsp = nPartIntorno(npi)
!
    if (Nsp > 0) then
!
!      mati = pg(npi)%imed
      p0i = Domain%prif
      pi = pg(npi)%pres
!
      sompW = zero
      Appunity = zero
!
!.. valore mediato delle differenze di pressione rispetto alla pressione della particella
!
      do j = 1,Nsp
!
        npartint = (npi-1)* NMAXPARTJ + j
        npj = PartIntorno(npartint)
        
!$$$
!        matj = pg(npj)%imed
!.. prova da fare  smoothing sullo stesso tipo di mezzo della particella corrente
!!        if (Med(mati)%tipo == Med(matj)%tipo) then            !!!!!!!!!!!!
          pesoj = pg(npj)%mass * PartKernel(4,npartint) / pg(npj)%dens
          Appunity = Appunity + pesoj
          sompW = sompW + (pg(npj)%pres - pi) * pesoj
!!        end if                                                 !!!!!!!!!!!!
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
      Ncbf = BoundaryDataPointer(1,npi)
      ibdt = BoundaryDataPointer(3,npi)
!
      if (Ncbf > 0) then
!
        do icbf = 1, Ncbf
          ibdp = ibdt + icbf - 1
          if (BoundaryDataTab(ibdp)%LocXYZ(3) > zero) then
            IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(2)
            AppUnity = AppUnity + IntWdV
          end if
        end do
!
      end if
!
      if (SmoothingNormalisation == "complete  ") then
        if (AppUnity < MinTotUnit) then
!AA504 comment            
          DiffP = sompW + (p0i - pi) * (one - AppUnity)   ! in case of non null reference pressure
        else
          DiffP = sompW / AppUnity
        end if
      else if (SmoothingNormalisation == "incomplete") then
        DiffP = sompW / AppUnity
      end if
!
!
!AA406
    endif !Domain%tipo
!
    end if  !(nsp > 0)
!
!AA406 sub
     if (Domain%tipo == "semi") pg(npi)%vpres = DiffP
!
    end if  !(koddens == 0)
!
  end do
!
!$omp end parallel do
!
!!!  call cpu_time(cpu_loop1b)
!!!  call cpu_time(cpu_loop2a)
!
!$omp parallel do default(none) &
!$omp private(npi,ii,smoothpi,TetaP1) &
!$omp shared(nag,pg,Med,Domain,dt,indarrayFlu,Array_Flu,esplosione)
!
!.. Ciclo di smoothing della pressione e della densita'

!!!!!!  do npi = 1, nag
!!!!!!!
!!!!!!    if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" .or. pg(npi)%state == "sol") cycle
!$$$$$
  do ii = 1,indarrayFlu
    npi = Array_Flu(ii)
!$$$$$
!.. esclusione particelle vicino alla faccia con condizione 'flow', 'velo' e 'sour'
    if (pg(npi)%koddens == 0) then 
!
!.. calcolo TetaP adeguato al passo temporale
      if (esplosione) then
        TetaP1 = Domain%TetaP * pg(npi)%Csound * dt / Domain%h
      else
        TetaP1 = Domain%TetaP * Med(pg(npi)%imed)%Celerita * dt / Domain%h
      end if
      smoothpi = pg(npi)%pres + TetaP1 * pg(npi)%vpres
      pg(npi)%Pres = smoothpi
      pg(npi)%dens = Med(pg(npi)%imed)%den0 * (one + (smoothpi - Domain%Prif) / Med(pg(npi)%imed)%eps)
!
    end if
!
  end do
!
!$omp end parallel do

!
!!!  call cpu_time(cpu_loop2b)
!

!AA501b start
!Deallocations
  if (n_bodies > 0) then
     deallocate(sompW_vec)
     deallocate(AppUnity_vec)
  endif
!AA501b end

return
end subroutine PressureSmoothing_3D
!---split

