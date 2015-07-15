!cfile Loop_Irre_2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : Loop_Irre_2D
!
! Last updating : April 08, 2014
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Amicarelli/Agate  14Jun11        Algorithm adaptation for RK1 scheme (not staggered)
! 04  Amicarelli/Agate  30Jun11        Algorithm adaptation for RK2 scheme 
! 05  Amicarelli/Agate  30Nov11        BSPH approach
!
!AA501b
! 06  Amicarelli/Agate  13nov12        Body dynamics 
!AA504
! 07  Amicarelli        08Apr14        (v5.04) Minor modifications related to: elapsed time estimation, erosion criterion, density thresholds, mixture viscosity,  
!                                      output for granular flows, hydrographs. 
! 08  Amicarelli        26Jan15        DBSPH-input (AA601). Several modifications (e.g.: calling DBSPH post-processing, removing fictitious air reservoirs,...)
!
!************************************************************************************
! Module purpose : Module Loop
!
! Calling routine: Gest_Trans
!
! Called routines: s_secon2
!                  memo_results
!                  inidt2
!                  NormFix
!                  subCalcPreIdro
!                  rundt2
!                  vellaw
!                  calcpre
!                  inter_SmoothPres
!                  inter_EqMoto
!                  AddBoundaryContributions_to_ME2D
!                  AddSuppBoundaryContribution_to_ME
!                  AddBoundaryExternalReaction
!                  contrmach
!                  inter_SmoothVelo_2D
!                  CancelOutgoneParticles_2D
!                  GenerateSourceParticles_2D
!                  OrdGrid1
!                  inter_EqCont_2D
!                  AddBoundaryContribution_to_CE2D
!                  aggdens
!                  viscapp
!                  CalcVarp
!                  CalcVarLength
!                  memo_ctl
!                  calc_pelo
!                  memo_restart
!                  writime2

!AA501b
!                  time_integration_body_dynamics
!AA504
!                  Shields
!
!************************************************************************************
!
  subroutine Loop_Irre_2D  
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
!include "omp_lib.h"
!
!integer(4), parameter :: iterazione = 2
!
!.. Local Scalars ..
!  integer(4)       :: id,nthreads
  integer(4)       :: Ncbs, IntNcbs, i, ii, num_out, npi
!
  integer(4)       :: it, it_print, it_memo, it_rest, ir
  integer(4)       :: OpCountot,SpCountot,EpCountot,EpOrdGridtot
!AA504
  integer(4)       :: ncel,aux,igridi,jgridi,kgridi,machine_Julian_day,machine_hour,machine_minute,machine_second 
!AA504 sub  
  double precision :: BCrodivV, dtvel, DtPreviousStep,TetaV1,xmax,ymax,appo1,appo2,appo3,pretot
  character(len=lencard)  :: nomsub = "Loop_Irre_2D"
!
!.. local Arrays ..
  double precision,dimension(1:SPACEDIM) :: tpres,tdiss,tvisc
  double precision,dimension(2)          :: ti, tf, ti_iter, tf_iter, tot_iter
  double precision,dimension(1:SPACEDIM) :: BoundReaction
!AA504 rm line

!.. External Routines ..
  integer(4),external :: ParticleCellNumber, CellIndices, CellNumber

!.. Executable Statements ..
!AA504 rm line
!  id = omp_get_thread_num()
!  if ( id == 0 ) then
!    nthreads = omp_get_num_threads()
!    write (*,*) 'There are ', nthreads, ' threads'
!  end if
!
!.. initializations
!
!................................... 2011 mar 08
!.. initialization of Specific Internal Energy
  if (esplosione) then
    do npi = 1,nag
!      if (pg(npi)%IntEn > 10000.0) then ! modificare per riconoscere il mezzo gas esplosione
      if (index(Med(pg(npi)%imed)%tipo,"gas") > 0) then
        pg(npi)%pres = (Med(pg(npi)%imed)%gamma - one) * pg(npi)%IntEn * pg(npi)%dens
        pg(npi)%Csound = Dsqrt(Med(pg(npi)%imed)%gamma * (Med(pg(npi)%imed)%gamma - one) * pg(npi)%IntEn)
      else
        pg(npi)%Csound = Med(pg(npi)%imed)%Celerita
        pg(npi)%IntEn = pg(npi)%pres / ((Med(pg(npi)%imed)%gamma - one) * pg(npi)%dens)
      end if
      pg(npi)%state = 'flu'
    end do
  end if
!...................................
!
  SpCount  = 0
  OpCount  = 0
  EpCount  = 0
  EpOrdGrid  = 0
  tot_iter = zero
  num_out  = 0
  SourceSide   = 0
  it_eff   = it_start
  it_print = it_start
  it_memo  = it_start
  it_rest  = it_start
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ DA VERIFICARE $$$$$$$$$$$$$$$$$$$$$$$$$$$
  call s_secon2 ( ti )
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!.. initialization of support array to count the particles not 'sol' in the case with 'granular'
!
!  indarraySol = 0
  indarrayFlu = 0
  do npi = 1,nag
    if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!    if (index(Med(pg(npi)%imed)%tipo,"granular") > 0) then
!      indarraySol = indarraySol + 1
!      Array_Sol(indarraySol) = npi
!    end if
    if (pg(npi)%state == "flu") then
      indarrayFlu = indarrayFlu + 1
!.. check on max dimension for the array controllo di non superare le dimensioni ipotizzate per il vettore e nel caso li ridefinisco opportunamente
      if (indarrayFlu > PARTICLEBUFFER) then
        call diagnostic (arg1=9,arg2=1,arg3=nomsub)
      end if
      Array_Flu(indarrayFlu) = npi
    end if
  end do
!
!.. Prepara condizioni per eventuali boundary di tipo sorgente
!.. Genera e immette nel campo nuove particelle di sorgente per simulare una portata fluida entrante
!
  call PreSourceParticles_2D
!
!Initializing the stage
!
  if (Domain%time_split==0) Domain%time_stage = 1
!
!.. evaluates the distances between the one particle and all the nearest ones, storing them 
!.. in the general storage array
!.. Computation and storage of the interacting particle index, relative position, distance and kernel gradient (Sphera_Tools.f90)
!
  call start_and_stop(2,10)
!AA601
  if ((it_corrente==it_start).and.(Domain%tipo=="bsph")) it_corrente = -2  
  call CalcVarLength
!AA601
  if ((it_corrente==-2).and.(Domain%tipo=="bsph")) it_corrente = it_start  
  call start_and_stop(3,10)
!AA601 start
! Removing fictitious reservoirs used for DBSPH initialization
  call start_and_stop(2,9)
  if ((it_corrente==it_start).and.(Domain%tipo=="bsph")) then
!$omp parallel do default(none) shared(nag,pg,OpCount) private(npi) 
     do npi=1,nag
        if (pg(npi)%imed/=pg(1)%imed) then
           OpCount(pg(npi)%imed) = OpCount(pg(npi)%imed) + 1    
           pg(npi)%cella = -1
        endif
     end do
!$omp end parallel do
     call OrdGrid1 (nout)
     indarrayFlu = nag
     call start_and_stop(3,9)
     call start_and_stop(2,10)
!This fictitious value avoid CalcVarlength computing Gamma, sigma and density for this very call      
     it_corrente = -1
     call CalcVarLength
! The correct value of it_corrente is restored
     it_corrente = it_start
     call start_and_stop(3,10)
  endif
!AA601 end 
!AA406 start
! subroutine for wall BC treatment (BSPH)
! density and pressure updates for wall elements: MUSCL + partial Riemann solver + state equation 
       call start_and_stop(2,18)
       if ((Domain%tipo == "bsph") .and. (nag > 0) .and. (DBSPH%n_w > 0)) then
          call Gradients_to_MUSCL   
          call BC_wall_elements
       endif
       call start_and_stop(3,18)
!AA406 end

!AA501b start
! pressure initialization for body particles
 if (n_bodies > 0) then
     call start_and_stop(2,19)
     call body_pressure_mirror
     call body_pressure_postpro
     call start_and_stop(3,19)
  endif
!AA501b end

!
!.. evaluates the close boundaries and the integrals
!.. for the current particle in every loop, storing them in the general storage array
!.. computation and storage of the intersections between the kernel support and the frontier and the corresponding boundary integrals 
!.. &(semi-analytic approach)
!
!AA406 sub start
  if (Domain%tipo == "semi") then
     call start_and_stop(2,11)
     call ComputeBoundaryDataTab
     call start_and_stop(3,11)
  endif
!AA406 sub end
!
!.. evaluates the properties that must be attributed to the fixed particles
!
  if ( Domain%NormFix )  call NormFix
!
  if ( nscr > 0 ) write (nscr,"(a,1x,a)") " Running case:",trim(nomecas2)
!
!.. evaluates the new pressure field
!
!§  call calcpre  
!
!§  smoothing sulla VF per ridurre i gradienti all'interfaccia
!
!  if (diffusione) then
!    do ii = 1,5
!      do npi = 1,nag
!        if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!        call inter_SmoothVF (npi,appo1)
!        pg(npi)%VolFra = pg(npi)%VolFra + appo1 / pg(npi)%uni
!      end do
!    end do
!  end if
!
!§ 
!
!.. initializes the output files
!
  if ( nout > 0) then
    it_print = it_eff
    call print_results ( it_eff, it_print, 'inizio' )
  end if
  if ( nres > 0 ) then
    it_memo = it_eff
    it_rest = it_eff
    call memo_results ( it_eff, it_memo, it_rest, dtvel, 'inizio' )
  end if
  if ( vtkconv ) then
    call result_converter ('inizio')
  end if
!
!.. evaluates the initial value of the time step 
!
  if ( it_start == 0 ) call inidt2 
!
!-----------------------------------------------
!.. starts the time loop
!-----------------------------------------------
!
!AA504 start
! Computing the elapsed time before the iterations
 if (exetype == "linux") then
    if (Domain%tmax>0.d0) then
       call system("date +%j%H%M%S > date_pre_iterations.txt")
       open (unit_time_elapsed,file='date_pre_iterations.txt',status="unknown",form="formatted")
       read (unit_time_elapsed,'(i3,i2,i2,i2)') machine_Julian_day,machine_hour,machine_minute,machine_second
       close (unit_time_elapsed)
       Domain%t_pre_iter = machine_Julian_day*24*60*60+machine_hour*60*60+machine_minute*60+machine_second
    endif
 endif
!AA504 end

  it = it_start
!
  ITERATION_LOOP: do while ( it < Domain%itmax )
!
!.. set the iteration counter
!
    it = it + 1
    it_corrente = it
    it_eff = it
!
!.. stores the old time step value and evaluates the new value of time step duration and the total time value
!
    DtPreviousStep = dt
!
!.. CFL condition
    if (nag > 0) call rundt2     !$$$$$$$$$$ VERIFICARE IL SENSO DEL TEST $$$$$$$$$$$$$ se nag == 0 fermare esecuzione 
!
    tempo = tempo + dt
!
    if ( nscr > 0 ) write (nscr,"(a,i8,a,g13.6,a,g12.4,a,i8,a,i5)") " it= ",it,"   time= ",tempo,"  dt= ",dt, &
                                "    npart= ",nag,"   out= ",num_out
!
!AA402
1100 continue
!
!AA401 start
    if (Domain%time_split == 0) then
!.. non staggered time integration scheme
!
!.. kernel computations and boundary contributions
!.. Computation and storage of the neighbouring indices, distances and kernel gradients
      call start_and_stop(2,10)
!AA402 sub start
!
!AA406 sub
       if ((Domain%tipo == "semi").and.(Domain%time_stage < 2)) call CalcVarLength

!
!AA402 sub end
      call start_and_stop(3,10)
      
!
!.. evaluates the close boundaries and the integrals
!.. for the current particle in every loop, storing them in the general storage array
!.. Computation and storage of the boundary integrals
!
!AA406
      if (Domain%tipo == "semi") then
!
         call start_and_stop(2,11)
         call ComputeBoundaryDataTab
         call start_and_stop(3,11)
!
!AA406
      endif
!
    end if
!
!AA401 end
!
!AA406 
    if (Domain%tipo == "semi") then
        
!.. verifies if particles belonging to a rigid block and moving with a fixed velocity are stopped
!
    do ir = 1,NPartZone !!! TRASFORMARE il DO sulle particelle per eliminare le zone
!
      if ( Partz(ir)%move /= "fix" ) cycle      
!
!.. applies the movement with a law of velocity (npointv points of velocity assigned in time)
!
      if ( Partz(ir)%npointv > 1 ) then
!
        call vellaw(Partz(ir)%vlaw,Partz(ir)%vel,Partz(ir)%npointv)
!
        write(nout,"(f12.4,a,i2,a,3f14.9)") tempo,"  zona",ir,"  vel.",Partz(ir)%vel
!
        do npi = Partz(ir)%limit(1),Partz(ir)%limit(2) !!! TRASFORMARE il DO sulle particelle per eliminare le zone
          if ( pg(npi)%cella == 0 ) cycle
          pg(npi)%var(:) = Partz(ir)%vel(:)
          if ( tempo >= pg(npi)%tstop ) then
            pg(npi)%vstart(:) = zero
            pg(npi)%vel(:)    = zero
          else
            pg(npi)%vstart(:) = Partz(ir)%vel(:)
            pg(npi)%vel(:) = Partz(ir)%vel(:)
          end if
        end do
!
!.. applies the movement with a fixed value of velocity (npointv = 1)
!
      else if ( Partz(ir)%npointv == 1 ) then
!
        do npi = Partz(ir)%limit(1),Partz(ir)%limit(2) !!! TRASFORMARE il DO sulle particelle per eliminare le zone
          if ( pg(npi)%cella == 0 ) cycle
          pg(npi)%var(:) = Partz(ir)%vel(:)
          if ( tempo >= pg(npi)%tstop ) then
            pg(npi)%vstart(:) = zero
            pg(npi)%vel(:)    = zero
          else
            pg(npi)%vel(:) = Partz(ir)%vel(:)
          end if
        end do
!
      end if

    end do
!
!AA406 
    endif
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ DA VERIFICARE $$$$$$$$$$$$$$$$$$$$$$$$$$$
    call s_secon2 ( ti_iter )
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!AA403 start
!=== Erosion test (time_split=0) ################
    if ((Domain%time_split == 0) .and. (Domain%time_stage == 1)) then 
!.. erosion criterium + continuity equation right hand side + deformation tensor norm 
      call start_and_stop(2,12)
!
      if (erosione .and. .not. esplosione) then  ! erosione
!
!.. Calcolo stato ('flu' o 'sol') delle particelle di tipo 'granulare'
!.. test modello erosione
!
!AA504 sub start
        if (Granular_flows_options%ID_erosion_criterion>1) then
! Erosion criterion
          if (index(modelloerosione,"shields") > 0) then
!.. Shields criterion
!$omp parallel do default(none) shared(pg,nag) private(npi) 
            do npi=1,nag
              call Shields(npi) 
            end do
!$omp end parallel do
          elseif (index(modelloerosione,"mohr") > 0) then
! Compute second invariant of the rate strain tensor and density derivatives
            call inter_EqCont_2D 
! Mohr-Coulomb criterion
            call MohrC
          endif
        endif
!AA504 sub end
!
!.. aggiornamento vettore appoggio per conteggio particelle non 'sol' caso con 'granular'
!
        indarrayFlu = 0
        do npi = 1,nag
          if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
          if (pg(npi)%state == "flu") then
            indarrayFlu = indarrayFlu + 1
!.. controllo di non superare le dimensioni ipotizzate per il vettore e nel caso li ridefinisco opportunamente
            if (indarrayFlu > PARTICLEBUFFER) then
              call diagnostic (arg1=9,arg2=2,arg3=nomsub)
            end if
            Array_Flu(indarrayFlu) = npi
          endif
        end do
! 
      else  ! erosione
!.. non sto eseguendo un calcolo con un modello di erosione 
!.. devo aggiornare vettore appoggio per conteggio particelle 'flu' a tutte nag
        indarrayFlu = 0
        do npi = 1,nag
          if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
          indarrayFlu = indarrayFlu + 1
!.. controllo di non superare le dimensioni ipotizzate per il vettore e nel caso li ridefinisco opportunamente
          if (indarrayFlu > PARTICLEBUFFER) then
            call diagnostic (arg1=9,arg2=2,arg3=nomsub)
          end if
          Array_Flu(indarrayFlu) = npi
        end do
!
!.. le particelle escluse eventuali con cella=0 (???) e non 'std'
      endif  ! erosione
!
    endif 
!AA403 time_split
!
!=== EQUAZIONE DEL MOTO ################
!
!..
!.. solve the equations for the motion field (velocity components)
!.. and time increment of Specific Internal Energy
!
    call start_and_stop(2,6)
!
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!if (it == iterazione) open (unit=99, file="File_eqmoto_v40x.txt", form="formatted", status="unknown")
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!
!$omp parallel do default(none) &
!$omp private(npi,ii,tpres,tdiss,tvisc,Ncbs,IntNcbs,BoundReaction) &
!$omp shared(nag,Pg,Domain,BoundaryDataPointer,indarrayFlu,Array_Flu,it)
!
!.. loops on all the particles
!
!!!!!!    do npi = 1,nag
!!!!!!      if ( pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" .or. pg(npi)%state == "sol") cycle
!$$$$$
!.. computational particle loop
    do ii = 1,indarrayFlu
      npi = Array_Flu(ii)
!$$$$$
!
!.. evaluates the terms of acceleration depending on the pressure, dissipation and viscosity, one components
!.. for each X,Y,Z direction (including the sign)
!
!.. evaluates the motion gradients and accelerations terms
!
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!if (it == iterazione .and. npi > 16900 .and. npi < 17200) then
!  write (99,*) ' it = ',it,' n. particella = ',npi
!  write (99,'(a,3d23.15)') '  pg(npi)%  coord ',pg(npi)%coord
!  write (99,'(a,9d23.15)') '  eq. moto   pg(npi)%  vel,var,dens,mass,pres ',pg(npi)%vel,pg(npi)%var,pg(npi)%dens,pg(npi)%mass,pg(npi)%pres
!end if
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!.. solving the momentum equation
      call inter_EqMoto (npi,tpres,tdiss,tvisc)
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!if (it == iterazione .and. npi > 16900 .and. npi < 17200) then
!  write (99,'(a,3d23.15)') ' tpres ',tpres
!  write (99,'(a,3d23.15)') ' tdiss ',tdiss
!  write (99,'(a,3d23.15)') ' tvisc ',tvisc
!end if
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!
!.. searches for the boundary sides nearest the npi-th current particle
!
!AA404 sub
      if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
        pg(npi)%kodvel = 0
        pg(npi)%velass = zero
      end if
!
!.. some nearest boundary has been detected
!.. selects the boundary sides that effectively contribute to the motion field terms
!AA406 sub start
!   semi-analytic approach
!AA406 sub
      if (Domain%tipo == "semi") then
         Ncbs = BoundaryDataPointer(1,npi)
         IntNcbs = BoundaryDataPointer(2,npi)
      endif
!AA406 sub end
!
!.. evaluates these contributes and adds them to the RHS
!
!AA406 sub
      if ((Ncbs > 0 .and. IntNcbs > 0).and.(Domain%tipo == "semi")) then
!
!.. computation of the boundary terms for the balance equations (using the geometrical integrals already computed)           
        call AddBoundaryContributions_to_ME2D (npi,IntNcbs,tpres,tdiss,tvisc)                   
!
        if (pg(npi)%kodvel == 0) then
          BoundReaction = zero
!.. additional repulsive force (it shouldn't be available)
          call AddElasticBoundaryReaction_2D (npi, IntNcbs, BoundReaction)
!.. save the acceleration and the assigned boundary velocities components in the particle array
          pg(npi)%acc(:) = tpres(:) + tdiss(:) + tvisc(:) + Domain%grav(:) + BoundReaction(:)
        else
          pg(npi)%acc(:) = zero
        end if
      else
!
!AA406 sub start
!AA406 sub
        if  (Domain%tipo == "semi") then
           pg(npi)%acc(:) = tpres(:) + tdiss(:) + tvisc(:) + Domain%grav(:)
        else
              pg(npi)%acc(:) = (tpres(:)+tdiss(:)+tvisc(:))/pg(npi)%Gamma + Domain%grav(:)
        endif
!AA406 sub end
! 
      end if
!
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!if (it == iterazione .and. npi > 16900 .and. npi < 17200) then
!  write (99,'(a,9d23.15)') ' accel.   pg(npi)%  vel,var,acc ',pg(npi)%vel,pg(npi)%var,pg(npi)%acc
!  write (99,'(a,3d23.15)') ' tpres ',tpres
!  write (99,'(a,3d23.15)') ' tdiss ',tdiss
!  write (99,'(a,3d23.15)') ' tvisc ',tvisc
!  write (99,*) ' '
!end if
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!
    end do
!
!$omp end parallel do
!
!.. evaluates the final motion field in terms of final velocity components 
!
!.. At first, it evaluates the coefficient of the series development of 1st order
!.. v(n+1/2) = v(n-1/2) + half * (dt(n)+dt(n+1) * (Dv/Dt)(n) 
!

!AA501b start
  if (n_bodies > 0) then
     call start_and_stop(3,6)
     call start_and_stop(2,19)
     call RHS_body_dynamics
     call start_and_stop(3,19)
     call start_and_stop(2,6)
  endif
!AA501b end

!
!AA401 
!.. time integration scheme for velocity and energy
!
!.. non staggered time integration scheme
    if (Domain%time_split == 0) then   !if time_split
      call start_and_stop(3,6)
!
!.. staggered Euler option for time integration (velocity+E), 
!.. velocity smoothing, trajectory equation, BC, neighboring parameters (start)
    else if (Domain%time_split == 1) then   !if time_split
!
!.. dt computation for the momentum and E equation (staggered Euler time integration scheme)
    dtvel = half * (dt + DtPreviousStep)         
!
!$omp parallel do default(none) private(npi,ii) shared(nag,Pg,dtvel,indarrayFlu,Array_Flu,it)
!
!.. loops on all the particles
!
!!!!!!     do npi = 1,nag
!!!!!!!
!!!!!!       if (pg(npi)%cella == 0 .or.  pg(npi)%vel_type /= "std") cycle
!!!!!!!.. the particle has a solid condition, i.e. all the components are equal to zero by default
!!!!!!       if (pg(npi)%state == "sol" ) then
!!!!!!         pg(npi)%vel(:) = zero   ! azzerato in crit_erosion
!!!!!!         pg(npi)%var(:) = zero   ! azzerato in crit_erosion
!!!!!!       else
!$$$$$
    do ii = 1,indarrayFlu
      npi = Array_Flu(ii)
!$$$$$
!
!.. kodvel = 0: the particle is internal to the domain. All the velocity components are set using the 
!..             serie development, being DV/Dt = acceleration components
      if (pg(npi)%kodvel == 0) then 
         pg(npi)%vel(:) = pg(npi)%vel(:) + dtvel * pg(npi)%acc(:)
        
!.. kodvel = 1: the particle has a critical flux condition. The horizontal velocity components are set using the 
!..             serie development, while the vertical one is assigned
      else if (pg(npi)%kodvel == 1) then            
        pg(npi)%vel(:) = pg(npi)%vel(:) + dtvel * pg(npi)%acc(:)      
        pg(npi)%vel(3) = pg(npi)%velass(3)           
!
!.. kodvel = 2: the particle has an assigned normal velocity or source condition. All the velocity components 
!..             are set to the assigned normal components values
      else if (pg(npi)%kodvel == 2) then  
        pg(npi)%vel(:) = pg(npi)%velass(:)            
!
      end if
!
!!!!!!       end if
!
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!if (it == iterazione .and. npi > 16900 .and. npi < 17200) then
!  write (99,*) ' it = ',it,' n. particella = ',npi
!  write (99,'(a,9d23.15)') ' aggiornamento velocity   pg(npi)%  vel,var,dens,mass,pres ',pg(npi)%vel,pg(npi)%var,pg(npi)%dens,pg(npi)%mass,pg(npi)%pres
!end if
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$ DA PENSARE $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      !controllo velocita'
      !call contrmach ( i,celmax )   ! controllo del numero di mach delle particelle
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
    end do
!
!$omp end parallel do
!
    call start_and_stop(3,6)
!
!=== FINE AGGIORNAMENTO VEL ############################
!
!...................................... 2011 mar 08
!.. update of Specific Internal Energy
!.. computation of the E equation (explosion cases)
    if (esplosione) then
      do ii = 1,indarrayFlu
        npi = Array_Flu(ii)
        pg(npi)%IntEn = pg(npi)%IntEn + dtvel * pg(npi)%dEdT
!.. con dt circa stesso risultato che con dtvel
!        pg(npi)%IntEn = pg(npi)%IntEn + dt * pg(npi)%dEdT
      end do
    end if
!.. end update of Specific Internal Energy
!.........................................
!

!AA501b
! Time integration for body dynamics
  if (n_bodies > 0) then
     call start_and_stop(2,19)
     call time_integration_body_dynamics(dtvel)
     call start_and_stop(3,19)
  endif
!AA501b end

!=== CORREZIONE VELOCITA' ################
!
    call start_and_stop(2,7)
!
!.. smoothing of the velocity field in order to reduce oscillations
!
!.. velocity smoothing
    call inter_SmoothVelo_2D
!
!.. smoothing of the velocity
!
!$omp parallel do default(none) &
!$omp private(npi,ii,TetaV1) &
!$omp shared(nag,Pg,Med,Domain,dt,indarrayFlu,Array_Flu,esplosione,it)
!
!.. loops on all the active particles
!
!!!!!!     do npi = 1,nag
!!!!!!!
!!!!!!       if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" .or. pg(npi)%state == "sol") cycle
!$$$$$
    do ii = 1,indarrayFlu
      npi = Array_Flu(ii)
!$$$$$
      if (esplosione) then
!.. calcolando TetaV1 con Csound si ottengono gli stessi risultati di Celerita
        TetaV1 = Domain%TetaV * pg(npi)%Csound * dt / Domain%h
!        TetaV1 = Domain%TetaV * pg(npi)%Csound * dtvel / Domain%h   !forse e' piu' corretto cosi'
      else
!.. update TetaV aggiorno TetaV adeguato al passo temporale
        TetaV1 = Domain%TetaV * Med(pg(npi)%imed)%Celerita * dt / Domain%h
!        TetaV1 = Domain%TetaV * Med(pg(npi)%imed)%Celerita * dtvel / Domain%h   !forse e' piu' corretto cosi'
!.. utilizzando TetaV fisso sia per IntEn sia per Vel e' instabile
!       TetaV1 = Domain%TetaV
      end if
!
!...................................... 2011 mar 08
!.. update of Specific Internal Energy
      if (esplosione) pg(npi)%IntEn = pg(npi)%IntEn + TetaV1 * pg(npi)%Envar
!.. end update of Specific Internal Energy
!.........................................
!
!.. the particle is inside the domain and far from the boundaries
!     
      if (pg(npi)%kodvel == 0) then        
!
!.. evaluates the average velocity components
!
!.. update tetav aggiorno tetav adeguato al passo temporale
!        TetaV1 = Domain%TetaV * Med(pg(npi)%imed)%Celerita * dt / Domain%h
!.. upgrade the velocity component array with the average velocity values previously evaluated
        pg(npi)%var(:) = pg(npi)%vel(:) + TetaV1 * pg(npi)%var(:)      
        pg(npi)%vel(:) = pg(npi)%var(:)                                        
!
!.. the particle is close to a "source", "level" or "normal velocity boundary (kodvel = 1 or = 2)
!.. the final velocity is kept unmodified
!        
      else  
        pg(npi)%var(:) = pg(npi)%vel(:) 
      end if
!
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!if (it == iterazione .and. npi > 16900 .and. npi < 17200) then
!  write (99,*) ' it = ',it,' n. particella = ',npi
!  write (99,'(a,9d23.15)') ' smooth velocity   pg(npi)%  vel,var,dens,mass,pres ',pg(npi)%vel,pg(npi)%var,pg(npi)%dens,pg(npi)%mass,pg(npi)%pres
!end if
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!
    end do
!
!$omp end parallel do
!
    call start_and_stop(3,7)
!
!=== FINE CORREZIONE VELOCITA' ################
!
!
!=== AGGIORNAMENTO COEFFICIENTE DIFFUSIVO ################
!
!.. if the diffusion model is set, the diffusion coefficient is upgraded
!
    if (diffusione) then
      call start_and_stop(2,15)
!
!$omp parallel do default(none) private(npi,ii,appo1,appo2,appo3) shared(nag,Pg,Med,indarrayFlu,Array_Flu)
!
!!!!!!       do npi = 1,nag
!!!!!!!
!!!!!!         if (pg(npi)%state == "sol") cycle
!!!!!!         if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!$$$$$
!.. update the diffusion coefficient (diffused solid-liquid intefrace)
      do ii = 1,indarrayFlu
        npi = Array_Flu(ii)
!$$$$$
!
! DA rivedere
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
!
!$omp end parallel do
!
      call start_and_stop(3,15)
    end if
!
!=== FINE AGGIORNAMENTO COEFFICIENTE DIFFUSIVO ################
!
!.. update the position of the particles
!.. evaluates the new particle positions and upgrades the coordinate array
!
    call start_and_stop(2,8)
!
!AA406 sub test
! "flu" and fixed particles particles
!$omp parallel do default(none) private(npi) shared(nag,Pg,dt,med)
!
!.. loops on the active particles
!.. trajectory computation
!
    do npi = 1,nag
!
      if (pg(npi)%cella == 0) cycle
!
!.. store the old coordinates of the particles
      pg(npi)%CoordOld(:) = pg(npi)%coord(:)
!
      if (pg(npi)%vel_type /= "std") then
!.. the type of the velocity condition is not "std", the assigned velocity condition is used
        pg(npi)%coord(:) = pg(npi)%coord(:) + dt * pg(npi)%vstart(:)
!
      else
!.. otherwise, the smoothed value V^ already calculated are used to define the particle path
        pg(npi)%coord(:) = pg(npi)%coord(:) + dt * pg(npi)%var(:) 
!
      end if
!
    end do
!
!$omp end parallel do
!
!AA406 start
! wall element trajectories
    if (Domain%tipo == "bsph") then
!$omp parallel do default(none) private(npi) shared(DBSPH,pg_w,dt)
       do npi = 1,DBSPH%n_w
          if (pg_w(npi)%cella == 0) cycle
          pg_w(npi)%coord(:) = pg_w(npi)%coord(:) + dt * pg_w(npi)%vel(:)
       end do
!$omp end parallel do
    endif
!AA406 end
!
    call start_and_stop(3,8)
!
!.. checks for the particles gone out of the domain throughout the opened sides
!
    call start_and_stop(2,9)
!
    if (NumOpenSides > 0) call CancelOutgoneParticles_2D
!.. adds the new particles from the source conditions
!
    if (SourceSide /= 0) call GenerateSourceParticles_2D 
!
!.. reorders all the particles (olds and news) with respect the storage array and the virtual cell grid
!
    call OrdGrid1 ( nout )
!
    call start_and_stop(3,9)
!
!..??? CALCOLO QUANTITA' PER PARTICELLE FISSE - 20050516 ELIMINATO!
!
!.. set the parameters for the fixed particles 
!
    if ( Domain%NormFix )  call NormFix
!
!.. evaluates the distances between the one particle and all the nearest ones, storing them 
!.. in the general storage array
!.. Computation and storage of the neighbouring indices, distances and kernel gradients
!
    call start_and_stop(2,10)
    call CalcVarLength
    call start_and_stop(3,10)

!
!AA406 start (RK1stag)
! subroutine for wall BC treatment (BSPH)
! density and pressure updates for wall elements: MUSCL + partial Riemann solver + state equation 
    call start_and_stop(2,18)
    if ((Domain%tipo == "bsph") .and. (nag > 0) .and. (DBSPH%n_w > 0)) then
       call Gradients_to_MUSCL   
       call BC_wall_elements
    endif
    call start_and_stop(3,18)
!AA406 end
!
!.. evaluates the close boundaries and the integrals
!.. for the current particle in every loop, storing them in the general storage array
!.. Computation and storage of the boundary integrals
!
!AA406 sub start
!AA406 sub
    if (Domain%tipo == "semi") then 
       call start_and_stop(2,11)
       call ComputeBoundaryDataTab
       call start_and_stop(3,11)
    endif
!AA406 sub end
!
!AA401
!.. staggered Euler option for time integration (velocity+E), 
!.. velocity smoothing, trajectory equation, BC, neighboring parameters (end)
    end if   !if time_split
!
!=== EQUAZIONE DI CONTINUITA' ################
!
!.. erosion criterium + continuity equation right hand side + deformation tensor norm 
    call start_and_stop(2,12)
!
    if (erosione .and. .not. esplosione) then  ! erosione
!
      if (Domain%time_split == 1) then 
!
!.. Calcolo stato ('flu' o 'sol') delle particelle di tipo 'granulare'
!.. test modello erosione

!AA504 sub start
! Calling the proper subroutine for the erosion criterion 
    select case (Granular_flows_options%ID_erosion_criterion)
       case(1)
!$omp parallel do default(none) shared(pg,nag) private(npi,ncel)
          do npi=1,nag
             pg(npi)%vel_old(:) = pg(npi)%vel(:)
             pg(npi)%normal_int_old(:) = pg(npi)%normal_int(:)
             call initialization_fixed_granular_particle(npi)             
          end do
!$omp end parallel do 
!$omp parallel do default(none) shared(pg,nag) private(npi) 
          do npi=1,nag
              call Shields(npi) 
          end do
!$omp end parallel do
! Inizializzo valori default di viscosità per particelle fisse
!$omp parallel do default(none) shared(pg,nag,Granular_flows_options,Med) private(npi,ncel,aux,igridi,jgridi,kgridi)
          do npi=1,nag
             ncel = ParticleCellNumber(pg(npi)%coord)
             aux = CellIndices(ncel,igridi,jgridi,kgridi)
             if (Granular_flows_options%ID_erosion_criterion==1) then
                 if (pg(npi)%state=="sol") then
                   pg(npi)%mu = Med(Granular_flows_options%ID_main_fluid)%visc*Med(Granular_flows_options%ID_main_fluid)%den0
                   pg(npi)%visc = Med(Granular_flows_options%ID_main_fluid)%visc
                endif
             endif
          end do
!$omp end parallel do  
       case(2)
!$omp parallel do default(none) shared(pg,nag) private(npi)
          do npi=1,nag
             call Shields(npi) 
          end do
!$omp end parallel do 
       case(3)
! Compute second invariant of the rate strain tensor and density derivatives
          call inter_EqCont_2D 
! Mohr-Coulomb criterion
          call MohrC
       case default
    end select
!AA504 sub end

!.. aggiornamento vettore appoggio per conteggio particelle non 'sol' caso con 'granular'
!
!        indarraySol = 0
        indarrayFlu = 0
        do npi = 1,nag
          if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!          if (index(Med(pg(npi)%imed)%tipo,"granular") > 0) then
!            indarraySol = indarraySol + 1
!            Array_Sol(indarraySol) = npi
!          end if
          if (pg(npi)%state == "flu") then
            indarrayFlu = indarrayFlu + 1
!.. controllo di non superare le dimensioni ipotizzate per il vettore e nel caso li ridefinisco opportunamente
            if (indarrayFlu > PARTICLEBUFFER) then
              call diagnostic (arg1=9,arg2=2,arg3=nomsub)
            end if
            Array_Flu(indarrayFlu) = npi
          end if
        end do
      end if
!!write (999,*) '  '
!!write (999,*) ' ==> tempo =',tempo,'  indarrayFlu = ',indarrayFlu
!
!.. Aggiorno Calcolo Invariante secondo e derivata delle densita'
      call inter_EqCont_2D
!
    else  ! erosione
!.. non sto eseguendo un calcolo con un modello di erosione 
!
!.. Calcolo Invariante secondo e derivata delle densita'
      call inter_EqCont_2D 
!
      if (Domain%time_split == 1) then 
!
!.. devo aggiornare vettore appoggio per conteggio particelle 'flu' a tutte nag
        indarrayFlu = 0
        do npi = 1,nag
          if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
          indarrayFlu = indarrayFlu + 1
!.. controllo di non superare le dimensioni ipotizzate per il vettore e nel caso li ridefinisco opportunamente
          if (indarrayFlu > PARTICLEBUFFER) then
            call diagnostic (arg1=9,arg2=2,arg3=nomsub)
          end if
          Array_Flu(indarrayFlu) = npi
        end do
      end if
!.. le particelle escluse eventuali con cella=0 (???) e non 'std'
    end if  ! erosione
!
!AA404 sub
   if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
!$omp parallel do default(none) private(npi) shared(nag,pg)
      do npi = 1,nag
!.. set the default
        pg(npi)%koddens = 0
        pg(npi)%densass = pg(npi)%dens
      end do
!$omp end parallel do
    end if
!
!AA406 sub start
!.. loops on all the active particles
!.. Contributo contorno alla derivata della densita' e aggiornamento della densita'
    if (Domain%tipo == "semi") then
!$omp parallel do default(none) &
!$omp private(npi,ii,BCrodivV,Ncbs,IntNcbs) &
!$omp shared(nag,pg,BoundaryDataPointer,indarrayFlu,Array_Flu,it)
!AA406 sub end
!
!!!!!!    do npi = 1,nag
!!!!!!!
!!!!!!      if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!!!!!!!
!!!!!!!.. set the default
!!!!!!      pg(npi)%koddens = 0
!!!!!!      pg(npi)%densass = pg(npi)%dens
!!!!!!!
!!!!!!      if (pg(npi)%state == "flu") then
!
!$$$$$
!.. boundary contributions to the Continuity Equation
    do ii = 1,indarrayFlu
      npi = Array_Flu(ii)
!$$$$$
!
!.. Ricerca i lati di contorno vicini alla particella "npi"
      BCrodivV = zero
      Ncbs = BoundaryDataPointer(1,npi)
      IntNcbs = BoundaryDataPointer(2,npi)
!.. Seleziona i lati che danno effettivo contributo
      if (Ncbs > 0 .and. IntNcbs > 0) then     !Ci sono lati di contorno che danno contributo
        call AddBoundaryContribution_to_CE2D (npi, IntNcbs, BCrodivV)  !, Ncbs
      end if
!
!.. boundary type is fixe or tapis or level(?)
      if (pg(npi)%koddens == 0) then
        pg(npi)%dden = pg(npi)%dden - BCrodivV
!.. boundary type .... level(?)
      else if (pg(npi)%koddens == 1) then
        pg(npi)%dden = zero  ! da assegnare
! utile per prossime versioni       pg(npi)%dens = pg(npi)%densass  ! Viene mantenuta costante la densita'
! utile per prossime versioni       pg(npi)%koddens = 0
!.. boundary type is velocity or source
      else if (pg(npi)%koddens == 2) then
        pg(npi)%dden = zero
      end if
!
!!!      end if
!
!  id = omp_get_thread_num()
!  write (*,*) 'Item ',npi,'executed from thread', id, '  time = '
!
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!if (it == iterazione .and. npi > 16900 .and. npi < 17200) then
!  write (99,*) ' it = ',it,' n. particella = ',npi
!  write (99,'(a,9d23.15)') ' de density   pg(npi)%  vel,var,dens,mass,pres ',pg(npi)%vel,pg(npi)%var,pg(npi)%dens,pg(npi)%mass,pg(npi)%pres
!end if
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!
    end do
!
!$omp end parallel do
!
!AA406
    endif !(Domain%tipo == "semi") 
!
!.. loops on all the active particles
!.. Contributo contorno alla derivata della densita' e aggiornamento della densita'
!
!!!!    do npi = 1,nag
!$$$$$
!AA401
!.. time integration of the continuity equation
!
!.. non staggered time integration scheme
    if (Domain%time_split == 0) then   !if time_split
      call start_and_stop(3,12)
!
!.. staggered Euler option for time integration of the continuity equation + state equation
    else if (Domain%time_split == 1) then   !if time_split
!
!$omp parallel do default(none) private(npi,ii) shared(nag,pg,dt,indarrayFlu,Array_Flu,it,Domain,Med)
    do ii = 1,indarrayFlu
      npi = Array_Flu(ii)
!$$$$$
!
!!!!!! non serve      if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!
!AA406
      if (Domain%tipo == "bsph") pg(npi)%dden=pg(npi)%dden/pg(npi)%uni
!
!.. boundary type is fixe or tapis or level(?)
      if (pg(npi)%koddens == 0) then
!
!AA406 sub 
        if (Domain%tipo == "semi") pg(npi)%dens = pg(npi)%dens + dt * pg(npi)%dden
!
        pg(npi)%densass = zero
!.. boundary type .... level(?)
!      else if (pg(npi)%koddens == 1) then
! utile per prossime versioni       pg(npi)%dens = pg(npi)%densass  ! Viene mantenuta costante la densita'
! utile per prossime versioni       pg(npi)%koddens = 0
!.. boundary type is velocity or source
      else if (pg(npi)%koddens == 2) then
        pg(npi)%dens = pg(npi)%densass  ! Viene mantenuta costante la densita'
      end if
      
!AA504 start 
      if (Domain%density_thresholds == 1) then       
         if (pg(npi)%dens<(0.9*med(pg(npi)%imed)%den0))  pg(npi)%dens = 0.9*med(pg(npi)%imed)%den0
         if (pg(npi)%dens>(1.1*med(pg(npi)%imed)%den0))  pg(npi)%dens = 1.1*med(pg(npi)%imed)%den0
      endif
!AA504 end       
      
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!if (it == iterazione .and. npi > 16900 .and. npi < 17200) then
!  write (99,*) ' it = ',it,' n. particella = ',npi
!  write (99,'(a,9d23.15)') ' agg. density   pg(npi)%  vel,var,dens,mass,pres ',pg(npi)%vel,pg(npi)%var,pg(npi)%dens,pg(npi)%mass,pg(npi)%pres
!end if
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!
    end do
!
!$omp end parallel do
!
    call start_and_stop(3,12)
!
!=== FINE EQUAZIONE DI CONTINUITA' ################
!write(*,*) it, tempo, dt, 'dopo eq. continuita'
!
!
!.. modello bifluido
!.. solid-liquid interface diffusion model
!
    if (diffusione) then
      call start_and_stop(2,16)
      call aggdens
      call start_and_stop(3,16)
    end if
!
!.. modello bifluido
!
!=== FINE AGGIORNAMENTO DENS PARTICELLE ################
!write(*,*) it, tempo, dt, 'dopo aggdens'
!
!=== EQUAZIONE DI STATO #########################################
!
!.. evaluates the new pressure field
!.. state equation according to the speed of sound definition and eventually to the liquid-solid diffusion model (Sphera_Tools.f90)
!
    call start_and_stop(2,13)
    call calcpre  
    call start_and_stop(3,13)
    
!AA501btest
    if (n_bodies > 0) then
       call start_and_stop(2,19)
       call body_pressure_mirror
       call start_and_stop(3,19)
    endif
!AA501b end    
    
!
!AA401
!.. staggered Euler option for time integration of the continuity equation + state equation (end)
    end if   !if time_split
!
!=== FINE EQUAZIONE DI STATO ####################################
!
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!    do ii = 1,indarrayFlu
!      npi = Array_Flu(ii)
!if (it == iterazione .and. npi > 16900 .and. npi < 17200) then
!  write (99,*) ' it = ',it,' n. particella = ',npi
!  write (99,'(a,9d23.15)') ' calc pressure   pg(npi)%  vel,var,dens,mass,pres ',pg(npi)%vel,pg(npi)%var,pg(npi)%dens,pg(npi)%mass,pg(npi)%pres
!end if
!    end do
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!
!AA401
!.. time integration for RK schemes
    if (Domain%time_split == 0) call time_integration
!
!=== CORREZIONE (SMOOTHING) PRESSIONE PARTICELLE ################
!
!.. performs the pressure smoothing if the smoothing parameter is assigned
!.. pressure smoothing and density update
!
    call start_and_stop(2,14)
!
    if (Domain%TetaP > zero) then
!
      if (Domain%Psurf == 's') then
!
!.. set the pressure smoothed with the weight of atmospheric pressure and set the density using selective method
!
        call inter_SmoothPres
!
      else if (Domain%Psurf == 'a') then
!
!.. set the pressure smoothed with the weight of atmospheric pressure and set the density
!

!AA501b
        call PressureSmoothing_2D

!
      end if
!
    end if
!
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!    do ii = 1,indarrayFlu
!      npi = Array_Flu(ii)
!if (it == iterazione .and. npi > 16900 .and. npi < 17200) then
!  write (99,*) ' it = ',it,' n. particella = ',npi
!  write (99,'(a,9d23.15)') ' smooth. pressure   pg(npi)%  vel,var,dens,mass,pres ',pg(npi)%vel,pg(npi)%var,pg(npi)%dens,pg(npi)%mass,pg(npi)%pres
!end if
!    end do
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!
    call start_and_stop(3,14)
!

!AA501b
    if (n_bodies > 0) then
       call start_and_stop(2,19)
       call body_pressure_postpro
       call start_and_stop(3,19)
    endif

!.. modello di diffusione
!§
    if (diffusione) then
      call start_and_stop(2,16)
!
!$omp parallel do default(none) private(npi,ii) shared(nag,Pg,Med,indarrayFlu,Array_Flu)
!
!!!       do npi = 1,nag
!!!!
!!!         if (pg(npi)%state == "sol" .or. pg(npi)%koddens /= 0) cycle
!$$$$$
      do ii = 1,indarrayFlu
        npi = Array_Flu(ii)
        if (pg(npi)%koddens /= 0) cycle
!$$$$$
!
        if (pg(npi)%imed == 1) then
          pg(npi)%dens = pg(npi)%pres / (Med(1)%celerita * Med(1)%celerita) & 
                     + (Med(2)%den0 * VFmn + Med(1)%den0 * (one-VFmn))
        else if (pg(npi)%imed == 2) then
          pg(npi)%dens = pg(npi)%pres / (Med(2)%celerita * Med(2)%celerita) & 
                     + (Med(2)%den0 * VFmx + Med(1)%den0 * (one-VFmx))
        end if
        Pg(npi)%rhoc = pg(npi)%pres / (med(2)%celerita * med(2)%celerita) + med(2)%den0
        Pg(npi)%rhow = pg(npi)%pres / (med(1)%celerita * med(1)%celerita) + med(1)%den0
!
      end do
!
!$omp end parallel do
!
      call start_and_stop(3,16)
    end if
!§
!
    
!AA504 start
  call start_and_stop(2,20)
  if (Granular_flows_options%ID_erosion_criterion==1) call mixture_viscosity 
  call start_and_stop(3,20)
!AA504 end
    
!=== CALCOLO VISCOSITA' PARTICELLE ################ 
!.. apparent viscosity (erosion model)
!
!AA406
    if (diffusione .or. esplosione) then
!
    if ((Domain%time_split == 1) .or. (Domain%time_stage == Domain%RKscheme)) then
      call start_and_stop(2,15)
      call viscapp
      call start_and_stop(3,15)
    end if
!AA406
    endif
!
!=== FINE CALCOLO VISCOSITA' PARTICELLE ################
!
!AA406
     if (Domain%tipo == "semi") then
!
!=== Boundary Conditions (start) ################
! BC (checks for the particles gone out of the domain throughout the opened sides)
    if ((Domain%time_split == 0) .and. (Domain%time_stage == Domain%RKscheme)) then
      call start_and_stop(2,9)
      if (NumOpenSides > 0) call CancelOutgoneParticles_2D
!.. adds the new particles from the source conditions
      if (SourceSide /= 0) call GenerateSourceParticles_2D 
!.. reorders all the particles (old and new) with respect to the storage array and the virtual cell grid
      call OrdGrid1 (nout)
      call start_and_stop(3,9)
!  Set the parameters for the fixed particles 
      if ( Domain%NormFix )  call NormFix
!  (devo aggiornare vettore appoggio per conteggio particelle 'flu' a tutte nag)
      indarrayFlu = 0
     do npi = 1,nag
        if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
          indarrayFlu = indarrayFlu + 1
!  (controllo di non superare le dimensioni ipotizzate per il vettore e nel caso li ridefinisco opportunamente)
        if (indarrayFlu > PARTICLEBUFFER) then
          call diagnostic (arg1=9,arg2=2,arg3=nomsub)
        end if
        Array_Flu(indarrayFlu) = npi
      enddo 
    endif
!=== Boundary Conditions (end) ################
!
!AA406
     endif
!
!AA406 start
! subroutine for wall BC treatment (BSPH)
! density and pressure updates for wall elements: MUSCL + partial Riemann solver + state equation 
       call start_and_stop(2,18)
       if ((Domain%tipo == "bsph") .and. (nag > 0) .and. (DBSPH%n_w > 0)) then
          call Gradients_to_MUSCL   
          call BC_wall_elements
       endif
       call start_and_stop(3,18)
!AA406 end
!
!AA402 start
!.. update stage parameter
!
    if ((Domain%RKscheme > 1) .and. (Domain%time_split == 0)) then
      Domain%time_stage = MODULO(Domain%time_stage,Domain%RKscheme)
      Domain%time_stage = Domain%time_stage + 1
!.. End intermediary (time stage) loop
      if (Domain%time_stage < 2) then
        continue
      else
        go to 1100
      end if
    end if
!AA402 end

!=== SCRITTURA SU FILE RISULTATI ######################
!
!.. output
!AA401 
    if (Domain%time_split == 0) dtvel = dt
!
    call s_secon2 ( tf_iter )
    tot_iter = tot_iter + tf_iter - ti_iter
!
    if ( nout>0) then
      call print_results ( it_eff, it_print, 'loop__' )
    end if
    if ( nres>0 ) then
      call memo_results ( it_eff, it_memo, it_rest, dtvel, 'loop__' )
    end if
!
    !if (  ( mod(it+1,Domain%icpoi_fr) == 0 ) .OR. &
    !      ( mod(it+1,Domain%ipllb_fr) == 0 ) .OR. &
    !      ( mod(tempo,Domain%cpoi_fr) <= dtvel ) .OR. &
    !      ( mod(tempo,Domain%pllb_fr) <= dtvel ) &
    !   )  then
    !Calcolo grandezze nei punti di controllo (SEMPRE?)
    call CalcVarp
    !end if
!
    if (Domain%icpoi_fr > 0) then 
       if ( ( mod(it,Domain%icpoi_fr) == 0 ) .AND. npointst > 0 ) then
          call memo_ctl
       endif   
!AA501b
      if (n_bodies > 0) then
         if (mod(it,Domain%icpoi_fr) == 0)  call Body_dynamics_output      
      endif
      
      else if (Domain%cpoi_fr > zero) then
         if ( ( mod(tempo,Domain%cpoi_fr) <= dtvel ) .AND. npointst > 0 ) then
            call memo_ctl
          endif   
!AA501b
         if (n_bodies > 0) then
            if (mod(tempo,Domain%cpoi_fr) <= dtvel)  call Body_dynamics_output      
         endif
!AA601 start
! Writing DBSPH post-processing and update wet variable in pg array
         if ((Domain%tipo=="bsph").and.(mod(tempo,Domain%cpoi_fr)<=dtvel)) then
            if ((DBSPH%n_monitor_points>0).or.(DBSPH%n_monitor_regions==1)) call wall_elements_pp
!$omp parallel do default(none) shared(DBSPH,pg_w) private(npi)
            do npi=1,DBSPH%n_w
               pg_w(npi)%wet = 0 
            end do    
!$omp end parallel do
         endif
!AA601 end
    end if
!
    if (Domain%ipllb_fr > 0) then
      if ( ( mod(it,Domain%ipllb_fr) == 0 ) .AND. nlines > 0 ) then
        call  calc_pelo
      end if
    else if (Domain%pllb_fr > zero) then
      if ( ( mod(tempo,Domain%pllb_fr) <= dtvel ) .AND. nlines > 0 ) then
        call  calc_pelo
      end if
    end if
    
!AA504 start
    if (Granular_flows_options%monitoring_lines>0) then
       if ((int(tempo/Granular_flows_options%dt_out)>Granular_flows_options%it_out_last).or.(it_corrente==1)) then
           Granular_flows_options%it_out_last = int(tempo/Granular_flows_options%dt_out)
           call write_Granular_flows_interfaces
       endif
    endif
    if (Q_sections%n_sect>0) then
       if ((int(tempo/Q_sections%dt_out)>Q_sections%it_out_last).or.(it_corrente==1)) then
           Q_sections%it_out_last = int(tempo/Q_sections%dt_out)
           call sub_Q_sections
       endif
    endif
!AA504 end 
    
!
!=== FINE SCRITTURA SU FILE RISULTATI ######################
!
!
!=== INIZIO SCRITTURA SU FILE AVANZAMENTO FRONTE ######################
!.. post-processing for water front
!
    if (Domain%imemo_fr > 0) then
      if ( mod(it,Domain%imemo_fr) == 0 ) then
        xmax = -1.0d30
        ymax = -1.0d30
        do npi = 1,nag
          if ( pg(npi)%vel_type /= "std" .or. pg(npi)%cella == 0) cycle
          xmax = max(xmax,pg(npi)%coord(1))
          ymax = max(ymax,pg(npi)%coord(3))
        end do
        write (nfro,'(2g14.7,13x,a,g14.7)') tempo, xmax ,'-', ymax
      end if
    else if (Domain%memo_fr > zero) then
      if ( it > 1 .AND. mod(tempo,Domain%memo_fr) <= dtvel ) then
        xmax = -1.0d30
        ymax = -1.0d30
        do npi = 1,nag
          if ( pg(npi)%vel_type /= "std" .or. pg(npi)%cella == 0) cycle
          xmax = max(xmax,pg(npi)%coord(1))
          ymax = max(ymax,pg(npi)%coord(3))
        end do
        write (nfro,'(2g14.7,13x,a,g14.7)') tempo, xmax ,'-', ymax
      end if
    end if
!
!=== FINE SCRITTURA SU FILE AVANZAMENTO FRONTE ######################
!
!
!=== SCRITTURA SU FILE RESTART ######################
!    if ( nsav > 0 ) then
!      call memo_restart ( it_eff, it_rest, dtvel, nout, nsav )
!    end if
!write(*,*) it, tempo, dt, 'dopo calcolo memo_restart'
!=== FINE SCRITTURA SU FILE RESTART ######################
!
!
!=== SCRITTURA SU FILE FORMATO VTK ######################
! 
    if ( vtkconv ) then
      call result_converter ('loop__')
    end if
!
!write(*,*) it, tempo, dt, 'dopo calcolo result_converter'
!=== FINE SCRITTURA SU FORMATO VTK ######################
!
!
! check if file kill exists. if it exists the execution is stopped saving the results.
    inquire (file=nomefilekill, EXIST=kill_flag)
    if (kill_flag) exit ITERATION_LOOP
!
    if ( tempo >= Domain%tmax ) exit ITERATION_LOOP
!
!write(*,*) it, tempo, dt, 'dopo calcolo ITERATION_LOOP'
!
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!if (it == iterazione) then
!  close (unit=99)
!  stop
!end if
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!
  end do  ITERATION_LOOP !! fine loop sul it di interazione

! END LOOP
!-----------------------------------------------------------
!
!Salvataggio risultati e restart ULTIMO STEP (sempre?)
  if ( it_eff /= it_print .AND. nout>0) then
    it_print = it_eff
    call print_results ( it_eff, it_print, 'fine__' )
  end if
  if ( it_eff /= it_memo .AND. nres > 0 ) then
!    it_memo = it_eff
!    it_rest = it_eff
    call memo_results ( it_eff, it_memo, it_rest, dtvel, 'fine__' )
  end if
!  if ( it_eff /= it_rest .AND. nsav > 0 ) then
!    call memo_restart ( it_eff, it_rest, dtvel, nout, nsav )
!  end if
  if ( vtkconv ) then
    call result_converter ('fine__')
  end if
!
  call s_secon2 ( tf )
!
  if ( nout > 0 ) then

    write (nout,*) " "
    write (nout,'(a)')   "----------------------------------------------------------------------------------------"
    write (nout,*) " "
    SpCountot = 0
    do i = 1,Nmedium
      SpCountot = SpCountot + SpCount(i)
      write (nout,'(a,i15,a,a)') "Number of source particles        :  SpCount = ",SpCount(i)," medium ",Med(i)%tipo
    end do
    write (nout,'(a,i15)') "Number of total source particles  :  SpCountot = ",SpCountot
    write (nout,*) " "
    OpCountot = 0
    do i = 1,Nmedium
      OpCountot = OpCountot + OpCount(i)
      write (nout,'(a,i15,a,a)') "Number of outgone particles       :  OpCount = ",OpCount(i)," medium ",Med(i)%tipo
    end do
    write (nout,'(a,i15)') "Number of total outgone particles :  OpCountot = ",OpCountot
    write (nout,*) " "
    EpCountot = 0
    do i = 1,Nmedium
      EpCountot = EpCountot + EpCount(i)
      write (nout,'(a,i15,a,a)') "Number of escaped particles       :  EpCount = ",EpCount(i)," medium ",Med(i)%tipo
    end do
    write (nout,'(a,i15)') "Number of total escaped particles :  EpCountot = ",EpCountot
    write (nout,*) " "
    EpOrdGridtot = 0
    do i = 1,Nmedium
      EpOrdGridtot = EpOrdGridtot + EpOrdGrid(i)
      write (nout,'(a,i15,a,a)') "Number of escaped particles (OrdGrid1)      :  EpOrdGrid = ",EpOrdGrid(i)," medium ",Med(i)%tipo
    end do
    write (nout,'(a,i15)') "Number of total escaped particles (OrdGrid1):  EpOrdGridtot = ",EpOrdGridtot
    write (nout,*) " "
    write (nout,*) " "
    write (nout,*) "Final number of particles:       NAG = ",nag
!
!AA406
    if (Domain%tipo == "bsph") write (nout,*) "Final number of wall particles:       DBSPH%n_w = ",DBSPH%n_w
!   
    write (nout,*) " "
    write (nout,'(a)')   "----------------------------------------------------------------------------------------"
!
    call writime2 ( ti,tf,nout )
!
    write (nout,*) " "
    write (nout,'(a)')   "----------------------------------------------------------------------------------------"
    write (nout,*) " "
    write (nout,*) "Number of steps:          ",it_eff
    write (nout,*) "Total CPU steps time:     ",tot_iter(1)
    write (nout,*) "Total Elapsed steps time: ",tot_iter(2)
    write (nout,*) "Average CPU time for step:",tot_iter(1)/it_eff
    write (nout,*) " "
    write (nout,'(a)')   "----------------------------------------------------------------------------------------"
    write (nout,*) " "
!
  end if
!
!
!AA406 
  if ((Domain%tipo == "bsph").and.(DBSPH%n_w > 0)) deallocate (pg_w)

  return
  end subroutine Loop_Irre_2D
!---split

