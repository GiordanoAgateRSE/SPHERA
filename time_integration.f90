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

!cfile Euler.f90
!AA401 all the subroutine
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : Euler
!
! Last updating : April 18, 2013
!
! Improvement traceback:
! 00  Amicarelli       30Jun11   First development
! 01  Amicarelli/Agate 30Jun11   Adaptation for the first stage of RK2
! 02  Amicarelli/Agate 18apr13   Adaptation for sloshing case
! 03  Amicarelli       26Jan15   DBSPH-input (AA601). DBSPH kinematics and minor DBSPH modifications.
!
!************************************************************************************
! Module purpose : Euler time integration
!
! Calling routine: time_integration
!
! Called routines: start_and_stop, inter_SmoothVelo_2D, inter_SmoothVelo_3D,
!AA601 sub
!                  sloshing_tank_kinematics,DBSPH_kinematics,DBSPH_inlet_outlet
!
!************************************************************************************
!
  subroutine Euler  
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
!Scalar declaration
!AA601 sub
  integer(4)       :: npi,ii,i,j
  double precision :: TetaV1
!
!AA406
  double precision, dimension(:), allocatable :: uni_old
!
!AA406test
  double precision :: beta 
!AA601
integer(4),external :: ParticleCellNumber
!.. Executable Statements ..
! 
!$omp parallel do default(none) private(npi,ii) shared(nag,Pg,Domain,dt,indarrayFlu,Array_Flu,ts0_pg,esplosione)
!
!.. velocity and energy integration
   do ii = 1,indarrayFlu
     npi = Array_Flu(ii)
!
!AA402 start
!.. Store the results of the previous step for time stage parameters "ts0" (during the first stage of Heun's scheme)
!.. The coordinates and the smoothed velocities are saved during a following loop, even involving solid particles
     if (Domain%RKscheme == 2) then
       ts0_pg(npi)%ts_vel(:)=pg(npi)%vel(:)  
       ts0_pg(npi)%ts_acc(:)=pg(npi)%acc(:)                        
       ts0_pg(npi)%ts_dden=pg(npi)%dden
       ts0_pg(npi)%ts_dens=pg(npi)%dens
       if (esplosione) then
         ts0_pg(npi)%ts_IntEn=pg(npi)%IntEn
         ts0_pg(npi)%ts_dEdT=pg(npi)%dEdT
       end if
     end if
!AA402 end
!.. kodvel = 0: the particle is internal to the domain. All the velocity components are set using the 
!..             series development, being dv/dt = acceleration components
     if (pg(npi)%kodvel == 0) then 
        pg(npi)%vel(:) = pg(npi)%vel(:) + dt * pg(npi)%acc(:)
!
!AA406test
!AA601 sub start
! Impose boundary conditions at inlet and outlet sections (DB-SPH)
      if (Domain%tipo == "bsph") call DBSPH_inlet_outlet(npi)
!AA601 sub end
!.. kodvel = 1: the particle has a critical flux condition. The horizontal velocity components are set using the 
!..             series development, while the vertical one is assigned
     else if (pg(npi)%kodvel == 1) then 
       pg(npi)%vel(:) = pg(npi)%vel(:) + dt * pg(npi)%acc(:)      
       pg(npi)%vel(3) = pg(npi)%velass(3)           
!.. kodvel = 2: the particle has an assigned normal velocity or source condition. All the velocity components 
!..             are set to the assigned normal components values
     else if (pg(npi)%kodvel == 2) then  
       pg(npi)%vel(:) = pg(npi)%velass(:)            
     end if
! energy time integration (explosion option)
     if (esplosione) then
       pg(npi)%IntEn = pg(npi)%IntEn + dt * pg(npi)%dEdT
     end if
   end do
!
!$omp end parallel do
!
   call start_and_stop(3,17)
!
! velocity and energy smoothing
   call start_and_stop(2,7)
   if (ncord==2) then 
     call inter_SmoothVelo_2D
   else
     call inter_SmoothVelo_3D
   end if
!
!$omp parallel do default(none) &
!$omp private(npi,ii,TetaV1) &
!$omp shared(nag,Pg,Med,Domain,dt,indarrayFlu,Array_Flu,esplosione)
!
   do ii = 1,indarrayFlu
     npi = Array_Flu(ii)
     if (esplosione) then
!.. calcolando TetaV1 con Csound si ottengono gli stessi risultati di Celerita
       TetaV1 = Domain%TetaV * pg(npi)%Csound * dt / Domain%h
     else
!.. update TetaV aggiorno TetaV adeguato al passo temporale
       TetaV1 = Domain%TetaV * Med(pg(npi)%imed)%Celerita * dt / Domain%h
     end if
!
!.. update of Specific Internal Energy
     if (esplosione) pg(npi)%IntEn = pg(npi)%IntEn + TetaV1 * pg(npi)%Envar
!
!.. the particle is inside the domain and far from the boundaries
     if (pg(npi)%kodvel == 0) then        
!.. evaluates the average velocity components
!.. upgrade the velocity component array with the average velocity values previously evaluated
       pg(npi)%var(:) = pg(npi)%vel(:) + TetaV1 * pg(npi)%var(:)      
       pg(npi)%vel(:) = pg(npi)%var(:)                                        
     else  
!.. the particle is close to a "source", "level" or "normal velocity boundary (kodvel = 1 or = 2)
!.. the final velocity is kept unmodified
       pg(npi)%var(:) = pg(npi)%vel(:) 
     end if
   end do
!
!$omp end parallel do
!
   call start_and_stop(3,7)
!
   call start_and_stop(2,17)
!
!$omp parallel do default(none) private(npi) shared(nag,Pg,Domain,dt,ts0_pg,it_corrente)
!
!.. time integration for position and density (trajectories and continuity equation)
   do npi = 1,nag
     if (pg(npi)%cella == 0) cycle
!
! storing the coordinates of the previous time step (RK2) for time integration
       if (Domain%RKscheme == 2) then
         ts0_pg(npi)%ts_coord(:) = pg(npi)%coord(:)
         ts0_pg(npi)%ts_var(:) = pg(npi)%var(:)
       endif

!.. store the old coordinates of the particles
     pg(npi)%CoordOld(:) = pg(npi)%coord(:)
!
     if (pg(npi)%vel_type /= "std") then
!.. the type of the velocity condition is not "std", the assigned velocity condition is used
       pg(npi)%coord(:) = pg(npi)%coord(:) + dt * pg(npi)%vstart(:)
     else
!.. otherwise, the smoothed value V^ already calculated are used to define the particle path
           pg(npi)%coord(:) = pg(npi)%coord(:) + dt * pg(npi)%var(:) 
     end if
   end do
!
!$omp end parallel do
!AA601 sub start
! wall element trajectories
    if (((DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet)>0).and.(Domain%tipo=="bsph")) then
!commented part (maybe useful for further developments)
!!omp parallel do default(none) private(npi) shared(n_w,pg_w,dt)
!       do npi = 1,n_w
!          if (pg_w(npi)%cella == 0) cycle
!          pg_w(npi)%coord(:) = pg_w(npi)%coord(:) + dt * pg_w(npi)%vel(:)
!       end do
!!omp end parallel do
       if (DBSPH%n_w>0) then
       call DBSPH_kinematics
! Time integration for the first surface element velocity    
! Loop over the surface element to copy and paste, all over the elements, the imposed kinematics of the first surface element 
       pg_w(1)%coord(:) = pg_w(1)%coord(:) + dt * pg_w(1)%vel(:)
!$omp parallel do default(none) shared(DBSPH,pg_w,dt) private(npi)
       do npi=2,DBSPH%n_w
          pg_w(npi)%vel(:) = pg_w(1)%vel(:) 
          pg_w(npi)%coord(:) = pg_w(npi)%coord(:) + dt * pg_w(1)%vel(:) 
       enddo
!$omp end parallel do
      endif
! In-built motion of control lines
       if (DBSPH%in_built_monitors==.true.) then
          do i = 1,nlines
!.. loop sui punti della linea
             do j = control_lines(i)%Icont(1),control_lines(i)%Icont(2) 
                control_points(j)%coord(:) = control_points(j)%coord(:) + dt * pg_w(1)%vel(:)
                control_points(j)%cella = ParticleCellNumber(control_points(j)%coord(:))
             end do
          end do
       endif
!Boundary Conditions (start) 
! BC (checks for the particles gone out of the domain throughout the opened sides)
       if ((Domain%time_split == 0) .and. (Domain%time_stage == Domain%RKscheme)) then
          call start_and_stop(2,9)
          if (ncord==2) then 
             if (NumOpenSides > 0) call CancelOutgoneParticles_2D
! Adds the new particles from the source conditions
             if (SourceSide /= 0) call GenerateSourceParticles_2D 
             else
                if (NumOpenFaces > 0) call CancelOutgoneParticles_3D
! Adds the new particles from the source conditions
                if (SourceFace /= 0) call GenerateSourceParticles_3D 
          endif
! Reorders all the particles (old and new) with respect to the storage array and the virtual cell grid
          call OrdGrid1 (nout)
          call start_and_stop(3,9)
! Set the parameters for the fixed particles 
          if ( Domain%NormFix )  call NormFix
          allocate (uni_old(nag))
! Update auxiliary vector to count "flu" particles in nag 
          indarrayFlu = 0
          do npi = 1,nag
             if (pg(npi)%rhoSPH_new /= zero) then
                uni_old(npi) = pg(npi)%uni
                else
                   uni_old(npi) = zero
             endif
             if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
             indarrayFlu = indarrayFlu + 1
             Array_Flu(indarrayFlu) = npi
          enddo 
       endif
!Boundary Conditions (end) 
       do npi=1,nag 
          pg(npi)%Gamma = pg(npi)%Gamma + pg(npi)%dShep * dt
          if (pg(npi)%Gamma >= one) pg(npi)%Gamma = one 
       end do
       call CalcVarLength
!$omp parallel do default(none) private(npi,ii) shared(nag,pg,med,dt,indarrayFlu,Array_Flu,Domain,uni_old,beta,nPartIntorno_fw,it_corrente)
!Time integration for density (continuity equation)
       do ii = 1,indarrayFlu
          npi = Array_Flu(ii)
! Gamma=1 for particles in the inner domain
          if (nPartIntorno_fw(npi) == 0) pg(npi)%Gamma = one         
!AA406test (Gamma, uni o sigma? o niente?)
          if (Domain%tipo == "bsph") pg(npi)%dden = pg(npi)%dden / pg(npi)%Gamma
! non serve     if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!Boundary type is fixe or tapis or level(?)
          if (pg(npi)%koddens == 0) then
!AA406test
!       pg(npi)%dens = pg(npi)%dens + dt * pg(npi)%dden
!AA501test
!        pg(npi)%DensShep = pg(npi)%DensShep + pg(npi)%rhoSPH_new - pg(npi)%rhoSPH_old
             if (pg(npi)%FS == 1) then
!AA501test
                pg(npi)%dens = pg(npi)%rhoSPH_new / pg(npi)%sigma
!           pg(npi)%dens = pg(npi)%DensShep / pg(npi)%sigma
                if (pg(npi)%dens < (0.98*med(1)%den0)) pg(npi)%dens = 0.98*med(1)%den0 
                else
!AA501test
                   pg(npi)%dens = pg(npi)%rhoSPH_new / pg(npi)%Gamma
!         pg(npi)%dens = pg(npi)%DensShep / pg(npi)%Gamma
             endif
!AA406!!!test
!       beta = exp(-30000.*(min((pg(npi)%sigma/pg(npi)%Gamma),1.)-1.)**2)
!       pg(npi)%dens =  pg(npi)%rhoSPH_new / (beta*pg(npi)%Gamma+(1.-beta)*pg(npi)%sigma)
             pg(npi)%densass = zero
             else if (pg(npi)%koddens == 2) then
                pg(npi)%dens = pg(npi)%densass  
          end if
! Impose boundary conditions at inlet and outlet sections (DB-SPH)
          call DBSPH_inlet_outlet(npi)
       end do
!$omp end parallel do
       deallocate (uni_old)
       else
!$omp parallel do default(none) private(npi,ii) shared(nag,pg,dt,indarrayFlu,Array_Flu)
!Time integration for density (continuity equation)
          do ii = 1,indarrayFlu
             npi = Array_Flu(ii)
! non serve     if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!Boundary type is fixe or tapis or level(?)
             if (pg(npi)%koddens == 0) then
                pg(npi)%dens = pg(npi)%dens + dt * pg(npi)%dden
                pg(npi)%densass = zero
                else if (pg(npi)%koddens == 2) then
                   pg(npi)%dens = pg(npi)%densass  
             end if
          end do
!$omp end parallel do
    endif
!AA601 sub end

  return
  end subroutine Euler
!---split

!cfile Heun.f90
!AA402 all the subroutine
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : Heun
!
! Creation date : 30 June 2011 (A. Amicarelli, G. Agate) (AA402)
!
! Last updating : September 20, 2012
!
! Improvement traceback:
!
!************************************************************************************
! Module purpose : Heun (RK2) time integration (second stage)
!
! Calling routine: time_integration
!
! Called routines: start_and_stop, inter_SmoothVelo_2D, inter_SmoothVelo_3D
!
!************************************************************************************
!
  subroutine Heun  
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
!Scalar declaration
  integer(4)       :: npi,ii
  double precision :: TetaV1
!
!.. Executable Statements ..
!
!$omp parallel do default(none) private(npi,ii) shared(nag,Pg,dt,indarrayFlu,Array_Flu,ts0_pg,esplosione)
!
! velocity and energy integration
   do ii = 1,indarrayFlu
     npi = Array_Flu(ii)
!.. kodvel = 0: the particle is internal to the domain. All the velocity components are set using the 
!..             series development, being dv/dt = acceleration components
     if (pg(npi)%kodvel == 0) then 
       pg(npi)%vel(:) = ts0_pg(npi)%ts_vel(:) + half*dt* (pg(npi)%acc(:)+ts0_pg(npi)%ts_acc(:))
!.. kodvel = 1: the particle has a critical flux condition. The horizontal velocity components are set using the 
!..             series development, while the vertical one is assigned
     else if (pg(npi)%kodvel == 1) then            
       pg(npi)%vel(:) = ts0_pg(npi)%ts_vel(:) + half*dt* (pg(npi)%acc(:)+ts0_pg(npi)%ts_acc(:))    
       pg(npi)%vel(3) = pg(npi)%velass(3)           
!.. kodvel = 2: the particle has an assigned normal velocity or source condition. All the velocity components 
!..             are set to the assigned normal components values
     else if (pg(npi)%kodvel == 2) then  
       pg(npi)%vel(:) = pg(npi)%velass(:)            
     end if
! energy time integration (explosion option)
     if (esplosione) then
       pg(npi)%IntEn = ts0_pg(npi)%ts_IntEn + half*dt* (pg(npi)%dEdT+ts0_pg(npi)%ts_dEdT)
     end if
   end do
!
!$omp end parallel do
!
   call start_and_stop(3,17)
!
! velocity and energy smoothing
   call start_and_stop(2,7)
   if (ncord==2) then 
     call inter_SmoothVelo_2D
   else
     call inter_SmoothVelo_3D
   end if
!
!$omp parallel do default(none) &
!$omp private(npi,ii,TetaV1) &
!$omp shared(nag,Pg,Med,Domain,dt,indarrayFlu,Array_Flu,esplosione)
!
   do ii = 1,indarrayFlu
     npi = Array_Flu(ii)
     if (esplosione) then
!.. calcolando TetaV1 con Csound si ottengono gli stessi risultati di Celerita
       TetaV1 = Domain%TetaV * pg(npi)%Csound * dt / Domain%h
     else
!.. update TetaV aggiorno TetaV adeguato al passo temporale
       TetaV1 = Domain%TetaV * Med(pg(npi)%imed)%Celerita * dt / Domain%h
     end if
!
!.. update of Specific Internal Energy
     if (esplosione) pg(npi)%IntEn = pg(npi)%IntEn + TetaV1 * pg(npi)%Envar
!
!.. the particle is inside the domain and far from the boundaries
     if (pg(npi)%kodvel == 0) then        
!.. evaluates the average velocity components
!.. upgrade the velocity component array with the average velocity values previously evaluated
       pg(npi)%var(:) = pg(npi)%vel(:) + TetaV1 * pg(npi)%var(:)      
       pg(npi)%vel(:) = pg(npi)%var(:)                                        
     else  
!.. the particle is close to a "source", "level" or "normal velocity boundary (kodvel = 1 or = 2)
!.. the final velocity is kept unmodified
       pg(npi)%var(:) = pg(npi)%vel(:) 
     end if
   end do
!
!$omp end parallel do
!
   call start_and_stop(3,7)
!
   call start_and_stop(2,17)
!
!$omp parallel do default(none) private(npi) shared(nag,Pg,dt,ts0_pg)
!
!.. time integration for position and density (trajectories and continuity equation)
   do npi = 1,nag
     if (pg(npi)%cella == 0) cycle
!.. store the old coordinates of the particles (not here anymore, done at the first stage only)
     if (pg(npi)%vel_type /= "std") then
!.. the type of the velocity condition is not "std", the assigned velocity condition is used (it is not necessary, but it formally helps)
       pg(npi)%coord(:) = ts0_pg(npi)%ts_coord(:) + half*dt* (pg(npi)%vstart(:)+pg(npi)%vstart(:))
     else
!.. otherwise, the smoothed value V^ already calculated are used to define the particle path
       pg(npi)%coord(:) = ts0_pg(npi)%ts_coord(:) + half*dt* (ts0_pg(npi)%ts_var(:)+pg(npi)%var(:))  
     end if
   end do
!
!$omp end parallel do
!
!$omp parallel do default(none) private(npi,ii) shared(nag,pg,dt,indarrayFlu,Array_Flu,ts0_pg)
!
!.. time integration for density (continuity equation)
   do ii = 1,indarrayFlu
     npi = Array_Flu(ii)
     if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!.. boundary type is fixed or tapis or level(?)
     if (pg(npi)%koddens == 0) then
       pg(npi)%dens = ts0_pg(npi)%ts_dens + half*dt* (ts0_pg(npi)%ts_dden+pg(npi)%dden)
       pg(npi)%densass = zero
     else if (pg(npi)%koddens == 2) then
       pg(npi)%dens = pg(npi)%densass  ! Viene mantenuta costante la densita'
     end if
   end do
!
!$omp end parallel do
!
  return
  end subroutine Heun
!---split

!cfile time_integration_body_dynamics.f90
!AA501b the whole subroutine 
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : time_integration_body_dynamics
!
! Creation      : 12Jul12 (Amicarelli-Agate)
!
! Last updating : September 20, 2012
!
!************************************************************************************
! Module purpose : Euler time integration for body dynamics
!
! Calling routines: Loop_Irre_2D,Loop_Irre_3D
!
! Called routines: /
!
!************************************************************************************

  subroutine time_integration_body_dynamics(dtvel)  

! Used modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE

! Declarations 
  implicit none
  double precision, intent(in) :: dtvel
  integer(4)       :: i,npi,j
  double precision :: mod_normal,aux_umax
  double precision :: vec_temp(3),vec2_temp(3),vec3_temp(3),domega_dt(3),aux_vel(3)
  double precision,allocatable,dimension(:,:) :: teta

!Allocations
  allocate(teta(n_bodies,3))

! Statements 

! Loop over bodies (body dynamics)
!$omp parallel do default(none) private(i,vec_temp,vec2_temp,vec3_temp,domega_dt,aux_umax) &
!$omp shared(n_bodies,body_arr,dt,ncord,teta,dtvel,tempo)
  do i=1,n_bodies
!staggered parameters (velocity and angular velocity)

     if (body_arr(i)%imposed_kinematics == 0 ) then
! Computed kinematics       
        if (ncord == 2) then
           body_arr(i)%Force(2) = zero
           body_arr(i)%Moment(1) = zero
           body_arr(i)%Moment(3) = zero
        endif
        body_arr(i)%u_CM(:) = body_arr(i)%u_CM(:) + body_arr(i)%Force(:) / body_arr(i)%mass * dtvel 
        if (ncord == 3) then
           vec_temp(1) = dot_product(body_arr(i)%Ic(1,:),body_arr(i)%omega)
           vec_temp(2) = dot_product(body_arr(i)%Ic(2,:),body_arr(i)%omega)
           vec_temp(3) = dot_product(body_arr(i)%Ic(3,:),body_arr(i)%omega)
           call Vector_Product(body_arr(i)%omega,vec_temp,vec3_temp,3)
           vec2_temp = body_arr(i)%Moment(:) - vec3_temp(:)
        endif
        if (ncord == 2) vec2_temp = body_arr(i)%Moment(:) 
        domega_dt(1) = dot_product(body_arr(i)%Ic_inv(1,:),vec2_temp)
        domega_dt(2) = dot_product(body_arr(i)%Ic_inv(2,:),vec2_temp)
        domega_dt(3) = dot_product(body_arr(i)%Ic_inv(3,:),vec2_temp)
        body_arr(i)%omega(:) = body_arr(i)%omega(:) + domega_dt(:) * dtvel
        if (ncord == 2) then
           body_arr(i)%x_CM(2) = zero 
           body_arr(i)%u_CM(2) = zero
           body_arr(i)%omega(1) = zero
           body_arr(i)%omega(3) = zero
        endif
        else
! Imposed kinematics
           do j=1,body_arr(i)%n_records
              if (body_arr(i)%body_kinematics(j,1)>=tempo) then
                 if (body_arr(i)%body_kinematics(j,1)==tempo) then
                    body_arr(i)%u_CM(:) = body_arr(i)%body_kinematics(j,2:4)
                    body_arr(i)%omega(:) = body_arr(i)%body_kinematics(j,5:7)
                    else
                    body_arr(i)%u_CM(:) = body_arr(i)%body_kinematics(j-1,2:4) + &
                                         (body_arr(i)%body_kinematics(j,2:4)-body_arr(i)%body_kinematics(j-1,2:4))/ &
                                         (body_arr(i)%body_kinematics(j,1)-body_arr(i)%body_kinematics(j-1,1)) * &
                                         (tempo-body_arr(i)%body_kinematics(j-1,1))
                    body_arr(i)%omega(:) = body_arr(i)%body_kinematics(j-1,5:7) + &
                                          (body_arr(i)%body_kinematics(j,5:7)-body_arr(i)%body_kinematics(j-1,5:7))/ &
                                         (body_arr(i)%body_kinematics(j,1)-body_arr(i)%body_kinematics(j-1,1)) * &
                                         (tempo-body_arr(i)%body_kinematics(j-1,1))                                      
                 endif
                 exit
              endif
           enddo   
     endif
!     
     
!non-staggered parameters     
     body_arr(i)%x_CM(:) = body_arr(i)%x_CM(:) + body_arr(i)%u_CM(:) * dt 
     teta(i,:) =  body_arr(i)%omega(:) * dt 
     body_arr(i)%alfa(:) = body_arr(i)%alfa(:) + teta(i,:) 
! Initializing umax (maximum particle velocity of the body)
     body_arr(i)%umax = zero 
  end do
!$omp end parallel do

! Loop over body particles (static kinematics)
!$omp parallel do default(none) private(npi,vec_temp,mod_normal,vec2_temp,aux_vel) &
!$omp shared(n_body_part,body_arr,bp_arr,dt,ncord,teta,dtvel)
  do npi=1,n_body_part
! staggered parameter  
     call Vector_Product(body_arr(bp_arr(npi)%body)%omega,bp_arr(npi)%rel_pos,vec_temp,3)
     aux_vel(:) = bp_arr(npi)%vel(:) 
     bp_arr(npi)%vel(:) = body_arr(bp_arr(npi)%body)%u_CM(:) + vec_temp(:)
     if (ncord == 2) bp_arr(npi)%vel(2) = zero 
     bp_arr(npi)%acc(:) = (bp_arr(npi)%vel(:)-aux_vel(:))/dtvel
! non-staggered parameters     
     vec2_temp(:) = teta(bp_arr(npi)%body,:)
     call vector_rotation(bp_arr(npi)%rel_pos,vec2_temp)
     if (ncord == 2) bp_arr(npi)%rel_pos(2) = zero     
     bp_arr(npi)%pos(:) = bp_arr(npi)%rel_pos(:) + body_arr(bp_arr(npi)%body)%x_CM(:)
     call vector_rotation(bp_arr(npi)%normal,vec2_temp)
     mod_normal = dsqrt(dot_product(bp_arr(npi)%normal,bp_arr(npi)%normal))
     if (mod_normal > one) bp_arr(npi)%normal(:) = bp_arr(npi)%normal(:) / mod_normal    
     if (ncord == 2) then
        bp_arr(npi)%rel_pos(2) = zero
        bp_arr(npi)%normal(2) = zero
     endif 
  end do
!$omp end parallel do
!

!Updating max velocity within every body
  do npi=1,n_body_part
     aux_umax = dsqrt(dot_product(bp_arr(npi)%vel,bp_arr(npi)%vel))        
     body_arr(bp_arr(npi)%body)%umax = max(body_arr(bp_arr(npi)%body)%umax,aux_umax)   
  end do

! Part good for RK1, not for RK1stag
!! Loop over body particles (kinematics)
!!$omp parallel do default(none) private(npi) shared(n_body_part,body_arr,bp_arr,dt,ncord)
!  do npi=1,n_body_part
!     bp_arr(npi)%pos(:) = bp_arr(npi)%pos(:) + bp_arr(npi)%vel(:) * dt
!     if (ncord == 2) bp_arr(npi)%pos(2) = zero 
!  end do
!!$omp end parallel do

!! Loop over bodies (body dynamics)
!!$omp parallel do default(none) private(i,vec_temp,vec2_temp,vec3_temp,domega_dt,aux_umax) shared(n_bodies,body_arr,dt,ncord,teta)
!  do i=1,n_bodies
!     if (ncord == 2) then
!        body_arr(i)%Force(2) = zero
!        body_arr(i)%Moment(1) = zero
!        body_arr(i)%Moment(3) = zero
!     endif
!     body_arr(i)%x_CM(:) = body_arr(i)%x_CM(:) + body_arr(i)%u_CM(:) * dt 
!     teta(i,:) =  body_arr(i)%omega(:) * dt 
!     body_arr(i)%alfa(:) = body_arr(i)%alfa(:) + teta(i,:) 
!     body_arr(i)%u_CM(:) = body_arr(i)%u_CM(:) + body_arr(i)%Force(:) / body_arr(i)%mass * dt 
!     if (ncord == 3) then
!        vec_temp(1) = dot_product(body_arr(i)%Ic(1,:),body_arr(i)%omega)
!        vec_temp(2) = dot_product(body_arr(i)%Ic(2,:),body_arr(i)%omega)
!        vec_temp(3) = dot_product(body_arr(i)%Ic(3,:),body_arr(i)%omega)
!        call Vector_Product(body_arr(i)%omega,vec_temp,vec3_temp,3)
!        vec2_temp = body_arr(i)%Moment(:) - vec3_temp(:)
!     endif
!     if (ncord == 2) vec2_temp = body_arr(i)%Moment(:) 
!     domega_dt(1) = dot_product(body_arr(i)%Ic_inv(1,:),vec2_temp)
!     domega_dt(2) = dot_product(body_arr(i)%Ic_inv(2,:),vec2_temp)
!     domega_dt(3) = dot_product(body_arr(i)%Ic_inv(3,:),vec2_temp)
!     body_arr(i)%omega(:) = body_arr(i)%omega(:) + domega_dt(:) * dt
!     if (ncord == 2) then
!        body_arr(i)%x_CM(2) = zero 
!        body_arr(i)%u_CM(2) = zero
!        body_arr(i)%omega(1) = zero
!        body_arr(i)%omega(3) = zero
!     endif
!! Initializing umax (maximum particle velocity of the body)
!     body_arr(i)%umax = zero 
!  end do
!!$omp end parallel do
!
!! Loop over body particles (static kinematics)
!!$omp parallel do default(none) private(npi,vec_temp,mod_normal,vec2_temp,aux_umax) shared(n_body_part,body_arr,bp_arr,dt,ncord,teta)
!  do npi=1,n_body_part
!     call Vector_Product(body_arr(bp_arr(npi)%body)%omega,bp_arr(npi)%rel_pos,vec_temp,3)
!     bp_arr(npi)%vel(:) = body_arr(bp_arr(npi)%body)%u_CM(:) + vec_temp(:)
!     bp_arr(npi)%rel_pos(:) = bp_arr(npi)%pos(:) - body_arr(bp_arr(npi)%body)%x_CM(:)
!     vec2_temp(:) = teta(bp_arr(npi)%body,:)
!     call vector_rotation(bp_arr(npi)%normal,vec2_temp)
!     mod_normal = dsqrt(dot_product(bp_arr(npi)%normal,bp_arr(npi)%normal))
!     if (mod_normal > 1.) bp_arr(npi)%normal(:) = bp_arr(npi)%normal(:) / mod_normal    
!     if (ncord == 2) then
!        bp_arr(npi)%vel(2) = zero 
!        bp_arr(npi)%rel_pos(2) = zero
!        bp_arr(npi)%normal(2) = zero 
!     endif 
!! Update of umax
!     aux_umax = dsqrt(dot_product(bp_arr(npi)%vel,bp_arr(npi)%vel))        
!     body_arr(bp_arr(npi)%body)%umax = max(body_arr(bp_arr(npi)%body)%umax,aux_umax)   
!  end do
!!$omp end parallel do

!Deallocations
  deallocate(teta)

  return
  end subroutine time_integration_body_dynamics
!---split