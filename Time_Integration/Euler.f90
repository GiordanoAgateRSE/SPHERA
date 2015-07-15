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

