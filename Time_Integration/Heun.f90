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

