!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-Ricerca di Sistema)



! SPHERA authors and email contact are provided on SPHERA documentation.

! This file is part of SPHERA v.8.0.
! SPHERA v.8.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with SPHERA. If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: Euler                                           
! Description: Explicit RK1 time integration scheme (Euler scheme).  
!----------------------------------------------------------------------------------------------------------------------------------

subroutine Euler  
!------------------------
! Modules
!------------------------ 
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: npi,ii,i,j
double precision :: TetaV1,beta
double precision,dimension(:),allocatable :: uni_old
integer(4),external :: ParticleCellNumber
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
!$omp parallel do default(none)                                                &
!$omp private(npi,ii)                                                          &
!$omp shared(nag,Pg,Domain,dt,indarrayFlu,Array_Flu,ts0_pg,esplosione)
! Time integration for momentum and energy equations
do ii = 1,indarrayFlu
   npi = Array_Flu(ii)
! To save the results of the previous step for the time stage parameters "ts0" 
! (only during the first stage of Heun's scheme)
! The coordinates and the smoothed velocities are saved during a following 
! loop, even involving solid particles
   if (Domain%RKscheme==2) then
      ts0_pg(npi)%ts_vel(:) = pg(npi)%vel(:)  
      ts0_pg(npi)%ts_acc(:) = pg(npi)%acc(:)                        
      ts0_pg(npi)%ts_dden = pg(npi)%dden
      ts0_pg(npi)%ts_dens = pg(npi)%dens
      if (esplosione) then
         ts0_pg(npi)%ts_IntEn = pg(npi)%IntEn
         ts0_pg(npi)%ts_dEdT = pg(npi)%dEdT
      endif
   endif
! kodvel = 0: the particle is internal to the domain. 
   if (pg(npi)%kodvel==0) then 
      pg(npi)%vel(:) = pg(npi)%vel(:) + dt * pg(npi)%acc(:)
! Impose boundary conditions at inlet and outlet sections (DB-SPH)
      if (Domain%tipo=="bsph") call DBSPH_inlet_outlet(npi)
! kodvel = 1: the particle has a critical flux condition. Vertical velocity is 
! assigned.
      elseif (pg(npi)%kodvel==1) then 
         pg(npi)%vel(:) = pg(npi)%vel(:) + dt * pg(npi)%acc(:)      
         pg(npi)%vel(3) = pg(npi)%velass(3)           
! kodvel = 2: the particle has an assigned normal velocity (inlet section). 
! Velocity is imposed. 
         elseif (pg(npi)%kodvel==2) then  
            pg(npi)%vel(:) = pg(npi)%velass(:)            
   endif
   if (esplosione) then
      pg(npi)%IntEn = pg(npi)%IntEn + dt * pg(npi)%dEdT
   endif
enddo
!$omp end parallel do
call start_and_stop(3,17)
! Velocity and energy partial smoothing
call start_and_stop(2,7)
if (ncord==2) then 
   call inter_SmoothVelo_2D
   else
      call inter_SmoothVelo_3D
endif
!$omp parallel do default(none)                                                &
!$omp private(npi,ii,TetaV1)                                                   &
!$omp shared(nag,Pg,Med,Domain,dt,indarrayFlu,Array_Flu,esplosione)
do ii=1,indarrayFlu
   npi = Array_Flu(ii)
   if (esplosione) then
      TetaV1 = Domain%TetaV * pg(npi)%Csound * dt / Domain%h
      else
! To update "TetaV", depending on dt.
         TetaV1 = Domain%TetaV * Med(pg(npi)%imed)%Celerita * dt / Domain%h
   endif
   if (esplosione) pg(npi)%IntEn = pg(npi)%IntEn + TetaV1 * pg(npi)%Envar
   if (pg(npi)%kodvel==0) then        
! Particle is inside the domain and far from the boundaries.
! Velocity is partially smoothed.
      pg(npi)%var(:) = pg(npi)%vel(:) + TetaV1 * pg(npi)%var(:)      
      pg(npi)%vel(:) = pg(npi)%var(:)                                        
      else  
! The particle is close to a boundary of type: "source", "level" or 
! "normal velocity boundary" (kodvel = 1 or = 2). No partial smoothing for 
! velocity.
         pg(npi)%var(:) = pg(npi)%vel(:) 
   endif
enddo
!$omp end parallel do
call start_and_stop(3,7)
call start_and_stop(2,17)
!$omp parallel do default(none)                                                &
!$omp private(npi)                                                             &
!$omp shared(nag,Pg,Domain,dt,ts0_pg,it_corrente)
! Time integration for position and density (trajectories and continuity
!  equation)
do npi=1,nag
   if (pg(npi)%cella==0) cycle
! To save the coordinates of the previous time step (RK2) for time integration
   if (Domain%RKscheme==2) then
      ts0_pg(npi)%ts_coord(:) = pg(npi)%coord(:)
      ts0_pg(npi)%ts_var(:) = pg(npi)%var(:)
   endif
! To save the old particle coordinates 
   pg(npi)%CoordOld(:) = pg(npi)%coord(:)
   if (pg(npi)%vel_type /= "std") then
! With imposed velocity (BC)
      pg(npi)%coord(:) = pg(npi)%coord(:) + dt * pg(npi)%vstart(:)
      else
! With computed velocity 
         pg(npi)%coord(:) = pg(npi)%coord(:) + dt * pg(npi)%var(:) 
   endif
enddo
!$omp end parallel do
! Wall element trajectories
if (((DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet)>0).and.(Domain%tipo=="bsph")) &
   then
! This commented part should be useful for further developments
!       do npi = 1,n_w
!          if (pg_w(npi)%cella==0) cycle
!          pg_w(npi)%coord(:) = pg_w(npi)%coord(:) + dt * pg_w(npi)%vel(:)
!       enddo
   if (DBSPH%n_w>0) then
      call DBSPH_kinematics
! Time integration for the first surface element velocity    
! Loop over the surface element to copy and paste, all over the elements,
! the imposed kinematics of the first surface element 
      pg_w(1)%coord(:) = pg_w(1)%coord(:) + dt * pg_w(1)%vel(:)
!$omp parallel do default(none) shared(DBSPH,pg_w,dt) private(npi)
      do npi=2,DBSPH%n_w
         pg_w(npi)%vel(:) = pg_w(1)%vel(:) 
         pg_w(npi)%coord(:) = pg_w(npi)%coord(:) + dt * pg_w(1)%vel(:) 
      enddo
!$omp end parallel do
   endif
! In-built motion of monitoring lines
   if (DBSPH%in_built_monitors.eqv..true.) then
      do i=1,nlines
! Loop over the line points
         do j=control_lines(i)%Icont(1),control_lines(i)%Icont(2) 
            control_points(j)%coord(:) = control_points(j)%coord(:) + dt *  &
               pg_w(1)%vel(:)
            control_points(j)%cella =                                       &
               ParticleCellNumber(control_points(j)%coord(:))
         enddo
      enddo
   endif
! Boundary conditions: start 
! BC (checks for the particles gone out of the domain throughout the opened
! sides)
   if ((Domain%time_split==0).and.(Domain%time_stage==Domain%RKscheme)) then
      call start_and_stop(2,9)
      if (ncord==2) then 
         if (NumOpenSides>0) call CancelOutgoneParticles_2D
! Adds new particles at the inlet sections
         if (SourceSide /= 0) call GenerateSourceParticles_2D 
         else
            if (NumOpenFaces>0) call CancelOutgoneParticles_3D
! Adds new particles at the inlet sections
            if (SourceFace/=0) call GenerateSourceParticles_3D 
      endif
! Particle reordering on the backround positioning grid 
      call OrdGrid1 (nout)
      call start_and_stop(3,9)
! To set the parameters for the fixed particles 
      if (Domain%NormFix)  call NormFix
      allocate (uni_old(nag))
! To update the auxiliary vector to count mobile "flu" particles (in "nag") 
      indarrayFlu = 0
      do npi=1,nag
         if (pg(npi)%rhoSPH_new /= zero) then
            uni_old(npi) = pg(npi)%uni
            else
               uni_old(npi) = zero
         endif
         if ((pg(npi)%cella==0).or.(pg(npi)%vel_type/="std")) cycle
         indarrayFlu = indarrayFlu + 1
         Array_Flu(indarrayFlu) = npi
      enddo 
   endif
! Boundary conditions: end
   do npi=1,nag 
      pg(npi)%Gamma = pg(npi)%Gamma + pg(npi)%dShep * dt
      if (pg(npi)%Gamma>=one) pg(npi)%Gamma = one 
   enddo
   call CalcVarLength
!$omp parallel do default(none)                                                &
!$omp private(npi,ii)                                                          &
!$omp shared(nag,pg,med,dt,indarrayFlu,Array_Flu,Domain,uni_old)               &
!$omp shared(beta,nPartIntorno_fw,it_corrente)
! Time integration of the continuity equation
   do ii=1,indarrayFlu
      npi = Array_Flu(ii)
! Gamma=1 for particles in the inner domain
      if (nPartIntorno_fw(npi)==0) pg(npi)%Gamma = one         
      if (Domain%tipo=="bsph") pg(npi)%dden = pg(npi)%dden / pg(npi)%Gamma
! Boundary type is fixe or tapis or level(?)
      if (pg(npi)%koddens==0) then
         if (pg(npi)%FS==1) then
            pg(npi)%dens = pg(npi)%rhoSPH_new / pg(npi)%sigma
            if (pg(npi)%dens<(0.98d0*med(1)%den0)) pg(npi)%dens = 0.98d0 *  &
               med(1)%den0 
            else
               pg(npi)%dens = pg(npi)%rhoSPH_new / pg(npi)%Gamma
         endif
! Interesting test, according to Ferrand et al. (2013)
! beta = exp(-30000.*(min((pg(npi)%sigma/pg(npi)%Gamma),1.)-1.)**2)
! pg(npi)%dens =  pg(npi)%rhoSPH_new / (beta*pg(npi)%Gamma+(1.-beta)
!    *pg(npi)%sigma)
         pg(npi)%densass = zero
         elseif (pg(npi)%koddens==2) then
            pg(npi)%dens = pg(npi)%densass  
      endif
! Impose boundary conditions at inlet and outlet sections (DB-SPH)
      call DBSPH_inlet_outlet(npi)
   enddo
!$omp end parallel do
   deallocate(uni_old)
   else
!$omp parallel do default(none)                                                &
!$omp private(npi,ii)                                                          &
!$omp shared(nag,pg,dt,indarrayFlu,Array_Flu)
! Time integration for density (continuity equation)
      do ii=1,indarrayFlu
         npi = Array_Flu(ii)
! Boundary type is fixe or tapis or level(?)
         if (pg(npi)%koddens==0) then
            pg(npi)%dens = pg(npi)%dens + dt * pg(npi)%dden
            pg(npi)%densass = zero
            elseif (pg(npi)%koddens==2) then
               pg(npi)%dens = pg(npi)%densass  
         endif
      enddo
!$omp end parallel do
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine Euler

