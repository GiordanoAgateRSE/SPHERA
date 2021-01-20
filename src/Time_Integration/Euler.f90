!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2021 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)

! SPHERA authors and email contact are provided on SPHERA documentation.

! This file is part of SPHERA v.9.0.0
! SPHERA v.9.0.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with SPHERA. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: Euler                                           
! Description: Explicit RK1 time integration scheme (Euler scheme).  
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine Euler(BC_zmax_flag)
#elif defined SPACE_2D
subroutine Euler
#endif
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
#ifdef SPACE_3D
logical,intent(in) :: BC_zmax_flag
#endif
integer(4) :: npi,ii,i,j
double precision :: TetaV1
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
!$omp shared(nag,Pg,Domain,dt,indarrayFlu,Array_Flu,ts0_pg)
! Time integration for momentum equations
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
enddo
!$omp end parallel do
call start_and_stop(3,17)
! Velocity partial smoothing
call start_and_stop(2,7)
#ifdef SPACE_3D
if (input_any_t%TetaV>0.d0) call velocity_smoothing(BC_zmax_flag)
#elif defined SPACE_2D
if (input_any_t%TetaV>0.d0) call velocity_smoothing
#endif
!$omp parallel do default(none)                                                &
!$omp private(npi,ii,TetaV1)                                                   &
!$omp shared(nag,Pg,Med,Domain,dt,indarrayFlu,Array_Flu,input_any_t)
do ii=1,indarrayFlu
   npi = Array_Flu(ii)
   if (input_any_t%TetaV>0.d0) then
! To update "TetaV", depending on dt.
      TetaV1 = input_any_t%TetaV * Med(pg(npi)%imed)%Celerita * dt / Domain%h
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
      else
         pg(npi)%var(:) = pg(npi)%vel(:)
   endif
enddo
!$omp end parallel do
call start_and_stop(3,7)
call start_and_stop(2,17)
!$omp parallel do default(none)                                                &
!$omp private(npi)                                                             &
!$omp shared(nag,pg,Domain,dt,ts0_pg)
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
   if (pg(npi)%vel_type/="std") then
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
   call DBSPH_kinematics
! Time integration for the surface element velocity      
!$omp parallel do default(none) shared(DBSPH,pg_w,dt) private(npi)
   do npi=1,DBSPH%n_w
      pg_w(npi)%coord(:) = pg_w(npi)%coord(:) + dt * pg_w(npi)%vel(:)
   enddo
!$omp end parallel do
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
#ifdef SPACE_2D
         if (NumOpenSides>0) call CancelOutgoneParticles_2D
! Adds new particles at the inlet sections
         if (SourceSide/=0) call GenerateSourceParticles
#elif defined SPACE_3D
            if (NumOpenFaces>0) call CancelOutgoneParticles_3D
! Adds new particles at the inlet sections
            if (SourceFace/=0) call GenerateSourceParticles 
#endif
! Particle reordering on the backround positioning grid 
      call OrdGrid1
      call start_and_stop(3,9)
! To set the parameters for the fixed particles 
      if (Domain%NormFix)  call NormFix
      allocate(uni_old(nag))
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
      if ((DBSPH%Gamma_limiter_flag.eqv..true.).and.(pg(npi)%Gamma>=one))      &
         pg(npi)%Gamma = one 
   enddo
   call CalcVarLength
!$omp parallel do default(none)                                                &
!$omp private(npi,ii)                                                          &
!$omp shared(nag,pg,med,dt,indarrayFlu,Array_Flu,Domain,uni_old)               &
!$omp shared(nPartIntorno_fw,NMedium,DBSPH)
! Time integration of the continuity equation
   do ii=1,indarrayFlu
      npi = Array_Flu(ii)
      if (DBSPH%Gamma_limiter_flag.eqv..true.) then
! Gamma=1 for particles in the inner domain
         if (nPartIntorno_fw(npi)==0) pg(npi)%Gamma = one 
      endif      
      pg(npi)%dden = pg(npi)%dden / pg(npi)%Gamma
! SPH approximation of density (alternative to the continuity equation)
     if (pg(npi)%FS==0) then
        pg(npi)%dens = pg(npi)%rhoSPH_new / pg(npi)%Gamma
        else
           pg(npi)%dens = pg(npi)%rhoSPH_new / pg(npi)%sigma 
           if (pg(npi)%dens<(0.98d0*med(1)%den0)) pg(npi)%dens = 0.98d0 *      &
                                                                 med(1)%den0 
     endif
! Possible interesting test, according to Ferrand et al. (2013):
! beta = exp(-30000.*(min((pg(npi)%sigma/pg(npi)%Gamma),1.)-1.)**2) ,
! pg(npi)%dens =  pg(npi)%rhoSPH_new / (beta*pg(npi)%Gamma+(1.-beta) 
! * pg(npi)%sigma)
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
         pg(npi)%dens = pg(npi)%dens + dt * pg(npi)%dden
      enddo
!$omp end parallel do
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine Euler
