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
! Program unit: Heun                                           
! Description: Heun scheme: explicit RK2 time integration scheme.  
!----------------------------------------------------------------------------------------------------------------------------------

subroutine Heun  
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
integer(4) :: npi,ii
double precision :: TetaV1
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
!$omp shared(nag,Pg,dt,indarrayFlu,Array_Flu,ts0_pg,esplosione)
! Time integration for for momentum and energy equations
do ii=1,indarrayFlu
   npi = Array_Flu(ii)
! kodvel = 0: the particle is internal to the domain. 
   if (pg(npi)%kodvel==0) then 
      pg(npi)%vel(:) = ts0_pg(npi)%ts_vel(:) + half * dt * (pg(npi)%acc(:)     &
         + ts0_pg(npi)%ts_acc(:))
! kodvel = 1: the particle has a critical flux condition. Vertical velocity is 
! assigned.
      elseif (pg(npi)%kodvel==1) then            
         pg(npi)%vel(:) = ts0_pg(npi)%ts_vel(:) + half * dt *                  &
            (pg(npi)%acc(:)+ts0_pg(npi)%ts_acc(:))    
         pg(npi)%vel(3) = pg(npi)%velass(3)           
! kodvel = 2: the particle has an assigned normal velocity (inlet section). 
! Velocity is imposed. 
         elseif (pg(npi)%kodvel==2) then  
            pg(npi)%vel(:) = pg(npi)%velass(:)            
   endif
   if (esplosione) then
      pg(npi)%IntEn = ts0_pg(npi)%ts_IntEn + half * dt * (pg(npi)%dEdT +       &
         ts0_pg(npi)%ts_dEdT)
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
! Time integration to get position and density (trajectories and continuity
!  equation)
!$omp parallel do default(none) private(npi) shared(nag,Pg,dt,ts0_pg)
do npi=1,nag
   if (pg(npi)%cella==0) cycle
   if (pg(npi)%vel_type/="std") then
      pg(npi)%coord(:) = ts0_pg(npi)%ts_coord(:) + half * dt *                 &
         (pg(npi)%vstart(:) + pg(npi)%vstart(:))
      else
         pg(npi)%coord(:) = ts0_pg(npi)%ts_coord(:) + half * dt *              &
            (ts0_pg(npi)%ts_var(:) + pg(npi)%var(:))  
   endif
enddo
!$omp end parallel do
!$omp parallel do default(none)                                                &
!$omp private(npi,ii)                                                          &
!$omp shared(nag,pg,dt,indarrayFlu,Array_Flu,ts0_pg)
! Time integration of the continuity equation
do ii=1,indarrayFlu
   npi = Array_Flu(ii)
   if ((pg(npi)%cella==0).or.(pg(npi)%vel_type/="std")) cycle
   if (pg(npi)%koddens==0) then
      pg(npi)%dens = ts0_pg(npi)%ts_dens + half * dt * (ts0_pg(npi)%ts_dden    &
         + pg(npi)%dden)
      pg(npi)%densass = zero
      elseif (pg(npi)%koddens==2) then
         pg(npi)%dens = pg(npi)%densass  
   endif
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine Heun

