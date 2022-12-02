!-------------------------------------------------------------------------------
! SPHERA v.10.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)

! SPHERA authors and email contact are provided on SPHERA documentation.

! This file is part of SPHERA v.10.0.0
! SPHERA v.10.0.0 is free software: you can redistribute it and/or modify
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
! Program unit: Heun                                           
! Description: Heun scheme: explicit RK2 time integration scheme.  
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine Heun(BC_zmax_flag)
#elif defined SPACE_2D
subroutine Heun
#endif
!------------------------
! Modules
!------------------------
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use Memory_I_O_interface_module
!------------------------
! Declarations
!------------------------
implicit none
#ifdef SPACE_3D
logical,intent(in) :: BC_zmax_flag
#endif
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
!$omp shared(nag,Pg,dt,indarrayFlu,Array_Flu,ts0_pg)
! Time integration for for momentum equations
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
   pg(npi)%dens = ts0_pg(npi)%ts_dens + half * dt * (ts0_pg(npi)%ts_dden +     &
                  pg(npi)%dden)
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine Heun
