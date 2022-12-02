!-------------------------------------------------------------------------------
! SPHERA v.10.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
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
! Program unit: time_step_post_processing
! Description: Post processing each time step
!-------------------------------------------------------------------------------
subroutine time_step_post_processing(it_print,it_memo,it_rest,it,dtvel)
!------------------------
! Modules
!------------------------
use Dynamic_allocation_module
use Static_allocation_module
use Hybrid_allocation_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: it
double precision,intent(in) :: dtvel
integer(4),intent(inout) :: it_print,it_memo,it_rest
integer(4) :: npi
double precision :: xmax,ymax,zmax
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
#ifdef SPACE_3D
if (allocated(Z_fluid_max)) then
   if ((int(simulation_time/input_any_t%depth_dt_out)>Domain%depth_it_out_last)&
      .or.(on_going_time_step==1)) then
      Domain%depth_it_out_last = int(simulation_time / input_any_t%depth_dt_out)
      call Update_Zmax_at_grid_vert_columns(1)
      else
         call Update_Zmax_at_grid_vert_columns(0)
   endif
endif
#endif
if (ulog>0) then
   call Print_Results(it_eff,it_print,'loop__')
endif
if (nres>0) then
   call Memo_Results(it_eff,it_memo,it_rest,dtvel,'loop__')
endif
! Post-processing time series for: monitoring points, monitoring lines (but the 
! free surface), body dynamics, DB-SPH elements.
call CalcVarp
if (Domain%icpoi_fr>0) then
   if ((mod(it,Domain%icpoi_fr)==0).and.npointst>0) then
      call Memo_Ctl
   endif
#ifdef SOLID_BODIES
      if (mod(it,Domain%icpoi_fr)==0) call Body_dynamics_output
#endif
   if (mod(it,Domain%icpoi_fr)==0) call fluid_global_quantities
   elseif (input_any_t%cpoi_fr>zero) then
      if ((dmod(simulation_time,input_any_t%cpoi_fr)<=dtvel).and.npointst>0)   &
         call Memo_Ctl
#ifdef SOLID_BODIES
         if (dmod(simulation_time,input_any_t%cpoi_fr)<=dtvel) then
            call Body_dynamics_output
         endif         
#endif
      if (dmod(simulation_time,input_any_t%cpoi_fr)<=dtvel) then
         call fluid_global_quantities
      endif
! DB-SPH post-processing and update of the variable "wet" in the array "pg"
      if ((Domain%tipo=="bsph").and.                                           &
         (dmod(simulation_time,input_any_t%cpoi_fr)<=dtvel)) then
         if ((DBSPH%n_monitor_points>0).or.(DBSPH%n_monitor_regions==1)) then
            call wall_elements_pp
         endif
!$omp parallel do default(none)                                                &
!$omp shared(DBSPH,pg_w)                                                       &
!$omp private(npi)
         do npi=1,DBSPH%n_w
            pg_w(npi)%wet = 0 
         enddo    
!$omp end parallel do
      endif
endif
! Post-processing for the time series of the free-surface
if (Domain%ipllb_fr>0) then
   if ((mod(it,Domain%ipllb_fr)==0).and.nlines>0) then
      call calc_pelo
   endif
   elseif (input_any_t%pllb_fr>zero) then
      if ((dmod(simulation_time,input_any_t%pllb_fr)<=dtvel).and.nlines>0) then
         call calc_pelo
      endif
endif
if (Granular_flows_options%monitoring_lines>0) then
   if ((int(simulation_time/Granular_flows_options%dt_out)>                    &
      Granular_flows_options%it_out_last).or.(on_going_time_step==1)) then
      Granular_flows_options%it_out_last = int(simulation_time /               &
         Granular_flows_options%dt_out)
      call interface_post_processing
   endif
endif
#ifdef SPACE_3D
if (Q_sections%n_sect>0) then
   if ((int(simulation_time/Q_sections%dt_out)>Q_sections%it_out_last).or.     &
      (on_going_time_step==1)) then
      Q_sections%it_out_last = int(simulation_time / Q_sections%dt_out)
      call sub_Q_sections
   endif
endif
if (substations%n_sub>0) then
   if ((int(simulation_time/substations%dt_out)>substations%it_out_last).or.   &
      (on_going_time_step==1)) then
      substations%it_out_last = int(simulation_time / substations%dt_out)
      call electrical_substations
   endif
endif
#endif
! Post-processing for the water front
if (Domain%imemo_fr>0) then
   if (mod(it,Domain%imemo_fr)==0) then
      xmax = -1.d30
      ymax = -1.d30
      zmax = -1.d30
      do npi=1,nag
         if ((pg(npi)%vel_type/="std").or.(pg(npi)%cella==0)) cycle
         xmax = max(xmax,pg(npi)%coord(1))
#ifdef SPACE_3D
            ymax = max(ymax,pg(npi)%coord(2))
            zmax = max(zmax,pg(npi)%coord(3))
#elif defined SPACE_2D
               ymax = max(ymax,pg(npi)%coord(3))
#endif
      enddo
#ifdef SPACE_3D
         write(nfro,'(4g14.7)') simulation_time,xmax,ymax,zmax
#elif defined SPACE_2D
            write(nfro,'(2g14.7,13x,a,g14.7)') simulation_time,xmax,'-',ymax
#endif
   endif
   elseif (input_any_t%memo_fr>zero) then
      if (it>1.and.dmod(simulation_time,input_any_t%memo_fr)<=dtvel) then
         xmax = - 1.d30
         ymax = - 1.d30
         zmax = - 1.d30
         do npi=1,nag
            if (pg(npi)%vel_type/="std".or.pg(npi)%cella==0) cycle
            xmax = max(xmax,pg(npi)%coord(1))
#ifdef SPACE_3D
               ymax = max(ymax,pg(npi)%coord(2))
               zmax = max(zmax,pg(npi)%coord(3))
#elif defined SPACE_2D
                  ymax = max(ymax,pg(npi)%coord(3))
#endif
         enddo
#ifdef SPACE_3D
            write(nfro,'(4g14.7)') simulation_time,xmax,ymax,zmax
#elif defined SPACE_2D
               write(nfro,'(2g14.7,13x,a,g14.7)') simulation_time,xmax,'-',ymax
#endif
      endif
endif
! Paraview output and .txt file concatenation
if (vtkconv) then
   call result_converter('loop__')
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine time_step_post_processing
