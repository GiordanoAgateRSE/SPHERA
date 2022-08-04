!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
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
! Program unit: cat_post_proc                  
! Description: To concatenate the ".txt" output files and remove the original 
!              ones.              
!-------------------------------------------------------------------------------
subroutine cat_post_proc
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
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
#ifdef SOLID_BODIES
   call system("touch Body_dynamics.txt draft_Body_dynamics_draft.txt")
   call system("cat Body_dynamics.txt *Body_dynamics_* > Body-dynamics.txt")
   call system("rm -f *Body_dynamics*")
   call system("mv Body-dynamics.txt Body_dynamics.txt")
   call system("touch Body_particles.txt draft_Body_particles_draft.txt")
   call system("cat Body_particles.txt *Body_particles_* > Body-particles.txt")
   call system("rm -f *Body_particles*")
   call system("mv Body-particles.txt Body_particles.txt")
#endif
#ifdef SPACE_3D
if (Q_sections%n_sect>0) then
   call system("touch Q_sections.txt draft_Q_sections_draft.txt")
   call system("cat Q_sections.txt *Q_sections_* > Q-sections.txt")
   call system("rm -f *Q_sections*")
   call system("mv Q-sections.txt Q_sections.txt")
endif
if (substations%n_sub>0) then
   call system("touch substations.txt draft_substations_draft.txt")
   call system("cat substations.txt *substations_* > sub-stations.txt")
   call system("rm -f *substations*")
   call system("mv sub-stations.txt substations.txt")
endif
#endif
if (nlines>0) then
    call system("touch monitoring_lines.txt draft.cln")
    call system("cat monitoring_lines.txt *.cln > monitoring-lines.txt")
    call system("rm -f monitoring_lines.txt *.cln")
    call system("mv monitoring-lines.txt monitoring_lines.txt")
endif
if (npointst>0) then
   call system("touch monitoring_points.txt draft.cpt")
   call system("cat monitoring_points.txt *.cpt > monitoring-points.txt")
   call system("rm -f monitoring_points.txt *.cpt")
   call system("mv monitoring-points.txt monitoring_points.txt")
endif
if (Granular_flows_options%monitoring_lines>0) then
   call system("touch blt_interfaces.txt draft_blt_interfaces_draft.txt")
   call system("cat blt_interfaces.txt *blt_interfaces_* > blt-interfaces.txt")
   call system("rm -f *blt_interfaces*")
   call system("mv blt-interfaces.txt blt_interfaces.txt")
endif
if (DBSPH%n_w>0) then
   if (DBSPH%n_monitor_points>0) then
      call system("touch wall_IDs.txt draft_wall_IDs_draft.txt")
      call system("cat wall_IDs.txt *wall_IDs_* > wall-IDs.txt")
      call system("rm -f *wall_IDs_*")
      call system("mv wall-IDs.txt wall_IDs.txt")
   endif
   if (DBSPH%n_monitor_regions>0) then
      call system("touch wall_regions.txt draft_wall_regions_draft.txt")
      call system("cat wall_regions.txt *wall_regions_* > wall-regions.txt")
      call system("rm -f *wall_regions_*")
      call system("mv wall-regions.txt wall_regions.txt")
      call system("touch wall_Fx.txt draft_wall_Fx_draft.txt")
      call system("cat wall_Fx.txt *wall_Fx_* > wall-Fx.txt")
      call system("rm -f *wall_Fx_*")
      call system("mv wall-Fx.txt wall_Fx.txt")
   endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine cat_post_proc
