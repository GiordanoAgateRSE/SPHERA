!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2019 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: body_boundary_for_sliding_friction_normal_reaction
! Description: Body-boundary interactions. Intermediate computations for the 
!              sliding friction force and the normal reaction force under 
!              sliding. Generic "body particle - frontier" interaction. 
!-------------------------------------------------------------------------------
subroutine body_boundary_for_sliding_friction_normal_reaction(i_bp,            &
   bp_bound_interactions,normal_plane,bp_pos,interface_sliding_vel_max,        &
   mean_bound_normal,mean_bound_pos,aux_gravity)
!------------------------
! Modules
!------------------------
use Dynamic_allocation_module
use Hybrid_allocation_module
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
! ID of the current body particle
integer(4),intent(in) :: i_bp
! ID of the current neighbouring frontier
integer(4),intent(inout) :: bp_bound_interactions
double precision,intent(in) :: normal_plane(3)
! Position of the body particle
double precision,intent(in) :: bp_pos(3)
double precision,intent(inout) :: interface_sliding_vel_max
double precision,intent(inout) :: mean_bound_normal(3)
double precision,intent(inout) :: mean_bound_pos(3)
double precision,intent(out) :: aux_gravity(3)
double precision :: aux_scalar
double precision :: aux_vec(3)
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
if ((friction_angle>-1.d-9).and.                                               &
   (simulation_time>time_max_no_body_gravity_force)) then
! Active body-frontier interactions with non-negative friction angle. 
! Computations for the sliding friction force and explicit normal reaction 
! force under sliding.
! Sum of the normal vectors of the neighbouring frontiers: update
   mean_bound_normal(:) = mean_bound_normal(:) + normal_plane(:)
! Sum of the positions of the particles of the body each interaction: update
   mean_bound_pos(:) = mean_bound_pos(:) + bp_pos(:)
! Update of the number of "body particle - frontier" interactions for the 
! current body
   bp_bound_interactions = bp_bound_interactions + 1
! To update the maximum tangential velocity
   aux_vec(:) = dot_product(bp_arr(i_bp)%vel,normal_plane) * normal_plane(:)
   aux_vec(:) = bp_arr(i_bp)%vel(:) - aux_vec(:)
   aux_scalar = dsqrt(dot_product(aux_vec,aux_vec))
   interface_sliding_vel_max = max(interface_sliding_vel_max,aux_scalar)
   elseif (body_arr(bp_arr(i_bp)%body)%pmax<1.d-5) then
! "inactive gravity" ("balanced gravity") or negative friction angle. Dry body.
! Contribution of the normal reaction force under sliding and of the sliding 
! friction force depending on the local slope angle (instead of the friction 
! angle)
      aux_gravity(:) = 0.d0
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine body_boundary_for_sliding_friction_normal_reaction
