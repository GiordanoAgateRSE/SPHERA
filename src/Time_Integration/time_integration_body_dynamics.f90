!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)

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
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: time_integration_body_dynamics                                           
! Description: Euler time integration for body transport in fluid flows.  
!-------------------------------------------------------------------------------
subroutine time_integration_body_dynamics(dtvel) 
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
double precision,intent(in) :: dtvel
integer(4) :: i,npi,j
double precision :: mod_normal,aux_umax
double precision :: vec_temp(3),vec2_temp(3),vec3_temp(3),domega_dt(3)
double precision :: aux_vel(3)
double precision,allocatable,dimension(:,:) :: teta
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
allocate(teta(n_bodies,3))
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
! Loop over bodies (body dynamics)
!$omp parallel do default(none)                                                &
!$omp private(i,vec_temp,vec2_temp,vec3_temp,domega_dt,aux_umax)               &
!$omp shared(n_bodies,body_arr,dt,ncord,teta,dtvel,simulation_time)
do i=1,n_bodies
! Staggered parameters (velocity and angular velocity)
   if (body_arr(i)%imposed_kinematics==0) then
! Computed kinematics       
      if (ncord==2) then
         body_arr(i)%Force(2) = zero
         body_arr(i)%Moment(1) = zero
         body_arr(i)%Moment(3) = zero
      endif
      body_arr(i)%u_CM(:) = body_arr(i)%u_CM(:) + body_arr(i)%Force(:) /       &
         body_arr(i)%mass * dtvel 
      if (ncord==3) then
         vec_temp(1) = dot_product(body_arr(i)%Ic(1,:),body_arr(i)%omega)
         vec_temp(2) = dot_product(body_arr(i)%Ic(2,:),body_arr(i)%omega)
         vec_temp(3) = dot_product(body_arr(i)%Ic(3,:),body_arr(i)%omega)
         call Vector_Product(body_arr(i)%omega,vec_temp,vec3_temp,3)
         vec2_temp = body_arr(i)%Moment(:) - vec3_temp(:)
      endif
      if (ncord==2) vec2_temp = body_arr(i)%Moment(:) 
      domega_dt(1) = dot_product(body_arr(i)%Ic_inv(1,:),vec2_temp)
      domega_dt(2) = dot_product(body_arr(i)%Ic_inv(2,:),vec2_temp)
      domega_dt(3) = dot_product(body_arr(i)%Ic_inv(3,:),vec2_temp)
      body_arr(i)%omega(:) = body_arr(i)%omega(:) + domega_dt(:) * dtvel
      if (ncord==2) then
         body_arr(i)%x_CM(2) = zero 
         body_arr(i)%u_CM(2) = zero
         body_arr(i)%omega(1) = zero
         body_arr(i)%omega(3) = zero
      endif
      else
! Imposed kinematics
         do j=1,body_arr(i)%n_records
            if (body_arr(i)%body_kinematics(j,1)>=simulation_time) then
               if (body_arr(i)%body_kinematics(j,1)==simulation_time) then
                  body_arr(i)%u_CM(:) = body_arr(i)%body_kinematics(j,2:4)
                  body_arr(i)%omega(:) = body_arr(i)%body_kinematics(j,5:7)
                  else
                     body_arr(i)%u_CM(:) = body_arr(i)%body_kinematics(j-1,2:4)&
                        + (body_arr(i)%body_kinematics(j,2:4) -                &
                        body_arr(i)%body_kinematics(j-1,2:4)) /                &
                        (body_arr(i)%body_kinematics(j,1) -                    &
                        body_arr(i)%body_kinematics(j-1,1)) *                  &
                        (simulation_time - body_arr(i)%body_kinematics(j-1,1))
                     body_arr(i)%omega(:) =                                    &
                        body_arr(i)%body_kinematics(j-1,5:7) +                 &
                        (body_arr(i)%body_kinematics(j,5:7) -                  &
                        body_arr(i)%body_kinematics(j-1,5:7)) /                &
                        (body_arr(i)%body_kinematics(j,1) -                    &
                        body_arr(i)%body_kinematics(j-1,1)) *                  &
                        (simulation_time-body_arr(i)%body_kinematics(j-1,1))                                      
               endif
               exit
            endif
         enddo   
   endif
! Non-staggered parameters     
   body_arr(i)%x_CM(:) = body_arr(i)%x_CM(:) + body_arr(i)%u_CM(:) * dt 
   teta(i,:) = body_arr(i)%omega(:) * dt
   body_arr(i)%alfa(:) = body_arr(i)%alfa(:) + teta(i,:)
! Initializing umax (maximum particle velocity of the body)
   body_arr(i)%umax = zero
enddo
!$omp end parallel do
! Loop over body particles (static kinematics)
!$omp parallel do default(none)                                                &
!$omp private(npi,vec_temp,mod_normal,vec2_temp,aux_vel)                       &
!$omp shared(n_body_part,body_arr,bp_arr,dt,ncord,teta,dtvel)
do npi=1,n_body_part
! Staggered parameters
   call Vector_Product(body_arr(bp_arr(npi)%body)%omega,bp_arr(npi)%rel_pos,   &
      vec_temp,3)
   aux_vel(:) = bp_arr(npi)%vel(:) 
   bp_arr(npi)%vel(:) = body_arr(bp_arr(npi)%body)%u_CM(:) + vec_temp(:)
   if (ncord==2) bp_arr(npi)%vel(2) = zero 
   bp_arr(npi)%acc(:) = (bp_arr(npi)%vel(:) - aux_vel(:)) / dtvel
! Non-staggered parameters
   vec2_temp(:) = teta(bp_arr(npi)%body,:)
   call vector_rotation_Euler_angles(bp_arr(npi)%rel_pos,vec2_temp)
   if (ncord==2) bp_arr(npi)%rel_pos(2) = zero
   bp_arr(npi)%pos(:) = bp_arr(npi)%rel_pos(:) +                               &
      body_arr(bp_arr(npi)%body)%x_CM(:)
   call vector_rotation_Euler_angles(bp_arr(npi)%normal,vec2_temp)
   mod_normal = dsqrt(dot_product(bp_arr(npi)%normal,bp_arr(npi)%normal))
   if (mod_normal>one) bp_arr(npi)%normal(:) = bp_arr(npi)%normal(:) /         &
      mod_normal
   if (ncord==2) then
      bp_arr(npi)%rel_pos(2) = zero
      bp_arr(npi)%normal(2) = zero
   endif
enddo
!$omp end parallel do
! Updating max velocity within every body
do npi=1,n_body_part
   aux_umax = dsqrt(dot_product(bp_arr(npi)%vel,bp_arr(npi)%vel))        
   body_arr(bp_arr(npi)%body)%umax = max(body_arr(bp_arr(npi)%body)%umax,      &
      aux_umax)   
enddo
! Interesting part for RK1, not for Leapfrog: start
! Loop over body particles (kinematics)
!  do npi=1,n_body_part
!     bp_arr(npi)%pos(:) = bp_arr(npi)%pos(:) + bp_arr(npi)%vel(:) * dt
!     if (ncord==2) bp_arr(npi)%pos(2) = zero 
!  enddo
! Loop over bodies (body dynamics)
!  do i=1,n_bodies
!     if (ncord==2) then
!        body_arr(i)%Force(2) = zero
!        body_arr(i)%Moment(1) = zero
!        body_arr(i)%Moment(3) = zero
!     endif
!     body_arr(i)%x_CM(:) = body_arr(i)%x_CM(:) + body_arr(i)%u_CM(:) * dt 
!     teta(i,:) =  body_arr(i)%omega(:) * dt 
!     body_arr(i)%alfa(:) = body_arr(i)%alfa(:) + teta(i,:) 
!     body_arr(i)%u_CM(:) = body_arr(i)%u_CM(:) + body_arr(i)%Force(:) / &
!        body_arr(i)%mass * dt 
!     if (ncord==3) then
!        vec_temp(1) = dot_product(body_arr(i)%Ic(1,:),body_arr(i)%omega)
!        vec_temp(2) = dot_product(body_arr(i)%Ic(2,:),body_arr(i)%omega)
!        vec_temp(3) = dot_product(body_arr(i)%Ic(3,:),body_arr(i)%omega)
!        call Vector_Product(body_arr(i)%omega,vec_temp,vec3_temp,3)
!        vec2_temp = body_arr(i)%Moment(:) - vec3_temp(:)
!     endif
!     if (ncord==2) vec2_temp = body_arr(i)%Moment(:) 
!     domega_dt(1) = dot_product(body_arr(i)%Ic_inv(1,:),vec2_temp)
!     domega_dt(2) = dot_product(body_arr(i)%Ic_inv(2,:),vec2_temp)
!     domega_dt(3) = dot_product(body_arr(i)%Ic_inv(3,:),vec2_temp)
!     body_arr(i)%omega(:) = body_arr(i)%omega(:) + domega_dt(:) * dt
!     if (ncord==2) then
!        body_arr(i)%x_CM(2) = zero 
!        body_arr(i)%u_CM(2) = zero
!        body_arr(i)%omega(1) = zero
!        body_arr(i)%omega(3) = zero
!     endif
! Initializing umax (maximum particle velocity of the body)
!     body_arr(i)%umax = zero 
!  enddo
!
! Loop over body particles (static kinematics)
!  do npi=1,n_body_part
!     call Vector_Product(body_arr(bp_arr(npi)%body)%omega,bp_arr(npi)%rel_pos,&
!        vec_temp,3)
!     bp_arr(npi)%vel(:) = body_arr(bp_arr(npi)%body)%u_CM(:) + vec_temp(:)
!     bp_arr(npi)%rel_pos(:) = bp_arr(npi)%pos(:) -       &
!        body_arr(bp_arr(npi)%body)%x_CM(:)
!     vec2_temp(:) = teta(bp_arr(npi)%body,:)
!     call vector_rotation_Euler_angles(bp_arr(npi)%normal,vec2_temp)
!     mod_normal = dsqrt(dot_product(bp_arr(npi)%normal,bp_arr(npi)%normal))
!     if (mod_normal>1.) bp_arr(npi)%normal(:) = bp_arr(npi)%normal(:) / &
!         mod_normal    
!     if (ncord==2) then
!        bp_arr(npi)%vel(2) = zero 
!        bp_arr(npi)%rel_pos(2) = zero
!        bp_arr(npi)%normal(2) = zero
!     endif 
! Update of umax
!     aux_umax = dsqrt(dot_product(bp_arr(npi)%vel,bp_arr(npi)%vel))        
!     body_arr(bp_arr(npi)%body)%umax = max(body_arr(bp_arr(npi)%body)%umax, &
!       aux_umax)   
!  enddo
! Interesting part for RK1, not for Leapfrog: end
!------------------------
! Deallocations
!------------------------
deallocate(teta)
return
end subroutine time_integration_body_dynamics

