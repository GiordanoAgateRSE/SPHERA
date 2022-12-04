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
! Program unit: calc_pelo
! Description: Post-processing to write the time evolution of the free 
!              surface height at the monitoring lines. Filtering atomization of 
!              the fluid domain. The “lower fluid top height” 
!              (“z_lower_fluid_top”) is suitable to estimate the liquid depth 
!              as a filter of the “maximum liquid height” (“pelolib”) to cut 
!              off atomization and wave breaking effects. Filter criterion: 
!              during the upward search (starting from the bottom) for SPH 
!              particles along the monitoring line, stop at the first empty 
!              background cell, if at least one cell below is filled with the 
!              searched fluid.
!-------------------------------------------------------------------------------
subroutine calc_pelo
!------------------------
! Modules
!------------------------
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
logical :: fluid_below_flag,zlft_flag
integer(4) :: i,j,ncel,npartcel,mm,jj,minpos1,minpos2,alloc_stat
double precision,dimension(3) :: ragtemp
integer(4),dimension(:),allocatable :: particle_ID
double precision,dimension(:),allocatable :: particle_point_dis2
double precision,dimension(3,nlines) :: pelolib,z_lower_fluid_top
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
if (.not.allocated(particle_ID)) then
   allocate(particle_ID(NMAXPARTJ),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(ulog,*) 'Subroutine "calc_pelo". Allocation of the array ',        &
         '"particle_ID" failed; the simulation stops here.'
      stop
   endif
endif
if (.not.allocated(particle_point_dis2)) then
   allocate(particle_point_dis2(NMAXPARTJ),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(ulog,*) 'Subroutine "calc_pelo". Allocation of the array ',        &
         '"particle_point_dis2" failed; the simulation stops here.'
      stop
   endif
endif
!------------------------
! Initializations
!------------------------
pelolib(:,:) = 0.d0
z_lower_fluid_top(:,:) = 0.d0
!------------------------
! Statements
!------------------------
! Loop over the lines
!$omp parallel do default(none)                                                &
!$omp shared(nlines,control_lines,control_points,Icont,NPartOrd,pg,pelolib)    &
!$omp shared(Domain,z_lower_fluid_top)                                         &
!$omp private(i,zlft_flag,fluid_below_flag,j,ncel,npartcel,particle_ID)        &
!$omp private(particle_point_dis2,mm,jj,ragtemp,minpos1,minpos2)
do i=1,nlines
! Initialization of the flag for the lower fluid top height
   zlft_flag = .true.
! Initialization of the flag for the presence/absence of fluid particles below 
! the current cell
   fluid_below_flag = .false.
! Loop over the line points
   do j=control_lines(i)%Icont(1)+1,control_lines(i)%Icont(2)
! Cell of the line point
      ncel = control_points(j)%cella
! Particle out of the numerical domain
      if (ncel==0) cycle
! No particle in the cell
      if (Icont(ncel+1)<=Icont(ncel)) then
         if (fluid_below_flag.eqv..true.) zlft_flag = .false.
         cycle
      endif
! Search for the two particles, which are the closest to the monitoring point 
! and in the same background cell: computation of their averaged position.
      npartcel = 0
      particle_ID(:) = 0
      particle_point_dis2(:) = 9.d12
! Loop over the cell particles 
      do mm = Icont(ncel),Icont(ncel+1)-1  
         jj = NPartOrd(mm)
         ragtemp(1:3) = control_points(j)%coord(1:3) - pg(jj)%coord(1:3)
         npartcel = npartcel + 1
         particle_ID(npartcel) = jj
         particle_point_dis2(npartcel) = ragtemp(1) * ragtemp(1) + ragtemp(2)  &
                                         * ragtemp(2) + ragtemp(3) * ragtemp(3)
      enddo
      minpos1 = minloc(particle_point_dis2,1)
      particle_point_dis2(minpos1) = 9.d12
! In case there is just one particle in the cell, this is re-selected as the 
! second particle
      minpos2 = minloc(particle_point_dis2,1)
      pelolib(1:3,i) = (pg(particle_ID(minpos1))%coord(1:3) +                  &
                       pg(particle_ID(minpos2))%coord(1:3)) * half +           &
                       Domain%dx / 2.d0
      if (zlft_flag.eqv..true.) z_lower_fluid_top(1:3,i) = pelolib(1:3,i)
      fluid_below_flag = .true.
   enddo
enddo
!$omp end parallel do
write(nplb,'(30g14.7)') simulation_time,pelolib
write(uzlft,'(30g14.7)') simulation_time,z_lower_fluid_top
!------------------------
! Deallocations
!------------------------
if (allocated(particle_ID)) then
   deallocate(particle_ID,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(ulog,*) 'Subroutine "calc_pelo". Deallocation of the array ',      &
         '"particle_ID" failed; the simulation stops here.'
      stop
   endif
endif
if (allocated(particle_point_dis2)) then
   deallocate(particle_point_dis2,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(ulog,*) 'Subroutine "calc_pelo". Deallocation of the array ',      &
         '"particle_point_dis2" failed; the simulation stops here.'
      stop
   endif
endif
return
end subroutine calc_pelo
