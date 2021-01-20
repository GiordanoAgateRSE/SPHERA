!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2021 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: BC_zmax_t0
! Description: Selection of the zone vertices at the beginning of the 
!              simulation for Dirichlet's BCs for the maximum fluid height. 
!              Only active in 3D with SASPH boundary treatment. The associated 
!              DEM-DTM or bottom has to be locally a uniform and isotropic 
!              Cartesian grid (in plan view). Assessment of the minimum local 
!              bottom height (influential in the presence of dry bottom).
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine BC_zmax_t0
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Dynamic_allocation_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
logical :: min_flag
integer(4) :: i_zone,i_vert,test_xy,n_vert_selected,n_aux,alloc_stat
integer(4),dimension(:),allocatable :: aux_array,select_vert_array
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine point_inout_convex_non_degenerate_polygon(point,n_sides,         &
      point_pol_1,point_pol_2,point_pol_3,point_pol_4,point_pol_5,point_pol_6, &
      test)
      implicit none
      integer(4),intent(in) :: n_sides
      double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
      double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
      double precision,intent(in) :: point_pol_5(2),point_pol_6(2)
      integer(4),intent(inout) :: test
   end subroutine point_inout_convex_non_degenerate_polygon
end interface
interface
   subroutine z_min_max_DEM_DTM_9p_stencil(min_flag,i_zone,i_vertex,z_aux)
      implicit none
      logical,intent(in) :: min_flag
      integer(4),intent(in) :: i_zone,i_vertex
      double precision,intent(inout) :: z_aux
   end subroutine z_min_max_DEM_DTM_9p_stencil
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
min_flag = .false.
!------------------------
! Statements
!------------------------
! Loop over the zones
do i_zone=1,NPartZone
! To select any zone of type "zmax"
   if (Partz(i_zone)%tipo/="zmax") cycle
! Initialization of the first dimension of the array of the selected vertices
   n_aux = 100
! Initialization of the number of vertices selected
   n_vert_selected = 0
! Allocate the array of the selected vertices
   if (.not.allocated(select_vert_array)) then
      allocate(select_vert_array(n_aux),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(uerr,*) 'Allocation of the array "select_vert_array"',          &
            ' failed; the program stops here (program unit "BC_zmax_t0"). '
         stop
      endif
   endif
! Loop over the vertices of the associated DEM-DTM or bottom
   do i_vert=Partz(i_zone)%ID_first_vertex_sel,Partz(i_zone)%ID_last_vertex_sel
! Initialize the test variable
      test_xy = 0
! Check if the vertex lies inside the BC zone
      call point_inout_convex_non_degenerate_polygon(                          &
         Vertice(1:2,i_vert),Partz(i_zone)%plan_reservoir_points,              &
         Partz(i_zone)%plan_reservoir_pos(1,1:2),                              &
         Partz(i_zone)%plan_reservoir_pos(2,1:2),                              &
         Partz(i_zone)%plan_reservoir_pos(3,1:2),                              &
         Partz(i_zone)%plan_reservoir_pos(4,1:2),                              &
         Partz(i_zone)%plan_reservoir_pos(4,1:2),                              &
         Partz(i_zone)%plan_reservoir_pos(4,1:2),test_xy)
      if (test_xy==1) then
! The vertex lies inside the BC zone
! Update the number of selected vertices
         n_vert_selected = n_vert_selected + 1
         if (n_vert_selected>n_aux) then
! The length of the array of the selected vertices is too short
! Allocation of an auxiliary array
            if (.not.allocated(aux_array)) then
               allocate(aux_array(n_aux),STAT=alloc_stat)
               if (alloc_stat/=0) then
                  write(uerr,*) 'Allocation of the auxiliary array ',          &
                     '"aux_array" failed; the program stops here ',            &
                     '(program unit "BC_zmax_t0"). '
                  stop
               endif
            endif
! Copy the array of the selected vertices (it has to be reallocated) in its 
! auxiliary array
            aux_array(:) = select_vert_array(:)
! Deallocate the array of the selected vertices (being too short)
            if(allocated(select_vert_array)) then
               deallocate(select_vert_array,STAT=alloc_stat)
               if (alloc_stat/=0) then
                  write(uerr,*) 'Deallocation of the array ',                  &
                     '"select_vert_array" in the program unit "BC_zmax_t0" ',  &
                     'failed; the execution terminates here. '
                  stop
               endif
            endif
! Update the dimension of the auxiliary array for reallocation 
            n_aux = n_aux + 100
! Reallocate the array of the selected vertices
            if (.not.allocated(select_vert_array)) then
               allocate(select_vert_array(n_aux),STAT=alloc_stat)
               if (alloc_stat/=0) then
                  write(uerr,*) 'Reallocation of the array ',                  &
                     '"select_vert_array" failed; the program stops here ',    &
                     '(program unit BC_zmax_t0"). '
                  stop
               endif
            endif
! Initialize the reallocated array
            select_vert_array(1:size(aux_array)) = aux_array(1:size(aux_array))
            select_vert_array(size(aux_array)+1:n_aux) = 0
! Deallocate the auxiliary array
            if(allocated(aux_array)) then
               deallocate(aux_array,STAT=alloc_stat)
               if (alloc_stat/=0) then
                  write(uerr,*) 'Deallocation of the auxiliary array ',        &
                     '"aux_array" in the program unit "BC_zmax_t0" ',          &
                     'failed; the execution terminates here. '
                  stop
               endif
            endif
         endif
! Save the current vertex ID in the array of the selected vertices  
         select_vert_array(n_vert_selected) = i_vert
         select_vert_array(n_vert_selected) = i_vert
      endif
   enddo
! Allocation of zone array "BC_zmax_vertices" (i.e., the vertex array of the 
! "zmax" zone). The zone arrays "BC_zmax_vertices" are implicitly 
! deallocated with the deallocation of the derived type they belong to ("Partz")
   if (.not.allocated(Partz(i_zone)%BC_zmax_vertices)) then
      allocate(Partz(i_zone)%BC_zmax_vertices(n_vert_selected,3),              &
         STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(uerr,*) 'Allocation of the zone array "BC_zmax_vertices" ',     &
            'failed; the program stops here (program unit "BC_zmax_t0"). '
         stop
         else
            write(ulog,'(3a,i4,a,i9,a)') 'The allocation of the zone array ',  &
               '"BC_zmax_vertices" has successfully completed in the program', &
               ' unit "BC_zmax_t0": izone = ',i_zone,'; n_vert_selected = ',   &
               n_vert_selected,'.'
      endif
   endif
! Loop over the selected vertices (here "i_vert" is a relative vertex ID)
   do i_vert=1,n_vert_selected
! Copy the horizontal coordinates of the active elements of the array of 
! the selected vertices into the reference zone array
      Partz(i_zone)%BC_zmax_vertices(i_vert,1) =                               &
         Vertice(1,select_vert_array(i_vert))
      Partz(i_zone)%BC_zmax_vertices(i_vert,2) =                               &
         Vertice(2,select_vert_array(i_vert))
! Initialize the minimum height to start the extrusion (it will be influential 
! only in the presence of dry bottom)
      call z_min_max_DEM_DTM_9p_stencil(min_flag,i_zone,                       &
         select_vert_array(i_vert),Partz(i_zone)%BC_zmax_vertices(i_vert,3))
   enddo
! Deallocate the array of the selected vertices
   if(allocated(select_vert_array)) then
      deallocate(select_vert_array,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(uerr,*) 'Deallocation of the array ',                           &
            '"select_vert_array" in the program unit "BC_zmax_t0" failed; ',   &
            'the execution stops here. '
         stop
      endif
   endif
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine BC_zmax_t0
#endif
