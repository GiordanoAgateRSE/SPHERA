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
! Program unit: Update_Zmax_at_grid_vert_columns                 
! Description: Updating the 2D arrays of the maximum values of the fluid 
!              particle height for each grid column (only in 3D).
!              Printing the current 2D fields of the free-surface height 
!              (filtered and unfiltered), water depth (filtered and 
!              unfiltered), specific flow rate components (filtered) and 
!              height-averaged velocity (filtered), according to the output 
!              frequency chosen in the input file (only in 3D).
!              The filter zeroes the effects associated with the presence of 
!              multiple fluid depths along the same vertical (atomization, 
!              breaking) and the liquid detachment from the topographic surface.
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine Update_Zmax_at_grid_vert_columns(print_flag)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_file_module
use Memory_I_O_interface_module
!------------------------
! Declarations
!------------------------
implicit none
logical :: fluid_below_flag
integer(4),intent(in) :: print_flag
integer(4) :: npi,GridColumn,i_zone,i_vertex,i_aux,i_grid,j_grid
integer(4) :: alloc_stat,dealloc_stat,i_cell,i_cp
double precision :: q_step,U_step,eps,h_filt_step,h_step,qx_step,qy_step
double precision :: pos(3)
logical,allocatable,dimension(:) :: compact_fluid
integer(4),allocatable,dimension(:) :: n_part_step
double precision,allocatable,dimension(:) :: qx_step_grid,qy_step_grid
double precision,allocatable,dimension(:) :: Z_fluid_step_min
character(100) :: array_name
character(255) :: nomefile_h_step
integer(4),external :: ParticleCellNumber
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
if (.not.allocated(qx_step_grid)) then
   allocate(qx_step_grid(Grid%ncd(1)*Grid%ncd(2)),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(ulog,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Allocation of the array "qx_step_grid" failed; the simulation ',     &
         'stops here.'
      stop
   endif
endif
if (.not.allocated(qy_step_grid)) then
   allocate(qy_step_grid(Grid%ncd(1)*Grid%ncd(2)),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(ulog,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Allocation of the array "qy_step_grid" failed; the simulation ',     &
         'stops here.'
      stop
   endif
endif
if (.not.allocated(n_part_step)) then
   allocate(n_part_step(Grid%ncd(1)*Grid%ncd(2)),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(ulog,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Allocation of the array "n_part_step" failed; the simulation stops ',&
         'here.'
      stop
   endif
endif
if (.not.allocated(compact_fluid)) then
   allocate(compact_fluid(Grid%ncd(1)*Grid%ncd(2)*Grid%ncd(3)),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(ulog,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Allocation of the array "compact_fluid" failed; the simulation ',    &
         'here.'
      stop
   endif
endif
array_name = "Z_fluid_step_min"
call allocate_de_dp_r1(.true.,Z_fluid_step_min,Grid%ncd(1)*Grid%ncd(2),        &
   array_name,.false.)
!------------------------
! Initializations
!------------------------
eps = 1.d-3
Z_fluid_step(:,:) = -999.d0
Z_fluid_step_min(:) = 999999.d0
qx_step_grid(:) = 0.d0
qy_step_grid(:) = 0.d0
n_part_step(:) = 0
compact_fluid(:) = .false.
if (on_going_time_step==1) z_topog_max(:) = max_negative_number
!------------------------
! Statements
!------------------------
! To assess the maximum height of the topographic surface per column of the 
! backgournd positioning grid
if (on_going_time_step==1) then
   do i_zone=1,NPartZone
      if (Partz(i_zone)%ID_first_vertex_sel>0) then
!$omp parallel do default(none)                                                &
!$omp shared(Partz,Vertice,Grid,i_zone,z_topog_max)                            &
!$omp private(i_vertex,GridColumn,pos)
         do i_vertex=Partz(i_zone)%ID_first_vertex_sel,                        &
            Partz(i_zone)%ID_last_vertex_sel
            pos(1) = Vertice(1,i_vertex)
            pos(2) = Vertice(2,i_vertex)
            pos(3) = Grid%extr(3,1) + 1.d-7
            GridColumn = ParticleCellNumber(pos)
!$omp critical (omp_z_topog_max)
            z_topog_max(GridColumn) = max(z_topog_max(GridColumn),             &
                                      Vertice(3,i_vertex))
!$omp end critical (omp_z_topog_max)
         enddo
!$omp end parallel do
         exit
      endif
   enddo
endif
! To detect the cells associated with the conditions of "compact fluid"  
! (suitable for detecting the fluid depth without neither atomization nor wave 
! breaking effects): start.
! Cycle over the background grid x-axis rows
!$omp parallel do default(none)                                                &
!$omp shared(Grid,Icont,compact_fluid,Z_fluid_step_min,pg,Domain,NPartOrd)     &
!$omp private(i_grid,j_grid,GridColumn,fluid_below_flag,i_cell,i_cp,npi)
do i_grid=1,Grid%ncd(1)
! Cycle over the background grid y-axis rows
   do j_grid=1,Grid%ncd(2)
! Column ID 
      GridColumn = (j_grid - 1) * Grid%ncd(1) + i_grid
! Initialization of the flag to detect the presence/absence of fluid particles 
! below a given cell
      fluid_below_flag = .false.
! Cycle along the cells of the grid column
      do i_cell=GridColumn,(Grid%ncd(3)-1)*(Grid%ncd(1)*Grid%ncd(2))+GridColumn,&
         Grid%ncd(1)*Grid%ncd(2)
         if (Icont(i_cell+1)<=Icont(i_cell)) then
! No particle in the cell "i_cell"
            if (fluid_below_flag.eqv..true.) exit
            else
! This cell contains the particle associated with the "lower fluid top height" 
! (along this background grid column) or contains fluid particles below the 
! cited height.
               compact_fluid(i_cell) = .true.
               if (fluid_below_flag.eqv..false.) then
! Loop over the cell particles
                  do i_cp=Icont(i_cell),Icont(i_cell+1)-1
                     npi = NPartOrd(i_cp)
                     Z_fluid_step_min(GridColumn) = min(                       &
                        Z_fluid_step_min(GridColumn),(pg(npi)%coord(3) -       &
                        0.5d0 * Domain%dx))
                  enddo
               endif
               fluid_below_flag = .true.
         endif
      enddo
   enddo
enddo
!$omp end parallel do
! To detect the cells associated with the conditions of "compact fluid": end.
!$omp parallel do default(none)                                                &
!$omp shared(nag,pg,Grid,Z_fluid_max,Domain,Z_fluid_step,n_part_step)          &
!$omp shared(qx_step_grid,qy_step_grid,compact_fluid,Z_fluid_step_min)         &
!$omp shared(z_topog_max)                                                      &
!$omp private(npi,GridColumn,pos)
do npi=1,nag
   pos(1) = pg(npi)%coord(1)
   pos(2) = pg(npi)%coord(2)
   pos(3) = Grid%extr(3,1) + 1.d-7
   GridColumn = ParticleCellNumber(pos)
!$omp critical (omp_update_2D_max)
   Z_fluid_step(GridColumn,1) = max(Z_fluid_step(GridColumn,1),                &
                                (pg(npi)%coord(3) + 0.5d0 * Domain%dx))
   Z_fluid_max(GridColumn,1) = max(Z_fluid_max(GridColumn,1),(pg(npi)%coord(3) &
                               + 0.5d0 * Domain%dx))
   if ((compact_fluid(pg(npi)%cella).eqv..true.).and.                          &
      (Z_fluid_step_min(GridColumn)<(z_topog_max(GridColumn)+0.5d0*Domain%dx)))&
      then
! In the columns where topographic vertices are absent (very coarse DTM/DEM 
! resolutions), the filter on the liquid detachment from the topographic 
! surface does not apply. However, there is no actual issue because the output 
! for all the 2D fields involved and the electrical substations only refers to 
! the topographic vertices.
! In the columns where SPH particles undergo the highest vertex of 
! the topographic surface, the filter on the liquid detachment from the 
! topographic surface does not apply.
      Z_fluid_step(GridColumn,2) = max(Z_fluid_step(GridColumn,2),             &
                                      (pg(npi)%coord(3) + 0.5d0 * Domain%dx))
      Z_fluid_max(GridColumn,2) = max(Z_fluid_max(GridColumn,2),               &
                                  (pg(npi)%coord(3) + 0.5d0 * Domain%dx))
      qx_step_grid(GridColumn) = qx_step_grid(GridColumn) + pg(npi)%vel(1)
      qy_step_grid(GridColumn) = qy_step_grid(GridColumn) + pg(npi)%vel(2)
      n_part_step(GridColumn) = n_part_step(GridColumn) + 1
   endif
!$omp end critical (omp_update_2D_max)
enddo
!$omp end parallel do
!$omp parallel do default(none)                                                &
!$omp shared(Grid,n_part_step,qx_step_grid,qy_step_grid)                       &
!$omp private(GridColumn,pos,i_grid,j_grid)
do i_grid=1,Grid%ncd(1)
   do j_grid=1,Grid%ncd(2)
! To compute "qx_step_grid" and "qy_step_grid" (at this stage they are 
! depth-averaged velocity components)
      pos(1) = (i_grid - 0.5) * Grid%dcd(1) + Grid%extr(1,1)
      pos(2) = (j_grid - 0.5) * Grid%dcd(2) + Grid%extr(2,1)
      pos(3) = Grid%extr(3,1) + 1.d-7
      GridColumn = ParticleCellNumber(pos)
      if (n_part_step(GridColumn)>0) then
         qx_step_grid(GridColumn) = qx_step_grid(GridColumn) /                 &
                                    n_part_step(GridColumn)
         qy_step_grid(GridColumn) = qy_step_grid(GridColumn) /                 &
                                    n_part_step(GridColumn)
         else
            qx_step_grid(GridColumn) = 0.d0
            qy_step_grid(GridColumn) = 0.d0
      endif
   enddo
enddo
!$omp end parallel do
! .txt file creation and heading
if (print_flag==1) then
   write(nomefile_h_step,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),       &
      '_hu_hf_qx_qy_step',on_going_time_step,".txt"
   open(u2DO,file=nomefile_h_step,status="unknown",form="formatted")
   write(u2DO,'(12(a))') "           x(m)","           y(m)","     hu_step(m)",&
                         "     hf_step(m)","  Zu_max_stp(m)","  Zf_max_stp(m)",&
                         "     z_topog(m)","     q_x(m^2/s)","     q_y(m^2/s)",&
                         "    U_step(m/s)","   Z_min_stp(m)"," z_topog_max(m)"
! Writing the 2D free surface field at the current time step
endif
do i_zone=1,NPartZone
   if (Partz(i_zone)%ID_first_vertex_sel>0) then
!$omp parallel do default(none)                                                &
!$omp shared(u2DO,Partz,Vertice,Grid,Z_fluid_step,i_zone)                      &
!$omp shared(qx_step_grid,qy_step_grid,q_max,U_max)                            &
!$omp shared(print_flag,eps,Z_fluid_step_min,z_topog_max)                      &
!$omp private(i_vertex,GridColumn,pos,i_aux,q_step,U_step,h_filt_step,h_step)  &
!$omp private(qx_step,qy_step)
      do i_vertex=Partz(i_zone)%ID_first_vertex_sel,                           &
         Partz(i_zone)%ID_last_vertex_sel
         i_aux = i_vertex - Partz(i_zone)%ID_first_vertex_sel + 1 
         pos(1) = Vertice(1,i_vertex)
         pos(2) = Vertice(2,i_vertex)
         pos(3) = Grid%extr(3,1) + 1.d-7
         GridColumn = ParticleCellNumber(pos)
         h_filt_step = max((Z_fluid_step(GridColumn,2) -                       &
                       Vertice(3,i_vertex)),0.d0)
         h_step = max((Z_fluid_step(GridColumn,1) - Vertice(3,i_vertex)),0.d0)
         qx_step = qx_step_grid(GridColumn)
         qy_step = qy_step_grid(GridColumn)
         U_step = 0.d0
         if (h_filt_step>=eps) then
! The topographic surface in the column is wet.
! Here "qx_step" and "qy_step" become specific flow rate components.
            qx_step = qx_step_grid(GridColumn) * h_filt_step
            qy_step = qy_step_grid(GridColumn) * h_filt_step
            q_step = dsqrt(qx_step ** 2 + qy_step ** 2)
            q_max(i_aux) = max(q_max(i_aux),q_step)
            U_step = q_step / h_filt_step
            U_max(i_aux) = max(U_max(i_aux),U_step)
            else
! The topographic surface in the column is dry
               qx_step = 0.d0
               qy_step = 0.d0
         endif
         if (print_flag==1) then
!$omp critical (omp_write_h_step)
            write(u2DO,'(12(f14.4,1x))') Vertice(1,i_vertex),                  &
               Vertice(2,i_vertex),h_step,h_filt_step,                         &
               Z_fluid_step(GridColumn,1),Z_fluid_step(GridColumn,2),          &
               Vertice(3,i_vertex),qx_step,qy_step,U_step,                     &
               Z_fluid_step_min(GridColumn),z_topog_max(GridColumn)
!$omp end critical (omp_write_h_step)
         endif
      enddo
!$omp end parallel do 
      exit
   endif
enddo
! Closing the file
if (print_flag==1) close(u2DO)
!------------------------
! Deallocations
!------------------------
if (allocated(qx_step_grid)) then
   deallocate(qx_step_grid,STAT=dealloc_stat)
   if (dealloc_stat/=0) then
      write(ulog,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Deallocation of the array "qx_step_grid" failed; the simulation ',   &
         'stops here.'
      stop
   endif
endif
if (allocated(qy_step_grid)) then
   deallocate(qy_step_grid,STAT=dealloc_stat)
   if (dealloc_stat/=0) then
      write(ulog,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Deallocation of the array "qy_step_grid" failed; the simulation ',   &
         'stops here.'
      stop
   endif
endif
if (allocated(n_part_step)) then
   deallocate(n_part_step,STAT=dealloc_stat)
   if (dealloc_stat/=0) then
      write(ulog,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Deallocation of the array "n_part_step" failed; the simulation ',    &
         'stops here.'
      stop
   endif
endif
if (allocated(compact_fluid)) then
   deallocate(compact_fluid,STAT=dealloc_stat)
   if (dealloc_stat/=0) then
      write(ulog,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Deallocation of the array "compact_fluid" failed; the simulation ',  &
         'stops here.'
      stop
   endif
endif
array_name = "Z_fluid_step_min"
call allocate_de_dp_r1(.false.,Z_fluid_step_min,array_name=array_name,         &
   ulog_flag=.false.)
return
end subroutine Update_Zmax_at_grid_vert_columns
#endif
