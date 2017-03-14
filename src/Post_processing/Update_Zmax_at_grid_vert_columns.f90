!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
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
! Program unit: Update_Zmax_at_grid_vert_columns                 
! Description: Updating the 2D array of the maximum values of the fluid particle
!              height, for each grid columns (only in 3D).
!              Printing the 2D field of the water depth (current time step), 
!              according to the output frequency chosen in the input file (only 
!              in 3D).
!              Printing the 2D fields of the specific flow rate components 
!              (current time step), at the same frequency of the water depth 
!              (only in 3D).            
!-------------------------------------------------------------------------------
subroutine Update_Zmax_at_grid_vert_columns(print_flag)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: print_flag
integer(4) :: npi,GridColumn,i_zone,i_vertex,i_aux,i_grid,j_grid,aux_integer
integer(4) :: alloc_stat,dealloc_stat
double precision :: pos(3)
double precision,allocatable,dimension(:) :: Z_fluid_step,h_step,qx_step
double precision,allocatable,dimension(:) :: qy_step,qx_step_grid,qy_step_grid
double precision,allocatable,dimension(:) :: n_part_step 
character(255) :: nomefile_h_step
integer(4),external :: ParticleCellNumber
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
do i_zone=1,NPartZone
   if (Partz(i_zone)%IC_source_type==2) then
      aux_integer = Partz(i_zone)%ID_last_vertex -                             &
                      Partz(i_zone)%ID_first_vertex + 1
      if (.not.allocated(h_step)) then
         allocate(h_step(aux_integer),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(nout,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',   &
               'Allocation of the array "h_step" failed; the simulation ',     &
               'stops here. '
            stop
         endif
      endif
      if (.not.allocated(qx_step)) then
         allocate(qx_step(aux_integer),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(nout,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',   &
               'Allocation of the array "qx_step" failed; the simulation ',    &
               'stops here.'
            stop
         endif
      endif
      if (.not.allocated(qy_step)) then
         allocate(qy_step(aux_integer),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(nout,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',   &
               'Allocation of the array "qy_step" failed; the simulation ',    &
               'stops here.'
            stop
         endif
      endif
      exit
   endif
enddo
if (.not.allocated(Z_fluid_step)) then
   allocate(Z_fluid_step(Grid%ncd(1)*Grid%ncd(2)),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(nout,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Allocation of the array "Z_fluid_step" failed; the simulation ',     &
         'stops here.'
      stop
   endif
endif
if (.not.allocated(qx_step_grid)) then
   allocate(qx_step_grid(Grid%ncd(1)*Grid%ncd(2)),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(nout,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Allocation of the array "qx_step_grid" failed; the simulation ',     &
         'stops here.'
      stop
   endif
endif
if (.not.allocated(qy_step_grid)) then
   allocate(qy_step_grid(Grid%ncd(1)*Grid%ncd(2)),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(nout,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Allocation of the array "qy_step_grid" failed; the simulation ',     &
         'stops here.'
      stop
   endif
endif
if (.not.allocated(n_part_step)) then
   allocate(n_part_step(Grid%ncd(1)*Grid%ncd(2)),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(nout,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Allocation of the array "n_part_step" failed; the simulation stops ',&
         'here.'
      stop
   endif
endif
!------------------------
! Initializations
!------------------------
Z_fluid_step = -999.d0
h_step = 0.d0
qx_step = 0.d0
qy_step = 0.d0
qx_step_grid = 0.d0
qy_step_grid = 0.d0 
n_part_step = 0
!------------------------
! Statements
!------------------------
!$omp parallel do default(none)                                                &
!$omp shared(nag,pg,Grid,Z_fluid_max,Domain,Z_fluid_step,n_part_step)          &
!$omp shared(qx_step_grid,qy_step_grid)                                        &
!$omp private(npi,GridColumn,pos)
do npi=1,nag
   pos(1) = pg(npi)%coord(1)
   pos(2) = pg(npi)%coord(2)
   pos(3) = Grid%extr(3,1) + 0.0000001d0
   GridColumn = ParticleCellNumber(pos)
   Z_fluid_step(GridColumn) = max(Z_fluid_step(GridColumn),(pg(npi)%coord(3) + &
                              0.5d0 * Domain%dx)) 
   qx_step_grid(GridColumn) = qx_step_grid(GridColumn) + pg(npi)%vel(1)
   qy_step_grid(GridColumn) = qy_step_grid(GridColumn) + pg(npi)%vel(2)
   n_part_step(GridColumn) = n_part_step(GridColumn) + 1
   Z_fluid_max(GridColumn) = max(Z_fluid_max(GridColumn),(pg(npi)%coord(3) +   &
                             0.5d0 * Domain%dx))
enddo
!$omp end parallel do
!$omp parallel do default(none)                                                &
!$omp shared(Grid,n_part_step,qx_step_grid,qy_step_grid)                       &
!$omp private(GridColumn,pos,i_grid,j_grid)
do i_grid=1,Grid%ncd(1)
   do j_grid=1,Grid%ncd(2)
      pos(1) = (i_grid - 0.5) * Grid%dcd(1) + Grid%extr(1,1)
      pos(2) = (j_grid - 0.5) * Grid%dcd(2) + Grid%extr(2,1)
      pos(3) = Grid%extr(3,1) + 0.0000001d0
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
! .txt file creation and heading (only at the beginning of the simulation)
if (on_going_time_step==1) then
   write(nomefile_h_step,"(a,a)") nomecaso(1:len_trim(nomecaso)),              &
      "_h_qx_qy_step.txt"
   open(ncpt,file=nomefile_h_step,status="unknown",form="formatted")
   write(ncpt,*) "Water depth (h), Specific flow rate components (q_x,q_y)"
   write(ncpt,'(7(a))') "           x(m)","           y(m)","  h_step(m)",     &
      " Z_fl_max_stp(m)","     z_topog(m)","     q_x(m^2/s)","     q_y(m^2/s)"
   flush(ncpt)
   else 
! Writing the 2D free surface field at the current time step
      if (print_flag==1) then
         write(nomefile_h_step,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)), &
            '_h_qx_qy_step',on_going_time_step,".txt"
         open(ncpt,file=nomefile_h_step,status="unknown",form="formatted")
      endif   
      do i_zone=1,NPartZone
         if (Partz(i_zone)%IC_source_type==2) then
!$omp parallel do default(none)                                                &
!$omp shared(ncpt,Partz,Vertice,Grid,h_step,Z_fluid_step,i_zone,qx_step)       &
!$omp shared(qy_step,qx_step_grid,qy_step_grid,n_part_step,q_max,print_flag)   &
!$omp private(i_vertex,GridColumn,pos,i_aux)
            do i_vertex=Partz(i_zone)%ID_first_vertex,                         &
               Partz(i_zone)%ID_last_vertex
               i_aux = i_vertex-Partz(i_zone)%ID_first_vertex + 1 
               pos(1) =  Vertice(1,i_vertex)
               pos(2) = Vertice(2,i_vertex)
               pos(3) = Grid%extr(3,1) + 0.0000001d0
               GridColumn = ParticleCellNumber(pos)
               h_step(i_aux) = max((Z_fluid_step(GridColumn) -                 &
                               Vertice(3,i_vertex)),0.d0)
               qx_step(i_aux) = qx_step_grid(GridColumn)
               qy_step(i_aux) = qy_step_grid(GridColumn)
               if (qx_step_grid(GridColumn)/=-999.d0) then
                  qx_step(i_aux) = qx_step_grid(GridColumn) * h_step(i_aux)
                  qy_step(i_aux) = qy_step_grid(GridColumn) * h_step(i_aux)
                  q_max(i_aux) = max(q_max(i_aux),dsqrt(qx_step(i_aux) ** 2 +  &
                                 qy_step(i_aux) ** 2))
               endif
               if (print_flag==1) then
!$omp critical (omp_write_h_step)
                  write (ncpt,'(7(f14.4,1x))') Vertice(1,i_vertex),            &
                     Vertice(2,i_vertex),h_step(i_aux),                        &
                     Z_fluid_step(GridColumn),Vertice(3,i_vertex),             &
                     qx_step(i_aux),qy_step(i_aux)
!$omp end critical (omp_write_h_step)
               endif
            enddo
!$omp end parallel do    
            exit
         endif
      enddo
endif
! Closing the file
if (print_flag==1) close(ncpt)
if (allocated(Z_fluid_step)) then
   deallocate(Z_fluid_step,STAT=dealloc_stat)
   if (dealloc_stat/=0) then
      write(nout,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Deallocation of the array "Z_fluid_step" failed; the simulation ',   &
         'stops here.'
      stop
   endif
endif
if (allocated(h_step)) then
   deallocate(h_step,STAT=dealloc_stat)
   if (dealloc_stat/=0) then
      write(nout,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Deallocation of the array "h_step" failed; the simulation ',         &
         'stops here.'
      stop
   endif
endif
if (allocated(qx_step)) then
   deallocate(qx_step,STAT=dealloc_stat)
   if (dealloc_stat/=0) then
      write(nout,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Deallocation of the array "qx_step" failed; the simulation ',        &
         'stops here.'
      stop
   endif
endif
if (allocated(qy_step)) then
   deallocate(qy_step,STAT=dealloc_stat)
   if (dealloc_stat/=0) then
      write(nout,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Deallocation of the array qy_step failed; the simulation ',          &
         'stops here.'
      stop
   endif
endif
if (allocated(qx_step_grid)) then
   deallocate(qx_step_grid,STAT=dealloc_stat)
   if (dealloc_stat/=0) then
      write(nout,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Deallocation of the array "qx_step_grid" failed; the simulation ',   &
         'stops here.'
      stop
   endif
endif
if (allocated(qy_step_grid)) then
   deallocate(qy_step_grid,STAT=dealloc_stat)
   if (dealloc_stat/=0) then
      write(nout,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Deallocation of the array "qy_step_grid" failed; the simulation ',   &
         'stops here.'
      stop
   endif
endif
if (allocated(n_part_step)) then
   deallocate(n_part_step,STAT=dealloc_stat)
   if (dealloc_stat/=0) then
      write(nout,*) 'Subroutine "Update_Zmax_at_grid_vert_columns". ',         &
         'Deallocation of the array "n_part_step" failed; the simulation ',    &
         'stops here.'
      stop
   endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine Update_Zmax_at_grid_vert_columns

