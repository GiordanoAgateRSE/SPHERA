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
! Program unit: write_h_max                  
! Description: To compute and write the 2D fields of the maximum values of the
!              water depth, at the nodes of the Cartesian topography, provided
!              as input data (only in 3D; two variants -hu, hf- based either on 
!              the unfiltered -Zu- or filtered -Zf- free surface height). Same  
!              task for the 2D field of the maximum (over time) specific flow 
!              rates.
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine write_h_max
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
integer(4) :: i_zone,i_vertex,GridColumn,alloc_stat,dealloc_stat,aux_integer
double precision :: pos(3)
double precision,dimension(:,:),allocatable :: h_max 
character(255) :: nomefile_h_max
integer(4),external :: ParticleCellNumber
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
! h_max .txt files: creation and headings 
write(nomefile_h_max,"(a,a)") nomecaso(1:len_trim(nomecaso)),                  &
   "_hmax_qmax_Umax.txt"
open(ncpt,file=nomefile_h_max,status="unknown",form="formatted")
write(ncpt,'(9(a))') "           x(m)","           y(m)","      hu_max(m)",    &
   "      hf_max(m)","  Zu_flu_max(m)","  Zf_flu_max(m)","     z_topog(m)",    &
   "   q_max(m^2/s)","     U_max(m/s)"
flush(ncpt) 
do i_zone=1,NPartZone
   if (Partz(i_zone)%ID_first_vertex_sel>0) then
! Allocating h_max
      if (.not.allocated(h_max)) then
         aux_integer = Partz(i_zone)%ID_last_vertex_sel -                      &
                       Partz(i_zone)%ID_first_vertex_sel + 1
         allocate(h_max(aux_integer,2),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(ulog,*) 'Subroutine "write_h_max". ',                        &
               'Allocation of the array "h_max" failed; the simulation ',      &
               'stops here. '
            stop
         endif
      endif
! Initializing h_max
      h_max(:,:) = 0.d0
!$omp parallel do default(none)                                                &
!$omp shared(Partz,Vertice,Grid,h_max,Z_fluid_max,ncpt,i_zone,q_max,U_max)     &
!$omp private(i_vertex,GridColumn,pos)
      do i_vertex=Partz(i_zone)%ID_first_vertex_sel,                           &
         Partz(i_zone)%ID_last_vertex_sel
         pos(1) = Vertice(1,i_vertex)
         pos(2) = Vertice(2,i_vertex)
         pos(3) = Grid%extr(3,1) + 1.d-7
         GridColumn = ParticleCellNumber(pos)
         h_max(i_vertex-Partz(i_zone)%ID_first_vertex_sel+1,1) =               &
            max((Z_fluid_max(GridColumn,1) - Vertice(3,i_vertex)),0.d0)
         h_max(i_vertex-Partz(i_zone)%ID_first_vertex_sel+1,2) =               &
            max((Z_fluid_max(GridColumn,2) - Vertice(3,i_vertex)),0.d0)
!$omp critical (omp_write_h_max)
         write(ncpt,'(9(f14.4,1x))')Vertice(1,i_vertex),Vertice(2,i_vertex),   &
            h_max(i_vertex-Partz(i_zone)%ID_first_vertex_sel+1,1),             &
            h_max(i_vertex-Partz(i_zone)%ID_first_vertex_sel+1,2),             &
            Z_fluid_max(GridColumn,1),Z_fluid_max(GridColumn,2),               &
            Vertice(3,i_vertex),                                               &
            q_max(i_vertex-Partz(i_zone)%ID_first_vertex_sel+1),               &
            U_max(i_vertex-Partz(i_zone)%ID_first_vertex_sel+1)
!$omp end critical (omp_write_h_max)
      enddo
!$omp end parallel do
      exit
   endif
enddo
! h_max .txt file: closing
close (ncpt)
if (allocated(h_max)) then
   deallocate(h_max,STAT=dealloc_stat)
   if (dealloc_stat/=0) then
      write(ulog,*) 'Subroutine "write_h_max". ',                              &
         'Deallocation of the array "h_max" failed; the simulation stops here.'
      stop
   endif
endif
if (allocated(q_max)) then
   deallocate(q_max,STAT=dealloc_stat)
   if (dealloc_stat/=0) then
      write(ulog,*) 'Subroutine "write_h_max". ',                              &
         'Deallocation of the array "q_max" failed; the simulation stops here.'
      stop
   endif
endif
if (allocated(U_max)) then
   deallocate(U_max,STAT=dealloc_stat)
   if (dealloc_stat/=0) then
      write(ulog,*) 'Subroutine "write_h_max". ',                              &
         'Deallocation of the array "U_max" failed; the simulation stops here.'
      stop
   endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine write_h_max
#endif
