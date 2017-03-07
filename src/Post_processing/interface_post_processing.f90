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
! Program unit: interface_post_processing                  
! Description: Post-processing the interfaces for bed-load transport phenomena.             
!-------------------------------------------------------------------------------
subroutine interface_post_processing
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
integer(4) :: i_grid,j_grid,i_cell,i_aux,i,k_grid,i2_grid,j2_grid
double precision :: pos_aux(3)
character(255) :: nomefile_blt_interfaces
integer(4),external :: CellIndices,ParticleCellNumber
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
! .txt file creation and heading
write(nomefile_blt_interfaces,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),  &
   '_blt_interfaces_',on_going_time_step,".txt"
open(ncpt,file=nomefile_blt_interfaces,status="unknown",form="formatted")
if (on_going_time_step==1) then
! First step
   write (ncpt,*) "Bed load transport interfaces "
   write (ncpt,'((7x,a),(5x,a),(5x,a),(7x,a),(7x,a),(6x,a),(6x,a),(2x,a))')    &
      " Time(s)"," x_grid(m)"," y_grid(m)"," z_FS(m)"," z_fm(m)"," z_bed(m)",  &
      " z_bot(m)"," z_sat_top(m)"
   flush(ncpt)
   else
! Other steps
! Loop over the monitoring lines
      do i=1,Granular_flows_options%monitoring_lines
         if (Granular_flows_options%lines(i,1)==-999.d0) then
            do i2_grid=1,Grid%ncd(1)
               pos_aux(1) = (i2_grid - 0.5) * Grid%dcd(1) + grid%extr(1,1)
               pos_aux(2) = Granular_flows_options%lines(i,2)
               pos_aux(3) = 0.5d0 * Grid%dcd(3) + grid%extr(3,1) 
! Note: pos_aux(3) = 0.d0
               i_cell = ParticleCellNumber(pos_aux) 
               i_aux = CellIndices(i_cell,i_grid,j_grid,k_grid)
               call write_Granular_flows_interfaces(i_grid,j_grid,pos_aux)
            enddo
         endif
         if (Granular_flows_options%lines(i,2)==-999.d0) then
            do j2_grid=1,Grid%ncd(2)
               pos_aux(1) = Granular_flows_options%lines(i,1)
               pos_aux(2) = (j2_grid - 0.5) * Grid%dcd(2) + grid%extr(2,1)
               pos_aux(3) = 0.5d0 * Grid%dcd(3) + grid%extr(3,1) 
! Note: pos_aux(3) = 0.d0
               i_cell = ParticleCellNumber(pos_aux) 
               i_aux = CellIndices(i_cell,i_grid,j_grid,k_grid)
               call write_Granular_flows_interfaces(i_grid,j_grid,pos_aux)
            enddo
         endif
         if ((Granular_flows_options%lines(i,1).ne.-999.d0).and.               &
            (Granular_flows_options%lines(i,2).ne.-999.d0)) then
            pos_aux(1) = Granular_flows_options%lines(i,1)
            pos_aux(2) = Granular_flows_options%lines(i,2)
            pos_aux(3) = 0.5d0 * Grid%dcd(3) + grid%extr(3,1) 
! Note: pos_aux(3) = 0.d0
            i_cell = ParticleCellNumber(pos_aux) 
            i_aux = CellIndices(i_cell,i_grid,j_grid,k_grid)
            call write_Granular_flows_interfaces(i_grid,j_grid,pos_aux)
         endif        
      enddo
endif
close(ncpt)
!------------------------
! Deallocations
!------------------------
return
end subroutine interface_post_processing

