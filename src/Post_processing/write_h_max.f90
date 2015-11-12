!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-Ricerca di Sistema)



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
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: write_h_max                  
! Description: To compute and write the 2D array of the maximum values of the water depth, at the nodes of the Cartesian 
!              topography, provided as input data (only in 3D). Same task for the 2D field of the maximum (over time) specific 
!              flow rates.              
!----------------------------------------------------------------------------------------------------------------------------------

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
integer(4) :: i_zone,i_vertex,GridColumn
double precision :: pos(3)
double precision,dimension(:),allocatable :: h_max 
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
! h_max .txt file: creation and heading 
write(nomefile_h_max,"(a,a)") nomecaso(1:len_trim(nomecaso)),"_h_max_q_max.txt"
open(ncpt,file=nomefile_h_max,status="unknown",form="formatted")
write(ncpt,*) "Maximum water depth (m) and specific flow rate (m^2/s)"
write(ncpt,'(6(a))') "           x(m)","           y(m)","       h_max(m)",    &
   " Z_fluid_max(m)","     z_topog(m)","   q_max(m^2/s)"
flush(ncpt) 
do i_zone=1,NPartZone
   if (Partz(i_zone)%IC_source_type==2) then
! Allocating h_max 
      allocate(h_max(Partz(i_zone)%npoints))
! Initializing h_max
      h_max = 0.d0
!$omp parallel do default(none)                                                &
!$omp shared(Partz,Vertice,Grid,h_max,Z_fluid_max,ncpt,i_zone,q_max)           &
!$omp private(i_vertex,GridColumn,pos)
      do i_vertex=Partz(i_zone)%ID_first_vertex,Partz(i_zone)%ID_last_vertex
         pos(1) = Vertice(1,i_vertex)
         pos(2) = Vertice(2,i_vertex)
         pos(3) = Grid%extr(3,1)+0.0000001d0
         GridColumn = ParticleCellNumber(pos)
         h_max(i_vertex-Partz(i_zone)%ID_first_vertex+1) =                     &
            max((Z_fluid_max(GridColumn) - Vertice(3,i_vertex)),0.d0)
!$omp critical (omp_write_h_max)
         write(ncpt,'(6(f14.4,1x))')Vertice(1,i_vertex),Vertice(2,i_vertex),   &
            h_max(i_vertex-Partz(i_zone)%ID_first_vertex+1),                   &
            Z_fluid_max(GridColumn),Vertice(3,i_vertex),                       &
            q_max(i_vertex-Partz(i_zone)%ID_first_vertex+1)        
!$omp end critical (omp_write_h_max)
      enddo
!$omp end parallel do    
      exit
   endif
enddo
! h_max .txt file: closing
close (ncpt)
if (allocated(h_max)) deallocate(h_max)
if (allocated(q_max)) deallocate(q_max)
!------------------------
! Deallocations
!------------------------
return
end subroutine write_h_max

