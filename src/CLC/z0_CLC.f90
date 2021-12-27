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
! Program unit: z0_CLC
! Description: Assessment of the 2D field of the CLC class. Assessment of the 
!              2D field of the roughness characteristic length (z0) as function 
!              of the CLC class. Writing of a unique 2D output file on the 
!              horizontal background grid with the following quantities: x, y, 
!              CLC class, z0.
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine z0_CLC
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_file_module
use Memory_I_O_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: i_face,ix,iy,i_hcel,i_pol,CLC_active_cells,io_stat,ier,test
double precision :: CLC_filling,z0_default
double precision :: pos_cell(3)
character(100) :: array_name,file_name
logical,external :: ReadCheck
integer(4),external :: CellNumber
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
   subroutine cell_indices2pos(ix,iy,iz,pos)
      implicit none
      integer(4),intent(in) :: ix,iy,iz
      double precision,intent(out) :: pos(3)
   end subroutine cell_indices2pos
end interface
!------------------------
! Allocations
!------------------------
! Allocate the CLC 2D array/grid, with the same sizes of the background grid 
! along the horizontal directions
array_name = "CLC%class_2D"
call allocate_de_int4_r2(.true.,CLC%class_2D,Grid%ncd(1),Grid%ncd(2),          &
   array_name)
! 2D array/grid of z0 with the same sizes of the background grid along the 
! horizontal directions
array_name = "CLC%z0"
call allocate_de_dp_r2(.true.,CLC%z0,Grid%ncd(1),Grid%ncd(2),array_name)
!------------------------
! Initializations
!------------------------
CLC_active_cells = 0
CLC%class_2D(:,:) = 0
!------------------------
! Statements
!------------------------
! Initialize the z0 2D array with "z0_default.txt"
file_name = "./input/20_CLC/z0_default.txt"
call open_close_file(.true.,max_file_unit_booked+1,file_name)
read(max_file_unit_booked+1,*,iostat=io_stat) z0_default
if (.not.ReadCheck(io_stat,ier,1,"./input/20_CLC/z0_default.txt","z0_default", &
   max_file_unit_booked+1,ulog)) return
file_name = "./input/20_CLC/z0_default.txt"
call open_close_file(.false.,max_file_unit_booked+1,file_name)
CLC%z0(:,:) = z0_default
! Open the 2D output file
file_name = "z0_CLC_output.txt"
call open_close_file(.true.,uCLC,file_name)
! 2D output file label
write(uCLC,'(3(a10),(a8))') "x(m)","y(m)","CLC_class","z0(m)" 
! Loop over the horizontal background-grid cells
!$omp parallel do default(none)                                                &
!$omp shared (Grid,n_neigh_hcell_CLCpol,CLC,CLC_active_cells)                  &
!$omp private(iy,ix,i_hcel,pos_cell,i_pol,i_face,test)
do iy=1,Grid%ncd(2)
   do_ix: do ix=1,Grid%ncd(1)
! Index of the horizontal grid cell
      i_hcel = CellNumber(ix,iy,1)
! Coordinates of the cell barycentre
      call cell_indices2pos(ix,iy,1,pos_cell(1:3))
! Loop over the neighbouring CLC polygons
      do i_pol=1,n_neigh_hcell_CLCpol(i_hcel)
! Loop over the triangles of the CLC polygon
         do i_face=1,CLC%polygons(i_pol)%n_faces
! Test if the cell barycentre lies within the current triangle
            call point_inout_convex_non_degenerate_polygon(pos_cell(1:2),3,    &
   CLC%polygons(i_pol)%vertices(CLC%polygons(i_pol)%faces(i_face,1),1:2),      &
   CLC%polygons(i_pol)%vertices(CLC%polygons(i_pol)%faces(i_face,2),1:2),      &
   CLC%polygons(i_pol)%vertices(CLC%polygons(i_pol)%faces(i_face,3),1:2),      &
   CLC%polygons(i_pol)%vertices(CLC%polygons(i_pol)%faces(i_face,3),1:2),      &
   CLC%polygons(i_pol)%vertices(CLC%polygons(i_pol)%faces(i_face,3),1:2),      &
   CLC%polygons(i_pol)%vertices(CLC%polygons(i_pol)%faces(i_face,3),1:2),test)
            if (test==1) then
! Assign the class of the CLC polygon to the current element of the 
! CLC%class_2D array
               CLC%class_2D(ix,iy) = CLC%polygons(i_pol)%CLC_class
! Assign z0
               call z0_CLC_table(ix,iy)
!$omp critical (omp_CLC_active_cells)
! Update the count of CLC active cells
               CLC_active_cells = CLC_active_cells + 1
!$omp end critical (omp_CLC_active_cells)
               cycle do_ix
            endif
         enddo
      enddo
   enddo do_ix
enddo
!$omp end parallel do
! Update the 2D output file
! Loop over the horizontal background-grid cells
do iy=1,Grid%ncd(2)
   do ix=1,Grid%ncd(1)
      write(uCLC,'(2(f10.3),(i10),(f8.4))') pos_cell(1:2),CLC%class_2D(ix,iy), &
         CLC%z0(ix,iy)
   enddo
enddo
! Close the 2D output file
call open_close_file(.false.,uCLC,file_name)
! Based on the number of cells with null CLC, write the percentage of the total 
! area occupied by the CLC polygons in the log file. In case it is below 95% of 
! the 2D domain area, suggest to possibly cut a larger CLC input file.
CLC_filling = dfloat(CLC_active_cells) / dfloat(Grid%ncd(1)) /                 &
              dfloat(Grid%ncd(2)) * 100.d0
write(ulog,'(a,f5.1,a)') "The CLC input data fill the ",CLC_filling,           &
   "% of the horizontal domain."
if (CLC_filling<95.d0) then
   write(ulog,*) "It is suggested to consider the possibility of feeding ",    &
      "SPHERA with a larger CLC input file."
endif
!------------------------
! Deallocations
!------------------------
array_name = "CLC%class_2D"
call allocate_de_int4_r2(.false.,CLC%class_2D,array_name=array_name)
array_name = "CLC%z0"
call allocate_de_dp_r2(.false.,CLC%z0,array_name=array_name)
array_name = "n_neigh_hcell_CLCpol"
call allocate_de_int4_r1(.false.,n_neigh_hcell_CLCpol,array_name=array_name)
array_name = "neigh_hcell_CLCpol"
call allocate_de_int4_r1(.false.,neigh_hcell_CLCpol,array_name=array_name)
do i_pol=1,CLC%n_polygons
   write(array_name,*) i_pol
   array_name = "CLC%polygons(" // trim(array_name) // ")%vertices"
   call allocate_de_dp_r2(.false.,CLC%polygons(i_pol)%vertices,                &
      array_name=array_name)
   write(array_name,*) i_pol
   array_name = "CLC%polygons(" // trim(array_name) // ")%faces"
   call allocate_de_int4_r2(.false.,CLC%polygons(i_pol)%faces,                 &
      array_name=array_name)
enddo
array_name = "CLC%polygons"
call allocate_de_CLCp_r1(.false.,CLC%polygons,array_name=array_name)
return
end subroutine z0_CLC
#endif
