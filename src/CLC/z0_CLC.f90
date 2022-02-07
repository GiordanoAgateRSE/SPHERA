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
use Memory_I_O_interface_module
use Neighbouring_Search_interface_module
!------------------------
! Declarations
!------------------------
implicit none
logical :: activation_flag
integer(4) :: i_face,ix,iy,i_hcel,i_pol,CLC_active_cells,io_stat,ier,test
integer(4) :: i_pol_nei,aux_int,F_VUC
double precision :: CLC_filling,z0_default
double precision :: pos_cell(3)
character(100) :: array_name,file_name
logical,external :: ReadCheck
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
!------------------------
! Allocations
!------------------------
! Allocate the CLC 2D array/grid, with the same sizes of the background grid 
! along the horizontal directions
array_name = "CLC%class_2D"
call allocate_de_int4_r2(.true.,CLC%class_2D,Grid%ncd(1),Grid%ncd(2),          &
   array_name,ulog_flag=.true.)
! 2D array/grid of z0 with the same sizes of the background grid along the 
! horizontal directions
array_name = "CLC%z0"
call allocate_de_dp_r2(.true.,CLC%z0,Grid%ncd(1),Grid%ncd(2),array_name,       &
   ulog_flag=.true.)
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
if (.not.ReadCheck(io_stat,ier,1,file_name,"z0_default",max_file_unit_booked+1,&
   ulog)) then
   write(uerr,*) "Error in reading the file ",file_name,". The execution ",    &
      "stops here."
   stop
endif
call open_close_file(.false.,max_file_unit_booked+1,file_name)
CLC%z0(:,:) = z0_default
! Open the 2D output file
file_name = "z0_CLC_output.txt"
call open_close_file(.true.,uCLC,file_name)
! 2D output file label
write(uCLC,'(3(a10),(a8))') "x(m)","y(m)","CLC_class","z0(m)" 
! Loop over the horizontal background-grid cells
!$omp parallel do default(none)                                                &
!$omp shared (Grid,n_neigh_hcell_CLCpol,CLC,CLC_active_cells,uCLC)             &
!$omp shared (neigh_hcell_CLCpol,NMAXPARTJ)                                    &
!$omp private(iy,ix,i_hcel,pos_cell,i_pol,i_pol_nei,i_face,test,aux_int,F_VUC) &
!$omp private(activation_flag)
do iy=1,Grid%ncd(2)
   do ix=1,Grid%ncd(1)
! Index of the horizontal grid cell
      i_hcel = CellNumber(ix,iy,1)
! Coordinates of the cell barycentre
      call cell_indices2pos(ix,iy,1,pos_cell(1:3))
! Initializions for the cell
      F_VUC = 1000000000
      activation_flag = .false.
! Loop over the neighbouring CLC polygons
      do_i_pol_nei: do i_pol_nei=1,n_neigh_hcell_CLCpol(i_hcel)
! Element index of the neighbouring CLC polygon
         aux_int = (i_hcel - 1) * NMAXPARTJ + i_pol_nei
         i_pol = neigh_hcell_CLCpol(aux_int)
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
               if (CLC%polygons(i_pol)%faces(i_face,4)<F_VUC) then
! So far, the current face is the only one selected or has the smallest vertex 
! usage code (in case of overlapping faces)
                  F_VUC = CLC%polygons(i_pol)%faces(i_face,4)
! Assign the class of the CLC polygon to the current element of the 
! CLC%class_2D array
                  CLC%class_2D(ix,iy) = CLC%polygons(i_pol)%CLC_class
! Assign z0
                  call z0_CLC_table(ix,iy)
                  if (activation_flag.eqv..false.) then
                     activation_flag = .true. 
!$omp critical (omp_CLC_active_cells)
! Update the count of CLC active cells
                     CLC_active_cells = CLC_active_cells + 1
!$omp end critical (omp_CLC_active_cells)
                  endif
               endif
            endif
         enddo
      enddo do_i_pol_nei
!$omp critical (omp_CLC_output)
! Update the 2D output file
      write(uCLC,'(2(f10.3),(i10),(f8.4))') pos_cell(1:2),CLC%class_2D(ix,iy), &
         CLC%z0(ix,iy)
!$omp end critical (omp_CLC_output)
   enddo
enddo
!$omp end parallel do
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
call allocate_de_int4_r2(.false.,CLC%class_2D,array_name=array_name,           &
   ulog_flag=.true.)
array_name = "n_neigh_hcell_CLCpol"
call allocate_de_int4_r1(.false.,n_neigh_hcell_CLCpol,array_name=array_name,   &
   ulog_flag=.true.)
array_name = "neigh_hcell_CLCpol"
call allocate_de_int4_r1(.false.,neigh_hcell_CLCpol,array_name=array_name,     &
   ulog_flag=.true.)
!$omp parallel do default(none)                                                &
!$omp shared(CLC)                                                              &
!$omp private(i_pol,array_name)
do i_pol=1,CLC%n_polygons
   write(array_name,*) i_pol
   array_name = "CLC%polygons(" // trim(adjustl(array_name)) // ")%vertices"
   call allocate_de_dp_r2(.false.,CLC%polygons(i_pol)%vertices,                &
      array_name=array_name,ulog_flag=.false.)
   write(array_name,*) i_pol
   array_name = "CLC%polygons(" // trim(adjustl(array_name)) // ")%faces"
   call allocate_de_int4_r2(.false.,CLC%polygons(i_pol)%faces,                 &
      array_name=array_name,ulog_flag=.false.)
   write(array_name,*) i_pol
   array_name = "CLC%polygons(" // trim(adjustl(array_name)) // ")%v_occurrence"
   call allocate_de_int4_r1(.false.,CLC%polygons(i_pol)%v_occurrence,          &
      array_name=array_name,ulog_flag=.false.)
enddo
!$omp end parallel do
array_name = "CLC%polygons"
call allocate_de_CLCp_r1(.false.,CLC%polygons,array_name=array_name,           &
   ulog_flag=.true.)
return
end subroutine z0_CLC
#endif
