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
! Program unit: CLC_pre_processing
! Description: Pre-processsing of the CLC (Corine Land Cover) input data
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine CLC_pre_processing
!------------------------
! Modules
!------------------------
use I_O_file_module
use Hybrid_allocation_module
use Static_allocation_module                            
use Dynamic_allocation_module
use Memory_I_O_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: io_stat,ier,i_rec,i_vert,i_face,n_sides,file_line,I_O_unit,i_pol
integer(4) :: n_hcells
character(100) :: file_name,array_name
logical,external :: ReadCheck
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine ply_headings(I_O_unit,file_name,n_vertices,n_faces)
      implicit none
      integer(4),intent(in) :: I_O_unit
      character(100),intent(in) :: file_name
      integer(4),intent(out) :: n_vertices,n_faces
   end subroutine ply_headings
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
! Extracting the archive of the «.ply» CLC files
call system ("tar zxvf ./input/20_CLC/CLC_*ply.tar.gz")
! Counting of the number of «.ply» CLC files
call system ("find CLC_*ply/* | wc -l > test_CLC.txt")
file_name = "test_CLC.txt"
call open_close_file(.true.,max_file_unit_booked+1,file_name)
read(max_file_unit_booked+1,*,iostat=io_stat) CLC%n_polygons
if (.not.ReadCheck(io_stat,ier,1,"test_CLC.txt","CLC%n_polygons",              &
   max_file_unit_booked+1,ulog)) return
call open_close_file(.false.,max_file_unit_booked+1,file_name)
call system ("rm -f test_CLC.txt")
if (CLC%n_polygons==0) then
! In case no ".ply" file is available, the execution stops here
   write(uerr,*) "There is no CLC polygon. The execution stops here."
   stop
endif
! Writing the input CLC data in the log file (first part)
write(ulog,*) "----------------------------------------------------------------"
write(ulog,*) "CLC polygons"
write(ulog,*) "Number of CLC polygons: ",CLC%n_polygons
file_name = "CLC polygons"
call check_max_file_unit_ID(CLC%n_polygons+max_file_unit_booked,file_name)
! Allocate and initialize the array of the number of the neighbouring CLC 
! polygons with respect to the reference backgorund-grid cell
n_hcells = Grid%ncd(1) * Grid%ncd(2)
array_name = "n_neigh_hcell_CLCpol"
call allocate_de_int4_r1(.true.,n_neigh_hcell_CLCpol,n_hcells,array_name)
n_neigh_hcell_CLCpol(:) = 0
! Allocate and initialize the neighbouring list linking the horizontal 
! background-grid cells with the CLC polygons
array_name = "neigh_hcell_CLCpol"
call allocate_de_int4_r1(.true.,neigh_hcell_CLCpol,n_hcells*NMAXPARTJ,         &
   array_name)
neigh_hcell_CLCpol(:) = 0
! Allocate the derived-type allocatable array of the CLC polygons
array_name = "CLC_polygons"
call allocate_de_CLCp_r1(.true.,CLC%polygons,CLC%n_polygons,array_name)
! Reading the CLC-polygon IDs from the «CLC*.txt» file into the associated 
! fields of the CLC derived-type array
file_name = "./input/20_CLC/CLC_*.txt"
call open_close_file(.true.,max_file_unit_booked+1,file_name)
read(max_file_unit_booked+1,*,iostat=io_stat)
if (.not.ReadCheck(io_stat,ier,1,file_name,"file label",                       &
   max_file_unit_booked+1,ulog)) return
i_rec = 1
read(max_file_unit_booked+1,*,iostat=io_stat) CLC%polygons(i_rec)%ID,          &
   CLC%polygons(i_rec)%CLC_class
if (.not.ReadCheck(io_stat,ier,i_rec+1,file_name,"CLC classes",                &
   max_file_unit_booked+1,ulog)) return
i_rec = 2
do
! If the rest of the input file only contains doubles, the "do" construct is 
! interrupted to prevent a memory error (array-element index larger of 1 than 
! the maximum element ID)
   if (i_rec>CLC%n_polygons) exit
   read(max_file_unit_booked+1,*,iostat=io_stat) CLC%polygons(i_rec)%ID,       &
      CLC%polygons(i_rec)%CLC_class
! End-of-file condition
   if (io_stat<0) exit
   if (.not.ReadCheck(io_stat,ier,i_rec,"CLC.txt","CLC classes",               &
      max_file_unit_booked+1,ulog)) return
   if (CLC%polygons(i_rec)%ID/=CLC%polygons(i_rec-1)%ID) i_rec = i_rec + 1
enddo
call open_close_file(.false.,max_file_unit_booked+1,file_name)
! Loop over the CLC polygons
!$omp parallel do default(none)                                                &
!$omp shared(CLC,ulog,max_file_unit_booked,uerr)                               &
!$omp private(i_pol,file_name,I_O_unit,array_name,i_vert,io_stat,ier,i_face)   &
!$omp private(n_sides,file_line)
do i_pol=1,CLC%n_polygons
! Reading the generic «.ply » CLC file into the associated fields of the CLC 
! derived-type array element
   write(file_name,*) CLC%polygons(i_pol)%ID
   file_name = "./CLC_*ply/" // trim(file_name) // ".ply"
! Read the headings of the on-going ".ply" file
   I_O_unit = max_file_unit_booked + i_pol
   call ply_headings(I_O_unit,file_name,CLC%polygons(i_pol)%n_vertices,        &
      CLC%polygons(i_pol)%n_faces)
! Allocation of the vertices of the CLC polygon
   write(array_name,*) i_pol
   array_name = "CLC%polygons(" // trim(array_name) // ")%vertices"
   call allocate_de_dp_r2(.true.,CLC%polygons(i_pol)%vertices,                 &
      CLC%polygons(i_pol)%n_vertices,2,array_name)
! Re-open the ".ply" file
   call open_close_file(.true.,I_O_unit,file_name)
! Read the vertices of the current CLC polygon
   do i_vert=1,CLC%polygons(i_pol)%n_vertices
      read(I_O_unit,*,iostat=io_stat) CLC%polygons(i_pol)%vertices(i_vert,1:2)
! No need for a critical section with "ulog" as it is used only in case of 
! errors
      if (.not.ReadCheck(io_stat,ier,11+i_vert,file_name,array_name,I_O_unit,  &
         ulog)) exit
   enddo
! Allocation of the faces of the CLC polygon
   write(array_name,*) i_pol
   array_name = "CLC%polygons(" // trim(array_name) // ")%faces"
   call allocate_de_int4_r2(.true.,CLC%polygons(i_pol)%faces,                  &
      CLC%polygons(i_pol)%n_faces,3,array_name)
! Read the faces of the current CLC polygon
   do i_face=1,CLC%polygons(i_pol)%n_faces
      read(I_O_unit,*,iostat=io_stat) n_sides,                                 &
         CLC%polygons(i_pol)%faces(i_face,1:3)
      file_line = 11 + CLC%polygons(i_pol)%n_vertices + i_face
! No need for a critical section with "ulog" as it is used only in case of 
! errors
      if (.not.ReadCheck(io_stat,ier,file_line,file_name,array_name,I_O_unit,  &
         ulog)) exit
! Check update: the generic CLC polygon has to be partitioned in triangles only
      if (n_sides>3) then
         write(uerr,*) "The CLC polygon ",CLC%polygons(i_pol)%ID,              &
            " is not partitioned in triangles only. File line: ",file_line,    &
            ". The execution stops here."
         stop
      endif
   enddo
! Re-close the ".ply" file
   call open_close_file(.false.,I_O_unit,file_name)
! Contribution to the neighbouring search
   call neighbouring_hcell_CLCp(i_pol)
enddo
!$omp end parallel do
! Loop over the CLC polygons
do i_pol=1,CLC%n_polygons
! Writing the input CLC data in the log file (second part)
   write(ulog,'(4(a12))') "Polygon_ID","CLC_class","n_vertices","n_faces"
   write(ulog,'(4(i12))') CLC%polygons(i_pol)%ID,                              &
      CLC%polygons(i_pol)%CLC_class,CLC%polygons(i_pol)%n_vertices,            &
      CLC%polygons(i_pol)%n_faces
   write(ulog,'(3(a10))') "ID_vertex","x(m)","y(m)"
   do i_vert=1,CLC%polygons(i_pol)%n_vertices
      write(ulog,'(i10,2(f10.3))') i_vert,                                     &
         CLC%polygons(i_pol)%vertices(i_vert,1),                               &
         CLC%polygons(i_pol)%vertices(i_vert,2)
   enddo
   write(ulog,'(4(a12))') "ID_face","ID_vertex_1","ID_vertex_2","ID_vertex_3"
   do i_face=1,CLC%polygons(i_pol)%n_faces
      write(ulog,'(4(i12))') i_face,                                           &
         CLC%polygons(i_pol)%faces(i_face,1),                                  &
         CLC%polygons(i_pol)%faces(i_face,2),                                  &
         CLC%polygons(i_pol)%faces(i_face,3)
   enddo
enddo
write(ulog,*) "----------------------------------------------------------------"
! Removing the temporary archive extraction
call system ("rm -rf CLC_*ply")
!------------------------
! Deallocations
!------------------------
return
end subroutine CLC_pre_processing
#endif
