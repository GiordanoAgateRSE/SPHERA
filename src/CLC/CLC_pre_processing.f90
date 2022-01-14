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
use Memory_I_O_interface_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: io_stat,ier,i_rec,i_vert,i_face,n_sides,file_line,I_O_unit,i_pol
integer(4) :: n_hcells,max_file_unit_ID,thread_id,n_threads,i_aux
character(100) :: file_name,array_name,aux_char,aux_char_2
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
   subroutine vertex_usage_code(V_0_1,V_0_2,V_0_3,F_VUC)
      implicit none
      integer(4),intent(in) :: V_0_1,V_0_2,V_0_3
      integer(4),intent(out) :: F_VUC
   end subroutine vertex_usage_code
   function OMP_get_num_threads()
      integer(4) :: OMP_get_num_threads
   end function OMP_get_num_threads
   function OMP_get_thread_num()
      integer(4) :: OMP_get_thread_num
   end function OMP_get_thread_num
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
thread_id = 1
n_threads = 1
!------------------------
! Statements
!------------------------
! Extracting the archive of the «.ply» CLC files
call system ("tar zxvf ./input/20_CLC/CLC_polygons_ply.tar.gz")
! Counting of the number of «.ply» CLC files
call system ("find CLC_polygons_ply/* | wc -l > test_CLC.txt")
file_name = "test_CLC.txt"
call open_close_file(.true.,max_file_unit_booked+1,file_name)
read(max_file_unit_booked+1,*,iostat=io_stat) CLC%n_polygons
if (.not.ReadCheck(io_stat,ier,1,file_name,"CLC%n_polygons",                   &
   max_file_unit_booked+1,ulog)) then
   write(uerr,*) "Error in reading the file ",file_name,". The execution ",    &
      "stops here."
   stop
endif
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
! Reading the CLC-polygon block numbers from the "CLC*.txt" file into the 
! associated fields of the CLC derived-type array, starting from the detection 
! of the exact file name
call system("ls input/20_CLC/CLC*.txt > CLC_txt_test.txt")
file_name = "CLC_txt_test.txt"
call open_close_file(.true.,max_file_unit_booked+1,file_name)
read(max_file_unit_booked+1,'(a13,a)',iostat=io_stat) aux_char,aux_char_2
if (.not.ReadCheck(io_stat,ier,1,file_name,file_name,                          &
   max_file_unit_booked+1,ulog)) then
   write(uerr,*) "Error in reading the file ",file_name,". The execution ",    &
      "stops here."
   stop
endif
call open_close_file(.false.,max_file_unit_booked+1,file_name)
call system("rm -f CLC_txt_test.txt")
file_name = "input/20_CLC/" // trim(adjustl(aux_char_2))
call open_close_file(.true.,max_file_unit_booked+1,file_name)
read(max_file_unit_booked+1,*,iostat=io_stat)
if (.not.ReadCheck(io_stat,ier,1,file_name,"file label",                       &
   max_file_unit_booked+1,ulog)) then
   write(uerr,*) "Error in reading the file ",file_name,". The execution ",    &
      "stops here."
   stop
endif
i_rec = 1
read(max_file_unit_booked+1,*,iostat=io_stat) CLC%polygons(i_rec)%ID,          &
   CLC%polygons(i_rec)%CLC_class
if (.not.ReadCheck(io_stat,ier,i_rec+1,file_name,"CLC classes",                &
   max_file_unit_booked+1,ulog)) then
   write(uerr,*) "Error in reading the file ",file_name,". The execution ",    &
      "stops here."
   stop
endif
i_rec = 2
do
! If the rest of the input file only contains doubles, the "do" construct is 
! interrupted to prevent a memory error (array-element index larger of 1 than 
! the number of elements)
   if (i_rec>CLC%n_polygons) exit
   read(max_file_unit_booked+1,*,iostat=io_stat) CLC%polygons(i_rec)%ID,       &
      CLC%polygons(i_rec)%CLC_class
! End-of-file condition
   if (io_stat<0) exit
   if (.not.ReadCheck(io_stat,ier,i_rec,file_name,"CLC classes",               &
      max_file_unit_booked+1,ulog)) then
      write(uerr,*) "Error in reading the file ",file_name,". The execution ", &
         "stops here."
      stop
   endif
   if (CLC%polygons(i_rec)%ID/=CLC%polygons(i_rec-1)%ID) i_rec = i_rec + 1
enddo
call open_close_file(.false.,max_file_unit_booked+1,file_name)
! Check on the maximum file unit (part 1 of 2)
file_name = "CLC polygons"
call check_max_file_unit_ID(max_file_unit_ID,max_file_unit_booked+1,file_name)
! Loop over the CLC polygons
!$omp parallel do default(none)                                                &
!$omp shared(CLC,ulog,max_file_unit_booked,uerr,n_threads)                     &
!$omp private(i_pol,file_name,I_O_unit,array_name,i_vert,io_stat,ier,i_face)   &
!$omp private(n_sides,file_line,thread_id)
do i_pol=1,CLC%n_polygons
! Check on the maximum file unit (part 2 of 2)
   if (i_pol==1) then
! Reading the number of threads
!$ n_threads = OMP_GET_NUM_THREADS()
      write(ulog,'(3a,i9)') "On the check for the file units of the CLC ",     &
         "polygons. The local number of OMP threads and of the file units ",   &
         "requested is: ",n_threads
   endif
! Reading the generic «.ply » CLC file into the associated fields of the CLC 
! derived-type array element
! First correct the CLC polygon ID, which was assigned as a Paraview block 
! number
   CLC%polygons(i_pol)%ID = CLC%polygons(i_pol)%ID - 1
   write(file_name,*) CLC%polygons(i_pol)%ID
   file_name = "./CLC_polygons_ply/CLC_" // trim(adjustl(file_name)) // ".ply"
! Read the headings of the on-going ".ply" file
! Read the OMP thread ID
!$ thread_id = OMP_GET_THREAD_NUM()
! Assign the file unit
   I_O_unit = max_file_unit_booked + thread_id
   call ply_headings(I_O_unit,file_name,CLC%polygons(i_pol)%n_vertices,        &
      CLC%polygons(i_pol)%n_faces)
! Allocation of the vertices of the CLC polygon
   write(array_name,*) i_pol
   array_name = "CLC%polygons(" // trim(adjustl(array_name)) // ")%vertices"
   call allocate_de_dp_r2(.true.,CLC%polygons(i_pol)%vertices,                 &
      CLC%polygons(i_pol)%n_vertices,2,array_name)
! Re-open the ".ply" file
   call open_close_file(.true.,I_O_unit,file_name)
! Skip the headings
   read(I_O_unit,'(10/)',iostat=io_stat)
! Read the vertices of the current CLC polygon
   do i_vert=1,CLC%polygons(i_pol)%n_vertices
      if (.not.ReadCheck(io_stat,ier,11+i_vert,file_name,array_name,I_O_unit,  &
         ulog)) then
         write(uerr,*) "Error in reading the file ",file_name,                 &
            ". The execution stops here."
         stop
      endif
      read(I_O_unit,*,iostat=io_stat) CLC%polygons(i_pol)%vertices(i_vert,1:2)
! No need for a critical section with "ulog" as it is used only in case of 
! errors
      if (.not.ReadCheck(io_stat,ier,11+i_vert,file_name,array_name,I_O_unit,  &
         ulog)) then
         write(uerr,*) "Error in reading the file ",file_name,                 &
            ". The execution stops here."
         stop
      endif
   enddo
! Allocation and initialization of the array of the vertex occurrence
   array_name = "CLC%polygons(" // trim(adjustl(array_name)) // ")%v_occurrence"
   call allocate_de_int4_r1(.true.,CLC%polygons(i_pol)%v_occurrence,           &
      CLC%polygons(i_pol)%n_vertices,array_name)
   CLC%polygons(i_pol)%v_occurrence(:) = 0
! Allocation of the faces of the CLC polygon
   write(array_name,*) i_pol
   array_name = "CLC%polygons(" // trim(adjustl(array_name)) // ")%faces"
   call allocate_de_int4_r2(.true.,CLC%polygons(i_pol)%faces,                  &
      CLC%polygons(i_pol)%n_faces,4,array_name)
! Read the faces of the current CLC polygon
   do i_face=1,CLC%polygons(i_pol)%n_faces
      read(I_O_unit,*,iostat=io_stat) n_sides,                                 &
         CLC%polygons(i_pol)%faces(i_face,1:3)
! ".ply" files count vertex ID starting from zero: "1" has to be added
      CLC%polygons(i_pol)%faces(i_face,1:3) =                                  &
         CLC%polygons(i_pol)%faces(i_face,1:3) + 1
      file_line = 11 + CLC%polygons(i_pol)%n_vertices + i_face
! No need for a critical section with "ulog" as it is used only in case of 
! errors
      if (.not.ReadCheck(io_stat,ier,file_line,file_name,array_name,I_O_unit,  &
         ulog)) then
         write(uerr,*) "Error in reading the file ",file_name,                 &
            ". The execution stops here."
         stop
      endif
! Check update: the generic CLC polygon has to be partitioned in triangles only
      if (n_sides>3) then
         write(uerr,*) "The CLC polygon ",CLC%polygons(i_pol)%ID,              &
            " is not partitioned in triangles only. File line: ",file_line,    &
            ". The execution stops here."
         stop
      endif
! Update the array of the vertex occurrence
      do i_aux=1,3
   CLC%polygons(i_pol)%v_occurrence(CLC%polygons(i_pol)%faces(i_face,i_aux)) = &
      CLC%polygons(i_pol)%v_occurrence(CLC%polygons(i_pol)%faces(i_face,i_aux))&
      + 1
      enddo
   enddo
! Re-close the ".ply" file
   call open_close_file(.false.,I_O_unit,file_name)
! Compute the "vertex usage code" of all the faces of the current CLC 
! triangulated polygon
   do i_face=1,CLC%polygons(i_pol)%n_faces
      call vertex_usage_code(                                                  &
         CLC%polygons(i_pol)%v_occurrence(CLC%polygons(i_pol)%faces(i_face,1)),&
         CLC%polygons(i_pol)%v_occurrence(CLC%polygons(i_pol)%faces(i_face,2)),&
         CLC%polygons(i_pol)%v_occurrence(CLC%polygons(i_pol)%faces(i_face,3)),&
         CLC%polygons(i_pol)%faces(i_face,4))
   enddo
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
   write(ulog,'(4(a10))') "ID_vertex","x(m)","y(m)","V_o"
   do i_vert=1,CLC%polygons(i_pol)%n_vertices
      write(ulog,'(i10,2(f10.3),i10)') i_vert,                                 &
         CLC%polygons(i_pol)%vertices(i_vert,1),                               &
         CLC%polygons(i_pol)%vertices(i_vert,2),                               &
         CLC%polygons(i_pol)%v_occurrence(i_vert)
   enddo
   write(ulog,'(5(a12))') "ID_face","ID_vertex_1","ID_vertex_2","ID_vertex_3", &
      "F_VUC"
   do i_face=1,CLC%polygons(i_pol)%n_faces
      write(ulog,'(5(i12))') i_face,                                           &
         CLC%polygons(i_pol)%faces(i_face,1),                                  &
         CLC%polygons(i_pol)%faces(i_face,2),                                  &
         CLC%polygons(i_pol)%faces(i_face,3),                                  &
         CLC%polygons(i_pol)%faces(i_face,4)
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
