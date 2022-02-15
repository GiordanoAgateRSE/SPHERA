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
! Program unit: Import_ply_surface_meshes
! Description: To import the surface meshes (generated by SnappyHexMesh
!              -OpenFOAM-), as converted by Paraview into .ply files.
!              This subroutine is mandatory and activated only for the DB-SPH
!              boundary treatment scheme. In 3D, SPHERA (DBSPH) works with
!              triangular faces, in 2D with quadrilateral faces. Input .ply
!              files must be triangular/quadrilateral/pentagonal/hexagonal
!              meshes in 3D or square meshes in 2D.               
!-------------------------------------------------------------------------------
subroutine Import_ply_surface_meshes
!------------------------
! Modules
!------------------------
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use Memory_I_O_interface_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: file_stat,n_vertices,old_size_vert,old_size_face,new_size_vert
integer(4) :: new_size_face,n_faces,face_vert_num,j,k,surface_mesh_file_ID
integer(4) :: n_faces_aux
integer(4) :: aux_face_vert(6)
character(100) :: file_name,array_name
integer(4),dimension(:),allocatable :: aux_surface_mesh_file_ID
type(vertex_der_type),dimension(:),allocatable :: aux_der_type_vert
type(face_der_type),dimension(:),allocatable :: aux_der_type_faces
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine area_triangle(P1,P2,P3,area,normal)
      implicit none
      double precision,intent(in)    :: P1(3),P2(3),P3(3)
      double precision,intent(out)   :: area
      double precision,intent(out)   :: normal(3)
   end subroutine area_triangle
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
new_size_face = 0
surface_mesh_file_ID = 0
!------------------------
! Statements
!------------------------
! Open the file name list
file_name = "surface_mesh_list.inp"
call open_close_file(.true.,unit_file_list,file_name)
read(unit_file_list,*,IOSTAT=file_stat)
if (file_stat/=0) then
   write(uerr,*) 'Error in reading surface_mesh_list.inp in ',                 &
      'Import_pl_surface_meshes; the program stops here'
   stop
endif
do
! To increment the file_number
   surface_mesh_file_ID = surface_mesh_file_ID + 1
! Read the file name    
   read(unit_file_list,'(a)',IOSTAT=file_stat) file_name
! Exit the cicle at the end of file
   if (file_stat/=0) exit 
   file_name = trim(adjustl(file_name))
! Read the headings of the on-going mesh file
   call ply_headings(unit_DBSPH_mesh,file_name,n_vertices,n_faces)
   call open_close_file(.true.,unit_DBSPH_mesh,file_name)
   if (.not.allocated(DBSPH%surf_mesh%vertices)) then
      array_name = "DBSPH%surf_mesh%vertices"
      call allocate_de_vertex_r1(.true.,DBSPH%surf_mesh%vertices,n_vertices,   &
         array_name,ulog_flag=.true.)
      array_name = "aux_der_type_vert"
      call allocate_de_vertex_r1(.true.,aux_der_type_vert,n_vertices,          &
         array_name,ulog_flag=.true.)
      old_size_vert = 0
      else
         old_size_vert = size(DBSPH%surf_mesh%vertices)
         new_size_vert = old_size_vert + n_vertices
         aux_der_type_vert(:) = DBSPH%surf_mesh%vertices(:)
         array_name = "DBSPH%surf_mesh%vertices"
         call allocate_de_vertex_r1(.false.,DBSPH%surf_mesh%vertices,          &
            array_name=array_name,ulog_flag=.true.)
         call allocate_de_vertex_r1(.true.,DBSPH%surf_mesh%vertices,           &
            new_size_vert,array_name,ulog_flag=.true.)
         DBSPH%surf_mesh%vertices(1:old_size_vert) = aux_der_type_vert(:)
         array_name = "aux_der_type_vert"
         call allocate_de_vertex_r1(.false.,aux_der_type_vert,                 &
            array_name=array_name,ulog_flag=.true.)
         call allocate_de_vertex_r1(.true.,aux_der_type_vert,new_size_vert,    &
            array_name,ulog_flag=.true.)
   endif
! Read the vertex coordinates: start      
   do j=(old_size_vert+1),(old_size_vert+n_vertices)
      read(unit_DBSPH_mesh,*) DBSPH%surf_mesh%vertices(j)%pos(:)
   enddo
! To allocate or resize DBSPH%surf_mesh%faces on the maximum number of faces
! (worst case with all hexagonal faces). To allocate                           &
! DBSPH%surf_mesh%surface_mesh_file_ID and auxiliary arrays, accordingly.
   if (.not.allocated(DBSPH%surf_mesh%faces)) then
#ifdef SPACE_3D
      n_faces_aux = (DBSPH%ply_n_face_vert - 2) * n_faces
#elif defined SPACE_2D
      n_faces_aux = n_faces
#endif
      array_name = "DBSPH%surf_mesh%faces"
      call allocate_de_face_r1(.true.,DBSPH%surf_mesh%faces,n_faces_aux,       &
         array_name,ulog_flag=.true.)
      array_name = "aux_der_type_faces"
      call allocate_de_face_r1(.true.,aux_der_type_faces,n_faces_aux,          &
         array_name,ulog_flag=.true.)
      array_name = "DBSPH%surf_mesh%surface_mesh_file_ID"
      call allocate_de_int4_r1(.true.,DBSPH%surf_mesh%surface_mesh_file_ID,    &
         n_faces_aux,array_name,ulog_flag=.true.)
      array_name = "aux_surface_mesh_file_ID"
      call allocate_de_int4_r1(.true.,aux_surface_mesh_file_ID,n_faces_aux,    &
         array_name,ulog_flag=.true.)
      old_size_face = 0
      else
         old_size_face = size(DBSPH%surf_mesh%faces)
#ifdef SPACE_3D
            new_size_face = old_size_face + (DBSPH%ply_n_face_vert - 2) *      &
                            n_faces
#elif defined SPACE_2D
               new_size_face = old_size_face + n_faces
#endif
         aux_der_type_faces(:) = DBSPH%surf_mesh%faces(:)
         array_name = "DBSPH%surf_mesh%faces"
         call allocate_de_face_r1(.false.,DBSPH%surf_mesh%faces,               &
            array_name=array_name,ulog_flag=.true.)
         call allocate_de_face_r1(.true.,DBSPH%surf_mesh%faces,new_size_face,  &
            array_name,ulog_flag=.true.)
         DBSPH%surf_mesh%faces(1:old_size_face) = aux_der_type_faces(:)
         array_name = "aux_der_type_faces"
         call allocate_de_face_r1(.false.,aux_der_type_faces,                  &
            array_name=array_name,ulog_flag=.true.)
         call allocate_de_face_r1(.true.,aux_der_type_faces,new_size_face,     &
            array_name,ulog_flag=.true.)
         aux_surface_mesh_file_ID(:) = DBSPH%surf_mesh%surface_mesh_file_ID(:)
         array_name = "DBSPH%surf_mesh%surface_mesh_file_ID"
         call allocate_de_int4_r1(.false.,DBSPH%surf_mesh%surface_mesh_file_ID,&
            array_name=array_name,ulog_flag=.true.)
         call allocate_de_int4_r1(.true.,DBSPH%surf_mesh%surface_mesh_file_ID, &
            new_size_face,array_name,ulog_flag=.true.)
         DBSPH%surf_mesh%surface_mesh_file_ID(1:old_size_face) =               &
            aux_surface_mesh_file_ID(:)
         array_name = "aux_surface_mesh_file_ID"
         call allocate_de_int4_r1(.false.,aux_surface_mesh_file_ID,            &
            array_name=array_name,ulog_flag=.true.)
         call allocate_de_int4_r1(.true.,aux_surface_mesh_file_ID,             &
            new_size_face,array_name,ulog_flag=.true.)
   endif
! Read the face vertices: start
   k = old_size_face + 1
   do j=1,n_faces
      read(unit_DBSPH_mesh,*) face_vert_num,aux_face_vert(1:face_vert_num)
! Assignation of vertices with eventual conversion of any 4/5/6-side face
! into 2/3/4 triangles; computation of area and normal
#ifdef SPACE_3D
         select case (face_vert_num)
            case(3)
! To import vertices of the triangular face
! Face 1: vertices 1,2,3
! ".ply" files count vertex ID starting from zero: "1" has to be added
               DBSPH%surf_mesh%faces(k)%vert_list(1:3) = old_size_vert +       &
                                                         aux_face_vert(1:3) + 1
               DBSPH%surf_mesh%faces(k)%vert_list(4) = 0
               DBSPH%surf_mesh%surface_mesh_file_ID(k) = surface_mesh_file_ID
               k = k+1
! Computation of area and normal of the face
               call area_triangle(                                             &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(1))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(2))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(3))%pos,         &
                  DBSPH%surf_mesh%faces(k-1)%area,                             &
                  DBSPH%surf_mesh%faces(k-1)%normal)
            case(4)
! To import vertices of the quadrilateral face, split in 2 triangular faces
! Face 1: vertices 1,2,3
            DBSPH%surf_mesh%faces(k)%vert_list(1:3) = old_size_vert +          &
                                                      aux_face_vert(1:3) + 1
            DBSPH%surf_mesh%faces(k)%vert_list(4) = 0
            DBSPH%surf_mesh%surface_mesh_file_ID(k) = surface_mesh_file_ID
! Face 2: vertices 1,3,4
            DBSPH%surf_mesh%faces(k+1)%vert_list(1) = old_size_vert +          &
                                                      aux_face_vert(1) + 1
            DBSPH%surf_mesh%faces(k+1)%vert_list(2:3) = old_size_vert +        &
                                                       aux_face_vert(3:4) + 1
            DBSPH%surf_mesh%faces(k+1)%vert_list(4) = 0
            DBSPH%surf_mesh%surface_mesh_file_ID(k+1) = surface_mesh_file_ID
            k = k+2
! Computation of area and normal of the 2 faces
            call area_triangle(                                                &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-2)%vert_list(1))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-2)%vert_list(2))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-2)%vert_list(3))%pos,         &
               DBSPH%surf_mesh%faces(k-2)%area,                                &
               DBSPH%surf_mesh%faces(k-2)%normal)
            call area_triangle(                                                &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(1))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(2))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(3))%pos,         &
               DBSPH%surf_mesh%faces(k-1)%area,                                &
               DBSPH%surf_mesh%faces(k-1)%normal)
            case(5)
! To import vertices of the pentagonal face, split in 3 triangular faces
! Face 1: vertices 1,2,5
            DBSPH%surf_mesh%faces(k)%vert_list(1:2) = old_size_vert +          &
                                                      aux_face_vert(1:2) + 1
            DBSPH%surf_mesh%faces(k)%vert_list(3) = old_size_vert +            &
                                                    aux_face_vert(5) + 1
            DBSPH%surf_mesh%faces(k)%vert_list(4) = 0
            DBSPH%surf_mesh%surface_mesh_file_ID(k) = surface_mesh_file_ID
! Face 2: vertices 2,3,5
            DBSPH%surf_mesh%faces(k+1)%vert_list(1:2) = old_size_vert +        &
                                                        aux_face_vert(2:3) + 1
            DBSPH%surf_mesh%faces(k+1)%vert_list(3)= old_size_vert +           &
                                                     aux_face_vert(5) + 1
            DBSPH%surf_mesh%faces(k+1)%vert_list(4) = 0
            DBSPH%surf_mesh%surface_mesh_file_ID(k+1) = surface_mesh_file_ID
! Face 3: vertices 3,4,5
            DBSPH%surf_mesh%faces(k+2)%vert_list(1:3) = old_size_vert +        &
                                                        aux_face_vert(3:5) + 1
            DBSPH%surf_mesh%faces(k+2)%vert_list(4) = 0
            DBSPH%surf_mesh%surface_mesh_file_ID(k+2) = surface_mesh_file_ID
            k = k+3
! Computation of area and normal of the 3 faces
            call area_triangle(                                                &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-3)%vert_list(1))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-3)%vert_list(2))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-3)%vert_list(3))%pos,         &
               DBSPH%surf_mesh%faces(k-3)%area,                                &
               DBSPH%surf_mesh%faces(k-3)%normal)
            call area_triangle(                                                &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-2)%vert_list(1))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-2)%vert_list(2))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-2)%vert_list(3))%pos,         &
               DBSPH%surf_mesh%faces(k-2)%area,                                &
               DBSPH%surf_mesh%faces(k-2)%normal)
            call area_triangle(                                                &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(1))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(2))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(3))%pos,         &
               DBSPH%surf_mesh%faces(k-1)%area,                                &
               DBSPH%surf_mesh%faces(k-1)%normal)
            case(6)
! To import vertices of the hexagonal face, split in 4 triangular faces
! Face 1: vertices 1,2,6
            DBSPH%surf_mesh%faces(k)%vert_list(1:2) = old_size_vert +          &
                                                      aux_face_vert(1:2) + 1
            DBSPH%surf_mesh%faces(k)%vert_list(3) = old_size_vert +            &
                                                    aux_face_vert(6) + 1
            DBSPH%surf_mesh%faces(k)%vert_list(4) = 0
            DBSPH%surf_mesh%surface_mesh_file_ID(k) = surface_mesh_file_ID
! Face 2: vertices 2,5,6
            DBSPH%surf_mesh%faces(k+1)%vert_list(1) = old_size_vert +          &
                                                      aux_face_vert(2) + 1
            DBSPH%surf_mesh%faces(k+1)%vert_list(2:3)= old_size_vert +         &
                                                       aux_face_vert(5:6) + 1
            DBSPH%surf_mesh%faces(k+1)%vert_list(4) = 0
            DBSPH%surf_mesh%surface_mesh_file_ID(k+1) = surface_mesh_file_ID
! Face 3: vertices 2,3,5
            DBSPH%surf_mesh%faces(k+2)%vert_list(1:2) = old_size_vert +        &
                                                        aux_face_vert(2:3) + 1
            DBSPH%surf_mesh%faces(k+2)%vert_list(3)= old_size_vert +           &
                                                     aux_face_vert(5) + 1
            DBSPH%surf_mesh%faces(k+2)%vert_list(4) = 0
            DBSPH%surf_mesh%surface_mesh_file_ID(k+2) = surface_mesh_file_ID
! Face 4: vertices 3,4,5
            DBSPH%surf_mesh%faces(k+3)%vert_list(1:3) = old_size_vert +        &
                                                        aux_face_vert(3:5) + 1
            DBSPH%surf_mesh%faces(k+3)%vert_list(4) = 0
            DBSPH%surf_mesh%surface_mesh_file_ID(k+3) = surface_mesh_file_ID
            k = k+4
! Computation of area and normal of the 4 faces
            call area_triangle(                                                &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-4)%vert_list(1))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-4)%vert_list(2))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-4)%vert_list(3))%pos,         &
               DBSPH%surf_mesh%faces(k-4)%area,                                &
               DBSPH%surf_mesh%faces(k-4)%normal)
            call area_triangle(                                                &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-3)%vert_list(1))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-3)%vert_list(2))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-3)%vert_list(3))%pos,         &
               DBSPH%surf_mesh%faces(k-3)%area,                                &
               DBSPH%surf_mesh%faces(k-3)%normal)
            call area_triangle(                                                &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-2)%vert_list(1))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-2)%vert_list(2))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-2)%vert_list(3))%pos,         &
               DBSPH%surf_mesh%faces(k-2)%area,                                &
               DBSPH%surf_mesh%faces(k-2)%normal)
            call area_triangle(                                                &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(1))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(2))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(3))%pos,         &
               DBSPH%surf_mesh%faces(k-1)%area,                                &
               DBSPH%surf_mesh%faces(k-1)%normal)
         endselect
#elif defined SPACE_2D
            DBSPH%surf_mesh%faces(k)%vert_list(1:4) = old_size_vert +          &
                                                      aux_face_vert(1:4) + 1
            DBSPH%surf_mesh%surface_mesh_file_ID(k) = surface_mesh_file_ID
            k = k + 1
! Computation of normal (area will be re-written)             
            call area_triangle(                                                &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(1))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(2))%pos,         &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(3))%pos,         &
               DBSPH%surf_mesh%faces(k-1)%area,                                &
               DBSPH%surf_mesh%faces(k-1)%normal)
! Computation of area in 2D (segment length)
            DBSPH%surf_mesh%faces(k-1)%area = dsqrt(dot_product(               &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(1))%pos          &
- DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(3))%pos,       &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(1))%pos          &
- DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(3))%pos)) /    &
                                              dsqrt(2.d0)
#endif
   enddo
   call open_close_file(.false.,unit_DBSPH_mesh,file_name)
! Read the face vertices: end   
! Resize DBSPH%surf_mesh%faces on the actual number of faces
! new_size_face = estimated_new_size_face - local overestimation
#ifdef SPACE_3D
      new_size_face = size(DBSPH%surf_mesh%faces) - ((DBSPH%ply_n_face_vert - &
                      2) * n_faces - (k - 1 - old_size_face))
#elif defined  SPACE_2D
         new_size_face = size(DBSPH%surf_mesh%faces) - (n_faces - (k - 1 -     &
                         old_size_face))
#endif
   old_size_face = size(DBSPH%surf_mesh%faces)
   if (new_size_face<old_size_face) then
      aux_der_type_faces(:) = DBSPH%surf_mesh%faces(:)
      array_name = "DBSPH%surf_mesh%faces"
      call allocate_de_face_r1(.false.,DBSPH%surf_mesh%faces,                  &
         array_name=array_name,ulog_flag=.true.)
      call allocate_de_face_r1(.true.,DBSPH%surf_mesh%faces,new_size_face,     &
         array_name,ulog_flag=.true.)
      DBSPH%surf_mesh%faces(:) = aux_der_type_faces(1:new_size_face)
      array_name = "aux_der_type_faces"
      call allocate_de_face_r1(.false.,aux_der_type_faces,                     &
         array_name=array_name,ulog_flag=.true.)
      call allocate_de_face_r1(.true.,aux_der_type_faces,new_size_face,        &
         array_name,ulog_flag=.true.)
      aux_surface_mesh_file_ID(:) = DBSPH%surf_mesh%surface_mesh_file_ID(:)
      array_name = "DBSPH%surf_mesh%surface_mesh_file_ID"
      call allocate_de_int4_r1(.false.,DBSPH%surf_mesh%surface_mesh_file_ID,   &
         array_name=array_name,ulog_flag=.true.)
      call allocate_de_int4_r1(.true.,DBSPH%surf_mesh%surface_mesh_file_ID, &
         new_size_face,array_name,ulog_flag=.true.)
      DBSPH%surf_mesh%surface_mesh_file_ID(:) =                                &
         aux_surface_mesh_file_ID(1:new_size_face)
      array_name = "aux_surface_mesh_file_ID"
      call allocate_de_int4_r1(.false.,aux_surface_mesh_file_ID,               &
         array_name=array_name,ulog_flag=.true.)
      call allocate_de_int4_r1(.true.,aux_surface_mesh_file_ID,new_size_face,  &
         array_name,ulog_flag=.true.)
   endif
enddo
file_name = "surface_mesh_list.inp"
call open_close_file(.false.,unit_file_list,file_name)
! Initializing the number of surface elements
DBSPH%n_w = new_size_face 
!------------------------
! Deallocations
!------------------------
call allocate_de_vertex_r1(.false.,aux_der_type_vert,array_name=array_name,    &
   ulog_flag=.true.)
call allocate_de_face_r1(.false.,aux_der_type_faces,array_name=array_name,     &
   ulog_flag=.true.)
call allocate_de_int4_r1(.false.,aux_surface_mesh_file_ID,                     &
   array_name=array_name,ulog_flag=.true.)
return
end subroutine Import_ply_surface_meshes
