!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2016 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: ply2SPHERA_perimeter
! Description: Format conversion from .ply to the format of VERTICES and FACES 
!              of SPHERA main input file.        
!-------------------------------------------------------------------------------
program ply2SPHERA_perimeter
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: open_stat,read_stat,close_stat,write_stat
integer(4) :: perimeter_first_vertex_ID,perimeter_first_face_ID,perimeter_ID
integer(4) :: face_ID,vertex_ID,n_vertices,n_faces,face_sides,i_vert,i_face,i
double precision :: vertex_pos(3),z_offset
integer(4) :: vertex_IDs_of_ply_face(6)
character(100) :: ply_file_name,char_aux_2,char_aux
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
! To read ply2SPHERA_perimeter.inp
open(14,file="ply2SPHERA_perimeter.inp",IOSTAT=open_stat)
if (open_stat/=0) then
   write(0,*) "Error in opening ply2SPHERA_perimeter.inp. ",                   &
      "The program stops here. "
   stop
endif
read(14,*,IOSTAT=read_stat) perimeter_first_vertex_ID,perimeter_first_face_ID, &
   perimeter_ID
if (read_stat/=0) then
   write(0,*) "Error in reading the first line of ply2SPHERA_perimeter.inp. ", &
      "The program stops here. "
   stop
endif
read(14,*,IOSTAT=read_stat) z_offset
if (read_stat/=0) then
   write(0,*) "Error in reading the second line of ply2SPHERA_perimeter.inp. ",&
      "The program stops here. "
   stop
endif
read(14,*,IOSTAT=read_stat) ply_file_name
if (read_stat/=0) then
   write(0,*) "Error in reading the second line of ply2SPHERA_perimeter.inp. ",&
      "The program stops here. "
   stop
endif
close(14,IOSTAT=close_stat)
if (close_stat/=0) then
   write(0,*) "Error in closing ply2SPHERA_perimeter.inp. ",                   &
      "The program stops here. "
   stop
endif
! To convert data format from .ply to .txt files for SPHERA VERTICES and FACES. 
open(11,file=trim(ply_file_name),IOSTAT=open_stat)
if (open_stat/=0) then
   write(0,*) "Error in opening .ply input file. The program stops here. "
   stop
endif
open(12,file="perimeter_vertices.txt",IOSTAT=open_stat)
if (open_stat/=0) then
   write(0,*) "Error in opening perimeter_vertices.txt. The program stops here."
   stop
endif
open(13,file="perimeter_faces.txt",IOSTAT=open_stat)
if (open_stat/=0) then
   write(0,*) "Error in opening perimeter_faces.txt. The program stops here. "
   stop
endif
read(11,'(/,/,/)',IOSTAT=read_stat)
if (read_stat/=0) then
   write(0,*) "Error in reading the lines 1-4 of .ply input file. ",           &
      "The program stops here. "
   stop
endif
read(11,*,IOSTAT=read_stat) char_aux,char_aux_2,n_vertices
if (read_stat/=0) then
   write(0,*) "Error in reading line 5 of .ply input file. ",                  &
      "The program stops here. "
   stop
endif
read(11,'(/,/)',IOSTAT=read_stat)
if (read_stat/=0) then
   write(0,*) "Error in reading lines 6-8 of .ply input file. ",               &
      "The program stops here. "
   stop
endif
read(11,*,IOSTAT=read_stat) char_aux,char_aux_2,n_faces
if (read_stat/=0) then
   write(0,*) "Error in reading line 9 of .ply input file. ",                  &
      "The program stops here. "
   stop
endif
read(11,'(/)',IOSTAT=read_stat)
if (read_stat/=0) then
   write(0,*) "Error in reading lines 10-11 of .ply input file. ",             &
      "The program stops here. "
   stop
endif
do i_vert=1,n_vertices
   read(11,*,IOSTAT=read_stat) vertex_pos(:)
   vertex_pos(3) = vertex_pos(3) + z_offset 
   if (read_stat/=0) then
      write(0,*) "Error in reading vertex_pos from .ply input file. ",         &
         "The program stops here. "
      stop
   endif
   vertex_ID = i_vert + perimeter_first_vertex_ID - 1
   write(12,'(i12,1x,3(g15.5,1x))',IOSTAT=write_stat) vertex_ID,vertex_pos(:)
   if (write_stat/=0) then
      write(0,*) "Error in writing vertex_pos on perimeter_vertices.txt. ",    &
         "The program stops here. "
      stop
   endif   
enddo
do i_face=1,n_faces
   read(11,'(i1)',IOSTAT=read_stat,ADVANCE="NO") face_sides
   if (read_stat/=0) then
      write(0,*) "Error in reading the face sides from .ply input file. ",     &
         "The program stops here. "
      stop
   endif      
   read(11,*,IOSTAT=read_stat) vertex_IDs_of_ply_face(1:face_sides)
   if (read_stat/=0) then
      write(0,*) "Error in reading face vertices from .ply input file. ",      &
         "The program stops here. "
      stop
   endif
   vertex_IDs_of_ply_face(1:face_sides) = vertex_IDs_of_ply_face(1:face_sides) &
                                          + perimeter_first_vertex_ID 
   face_ID = i_face + perimeter_first_face_ID - 1
   if (face_sides<6) vertex_IDs_of_ply_face((face_sides+1):6) = -1. 
   write(13,'(8(i8,1x))',IOSTAT=write_stat) face_ID,                           &
      vertex_IDs_of_ply_face(:),perimeter_ID
   if (write_stat/=0) then
      write(0,*) "Error in writing faces on perimeter_faces.txt. ",            &
         "The program stops here. "
      stop
   endif  
enddo
close(11)
close(12)
close(13)
!------------------------
! Deallocations
!------------------------
return
end program ply2SPHERA_perimeter

