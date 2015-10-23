!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-) 
!      
!     
!   
!      
!  

! This file is part of SPHERA.
!  
!  
!  
!  
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
!  
!  
!  
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: ply2SPHERA_perimeter
! Description: Format conversion from .ply to the format of VERTICES and FACES of SPHERA main input file.        
!----------------------------------------------------------------------------------------------------------------------------------
program ply2SPHERA_perimeter
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: open_stat,read_stat,close_stat,perimeter_first_vertex_ID
integer(4) :: perimeter_first_face_ID,perimeter_ID,vertex_ID
integer(4) :: perimeter_first_vertex_ID,face_ID,perimeter_first_face_ID
character(100) :: ply_file_name
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
   write(*,*) "Error in opening ply2SPHERA_perimeter.inp. ",                   &
      "The program stops here. "
   stop
endif
read(14,*,IOSTAT=read_stat) perimeter_first_vertex_ID,perimeter_first_face_ID, &
   perimeter_ID
if (read_stat/=0) then
   write(*,*) "Error in reading the first line of ply2SPHERA_perimeter.inp. ", &
      "The program stops here. "
   stop
endif   
read(14,*,IOSTAT=read_stat) ply_file_name
if (read_stat/=0) then
   write(*,*) "Error in reading the second line of ply2SPHERA_perimeter.inp. ",&
      "The program stops here. "
   stop
endif
close(14,IOSTAT=close_stat)
if (close_stat/=0) then
   write(*,*) "Error in closing ply2SPHERA_perimeter.inp. ",                   &
      "The program stops here. "
   stop
endif
vertex_ID = perimeter_first_vertex_ID
face_ID = perimeter_first_face_ID
! To convert data format from .ply to .txt files for SPHERA VERTICES and FACES. 
open(11,file=trim(plyfile),IOSTAT=open_stat)
if (open_stat/=0) then
   write(*,*) "Error in opening .ply input file. The program stops here. "
   stop
endif
open(12,file="perimeter_vertices.txt",IOSTAT=open_stat)
if (open_stat/=0) then
   write(*,*) "Error in opening perimeter_vertices.txt. The program stops here."
   stop
endif
open(13,file="perimeter_faces.txt",IOSTAT=open_stat)
if (open_stat/=0) then
   write(*,*) "Error in opening perimeter_faces.txt. The program stops here. "
   stop
endif
read(11,'(/,/,/)',IOSTAT=read_stat)
if (read_stat/=0) then
   write(*,*) "Error in reading the lines 1-4 of .ply input file. ",           &
      "The program stops here. "
   stop
endif
read(11,'(2a,i)',IOSTAT=read_stat) char_aux,n_vertices
if (read_stat/=0) then
   write(*,*) "Error in reading line 5 of .ply input file. ",                  &
      "The program stops here. "
   stop
endif
read(11,'(/,/)',IOSTAT=read_stat)
if (read_stat/=0) then
   write(*,*) "Error in reading lines 6-8 of .ply input file. ",               &
      "The program stops here. "
   stop
endif
read(11,'(2a,i)',IOSTAT=read_stat) char_aux,n_faces
if (read_stat/=0) then
   write(*,*) "Error in reading line 9 of .ply input file. ",                  &
      "The program stops here. "
   stop
endif
read(11,'(/)',IOSTAT=read_stat)
if (read_stat/=0) then
   write(*,*) "Error in reading lines 10-11 of .ply input file. ",              &
      "The program stops here. "
   stop
endif
do i_vert=1,n_vertices
   read(11,*,!AA!!!) vertex_pos(:)
   vertex_ID = vertex_ID + 1
   write(12,'(i12,3(g12.5))') vertex_ID,vertex_pos(:)
enddo
do i_face=1,n_faces
   read(11,*) face_sides,vertex_IDs_of_ply_face(:)
   face_ID = face_ID + 1
   if (face_sides<6) vertex_IDs_of_ply_face(face_sides+1:6) = -1.d0 
   write(13,'(i12,6(g12.5),i12)') face_ID,vertex_IDs_of_ply_face(:),perimeter_ID
enddo
close(11)
close(12)
close(13)
!------------------------
! Deallocations
!------------------------
return
end subroutine ply2SPHERA_perimeter

