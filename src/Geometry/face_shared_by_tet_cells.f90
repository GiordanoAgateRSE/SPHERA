!-------------------------------------------------------------------------------
! SPHERA v.10.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.10.0.0
! SPHERA v.10.0.0 is free software: you can redistribute it and/or modify
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
! Program unit: face_shared_by_tet_cells
! Description: Check if 2 tetrahedrcal cells have 1 face in common. 
!              “2 cells with 3 points in common” is equivalent to “2 cells with 
!              1 face in common”.
!-------------------------------------------------------------------------------
#if (defined SPACE_3D) && (defined SOLID_BODIES)
subroutine face_shared_by_tet_cells(cell_1,cell_2,test,shared_face_cell_1,     &
   shared_face_cell_2)
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
integer(4),dimension(4),intent(in) :: cell_1,cell_2
logical,intent(inout) :: test
! ID (1,2,3 or 4) of the face shared for both cell_1 and cell_2 
integer(4),intent(out) :: shared_face_cell_1,shared_face_cell_2
! Flags of the shared points defining the faces (point order according to 
! tetmesh ".vtu" convention): ".true." means shared point
logical,dimension(4) :: shared_points_cell_1,shared_points_cell_2
integer(4) :: shared_points_count,ii,i_point,j_point
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
shared_points_cell_1(1:4) = .false.
shared_points_cell_2(1:4) = .false.
shared_points_count = 0
test = .false.
!------------------------
! Statements
!------------------------
! Loops over the points of the two cells to cover all the point-point 
! comparisons
do i_point=1,4
   do j_point=1,4
      if (cell_1(i_point)==cell_2(j_point)) then 
! Shared point update in the definition of the faces
         shared_points_cell_1(i_point) = .true.
         shared_points_cell_2(j_point) = .true.
! Counter update
         shared_points_count = shared_points_count + 1
      endif
   enddo
enddo
if (shared_points_count==3) then
! Positive test: 1 face is shared
   test = .true.
! The position of the non-shared point in the point order of the face defines 
! the ID of the shared face
   do ii=1,4
! Shared face of cell_1
      if (.not.shared_points_cell_1(ii)) then
         select case (ii)
            case(1)
               shared_face_cell_1 = 4
            case(2)
               shared_face_cell_1 = 3
            case(3)
               shared_face_cell_1 = 2
            case(4)
               shared_face_cell_1 = 1
         endselect
      endif
! Shared face of cell_2
      if (.not.shared_points_cell_2(ii)) then
         select case (ii)
         case(1)
            shared_face_cell_2 = 4
         case(2)
            shared_face_cell_2 = 3
         case(3)
            shared_face_cell_2 = 2
         case(4)
            shared_face_cell_2 = 1
         endselect
      endif
   enddo
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine face_shared_by_tet_cells
#endif
