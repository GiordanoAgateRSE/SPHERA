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
! Program unit: adjacent_faces_isolated_points
! Description: Provided 2 adjacent triangular/quadrilateral faces, it finds at 
!              least 2 vertices not in common, at least one per face. They are 
!              ID_face1_iso and ID_face2_iso. In case the faces are not 
!              adjacent, then false_hyp=.true.         
!-------------------------------------------------------------------------------
subroutine adjacent_faces_isolated_points(face1,face2,ID_face1_iso,            &
                                          ID_face2_iso,false_hyp)
!------------------------
! Modules
!------------------------
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
double precision, dimension(4,3), intent(in) :: face1,face2
logical,intent(out) :: false_hyp
integer(4),intent(out) :: ID_face1_iso,ID_face2_iso
integer(4) :: test_face1,test_face2,i,j,n_vert
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
test_face1 = 0
test_face2 = 0
false_hyp = .false.
#ifdef SPACE_3D
   n_vert = 3
#elif defined SPACE_2D
      n_vert = 4
#endif
!------------------------
! Statements
!------------------------
! do over the 3/4 vertices of the first face
do_vertices_face1: do i=1,n_vert
! do over the 3/4 vertices of the second face   
   do j=1,n_vert
      if ((face1(i,1)==face2(j,1)).and.(face1(i,2)==face2(j,2)).and.           &
         (face1(i,3)==face2(j,3))) then
! In case the vertex is in common, the vertex of the first face is not anymore 
! isolated and it updates test_face1/2 (the sum 
! of the 2 non-isolated vertex IDs of face 1/2; in 2D the ID squares are 
! considered to avoid ambiguities in the following computations)
#ifdef SPACE_3D
            test_face1 = test_face1 + i 
            test_face2 = test_face2 + j 
#elif defined SPACE_2D
               test_face1 = test_face1 + i ** 2 
               test_face2 = test_face2 + j ** 2
#endif
         cycle do_vertices_face1
      endif
   enddo
enddo do_vertices_face1
! The only/first vertex (in 3D/2D) of the 1st face, not contributing to 
! test_face1, is finally found
#ifdef SPACE_3D
   select case (test_face1)
      case(5)
         ID_face1_iso = 1
      case(4)
         ID_face1_iso = 2
      case(3)
         ID_face1_iso = 3
      case default
         ID_face1_iso = 0
   endselect
#elif defined SPACE_2D
      select case (test_face1)
         case(13,20,25)
            ID_face1_iso = 1
         case(10,17)
            ID_face1_iso = 2
         case(5)
            ID_face1_iso = 3
         case default
            ID_face1_iso = 0
      endselect   
#endif
! The only/first vertex (in 3D/2D) of the 2nd face, not contributing to   
! test_face2, is finally found
#ifdef SPACE_3D
   select case (test_face2)
      case(5)
         ID_face2_iso = 1
      case(4)
         ID_face2_iso = 2
      case(3)
         ID_face2_iso = 3
      case default
         ID_face2_iso = 0
   endselect  
#elif defined SPACE_2D
      select case (test_face2)
         case(13,20,25)
            ID_face2_iso = 1
         case(10,17)
            ID_face2_iso = 2
         case(5)
            ID_face2_iso = 3
         case default
            ID_face2_iso = 0
      endselect 
#endif
if ((ID_face1_iso==0).or.(ID_face2_iso==0)) false_hyp = .true.
!------------------------
! Deallocations
!------------------------
return
end subroutine adjacent_faces_isolated_points
