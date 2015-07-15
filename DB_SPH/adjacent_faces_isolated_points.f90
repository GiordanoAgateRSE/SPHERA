!AA601 the whole subroutine
!cfile BC_wall_elements.f90
!************************************************************************************
!                          S P H E R A 5.0.3 
!
!                Smoothed Particle Hydrodinamics Code
!
!************************************************************************************
!
! File name     : adjacent_faces_isolated_points
! Creation      : Amicarelli A., 26Jan15
!
!************************************************************************************
! Module purpose : Provided 2 adjacent triangular/quadrilateral faces, find at least 2 vertices not in common, at least one per face.
!                  They are ID_face1_iso and ID_face2_iso. In case the faces are not adjacent, then false_hyp=.true.
!
! Calling routine: semi_particle_volumes
!
! Called routines: /
!
!************************************************************************************

subroutine adjacent_faces_isolated_points(face1,face2,ID_face1_iso,ID_face2_iso,false_hyp)

! Modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE

! Declarations
implicit none
double precision, dimension(4,3), intent(in) :: face1,face2
integer(4),intent(out) :: ID_face1_iso,ID_face2_iso
logical,intent(out) :: false_hyp 
integer(4) :: test_face1,test_face2,i,j,n_vert

! Interface blocks

! Allocations

! Initializations
test_face1 = 0
test_face2 = 0
false_hyp = .false.

! Statements
if (ncord==3) then
   n_vert = 3
   else
      n_vert = 4
endif
! do over the 3/4 vertices of the first face
do_vertices_face1: do i=1,n_vert
! do over the 3/4 vertices of the second face   
    do j=1,n_vert
       if ( (face1(i,1)==face2(j,1)) .and. (face1(i,2)==face2(j,2)) .and. (face1(i,3)==face2(j,3)) ) then
! In case the vertex is in common, the vertex of the first face is not anymore isolated and update test_face1/2 (the sum 
! of the 2 non-isolated vertex IDs of face 1/2; in 2D the ID squares are considered to avoid ambiguities in the following computations)
       if (ncord==3) then
          test_face1 = test_face1 + i 
          test_face2 = test_face2 + j 
          else
             test_face1 = test_face1 + i**2 
             test_face2 = test_face2 + j**2    
       endif
          cycle do_vertices_face1
       endif
   enddo
! end do over the vertices of the second face      
end do do_vertices_face1
! end do over the vertices of the first face
! The only/first vertex (in 3D/2D) of the 1st face, not contributing to test_face1, is finally found
if (ncord==3) then
   select case(test_face1)
      case(5)
      ID_face1_iso = 1
      case(4)
      ID_face1_iso = 2
      case(3)
      ID_face1_iso = 3
      case default
      ID_face1_iso = 0
   end select
   else
      select case(test_face1)
      case(13,20,25)
      ID_face1_iso = 1
      case(10,17)
      ID_face1_iso = 2
      case(5)
      ID_face1_iso = 3
      case default
      ID_face1_iso = 0
      end select   
endif
! The only/first vertex (in 3D/2D) of the 2nd face, not contributing to test_face2, is finally found
if (ncord==3) then
   select case(test_face2)
      case(5)
      ID_face2_iso = 1
      case(4)
      ID_face2_iso = 2
      case(3)
      ID_face2_iso = 3
      case default
      ID_face2_iso = 0
   end select  
   else
      select case(test_face2)
         case(13,20,25)
         ID_face2_iso = 1
         case(10,17)
         ID_face2_iso = 2
         case(5)
         ID_face2_iso = 3
         case default
         ID_face2_iso = 0
      end select 
endif
if ((ID_face1_iso==0).or.(ID_face2_iso==0)) false_hyp = .true.

! Deallocations

return
end subroutine adjacent_faces_isolated_points
!---split

