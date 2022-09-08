!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: IC_CAE_body_particles
! Description: Initial Conditions for the CAE-made solid body particles      
!-------------------------------------------------------------------------------
#if (defined SPACE_3D) && (defined SOLID_BODIES)
subroutine IC_CAE_body_particles
!------------------------
! Modules
!------------------------
use Hybrid_allocation_module
use Static_allocation_module
use Dynamic_allocation_module
use Memory_I_O_interface_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
logical :: test
integer(4) :: i_vtu_grid,i_vert,i_vert2,i_bp,ib,i_vtu_cell,aux_int,aux_int_2
integer(4) :: ib_aux,j_vtu_cell,shared_face_cell_1,shared_face_cell_2,i_face
double precision :: aux_scal,face_area
double precision,dimension(3) :: P1,P2,P3,P4,vec_P1_P2,vec_P1_P3
double precision,dimension(3) :: vec_P1_P4,aux_vec,face_normal
integer(4),dimension(4) :: cell_1,cell_2
character(100) :: array_name
logical,dimension(:,:),allocatable :: test_surface_faces_1,test_surface_faces_2
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine area_triangle(P1,P2,P3,area,normal)
      implicit none
      double precision,intent(in) :: P1(3),P2(3),P3(3)
      double precision,intent(out) :: area
      double precision,intent(out) :: normal(3)
   end subroutine area_triangle
   subroutine Vector_Product(uu,VV,ww,SPACEDIM)
      implicit none
      integer(4),intent(in) :: SPACEDIM
      double precision,intent(in),dimension(SPACEDIM) :: uu,VV
      double precision,intent(inout),dimension(SPACEDIM) :: ww
   end subroutine Vector_Product
   integer(4) function ParticleCellNumber(coord)
      implicit none
      double precision,intent(in) :: coord(3)
   end function ParticleCellNumber
   subroutine face_shared_by_tet_cells(cell_1,cell_2,test,shared_face_cell_1,  &
      shared_face_cell_2)
      implicit none
      integer(4),dimension(4),intent(in) :: cell_1,cell_2
      logical,intent(inout) :: test
      integer(4),intent(out) :: shared_face_cell_1,shared_face_cell_2
   end subroutine face_shared_by_tet_cells
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
! Loop over the CAE-made solid bodies (i.e., the ".vtu" grids)
do i_vtu_grid=1,n_bodies_CAE
! ID of the CAE-made body in the array of solid bodies
   ib = i_vtu_grid + (n_bodies - n_bodies_CAE)
! Number of body particles until the previous body (for later assessments of 
! the body particle IDs)
   aux_int = 0
   do ib_aux=1,(ib-1)
      aux_int = aux_int + body_arr(ib_aux)%npart
   enddo
! Allocation and initialization of the temporary arrays to test the faces of 
! the ".vtu" cells to define them as surface or internal.
if ((body_arr(ib)%surface_detection==1).or.                                    &
   (body_arr(ib)%surface_detection==3)) then
   array_name = "test_surface_faces_1"
   call allocate_de_log_r2(.true.,test_surface_faces_1,body_arr(ib)%npart,4,   &
      array_name,ulog_flag=.true.)
   test_surface_faces_1(1:body_arr(ib)%npart,1:4) = .true.
endif
if ((body_arr(ib)%surface_detection==2).or.                                    &
   (body_arr(ib)%surface_detection==3)) then
   array_name = "test_surface_faces_2"
   call allocate_de_log_r2(.true.,test_surface_faces_2,body_arr(ib)%npart,4,   &
      array_name,ulog_flag=.true.)
   test_surface_faces_2(1:body_arr(ib)%npart,1:4) = .true.
endif
! First loop over the body particles of the current body; they are bijectively 
! associated with the ".vtu" cells. 
! Particle variables assessed: volume, position, cell, surface flag, 
! velocity.
! Body variables assessed: volume, maximum velocity.
!$omp parallel do default(none)                                                &
!$omp shared(vtu_grids,bp_arr,ib,body_arr,i_vtu_grid,aux_int)                  &
!$omp shared(test_surface_faces_1,test_surface_faces_2)                        &
!$omp private(i_vtu_cell,i_bp,i_vert2,P1,P2,P3,P4,vec_P1_P2)                   &
!$omp private(vec_P1_P3,vec_P1_P4,aux_vec,aux_scal,i_vert,j_vtu_cell,test)     &
!$omp private(shared_face_cell_1,shared_face_cell_2,cell_1,cell_2,aux_int_2)
   do i_vtu_cell=1,body_arr(ib)%npart
! ID of the current body particle
      i_bp = i_vtu_cell + aux_int
! Initializations
      bp_arr(i_bp)%body = ib
! For surface detection method "1", this initialization has no effect, for the 
! methods "2" and "3", it is mandatory
      bp_arr(i_bp)%surface = .true.
      bp_arr(i_bp)%pos(1:3) = 0.d0
      bp_arr(i_bp)%area = 0.d0
      bp_arr(i_bp)%normal(1:3) = 0.d0
      bp_arr(i_bp)%pres = 0.d0
      bp_arr(i_bp)%acc(1:3) = 0.d0
      bp_arr(i_bp)%vel_mir(1:3) = 0.d0
! Volume of the body particle (tetrahedron)
      i_vert2 = vtu_grids(i_vtu_grid)%cells%points_IDs(i_vtu_cell,1)
      P1(1:3) = vtu_grids(i_vtu_grid)%points%vertex(i_vert2)%pos(1:3)
      i_vert2 = vtu_grids(i_vtu_grid)%cells%points_IDs(i_vtu_cell,2)
      P2(1:3) = vtu_grids(i_vtu_grid)%points%vertex(i_vert2)%pos(1:3)
      i_vert2 = vtu_grids(i_vtu_grid)%cells%points_IDs(i_vtu_cell,3)
      P3(1:3) = vtu_grids(i_vtu_grid)%points%vertex(i_vert2)%pos(1:3)
      i_vert2 = vtu_grids(i_vtu_grid)%cells%points_IDs(i_vtu_cell,4)
      P4(1:3) = vtu_grids(i_vtu_grid)%points%vertex(i_vert2)%pos(1:3)
      vec_P1_P2(1:3) = P2(1:3) - P1(1:3)
      vec_P1_P3(1:3) = P3(1:3) - P1(1:3)
      vec_P1_P4(1:3) = P4(1:3) - P1(1:3)
      call Vector_Product(vec_P1_P2,vec_P1_P3,aux_vec,3)
      aux_scal = dot_product(vec_P1_P4,aux_vec)
      bp_arr(i_bp)%volume = dabs(aux_scal) / 6.d0
!$omp critical (omp_IC_CAE_body_particles_volume)
! Update the body volume
      body_arr(ib)%volume = body_arr(ib)%volume + bp_arr(i_bp)%volume
!$omp end critical (omp_IC_CAE_body_particles_volume)
      aux_int_2 = 0
! Loop over the vertex relative IDs "i_vert" in the local list of the current 
! ".vtu" cell
      do i_vert=1,4
! "i_vert2" is the global vertex ID in the list of ".vtu" points
         i_vert2 = vtu_grids(i_vtu_grid)%cells%points_IDs(i_vtu_cell,i_vert)
! Particle position: cell volume barycentre (sums)
         bp_arr(i_bp)%pos(1:3) = bp_arr(i_bp)%pos(1:3) +                       &
            vtu_grids(i_vtu_grid)%points%vertex(i_vert2)%pos(1:3)
! Surface detection
         if ((body_arr(ib)%surface_detection==2).or.                           &
             (body_arr(ib)%surface_detection==3)) then
! Possibly invalidate the surface faces (of the cell) associated with an inner 
! cell point
            if ((bp_arr(i_bp)%surface).and.                                    &
               (.not.vtu_grids(i_vtu_grid)%points%surface(i_vert2))) then
! Body particles with no surface faces are inner body particles. Body particles 
! with at least a surface face are surface body particles.            
! Any body particle with at least two inner points is an inner body particle.
               aux_int_2 = aux_int_2 + 1
               if (aux_int_2>1) bp_arr(i_bp)%surface = .false.
               select case (i_vert)
! Invalidate the cell faces associated with the inner point as surface faces 
! (from the point to the faces)
                  case(1)
                     test_surface_faces_2(i_vtu_cell,1) = .false.
                     test_surface_faces_2(i_vtu_cell,2) = .false.
                     test_surface_faces_2(i_vtu_cell,3) = .false.
                  case(2)
                     test_surface_faces_2(i_vtu_cell,1) = .false.
                     test_surface_faces_2(i_vtu_cell,2) = .false.
                     test_surface_faces_2(i_vtu_cell,4) = .false.
                  case(3)
                     test_surface_faces_2(i_vtu_cell,1) = .false.
                     test_surface_faces_2(i_vtu_cell,3) = .false.
                     test_surface_faces_2(i_vtu_cell,4) = .false.
                  case(4)
                     test_surface_faces_2(i_vtu_cell,2) = .false.
                     test_surface_faces_2(i_vtu_cell,3) = .false.
                     test_surface_faces_2(i_vtu_cell,4) = .false.
                  case default
! This case cannot occur
               endselect
            endif
         endif
      enddo
! Surface detection
      if ((body_arr(ib)%surface_detection==1).or.                              &
         (body_arr(ib)%surface_detection==3)) then
! Tetmesh cell. 1 face shared: internal face; 0 faces shared: surface face; 
! >1 face shared: error.
! Loop over the foolowing ".vtu" cells for point-to-point (face-to-face) 
! comparisons
         do j_vtu_cell=i_vtu_cell+1,body_arr(ib)%npart
! Test on a possible sharing face
            cell_1(1:4) = vtu_grids(i_vtu_grid)%cells%points_IDs(i_vtu_cell,1:4)
            cell_2(1:4) = vtu_grids(i_vtu_grid)%cells%points_IDs(j_vtu_cell,1:4)
            call face_shared_by_tet_cells(cell_1,cell_2,test,                  &
               shared_face_cell_1,shared_face_cell_2)
            if (test) then
! Invalidate surface flag for shared faces
               test_surface_faces_1(i_vtu_cell,shared_face_cell_1) = .false.
               test_surface_faces_1(j_vtu_cell,shared_face_cell_2) = .false.
            endif
         enddo
      endif
! Body particle position (barycentre of the input tetrahedron before possible 
! rotation)
      bp_arr(i_bp)%pos(1:3) = bp_arr(i_bp)%pos(1:3) / 4.d0
#ifdef SPACE_3D
! Translation applied to the body particles of the current body
      bp_arr(i_bp)%pos(1:3) = bp_arr(i_bp)%pos(1:3) +                          &
         body_arr(ib)%vec_bp_CAE_trans(1:3)
#endif 
! Rotation around the centre of rotation provided in input: position
      bp_arr(i_bp)%rel_pos(1:3) = bp_arr(i_bp)%pos(1:3) -                      &
                                  body_arr(ib)%x_rotC(1:3)
      call vector_rotation_Rodrigues(body_arr(ib)%n_R_IO(:),                   &
         body_arr(ib)%teta_R_IO,bp_arr(i_bp)%rel_pos(:))
      bp_arr(i_bp)%pos(1:3) = bp_arr(i_bp)%rel_pos(1:3) +                      &
                              body_arr(ib)%x_rotC(1:3)
! Relative position
      bp_arr(i_bp)%rel_pos(1:3) = bp_arr(i_bp)%pos(1:3) - body_arr(ib)%x_CM(1:3)
! Particle cell
      bp_arr(i_bp)%cell = ParticleCellNumber(bp_arr(i_bp)%pos)
! Velocity
      call Vector_Product(body_arr(ib)%omega,bp_arr(i_bp)%rel_pos,aux_vec,3)
      bp_arr(i_bp)%vel(1:3) = body_arr(ib)%u_CM(1:3) + aux_vec(1:3)
! Update of umax
      aux_scal = dsqrt(dot_product(bp_arr(i_bp)%vel,bp_arr(i_bp)%vel))        
!$omp critical (omp_IC_CAE_body_particles_umax)
      body_arr(ib)%umax = max(body_arr(ib)%umax,aux_scal)
!$omp end critical (omp_IC_CAE_body_particles_umax)
   enddo
!$omp end parallel do
! Second loop over the ".vtu" cells of the current body.
! Generic body particle features: mass.
! Surface body particle features: area; normal (with possible IC rotation).
!$omp parallel do default(none)                                                &
!$omp shared(vtu_grids,i_vtu_grid,bp_arr,body_arr,ib,aux_int)                  &
!$omp shared(test_surface_faces_1,test_surface_faces_2,uerr)                   &
!$omp private(i_vtu_cell,i_bp,i_vert2,P1,P2,P3,face_area,face_normal,i_face)
   do i_vtu_cell=1,body_arr(ib)%npart
      i_bp = i_vtu_cell + aux_int
! Surface assignment
! Re-inizialization from scratch, even if for the "method 2" the values have 
! been already computed. This way is faster for treating any method, than 
! avoiding redundancy for the method "2".
      bp_arr(i_bp)%surface = .false.
! Loop over the cell faces
      do i_face=1,4
         select case (body_arr(ib)%surface_detection)
            case(1)
! Activate surface flag for faces not shared
               vtu_grids(i_vtu_grid)%cells%surface(i_vtu_cell,i_face) =        &
                  test_surface_faces_1(i_vtu_cell,i_face)
            case(2)
               vtu_grids(i_vtu_grid)%cells%surface(i_vtu_cell,i_face) =        &
                  test_surface_faces_2(i_vtu_cell,i_face)
            case(3)
               if ((test_surface_faces_1(i_vtu_cell,i_face)).or.               &
                  (test_surface_faces_2(i_vtu_cell,i_face))) then
                  vtu_grids(i_vtu_grid)%cells%surface(i_vtu_cell,i_face) =     &
                     .true.
               endif
            case default
               write(uerr,*) "Error. The surface detection method for ",       &
                  "CAE-made solid bodies has not benn properly assigned. ",    &
                  "The program stops here."
               stop
         endselect
         if (vtu_grids(i_vtu_grid)%cells%surface(i_vtu_cell,i_face)) then
! Activate surface flag for body particle in case of at least 1 surface face of 
! the associated vtu cell.
            bp_arr(i_bp)%surface = .true.
         endif
      enddo
      if (bp_arr(i_bp)%surface) then
! Surface body particles
! Particle area (internal faces are not considered)
         if (vtu_grids(i_vtu_grid)%cells%surface(i_vtu_cell,1)) then
! The first face of the current cell/particle is a surface face
! Find the vertices of the first surface face (from the face to the points)
            i_vert2 = vtu_grids(i_vtu_grid)%cells%points_IDs(i_vtu_cell,1)
            P1(1:3) = vtu_grids(i_vtu_grid)%points%vertex(i_vert2)%pos(1:3)
            i_vert2 = vtu_grids(i_vtu_grid)%cells%points_IDs(i_vtu_cell,2)
            P2(1:3) = vtu_grids(i_vtu_grid)%points%vertex(i_vert2)%pos(1:3)
            i_vert2 = vtu_grids(i_vtu_grid)%cells%points_IDs(i_vtu_cell,3)
            P3(1:3) = vtu_grids(i_vtu_grid)%points%vertex(i_vert2)%pos(1:3)
! Assess the area of the surface face
            call area_triangle(P1,P2,P3,face_area,face_normal)
! Update the normal of the surface body particle
            bp_arr(i_bp)%normal(1:3) = bp_arr(i_bp)%normal(1:3) + face_area *  &
                                       face_normal(1:3)
         endif
         if (vtu_grids(i_vtu_grid)%cells%surface(i_vtu_cell,2)) then
! The second face of the current cell/particle is a surface face
            i_vert2 = vtu_grids(i_vtu_grid)%cells%points_IDs(i_vtu_cell,4)
            P1(1:3) = vtu_grids(i_vtu_grid)%points%vertex(i_vert2)%pos(1:3)
            i_vert2 = vtu_grids(i_vtu_grid)%cells%points_IDs(i_vtu_cell,2)
            P2(1:3) = vtu_grids(i_vtu_grid)%points%vertex(i_vert2)%pos(1:3)
            i_vert2 = vtu_grids(i_vtu_grid)%cells%points_IDs(i_vtu_cell,1)
            P3(1:3) = vtu_grids(i_vtu_grid)%points%vertex(i_vert2)%pos(1:3)
            call area_triangle(P1,P2,P3,face_area,face_normal)
            bp_arr(i_bp)%normal(:) = bp_arr(i_bp)%normal(1:3) + face_area *    &
                                     face_normal(1:3)
         endif
         if (vtu_grids(i_vtu_grid)%cells%surface(i_vtu_cell,3)) then
! The third face of the current cell/particle is a surface face
            i_vert2 = vtu_grids(i_vtu_grid)%cells%points_IDs(i_vtu_cell,1)
            P1(1:3) = vtu_grids(i_vtu_grid)%points%vertex(i_vert2)%pos(1:3)
            i_vert2 = vtu_grids(i_vtu_grid)%cells%points_IDs(i_vtu_cell,3)
            P2(1:3) = vtu_grids(i_vtu_grid)%points%vertex(i_vert2)%pos(1:3)
            i_vert2 = vtu_grids(i_vtu_grid)%cells%points_IDs(i_vtu_cell,4)
            P3(1:3) = vtu_grids(i_vtu_grid)%points%vertex(i_vert2)%pos(1:3)
            call area_triangle(P1,P2,P3,face_area,face_normal)
            bp_arr(i_bp)%normal(:) = bp_arr(i_bp)%normal(1:3) + face_area *    &
                                     face_normal(1:3)
         endif
         if (vtu_grids(i_vtu_grid)%cells%surface(i_vtu_cell,4)) then
! The fourth face of the current cell/particle is a surface face
            i_vert2 = vtu_grids(i_vtu_grid)%cells%points_IDs(i_vtu_cell,4)
            P1(1:3) = vtu_grids(i_vtu_grid)%points%vertex(i_vert2)%pos(1:3)
            i_vert2 = vtu_grids(i_vtu_grid)%cells%points_IDs(i_vtu_cell,3)
            P2(1:3) = vtu_grids(i_vtu_grid)%points%vertex(i_vert2)%pos(1:3)
            i_vert2 = vtu_grids(i_vtu_grid)%cells%points_IDs(i_vtu_cell,2)
            P3(1:3) = vtu_grids(i_vtu_grid)%points%vertex(i_vert2)%pos(1:3)
            call area_triangle(P1,P2,P3,face_area,face_normal)
            bp_arr(i_bp)%normal(:) = bp_arr(i_bp)%normal(1:3) + face_area *    &
                                     face_normal(1:3)
         endif
! Area of the surface body particle
         bp_arr(i_bp)%area =                                                   &
                  dsqrt(dot_product(bp_arr(i_bp)%normal,bp_arr(i_bp)%normal))
! Normal of the surface body particle
         bp_arr(i_bp)%normal(1:3) = bp_arr(i_bp)%normal(1:3) / bp_arr(i_bp)%area
! Rotation around the centre of rotation provided in input: normal
         call vector_rotation_Rodrigues(body_arr(ib)%n_R_IO(:),                &
            body_arr(ib)%teta_R_IO,bp_arr(i_bp)%normal(:))
      endif
!$omp critical (omp_IC_CAE_body_particles_area_int_ndA)
! Update the body area
      body_arr(ib)%area = body_arr(ib)%area + bp_arr(i_bp)%area
! Update the integral of the normal over the external body surface
      body_arr(ib)%int_ndA(1:3) = body_arr(ib)%int_ndA(1:3) +                  &
                                  bp_arr(i_bp)%area * bp_arr(i_bp)%normal(1:3)
!$omp end critical (omp_IC_CAE_body_particles_area_int_ndA)
! Mass of the body particle
      bp_arr(i_bp)%mass = bp_arr(i_bp)%volume * body_arr(ib)%mass /            &
                          body_arr(ib)%volume
   enddo
!$omp end parallel do
enddo
!------------------------
! Deallocations
!------------------------
! Deallocation of the temporary arrays to test the faces of the 
! ".vtu" cells to define them as surface or internal.
if ((body_arr(ib)%surface_detection==1).or.                                    &
   (body_arr(ib)%surface_detection==3)) then
   array_name = "test_surface_faces_1"
   call allocate_de_log_r2(.false.,test_surface_faces_1,array_name=array_name, &
      ulog_flag=.true.)
endif
if ((body_arr(ib)%surface_detection==2).or.                                    &
   (body_arr(ib)%surface_detection==3)) then
   array_name = "test_surface_faces_2"
   call allocate_de_log_r2(.false.,test_surface_faces_2,array_name=array_name, &
      ulog_flag=.true.)
endif
return
end subroutine IC_CAE_body_particles
#endif
