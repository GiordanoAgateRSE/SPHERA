!AA601 the whole subroutine
!cfile BC_wall_elements.f90
!************************************************************************************
!                          S P H E R A 5.0.3 
!
!                Smoothed Particle Hydrodinamics Code
!
!************************************************************************************
!
! Subroutine : DBSPH_find_close_faces
!
! Creation   : Amicarelli A., 26Jan15; Main features (see Purpose)
!
!************************************************************************************
! Purpose            : Find the adjacent surface elements of a generic 
!                      computational surface element, both in 3D -triangles- and 2D 
!                      -quadrilaterals- (DBSPH)
! Calling procedures : gest_input
! Called procedures  : /
!************************************************************************************

subroutine DBSPH_find_close_faces

! Modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE

! Declarations
implicit none
integer(4) :: npi,irestocell,i_cell_comp,j_cell_comp,k_cell_comp,j_cell_min,j_cell_max,ID_cell,i_cell,j_cell,k_cell,npj,i_vert,j_vert,aux_adjacent_faces,ww,n_vert
integer(4) :: aux_common_vertices
integer(4),external :: CellIndices,CellNumber

! Interface blocks

! Allocations

! Initializations

! Statements
!Loop over the DBSPH surface elements
!$omp parallel do default(none) &
!$omp shared(DBSPH,pg_w,Icont_w,ncord,NPartOrd_w) &
!$omp private(npi,irestocell,i_cell_comp,j_cell_comp,k_cell_comp,j_cell_min,j_cell_max,ID_cell,i_cell,j_cell,k_cell,npj,i_vert,j_vert,aux_adjacent_faces,ww,aux_common_vertices,n_vert)
loop_surface_elements: do npi = 1,DBSPH%n_w
   aux_adjacent_faces = 0
!Loop over the adjacent cells (background grid) 
   if (pg_w(npi)%cella == 0) cycle
   irestocell = CellIndices (pg_w(npi)%cella,i_cell_comp,j_cell_comp,k_cell_comp)
! Adjacent cell indices
   j_cell_min = j_cell_comp - (ncord-2)
   j_cell_max = j_cell_comp + (ncord-2)
   do j_cell = j_cell_min,j_cell_max    
      do i_cell = i_cell_comp-1,i_cell_comp+1    
         do k_cell = k_cell_comp-1,k_cell_comp+1  
            ID_cell = CellNumber(i_cell,j_cell,k_cell)
            if (ID_cell==0) cycle    
            if (Icont_w(ID_cell+1) <= Icont_w(ID_cell)) cycle
! Loop over the neighbouring wall particles in the cell
            loop_neighbour_surface_elements: do ww = Icont_w(ID_cell),Icont_w(ID_cell+1)-1
               npj = NPartOrd_w(ww)
               if (npi==npj) cycle
! Avoid considering inlet and outlet DBSPH elements    
               if (npj>DBSPH%n_w) cycle
               aux_common_vertices = 0
               if (ncord==3) then
                  n_vert = 3
                  else
                     n_vert = 4
               endif      
! Loop over the vertices of the computational element  
               do i_vert=1,n_vert
! Loop over the vertices of the neighbouring element
                  do j_vert=1,n_vert
                     if ( (DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(npi)%vert_list(i_vert))%pos(1) ==      &
                           DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(npj)%vert_list(j_vert))%pos(1) ) .and. & 
                          (DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(npi)%vert_list(i_vert))%pos(2) == &
                           DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(npj)%vert_list(j_vert))%pos(2) ) .and. & 
                          (DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(npi)%vert_list(i_vert))%pos(3) == &
                           DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(npj)%vert_list(j_vert))%pos(3) )  ) then
                        aux_common_vertices = aux_common_vertices + 1
                        if (aux_common_vertices==2) then
                           pg_w(npi)%adjacent_faces(aux_adjacent_faces+1) = npj 
                           aux_adjacent_faces = aux_adjacent_faces + 1
                           if (aux_adjacent_faces==ncord) cycle loop_surface_elements
                           cycle loop_neighbour_surface_elements
                        endif
                     endif
                  enddo
               enddo
            end do loop_neighbour_surface_elements
         end do    
      end do 
   end do 
end do loop_surface_elements
!$omp end parallel do

! Deallocations

return
end subroutine DBSPH_find_close_faces
!end of the subroutine

