!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: DBSPH_find_close_faces
! Description: Finding the adjacent surface elements of a given surface element,
!              both using 3D -triangular elements- and 2D -quadrilateral raw
!              elements- configurations (DB-SPH).          
!-------------------------------------------------------------------------------
subroutine DBSPH_find_close_faces
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: npi,irestocell,i_cell_comp,j_cell_comp,k_cell_comp,j_cell_min   
integer(4) :: j_cell_max,ID_cell,i_cell,j_cell,k_cell,npj,i_vert,j_vert       
integer(4) :: aux_adjacent_faces,ww,n_vert,aux_common_vertices
integer(4),external :: CellIndices,CellNumber
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
! Loop over the DBSPH surface elements
!$omp parallel do default(none) shared(DBSPH,pg_w,Icont_w,ncord,NPartOrd_w)    &
!$omp private(npi,irestocell,i_cell_comp,j_cell_comp,k_cell_comp,j_cell_min)   &
!$omp private(j_cell_max,ID_cell,i_cell,j_cell,k_cell,npj,i_vert,j_vert)       &
!$omp private(aux_adjacent_faces,ww,aux_common_vertices,n_vert)
loop_surface_elements: do npi = 1,DBSPH%n_w
   aux_adjacent_faces = 0
! Loop over the adjacent cells (background grid) 
   if (pg_w(npi)%cella==0) cycle
   irestocell = CellIndices (pg_w(npi)%cella,i_cell_comp,j_cell_comp,          &
                             k_cell_comp)
! Adjacent cell indices
   j_cell_min = j_cell_comp - (ncord - 2)
   j_cell_max = j_cell_comp + (ncord - 2)
   do j_cell=j_cell_min,j_cell_max    
      do i_cell=i_cell_comp-1,i_cell_comp+1    
         do k_cell=k_cell_comp-1,k_cell_comp+1  
            ID_cell = CellNumber(i_cell,j_cell,k_cell)
            if (ID_cell==0) cycle    
            if (Icont_w(ID_cell+1)<=Icont_w(ID_cell)) cycle
! Loop over the neighbouring wall particles in the cell
            loop_neighbour_surface_elements: do ww=Icont_w(ID_cell),           &
               (Icont_w(ID_cell+1)-1)
               npj = NPartOrd_w(ww)
               if (npi==npj) cycle
! Avoid considering inlet and outlet DB-SPH elements    
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
if ((DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(npi)%vert_list(i_vert))%pos(1)==&
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(npj)%vert_list(j_vert))%pos(1)).and. & 
(DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(npi)%vert_list(i_vert))%pos(2)==    &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(npj)%vert_list(j_vert))%pos(2)).and. & 
(DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(npi)%vert_list(i_vert))%pos(3) ==   &
DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(npj)%vert_list(j_vert))%pos(3)) ) then
                        aux_common_vertices = aux_common_vertices + 1
                        if (aux_common_vertices==2) then
                           pg_w(npi)%adjacent_faces(aux_adjacent_faces+1) = npj 
                           aux_adjacent_faces = aux_adjacent_faces + 1
                           if (aux_adjacent_faces==ncord) cycle                &
                                                          loop_surface_elements
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
!------------------------
! Deallocations
!------------------------
return
end subroutine DBSPH_find_close_faces

