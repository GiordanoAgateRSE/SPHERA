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
! Program unit: neighbouring_hcell_CLCp                               
! Description: Contribution of the current polygon to the neighbouring search 
!              linking the cells of the horizontal background grid with the 
!              neighbouring CLC polygons
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine neighbouring_hcell_CLCp(i_pol)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_file_module
use Neighbouring_Search_interface_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: i_pol
integer(4) :: i_vert,i_hcel,aux_int,aux_int_2,ix,iy
double precision :: x_min,x_max,y_min,y_max,z_1
double precision :: pos(3)
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
! Identify the lower-left and the upper-right cell IDs of the rectangle 
! circumscribing the current CLC polygon in the horizontal projection of the 
! background grid
x_min = max_positive_number
x_max = max_negative_number
y_min = max_positive_number
y_max = max_negative_number
do i_vert=1,CLC%polygons(i_pol)%n_vertices
   x_min=min(x_min,CLC%polygons(i_pol)%vertices(i_vert,1))
   x_max=max(x_max,CLC%polygons(i_pol)%vertices(i_vert,1))
   y_min=min(y_min,CLC%polygons(i_pol)%vertices(i_vert,2))
   y_max=max(y_max,CLC%polygons(i_pol)%vertices(i_vert,2))
enddo
! Horizontal cell indices of the lower-left corner of the circumscribing 
! rectangle
z_1 = Grid%extr(3,1) + Grid%dcd(3) / 2.d0
pos(1) = x_min
pos(2) = y_min
pos(3) = z_1
i_hcel = ParticleCellNumber(pos)
! Horizontal cell indices of the lower-left corner of the circumscribing 
! rectangle
aux_int = CellIndices(i_hcel,CLC%polygons(i_pol)%ix_cel_ll,                    &
          CLC%polygons(i_pol)%iy_cel_ll,aux_int_2)
! Horizontal cell indices of the upper-right corner of the circumscribing 
! rectangle
pos(1) = x_max
pos(2) = y_max
pos(3) = z_1
i_hcel = ParticleCellNumber(pos)
! Horizontal cell indices of the upper-right corner of the circumscribing 
! rectangle
aux_int = CellIndices(i_hcel,CLC%polygons(i_pol)%ix_cel_ur,                    &
          CLC%polygons(i_pol)%iy_cel_ur,aux_int_2)
! In case a cell index lies outside the numerical domain, either it is replaced 
! by the associated edge domain index or the rectangle is discarded as it 
! completely lies outside the domain.
if (CLC%polygons(i_pol)%ix_cel_ll<1) then
   CLC%polygons(i_pol)%ix_cel_ll = 1
   else
      if (CLC%polygons(i_pol)%ix_cel_ll>Grid%ncd(1)) return
endif
if (CLC%polygons(i_pol)%iy_cel_ll<1) then
   CLC%polygons(i_pol)%iy_cel_ll = 1
   else
      if (CLC%polygons(i_pol)%iy_cel_ll>Grid%ncd(2)) return
endif
if (CLC%polygons(i_pol)%ix_cel_ur<1) then
   return
   else
      if (CLC%polygons(i_pol)%ix_cel_ur>Grid%ncd(1)) then
         CLC%polygons(i_pol)%ix_cel_ur = Grid%ncd(1)
      endif
endif
if (CLC%polygons(i_pol)%iy_cel_ur<1) then
   return
   else
      if (CLC%polygons(i_pol)%iy_cel_ur>Grid%ncd(2)) then
         CLC%polygons(i_pol)%iy_cel_ur = Grid%ncd(2)
      endif
endif
! Assign the CLC polygon ID to the neighbouring CLC-polygon list of all the 
! background-grid cells occupied by the circumscribing rectangle
do iy=CLC%polygons(i_pol)%iy_cel_ll,CLC%polygons(i_pol)%iy_cel_ur
   do ix=CLC%polygons(i_pol)%ix_cel_ll,CLC%polygons(i_pol)%ix_cel_ur
      i_hcel = CellNumber(ix,iy,1)
!$omp critical (omp_CLC_neighbouring)
      n_neigh_hcell_CLCpol(i_hcel) = n_neigh_hcell_CLCpol(i_hcel) + 1
      aux_int = (i_hcel - 1) * NMAXPARTJ + n_neigh_hcell_CLCpol(i_hcel)
      neigh_hcell_CLCpol(aux_int) = i_pol
!$omp end critical (omp_CLC_neighbouring)
   enddo
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine neighbouring_hcell_CLCp
#endif
