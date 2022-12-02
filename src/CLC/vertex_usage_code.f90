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
! Program unit: vertex_usage_code
! Description: Computation of the "vertex usage code" of a triangular face, 
!              provided the occurrences of its 3 vertices in the triangulated 
!              polygon the current face belongs to.
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine vertex_usage_code(V_0_1,V_0_2,V_0_3,F_VUC)
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
! Vertex occurrences in the triangulated polygon the current face belongs to
integer(4),intent(in) :: V_0_1,V_0_2,V_0_3
! Vertex usage code of the face
integer(4),intent(out) :: F_VUC
! Combination of the possible vertex usage codes
integer(4),dimension(6) :: VUC
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
VUC(1) = V_0_1 * 1000000 + V_0_2 * 1000 + V_0_3
VUC(2) = V_0_1 * 1000000 + V_0_3 * 1000 + V_0_2
VUC(3) = V_0_2 * 1000000 + V_0_1 * 1000 + V_0_3
VUC(4) = V_0_2 * 1000000 + V_0_3 * 1000 + V_0_1
VUC(5) = V_0_3 * 1000000 + V_0_1 * 1000 + V_0_2
VUC(6) = V_0_3 * 1000000 + V_0_2 * 1000 + V_0_1
F_VUC = minval(VUC)
!------------------------
! Deallocations
!------------------------
return
end subroutine vertex_usage_code
#endif
