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
! Program unit: Neighbouring_Search_interface_module            
! Description: Interfaces to the program units of the folder 
!              Neighbouring_Search
!-------------------------------------------------------------------------------
module Neighbouring_Search_interface_module
interface
#ifdef SPACE_3D
   subroutine cell_indices2pos(ix,iy,iz,pos)
      implicit none
      integer(4),intent(in) :: ix,iy,iz
      double precision,intent(out) :: pos(3)
   end subroutine cell_indices2pos
#endif
   integer(4) function CellIndices(nc,i,j,k)
      implicit none
      integer(4),intent(in) :: nc
      integer(4),intent(out) :: i,j,k
   end function CellIndices
   integer(4) function CellNumber(i,j,k)
      implicit none
      integer(4),intent(in) :: i,j,k
   end function CellNumber
   subroutine InterFix(npi,appo,unity)
      implicit none
      integer(4),intent(in) :: npi
      double precision,intent(inout) :: unity
      double precision,intent(inout),dimension(3) :: appo
   end subroutine InterFix
   integer(4) function ParticleCellNumber(coord)
      implicit none
      double precision,intent(in) :: coord(3)
   end function ParticleCellNumber
end interface 
end module Neighbouring_Search_interface_module
