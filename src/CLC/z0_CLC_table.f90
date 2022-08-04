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
! Program unit: z0_CLC_table
! Description: Assignment of the z0 value of the current horizontal background 
!              cell as function of the CLC class 
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine z0_CLC_table(ix,iy)
!------------------------
! Modules
!------------------------
use Hybrid_allocation_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: ix,iy
integer(4) :: CLC_class_available
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
! Assignment of z0 as function of the CLC class
CLC_class_available = int(CLC%class_2D(ix,iy) / 1000) * 1000
select case (CLC_class_available)
   case(11000)
      CLC%z0(ix,iy) = 1.0d0
   case(12000)
      CLC%z0(ix,iy) = 1.0d0
   case(13000)
      CLC%z0(ix,iy) = 0.15d0
   case(14000)
      CLC%z0(ix,iy) = 0.3d0
   case(21000)
      CLC%z0(ix,iy) = 0.03d0
   case(22000)
      CLC%z0(ix,iy) = 0.3d0
   case(23000)
      CLC%z0(ix,iy) = 0.03d0
   case(24000)
      CLC%z0(ix,iy) = 0.2d0
   case(31000)
      CLC%z0(ix,iy) = 1.4d0
   case(32000)
      CLC%z0(ix,iy) = 0.2d0
   case(33000)
      CLC%z0(ix,iy) = 0.02d0
   case(41000)
      CLC%z0(ix,iy) = 0.005d0
   case(42000)
      CLC%z0(ix,iy) = 0.005d0
   case(51000)
      CLC%z0(ix,iy) = 0.0002d0
   case(52000)
      CLC%z0(ix,iy) = 0.0002d0
   case default
      write(uerr,*) "The cell (",ix,",",iy,") of the horizontal background ",  &
         " grid is assigned a CLC class unrecognized. The execution stops here."
      stop 
end select
!------------------------
! Deallocations
!------------------------
return
end subroutine z0_CLC_table
#endif
