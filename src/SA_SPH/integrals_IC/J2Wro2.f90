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
! Program unit: J2Wro2
! Description: To compute the integral J2Wro2 (Di Monaco et al., 2011, EACFM)                        
!-------------------------------------------------------------------------------
#ifdef SPACE_2D
double precision function J2Wro2(ro)
!------------------------
! Modules
!------------------------
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
! 1/120
double precision,parameter :: a1 = 0.00833333333333333d0
! 8/30
double precision,parameter :: a2 = 0.26666666666666667d0
double precision,intent(in) :: ro
double precision :: ro2,ro3
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
ro2 = ro * ro
ro3 = ro2 * ro
!------------------------
! Statements
!------------------------
if ((zero<=ro).and.(ro<one)) then
   J2Wro2 = KERNELCONST2D * (0.25d0 - (a1 * (40.d0 - 36.d0 * ro2 + 15.d0 * ro3 &
            ) * ro3))
   elseif ((one<=ro).and.(ro<two)) then
      J2Wro2 = KERNELCONST2D * (a2 - (a1 * (80.d0 - 90.d0 * ro + 36.d0 * ro2 - &
               5.d0 * ro3) * ro3))
      else
         J2Wro2 = zero
endif
!------------------------
! Deallocations
!------------------------
return
end function J2Wro2
#endif
