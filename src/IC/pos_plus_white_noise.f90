!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2020 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: pos_plus_white_noise
! Description: Add a white noise to a particle position (initial conditions).
!-------------------------------------------------------------------------------
subroutine pos_plus_white_noise(max_rnd_eps,pos)
!------------------------
! Modules
!------------------------
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
! maximum random displacement (ratio with respect to dx)
double precision,intent(in) :: max_rnd_eps
double precision,intent(inout) :: pos(3)
! random variable (uniform distribution between 0.d0 and 1.d0)
double precision :: rnd(3)
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
call random_seed()
!------------------------
! Statements
!------------------------
call random_number(rnd(1))
call random_number(rnd(2))
call random_number(rnd(3))
pos(:) = pos(:) + (2.d0 * rnd(:) - 1.d0) * max_rnd_eps * Domain%dx
!------------------------
! Deallocations
!------------------------
return
end subroutine pos_plus_white_noise
