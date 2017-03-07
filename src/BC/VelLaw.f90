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
! Program unit: VelLaw
! Description:  To impose an input kinematics to particles.   
!-------------------------------------------------------------------------------
subroutine VelLaw (vlaw,vel,np)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(IN) :: np
double precision,intent(IN),dimension(0:3,MAXPOINTSVLAW) :: vlaw
double precision,intent(OUT),dimension(3) :: vel
integer(4) :: n
double precision :: fra
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
if (np<=1) return
do n=2,np
   if (simulation_time>vlaw(0,n)) cycle
   fra = (simulation_time - vlaw(0,n-1)) / (vlaw(0,n) - vlaw(0,n-1))   
   vel(1:3) = vlaw(1:3,n-1) +  (vlaw(1:3,n) - vlaw(1:3,n-1)) * fra
   return
end do
vel(1:3) = vlaw(1:3,np)
!------------------------
! Deallocations
!------------------------
return
end subroutine VelLaw

