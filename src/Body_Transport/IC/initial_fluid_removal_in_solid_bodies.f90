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
! Program unit: initial_fluid_removal_in_solid_bodies
! Description: Initial removal of possible SPH fluid particles within solid 
!              bodies (due to design) at the beginning of the simulation 
!              (initial conditions for fluid particle positions in case of 
!              Fluid - Structure interactions).
!-------------------------------------------------------------------------------
#ifdef SOLID_BODIES
subroutine initial_fluid_removal_in_solid_bodies
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
integer(4) :: npi,j,npartint,npj
double precision :: dis_ref,dis_bp_f
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
! Distance threshold to remove fluid particles
#ifdef SPACE_3D
dis_ref = c_ini_rem_fp_sb * dsqrt(3.d0) / 4.d0 * Domain%dx 
#elif defined SPACE_2D
dis_ref = c_ini_rem_fp_sb * dsqrt(2.d0) / 4.d0 * Domain%dx
#endif
! Loop over the body particles 
!$omp parallel do default(none)                                                &
!$omp shared(n_body_part,nPartIntorno_bp_f,NMAXPARTJ,PartIntorno_bp_f,pg)      &
!$omp shared(rag_bp_f,OpCount,dis_ref,bp_arr,Domain)                           &
!$omp private(npi,j,npartint,npj,dis_bp_f)
! Loop over solid particles
do npi=1,n_body_part
! Loop over the neighbouring fluid particles
   do j=1,nPartIntorno_bp_f(npi)
      npartint = (npi - 1) * NMAXPARTJ + j
      npj = PartIntorno_bp_f(npartint)
      dis_bp_f = dsqrt(dot_product(rag_bp_f(1:3,npartint),                     &
                 rag_bp_f(1:3,npartint)))
      if (dis_bp_f<=dis_ref) then
!$omp critical (omp_initial_fluid_removal_in_solid_bodies)
         if (pg(npj)%cella>-1) then
            pg(npj)%cella = -1
            OpCount(pg(npj)%imed) = OpCount(pg(npj)%imed) + 1
         endif
!$omp end critical (omp_initial_fluid_removal_in_solid_bodies)
      endif
   enddo
enddo
!$omp end parallel do
fluid_in_body_count = 0
!------------------------
! Deallocations
!------------------------
return
end subroutine initial_fluid_removal_in_solid_bodies
#endif
