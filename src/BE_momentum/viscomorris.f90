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
! Program unit: viscomorris
! Description: To compute the volume inter-particle contributions to the shear
! viscosity term in the momentum equation (Morris, 1997, JCP). Both interactions
! "fluid particle - fluid particle" and "fluid particle - semi-particle" (DBSPH) 
! are considered.
!-------------------------------------------------------------------------------
subroutine viscomorris(npi,npj,npartint,mass_comput_part,dens_comput_part,     &
kin_visc_comput_part,mass_neighbour,dens_neighbour,kin_visc_neighbour,         &
kernel_der,vel_type,rel_dis,dervel,rvw)
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
integer(4),intent(in) :: npi,npj,npartint
double precision,intent(in) :: mass_comput_part,dens_comput_part
double precision,intent(in) :: kin_visc_comput_part,mass_neighbour
double precision,intent(in) :: dens_neighbour,kin_visc_neighbour,kernel_der
double precision,intent(in) :: rel_dis(3),dervel(3)
character(3),intent(in) :: vel_type
double precision,intent(out) :: rvw(3)
double precision :: amassj,rhotilde,anuitilde,factivis
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
if (vel_type/="std") then
   amassj = mass_comput_part
   rhotilde = dens_comput_part
   anuitilde = two * kin_visc_comput_part
   else
      amassj = mass_neighbour
      rhotilde = (kin_visc_comput_part * dens_comput_part + kin_visc_neighbour &
                  * dens_neighbour + 0.001d0)
      anuitilde = 4.0d0 * (kin_visc_comput_part * kin_visc_neighbour)
endif
factivis = amassj * anuitilde / rhotilde
rvw(1:3) = factivis * (-dervel(1:3) * kernel_der * dot_product(rel_dis,rel_dis))
!------------------------
! Deallocations
!------------------------
return
end subroutine viscomorris

