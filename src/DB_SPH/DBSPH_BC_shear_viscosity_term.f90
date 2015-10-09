!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-; SPHERA has been authored for RSE SpA by 
!    Andrea Amicarelli, Antonio Di Monaco, Sauro Manenti, Elia Bon, Daria Gatti, Giordano Agate, Stefano Falappi, 
!    Barbara Flamini, Roberto Guandalini, David Zuccalà).
! Main numerical developments of SPHERA: 
!    Amicarelli et al. (2015,CAF), Amicarelli et al. (2013,IJNME), Manenti et al. (2012,JHE), Di Monaco et al. (2011,EACFM). 
! Email contact: andrea.amicarelli@rse-web.it

! This file is part of SPHERA.
! SPHERA is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with SPHERA. If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Program unit: DBSPH_BC_shear_viscosity_term
! Description: Computation of the contributions to the numerator of the boundary
!              shear viscosity term in DB-SPH-NS.          
!-------------------------------------------------------------------------------
subroutine DBSPH_BC_shear_viscosity_term(i_0,i_a,npartint,DBSPH_BC_she_vis_term)
!------------------------
! Modules
!------------------------ 
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: i_0,i_a,npartint
double precision,intent(inout) :: DBSPH_BC_she_vis_term(3)
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
! Semi-particle velocity gradient (times shear viscosity) in the Viscous  
! Sub-Layer of the Surface Neutral Boundary Layer
pg_w(i_a)%grad_vel_VSL_times_mu(:) = pg_w(i_a)%grad_vel_VSL_times_mu(:) /      &
                                     pg_w(i_a)%sigma
! Contribution to the numerator of the boundary shear viscosity term 
DBSPH_BC_she_vis_term(:) = DBSPH_BC_she_vis_term(:) -                          &
                           kernel_fw(1,npartint) * pg_w(i_a)%weight *          &
                           (pg_w(i_a)%grad_vel_VSL_times_mu(:) +               &
                           grad_vel_VSL_fw(:,npartint) * pg(i_0)%mu)
!------------------------
! Deallocations
!------------------------
return
end subroutine DBSPH_BC_shear_viscosity_term

