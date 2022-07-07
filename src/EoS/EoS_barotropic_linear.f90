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
! Program unit: EoS_barotropic_linear 
! Description: Barotropic linear Equation of State (EoS) for pressure 
!              estimation and its inverse for density estimation
!-------------------------------------------------------------------------------
subroutine EoS_barotropic_linear(k_bulk,rho_ref,p_ref,rho_in,p_in,rho_out,p_out)
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(in) :: k_bulk,rho_ref,p_ref
double precision,intent(in),optional :: rho_in,p_in
double precision,intent(out),optional :: rho_out,p_out
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
if (present(rho_in)) then
! EoS
   p_out = k_bulk / rho_in * (rho_in - rho_ref) + p_ref
   elseif (present(p_in)) then
! inverse of EoS
      rho_out = rho_ref / (1.d0 - (p_in - p_ref) / k_bulk)
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine EoS_barotropic_linear
