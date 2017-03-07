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
! Program unit: DBSPH_inlet_outlet
! Description: Impose boundary conditions at the inlet and outlet sections
!              (DB-SPH boundary treatment scheme).            
!-------------------------------------------------------------------------------
subroutine DBSPH_inlet_outlet(npi)
!------------------------
! Modules
!------------------------ 
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: npi
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
! By hypothesis a fluid particle cannot be close to two or more inlet/outlet 
! sections (otherwise it randomly takes an inlet velocity value)
! Inlet boundary conditions
if (pg(npi)%DBSPH_inlet_ID>0) pg(npi)%vel(:) =                                 &
   pg_w(pg(npi)%DBSPH_inlet_ID)%vel(:)
! Outlet boundary conditions
if (pg(npi)%DBSPH_outlet_ID>0) then
   pg(npi)%dens = Med(1)%den0
   pg(npi)%pres = pg_w(pg(npi)%DBSPH_outlet_ID)%pres
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine DBSPH_inlet_outlet

