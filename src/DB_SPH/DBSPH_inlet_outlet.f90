!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-) 
!      
!     
!   
!      
!  

! This file is part of SPHERA.
!  
!  
!  
!  
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
!  
!  
!  
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: DBSPH_inlet_outlet
! Description: Impose boundary conditions at the inlet and outlet sections (DB-SPH boundary treatment scheme).            
!----------------------------------------------------------------------------------------------------------------------------------

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

