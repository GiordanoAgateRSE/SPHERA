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
! Program unit: write_Granular_flows_interfaces                  
! Description: To print the interfaces for bed-load transport phenomena.             
!-------------------------------------------------------------------------------
subroutine write_Granular_flows_interfaces(i_grid,j_grid,pos_aux)
!------------------------
! Modules
!------------------------
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: i_grid,j_grid
double precision,intent(in) :: pos_aux(3)
double precision :: z_free_surface,z_BedLoad_PureFluid,z_bed,z_soil_bottom
double precision :: z_sat_top
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
if (ind_interfaces(i_grid,j_grid,1)>0) then
   z_free_surface = pg(ind_interfaces(i_grid,j_grid,1))%coord(3) + Domain%dx / &
                    2.d0
   else
      z_free_surface = -999.d0
endif
if (ind_interfaces(i_grid,j_grid,4)>0) then
   z_bed = pg(ind_interfaces(i_grid,j_grid,4))%coord(3) + Domain%dx / 2.d0
   else
      z_bed = -999.d0
endif
if (ind_interfaces(i_grid,j_grid,3)>0) then
   z_BedLoad_PureFluid = pg(ind_interfaces(i_grid,j_grid,3))%coord(3) +        &
                         Domain%dx / 2.d0
   else
      z_BedLoad_PureFluid = z_bed 
endif
if (ind_interfaces(i_grid,j_grid,5)>0) then
   z_soil_bottom = pg(ind_interfaces(i_grid,j_grid,5))%coord(3) + Domain%dx /  &
                   2.d0
   else
      z_soil_bottom = -999.d0 
endif
if (ind_interfaces(i_grid,j_grid,6)>0) then
   z_sat_top = pg(ind_interfaces(i_grid,j_grid,6))%coord(3) + Domain%dx / 2.d0
   else
      z_sat_top = z_soil_bottom
endif
write(ncpt,'(8(g14.6,1x))') simulation_time,pos_aux(1),pos_aux(2),             &
                            z_free_surface,z_BedLoad_PureFluid,z_bed,          &
                            z_soil_bottom,z_sat_top
!------------------------
! Deallocations
!------------------------
return
end subroutine write_Granular_flows_interfaces

