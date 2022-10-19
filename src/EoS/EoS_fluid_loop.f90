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
! Program unit: EoS_fluid_loop 
! Description: Loop over the fluid particles to apply the Equation of State 
!              (EoS) to assess pressure
!-------------------------------------------------------------------------------
subroutine EoS_fluid_loop
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
integer(4) :: npi
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine EoS_barotropic_linear(k_bulk,rho_ref,p_ref,rho_in,p_in,rho_out,  &
      p_out)
      implicit none
      double precision,intent(in) :: k_bulk,rho_ref,p_ref
      double precision,intent(in),optional :: rho_in,p_in
      double precision,intent(out),optional :: rho_out,p_out
   end subroutine EoS_barotropic_linear
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
if (on_going_time_step>0) then  
! Loop over all the particles
!$omp parallel do default(none) shared(nag,pg,Domain,Med) private(npi)
   do npi=1,nag
! Skip the outgone particles and the particles with velocity type different 
! from "std"
      if ((pg(npi)%cella<=0).or.(pg(npi)%vel_type/="std")) cycle
      call EoS_barotropic_linear(Med(pg(npi)%imed)%eps,Med(pg(npi)%imed)%den0, &
         Domain%prif,rho_in=pg(npi)%dens,p_out=pg(npi)%pres)
   enddo
!$omp end parallel do
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine EoS_fluid_loop
