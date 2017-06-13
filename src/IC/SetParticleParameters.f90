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
! Program unit: SetParticleParameters     
! Description: Setting initial particle parameters.               
!-------------------------------------------------------------------------------
subroutine SetParticleParameters(npi,Nz,Mate)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use I_O_file_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(IN) :: npi,Nz,Mate
double precision :: tstop
integer(4), external :: ParticleCellNumber
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
if (Domain%RKscheme>1) ts0_pg(npi) = ts_pgZero
pg(npi)%izona = Nz     
if (ncord == 2) then
   pg(npi)%mass = Domain%PVolume * Med(Mate)%den0 
   pg(npi)%coord(2) = zero              
   pg(npi)%CoordOld(2) = zero              
   else 
      pg(npi)%mass = Domain%PVolume * Med(Mate)%den0
endif
! Current velocity and initial velocity
pg(npi)%vel = partz(Nz)%vel
pg(npi)%vstart = partz(Nz)%vel
! To compute time stop for particle of type "law"
call stoptime (partz(Nz),tstop)
! TO compute velocity for particle of type "law"
call vellaw (partz(Nz)%vlaw,partz(Nz)%vel,partz(Nz)%npointv)
! Stopping time for blocks in movement
pg(npi)%tstop = tstop                    
! Material ID
pg(npi)%imed = Mate
! Viscosities
pg(npi)%visc = Med(Mate)%visc
pg(npi)%mu   = Med(Mate)%visc * Med(Mate)%den0
! DB-SPH parameters
if (Domain%tipo=="bsph") then
   pg(npi)%Gamma = one 
   pg(npi)%rhoSPH_new = zero
   pg(npi)%uni = zero
   pg(npi)%sigma = zero
   pg(npi)%sigma_same_fluid = zero
   pg(npi)%dShep = zero 
   pg(npi)%FS = 0 
   pg(npi)%Gamma_last_active = zero
   pg(npi)%DBSPH_inlet_ID = 0
   pg(npi)%DBSPH_outlet_ID = 0
endif
! Mixture density for granular SPH particles (bed-load transport) 
if (Granular_flows_options%ID_erosion_criterion==1) then
   if (Med(pg(npi)%imed)%tipo=="granular") then
      pg(npi)%dens = Med(pg(npi)%imed)%den0_s 
! IC viscosity from pure fluid, as I cannot calculate the right one at this 
! stage.         
      pg(npi)%mu = Med(Granular_flows_options%ID_main_fluid)%visc *            &
                   Med(Granular_flows_options%ID_main_fluid)%den0
      pg(npi)%visc = Med(Granular_flows_options%ID_main_fluid)%visc
   endif
   call initialization_fixed_granular_particle(npi)
   pg(npi)%sigma_prime_m = 0.0d0
   pg(npi)%pres_fluid = 0.0d0
endif
! Particle status, depending on the velocity components (fluid or solid).
if ((index(Med(Mate)%tipo,"liquid")>0).or.                                     &
   (index(Med(Mate)%tipo,"smagorin")>0)) then
   pg(npi)%state = "flu" 
   elseif ((index(Med(Mate)%tipo,"granular")>0).or.                           &
      (index(Med(Mate)%tipo,"general")>0)) then
      pg(npi)%state = "sol"
      elseif (index(Med(Mate)%tipo,"gas")>0) then
         pg(npi)%state = "flu" 
endif
if ((index(Med(Mate)%tipo,"granular")>0).or.                                   &
   (index(Med(Mate)%tipo,"general")>0)) then
   if ((pg(npi)%vel(1)/=zero).or.(pg(npi)%vel(2)/=zero).or.                    &
      (pg(npi)%vel(3)/=zero)) pg(npi)%state = "flu"     
endif
! Motion index
pg(npi)%vel_type = partz(Nz)%move        
if (partz(Nz)%move/="std") pg(npi)%visc = zero
! Boundary slip condition
pg(npi)%slip = partz(Nz)%slip  
! Grid cell 
pg(npi)%cella = ParticleCellNumber(pg(npi)%coord)
! Particle color definition, as defined in the input file.
call defcolpartzero (Nz,partz,pg(npi))
! Modulo diffusione
if (diffusione) then
   if (pg(npi)%imed==1) then
      pg(npi)%VolFra = VFmn
   end if
   if (pg(npi)%imed==2) then          
      pg(npi)%VolFra = VFmx
   end if
   else
      pg(npi)%VolFra = one
end if
if (esplosione) then
   pg(npi)%IntEn  = Med(pg(npi)%imed)%InitialIntEn
end if
!------------------------
! Deallocations
!------------------------
return
end subroutine SetParticleParameters

