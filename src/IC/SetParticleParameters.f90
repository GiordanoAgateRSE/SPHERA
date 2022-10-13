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
use Memory_I_O_interface_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: npi,Nz,Mate
double precision :: tstop
double precision,dimension(3) :: aux_vec
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
pg(npi)%volume = Domain%PVolume
pg(npi)%mass = pg(npi)%volume * Med(Mate)%den0
pg(npi)%dmass_dt = 0.d0
pg(npi)%dden_ALE12 = 0.d0
pg(npi)%p0_neg_ALE = .false.
#ifdef SPACE_2D
   pg(npi)%coord(2) = zero              
   pg(npi)%CoordOld(2) = zero
#endif
! Current velocity
call Vector_Product(partz(Nz)%omega,pg(npi)%coord,aux_vec,3)
pg(npi)%vel(1:3) = partz(Nz)%vel(1:3) + aux_vec(1:3)
pg(npi)%mom(1:3) = pg(npi)%vel(1:3) * pg(npi)%mass
! The initial velocity "vstart" is only influential in case of non-standard 
! motion
pg(npi)%vstart = partz(Nz)%vel
! To compute time stop for particle of type "law"
call stoptime(partz(Nz),tstop)
! To compute velocity for particle of type "law"
call vellaw(partz(Nz)%vlaw,partz(Nz)%vel,partz(Nz)%npointv)
! Initial ALE velocity increments and fluid velocity
pg(npi)%dvel_ALE1(1:3) = 0.d0
pg(npi)%dvel_ALE3(1:3) = 0.d0
pg(npi)%vel_fluid(1:3) = pg(npi)%vel(1:3)
! Stopping time for blocks in movement
pg(npi)%tstop = tstop
! Material ID
pg(npi)%imed = Mate
! Viscosities
pg(npi)%kin_visc = Med(Mate)%kin_visc
pg(npi)%mu = Med(Mate)%kin_visc * Med(Mate)%den0
! DB-SPH parameters
if (Domain%tipo=="bsph") then
   pg(npi)%Gamma = one 
   pg(npi)%rhoSPH_new = zero
   pg(npi)%uni = zero
   pg(npi)%sigma = zero
   pg(npi)%dShep = zero 
   pg(npi)%FS = 0 
   pg(npi)%Gamma_last_active = zero
   pg(npi)%DBSPH_inlet_ID = 0
   pg(npi)%DBSPH_outlet_ID = 0
endif
! Mixture density for granular SPH particles (bed-load transport) 
if (Granular_flows_options%KTGF_config>0) then
   if (Med(pg(npi)%imed)%tipo=="granular") then
      pg(npi)%dens = Med(pg(npi)%imed)%den0_s 
! IC viscosity from pure fluid, as I cannot calculate the right one at this 
! stage.         
      pg(npi)%mu = Med(Granular_flows_options%ID_main_fluid)%kin_visc *        &
                   Med(Granular_flows_options%ID_main_fluid)%den0
      pg(npi)%kin_visc = Med(Granular_flows_options%ID_main_fluid)%kin_visc
   endif
   call initialization_fixed_granular_particle(npi)
   pg(npi)%sigma_prime_m = 0.d0
   pg(npi)%pres_fluid = 0.d0
endif
! Particle status, depending on the velocity components (fluid or solid).
if (index(Med(Mate)%tipo,"liquid")>0) then
   pg(npi)%state = "flu"
   elseif (index(Med(Mate)%tipo,"granular")>0) then
      pg(npi)%state = "sol"
endif
if (index(Med(Mate)%tipo,"granular")>0) then
   if ((pg(npi)%vel(1)/=zero).or.(pg(npi)%vel(2)/=zero).or.                    &
      (pg(npi)%vel(3)/=zero)) pg(npi)%state = "flu"     
endif
! Motion index
pg(npi)%vel_type = partz(Nz)%move        
if (partz(Nz)%move/="std") pg(npi)%kin_visc = zero
! Grid cell 
pg(npi)%cella = ParticleCellNumber(pg(npi)%coord)
! Particle color definition, as defined in the input file.
call defcolpartzero(Nz,partz,pg(npi))
!------------------------
! Deallocations
!------------------------
return
end subroutine SetParticleParameters
