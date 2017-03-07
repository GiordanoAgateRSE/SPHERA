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
! Program unit: fixed_bed_slope_limited
! Description: Forced deposition (or no erosion) for particles at least 2h below
!              the fixed bed (as it is defined in the associated column) during
!              the same time step: i.e. the maximum slope of the fixed bed is 
!              2h/2h. This avoids eventual too fast propagation of erosion along
!              the vertical (erosion is an interface phenomenon).                    
!-------------------------------------------------------------------------------
subroutine fixed_bed_slope_limited(npi,igridi,jgridi,test)
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
integer(4),intent(in) :: npi,igridi,jgridi
logical,intent(inout) :: test
integer(4) :: aux_ID
double precision :: Velocity2,fixed_bed_tolerance,pretot
double precision :: aux_vec(3)
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
if (ncord==3) then
   fixed_bed_tolerance = 4.d0 * Domain%h
   elseif (ncord==2) then
      fixed_bed_tolerance = 2.d0 * Domain%h
endif    
if (ind_interfaces(igridi,jgridi,4)>0) then
   if (pg(npi)%coord(3)<(pg(ind_interfaces(igridi,jgridi,4))%coord(3)-         &
      fixed_bed_tolerance)) then
      if (pg(npi)%state=="flu") then
         pg(npi)%state = "sol"   
         pg(npi)%vel = 0.d0
         pg(npi)%var = 0.d0
         pg(npi)%sigma_prime_m = 0.0d0
         pg(npi)%pres_fluid = 0.d0  
      endif
      if (pg(npi)%indneighliqsol.ne.0) then
         aux_ID = pg(npi)%indneighliqsol
         elseif (pg(npi)%state=="flu") then 
            aux_ID = pg(npi)%ind_neigh_mob_for_granmob    
            else
               aux_ID = pg(npi)%ind_neigh_mix_bed 
      endif
      if (aux_ID.ne.0) then
         aux_vec(:) = pg(aux_ID)%vel_old(:) - pg(npi)%vel_old(:)
         Velocity2 = dot_product(aux_vec,aux_vec) 
         pretot = pg(aux_ID)%pres  + (pg(aux_ID)%coord(3) - pg(npi)%coord(3))  &
                  * ( - Domain%grav(3)) * med(pg(aux_ID)%imed)%den0 +          &
                  Velocity2 * half * Med(pg(aux_ID)%imed)%den0
         pg(npi)%dens = med(pg(npi)%imed)%den0 + (pretot /                     &
                        (Med(pg(npi)%imed)%celerita *                          &
                        Med(pg(npi)%imed)%celerita))   
         else
            if (ind_interfaces(igridi,jgridi,3)>0) then
               pretot = pg(ind_interfaces(igridi,jgridi,3))%pres +             &
                        (pg(ind_interfaces(igridi,jgridi,3))%coord(3) -        &
                        pg(npi)%coord(3)) * ( - Domain%grav(3)) *              &
                        med(pg(npi)%imed)%den0
               pg(npi)%dens = med(pg(npi)%imed)%den0 + (pretot /               &
                              (Med(pg(npi)%imed)%celerita *                    &
                              Med(pg(npi)%imed)%celerita))
            endif
      endif 
      test = .false.
   endif 
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine fixed_bed_slope_limited

