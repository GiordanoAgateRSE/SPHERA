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
! Program unit: fluid_particle_imposed_kinematics
! Description: Possible imposed velocities for the fluid particles.
!-------------------------------------------------------------------------------
subroutine fluid_particle_imposed_kinematics
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Dynamic_allocation_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: i_zone,npi
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
! One might loop over particles instead of over zones
do i_zone=1,NPartZone
   if (Partz(i_zone)%move/="fix") cycle
! It assigns the movement with an imposed kinematics ("npointv" velocity data)
   if (Partz(i_zone)%npointv>1) then
      call vellaw(Partz(i_zone)%vlaw,Partz(i_zone)%vel,Partz(i_zone)%npointv)
      write(ulog,"(f12.4,a,i2,a,3e12.4)") simulation_time,"  zona",i_zone,     &
         "  vel.",Partz(i_zone)%vel
! As an alternative, one may loop over particles 
      do npi=Partz(i_zone)%limit(1),Partz(i_zone)%limit(2) 
         if (pg(npi)%cella==0) cycle
         pg(npi)%var(:) = Partz(i_zone)%vel(:)
         if (simulation_time>=pg(npi)%tstop) then
            pg(npi)%vstart(:) = zero
            pg(npi)%vel(:) = zero
            else
               pg(npi)%vstart(:) = Partz(i_zone)%vel(:)
               pg(npi)%vel(:) = Partz(i_zone)%vel(:)
         endif
      enddo
! Fixed value of velocity (npointv=1)
      elseif (Partz(i_zone)%npointv==1) then
! As an alternative, one may loop over particles 
         do npi=Partz(i_zone)%limit(1),Partz(i_zone)%limit(2)
            if (pg(npi)%cella==0) cycle
            pg(npi)%var(:) = Partz(i_zone)%vel(:)
            if (simulation_time>=pg(npi)%tstop) then
               pg(npi)%vstart(:) = zero
               pg(npi)%vel(:) = zero
               else
                  pg(npi)%vel(:) = Partz(i_zone)%vel(:)
            endif
         enddo
   endif
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine fluid_particle_imposed_kinematics
