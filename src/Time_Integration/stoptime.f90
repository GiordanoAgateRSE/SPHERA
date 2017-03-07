!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)

! SPHERA authors and email contact are provided on SPHERA documentation.

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
! Program unit: stoptime                                          
! Description: Stopping time.
!-------------------------------------------------------------------------------

subroutine stoptime(partzlocal,tstop)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(INOUT) :: tstop
type(TyZone),intent(INOUT) :: partzlocal
logical :: out
integer(4) :: k,n,icord
double precision :: tstopc,acc,deltat,spo,dspo,rad
double precision,dimension(3) :: dxyz
double precision,dimension(3,2) :: vlimits,tlimits
character(len=lencard) :: nomsub = "stoptime"
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
tstop = max_positive_number
!------------------------
! Statements
!------------------------
! To check on fixed particles
if (partzlocal%move=="fix") then
   do n=1,ncord
      icord = icoordp(n,ncord-1)
      if (partzlocal%vel(icord)>zero) then
         tstopc = (Domain%coord(icord,2) - partzlocal%coordMM(icord,2) - one * &
            Domain%dx) / partzlocal%vel(icord)
         elseif (partzlocal%vel(icord)<zero) then
            tstopc = (Domain%coord(icord,1) - partzlocal%coordMM(icord,1) +    & 
               one * Domain%dx) / partzlocal%vel(icord)
            else
               tstopc = Domain%tmax
      endif
      tstop = min(tstop,tstopc)
   enddo
! To check on particles with imposed velocity 
   elseif (partzlocal%move=="law") then
! Trajectories
      vlimits = zero
      do n=1,ncord
         icord = icoordp(n,ncord-1)
         vlimits(icord,1) = Domain%coord(icord,1)
         vlimits(icord,2) = Domain%coord(icord,2)
      enddo
      tlimits = zero
      dxyz = zero
      out = .FALSE.
      LAW_ZONE_LOOP: do k=2,partzlocal%npointv
         COORDS_LOOP: do n=1,ncord
            icord = icoordp(n,ncord-1)
! Accelerations
            deltat = partzlocal%vlaw(0,k) - partzlocal%vlaw(0,k-1)
            acc = (partzlocal%vlaw(icord,k) - partzlocal%vlaw(icord,k-1)) /    &
               deltat
! Trajectories
            dspo = partzlocal%vlaw(icord,k-1) * deltat + acc * deltat * deltat &
               * half
            spo  = dxyz(icord) + dspo
! To check on the minimum limit has been overridden and how much time has 
! been required.
            if ((partzlocal%coordMM(icord,1)+spo)<vlimits(icord,1)) then
               out = .TRUE.
               dspo = vlimits(icord,1) - (partzlocal%coordMM(icord,1) +        &
                  dxyz(icord))
               if (acc==zero) then
                  if (partzlocal%vlaw(icord,k-1)==zero) then
                     deltat = max_positive_number
                     else
                        deltat = dspo / partzlocal%vlaw(icord,k-1)
                  endif
                  else
                     rad  = partzlocal%vlaw(icord,k-1) *                       &
                        partzlocal%vlaw(icord,k-1) - 4.0d0 * 0.5d0 * acc * dspo
                     if (rad>=zero) then
                        rad = Dsqrt(rad)  
                        deltat = (partzlocal%vlaw(icord,k-1) + rad) / (2.d0 *  &
                           0.5d0 * acc)
                        else
                           call diagnostic(arg1=10,arg2=88,arg3=nomsub)       
                     endif
               endif
            endif
! To add the interval time to the total time
            tlimits(icord,1) = tlimits(icord,1) + deltat  
! To check on the minimum limit has been overridden and how much time has 
! been required.
            if ((partzlocal%coordMM(icord,2)+spo)>vlimits(icord,2)) then 
               out = .TRUE.
               dspo = vlimits(icord,2) - (partzlocal%coordMM(icord,2) +       &
                  dxyz(icord))
               if (acc==zero) then
                  if (partzlocal%vlaw(icord,k-1)==zero) then
                     deltat = max_positive_number
                     else
                        deltat = dspo / partzlocal%vlaw(icord,k-1)
                  endif
                  else
                  rad = partzlocal%vlaw(icord,k-1) * partzlocal%vlaw(icord,k-1)&
                     - 4.0d0 * 0.5d0 * acc * dspo
                  if (rad>=zero) then
                     rad = Dsqrt(rad)  
                     deltat = (-partzlocal%vlaw(icord,k-1) + rad) / (2.d0 *    &
                        0.5d0 * acc)
                     else
                        call diagnostic(arg1=10,arg2=88,arg3=nomsub)       
                  endif
               endif
            endif
! To add the interval time to the total time
            tlimits(icord,2) = tlimits(icord,2) + deltat  
! To save the evaluated displacement along the trajectory
            dxyz(icord) = dxyz(icord) + dspo
         enddo COORDS_LOOP
! It ends if particle has gone out of the domain
         if (out) exit LAW_ZONE_LOOP
      enddo LAW_ZONE_LOOP
! To evaluate the minimum time
      do n=1,ncord
         icord = icoordp(n,ncord-1)
         tstop = min(tstop,tlimits(icord,1),tlimits(icord,2))
      enddo
      partzlocal%move = "fix"
      else
         tstop = Domain%tmax
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine stoptime

