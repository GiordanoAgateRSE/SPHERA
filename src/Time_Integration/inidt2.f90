!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-Ricerca di Sistema)



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
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: inidt2                                           
! Description: Initial time step. 
!----------------------------------------------------------------------------------------------------------------------------------

subroutine inidt2 
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
integer(4) :: j,npi,mate,ii
double precision :: dtmin,dt_CFL,dt_dif,dt_vis,diffmax,U
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
dtmin   = 1.0d+30                                                        
diffmax = zero
!------------------------
! Statements
!------------------------
! Time step depends on 3 conditions:
! 1) the CFL condition: dt_CFL=min(2h/(c+U))
! 2) viscous stability condition dt_vis=min(rho*h**2/(0.5*mu)) 
! 3) interface diffusion condition dt_diff=(h**2/2*teta)
do ii=1,indarrayFlu
   npi = Array_Flu(ii)
   mate = pg(npi)%imed
   U = sqrt(pg(npi)%vel(1) ** 2 + pg(npi)%vel(2) ** 2 + pg(npi)%vel(3) ** 2)
   dt_CFL = 2.d0 * Domain%h / (Med(mate)%celerita+U)
   dt_vis = pg(npi)%dens * squareh / (half * pg(npi)%mu)
   dtmin  = min(dtmin,dt_CFL,dt_vis)
enddo
do j=1,nmedium
   diffmax = max(diffmax,Med(j)%codif)
enddo
dt_dif = half * squareh / (diffmax+0.000000001d0)
dtmin = min (dtmin,dt_dif)
! Initial dt for a jet
if (indarrayFlu==0) dtmin = 2.*Domain%h/(Med(1)%celerita)
! CFL is used as a constant for every condition (the CFL, the viscous and 
! the diffusive one)
dt = Domain%CFL * dtmin
! To initialize the averaged time step 
dt_average = dt
!------------------------
! Deallocations
!------------------------
return
end subroutine inidt2

