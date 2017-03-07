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
! Program unit: inidt2                                           
! Description: Initial time step. 
!-------------------------------------------------------------------------------
subroutine inidt2 
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
use I_O_file_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: j,npi,ii
double precision :: dtmin,dt_CFL,dt_dif,dt_vis,diffmax,U,vsc_coeff
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
dtmin = 1.0d+30                                                        
diffmax = zero
!------------------------
! Statements
!------------------------
if (indarrayFlu==0) then
! In case there is no fluid (or mobile) particle
   dt = 0.d0
   else
! Time step depends on 3 stability criteria:
! 1) CFL criterion: dt_CFL=CFL*(2h/(c+U))
! 2) viscous term criterion: dt_vis=vsc_coeff*(rho*h**2/(0.5*mu)) 
! 3) interface diffusion criterion: dt_diff=(h**2/2*teta)
      do ii=1,indarrayFlu
         npi = Array_Flu(ii)
         U = sqrt(pg(npi)%vel(1) ** 2 + pg(npi)%vel(2) ** 2 + pg(npi)%vel(3)   & 
             ** 2)
         dt_CFL = Domain%CFL * 2.d0 * Domain%h / (Med(pg(npi)%imed)%celerita + &
                  U)
         dt_vis = Domain%vsc_coeff * pg(npi)%dens * squareh / (half *          &
                  pg(npi)%mu)
         dtmin = min(dtmin,dt_CFL,dt_vis)
      enddo
      do j=1,NMedium
         diffmax = max(diffmax,Med(j)%codif)
      enddo
      if (diffusione) then
         dt_dif = half * squareh / (diffmax+0.000000001d0)
         dtmin = min(dtmin,dt_dif)
      endif
      dt = dtmin
endif
! To initialize the averaged time step 
dt_average = dt
!------------------------
! Deallocations
!------------------------
return
end subroutine inidt2

