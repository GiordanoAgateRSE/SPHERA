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
! Program unit: rundt2                                          
! Description: Time step computation according to stability constraints (inertia terms, visosity term, interface diffusion -not 
!              recommended-). Plus. a special treatment for Monaghan artificial viscosity term and management of low-velocity 
!              SPH mixture particles for bed-laod transport phenomena.
!----------------------------------------------------------------------------------------------------------------------------------

subroutine rundt2
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
double precision :: dtmin,dt_CFL,dt_dif,dt_vis,diffmax,U,celiq
double precision :: visc_Mon_j,dt_Mon_j,dt_Mon
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
dt_Mon = max_positive_number 
!------------------------
! Statements
!------------------------
! Time step depends on 3 conditions:
! 1) the CFL condition: dt_CFL=min(2h/(c+U))
! 2) viscous stability condition dt_vis=min(rho*h**2/(0.5*mu)) 
! 3) interface diffusion condition dt_diff=(h**2/2*teta)
do ii=1,indarrayFlu
   npi = Array_Flu(ii)
   if (Granular_flows_options%ID_erosion_criterion==1) then
      if (pg(npi)%state=="sol") cycle
      if ((pg(npi)%coord(1)<Granular_flows_options%x_min_dt).or.               &
         (pg(npi)%coord(1)>Granular_flows_options%x_max_dt).or.                &
         (pg(npi)%coord(2)<Granular_flows_options%y_min_dt).or.                &
         (pg(npi)%coord(2)>Granular_flows_options%y_max_dt).or.                &
         (pg(npi)%coord(3)<Granular_flows_options%z_min_dt).or.                &
         (pg(npi)%coord(3)>Granular_flows_options%z_max_dt)) then
         cycle
      endif    
   endif    
   mate = pg(npi)%imed
   U = sqrt(pg(npi)%vel(1) ** 2 + pg(npi)%vel(2) ** 2 + pg(npi)%vel(3) ** 2)
   dt_CFL = 2.d0 * Domain%h / (Med(mate)%celerita + U)
   dt_vis = pg(npi)%dens * squareh / (half * pg(npi)%mu)
   dtmin = min(dtmin,dt_CFL,dt_vis)
   celiq = Med(mate)%eps / pg(npi)%dens
   if (celiq>=zero) pg(npi)%Csound = Dsqrt(celiq)
enddo
do j=1,nmedium
   diffmax = max(diffmax,Med(j)%codif)
   if (dt_alfa_Mon.eqv..true.) then
      visc_Mon_j = Med(j)%alfaMon * Med(j)%celerita * Domain%h /               &
         Med(j)%den0
      dt_Mon_j = squareh / (half * visc_Mon_j)
      dt_Mon = min(dt_Mon,dt_Mon_j)
   endif
enddo
dt_dif = half * squareh / (diffmax+0.000000001d0)
dtmin = min(dtmin,dt_dif)
if (dt_alfa_Mon.eqv..true.) dtmin = min(dtmin,dt_Mon)
! CFL is used as a constant for each of the 3 stability conditions 
dt = (one - pesodt) * Domain%CFL * dtmin + pesodt * dt_average
dt_average = (dt_average * (it_corrente - 1) + dt) / it_corrente
!------------------------
! Deallocations
!------------------------
return
end subroutine rundt2

