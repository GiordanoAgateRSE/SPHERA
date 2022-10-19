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
! Program unit: ALE3_term_momentum
! Description: ALE3 term: computation and contribution to the Momentum Equation.
!-------------------------------------------------------------------------------
subroutine ALE3_term_momentum(dt_previous_step,dtvel)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(in) :: dt_previous_step
double precision,intent(inout) :: dtvel
integer(4) :: ii,npi,contj,npartint,npj
double precision :: d_rho_dvelALE1(3),dvel(3)
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
! dt computation
dtvel = 0.5d0 * (dt + dt_previous_step)
!$omp parallel do default(none)                                                &
!$omp shared(pg,indarrayFlu,Array_Flu,Med,nPartIntorno,NMAXPARTJ,PartIntorno)  &
!$omp shared(PartKernel,rag,input_any_t,dtvel)                                 &
!$omp private(npi,ii,contj,npartint,npj,dvel,d_rho_dvelALE1)
! Loop over particles
do ii = 1,indarrayFlu
   npi = Array_Flu(ii)
! Initialization of the ALE3 velocity increment
   pg(npi)%dvel_ALE3(1:3) = 0.d0
! The mixture particles in the elastic-plastic strain regime are held fixed
   if (pg(npi)%mu>(Med(pg(npi)%imed)%mumx*(1.d0-1.d-9))) then
      cycle
   endif
   if (.not.(pg(npi)%p0_neg_ALE)) then
      do contj=1,nPartIntorno(npi)
         npartint = (npi - 1) * NMAXPARTJ + contj
         npj = PartIntorno(npartint)
         if (npi==npj) cycle
! Update of the ALE3 velocity increment (here expressed as acceleration). BC 
! ALE3 terms for ME+CV are null. At this stage ALE1 velocities are still 
! saved as accelerations.
         dvel(1:3) = pg(npj)%vel(1:3) - pg(npi)%vel(1:3)
         d_rho_dvelALE1(1:3) = (pg(npj)%dens * pg(npj)%dvel_ALE1(1:3) +        &
                               pg(npi)%dens * pg(npi)%dvel_ALE1(1:3)) * dtvel
         pg(npi)%dvel_ALE3(1:3) = pg(npi)%dvel_ALE3(1:3) + dvel(1:3) *         &
                                  pg(npj)%volume * PartKernel(1,npartint) *    &
                                  dot_product(d_rho_dvelALE1(1:3),             &
                                  rag(1:3,npartint)) / (2.d0 * pg(npi)%dens)
      enddo
! Contribution to the ALE3 term in the ME-VC
      pg(npi)%acc(1:3) = pg(npi)%acc(1:3) + pg(npi)%dvel_ALE3(1:3)
   endif
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine ALE3_term_momentum
