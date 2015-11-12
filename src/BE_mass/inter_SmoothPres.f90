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
! Program unit: inter_SmoothPres
! Description: To calculate a corrective term for pressure. 
!----------------------------------------------------------------------------------------------------------------------------------

subroutine inter_SmoothPres
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
integer(4) :: npi,npj,contj,npartint,ii
double precision :: unity,presi,presj,rhoj,amassj,pesoj,appo1,appo2,TetaP1
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
!$omp parallel do default(none) &
!$omp private(npi,ii,unity,appo1,appo2,contj,npartint,npj,presi,rhoj,presj)    &
!$omp private(amassj,pesoj)                                                    &
!$omp shared(nag,pg,Med,Domain,nPartIntorno,NMAXPARTJ,PartIntorno,PartKernel)  &
!$omp shared(indarrayFlu,Array_Flu)
! Loop over all the particles
do ii=1,indarrayFlu
   npi = Array_Flu(ii)
   unity = zero
   appo1 = zero
   appo2 = zero
   do contj=1,nPartIntorno(npi)
      npartint = (npi-1)* NMAXPARTJ + contj
      npj = PartIntorno(npartint)
      if (pg(npj)%vel_type/="std") cycle
      presi = pg(npi)%pres
      rhoj = pg(npj)%dens    
      presj = pg(npj)%pres    
      amassj = pg(npj)%mass
      pesoj = amassj * PartKernel(4,npartint) / rhoj
      unity = unity + pesoj  
      appo1 = appo1 + (presj - presi) * pesoj  
      appo2 = appo2 - Domain%grav(3) * Med(pg(npi)%imed)%den0 *                &
              (pg(npj)%coord(3) - pg(npi)%coord(3)) * pesoj     
   end do
   if (Domain%Psurf/='s') then
      pg(npi)%vpres = appo1
      else if (unity>0.8d0) then
         pg(npi)%vpres=appo1
         else
            pg(npi)%vpres = appo1 + appo2
   end if
   pg(npi)%uni = unity
end do
!$omp end parallel do
!$omp parallel do default(none)                                                &
!$omp private(npi,ii,TetaP1)                                                   &
!$omp shared(nag,Pg,Med,Domain,dt,indarrayFlu,Array_Flu,esplosione)
do ii = 1,indarrayFlu
   npi = Array_Flu(ii)
   if (esplosione) then
      TetaP1 = Domain%TetaP * pg(npi)%Csound * dt / Domain%h
      else
! Computation of TetaP depending on the time step
         TetaP1 = Domain%TetaP * Med(pg(npi)%imed)%Celerita * dt / Domain%h
   end if
   pg(npi)%pres = pg(npi)%pres + TetaP1 * pg(npi)%vpres / pg(npi)%uni
! To update density depending on pressure
   pg(npi)%dens = Med(pg(npi)%imed)%den0 * (one + pg(npi)%pres /               &
                  Med(pg(npi)%imed)%eps)
end do
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine inter_SmoothPres

