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
! Program unit: Gamma_boun 
! Description: Interpolative function defined by Monaghan (2005) for boundary force particles (Amicarelli et al.,2015,CAF).     
!----------------------------------------------------------------------------------------------------------------------------------

double precision function Gamma_boun(r,h)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
double precision :: r,h,q
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
q = r / h
!------------------------
! Statements
!------------------------
if (q<=(2.d0/3.d0)) then
   Gamma_boun = 2.d0 / 3.d0
   else if (q<=1.0d0) then
      Gamma_boun = 2.d0 * q - (3.d0 / 2.d0) * q * q 
      else if (q<=2.0d0) then
         Gamma_boun = 0.5d0 * ((2.d0 - q) ** 2)
         else
            Gamma_boun = 0.d0
end if
Gamma_boun = dabs(Gamma_boun)
!------------------------
! Deallocations
!------------------------
return
end function Gamma_boun 

