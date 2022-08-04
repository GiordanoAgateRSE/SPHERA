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
! Program unit: w
! Description: Kernel function. Cubic beta-spline function as defined in 
!              the paper of Monaghan & Lattanzio (1985).
!-------------------------------------------------------------------------------
double precision function w(r,h,coef)
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
double precision,parameter :: a1 = 0.666666667d0
double precision :: r,h,s,q,dms,coef
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
s = r / h
if (s<=1.d0) then
   q = a1 + s * s * (s * 0.5d0 - 1.d0)
   elseif (s>=2.d0) then
      q = 0.d0
      else 
         dms = 2.d0 - s
         q = dms * dms * dms / 6.d0
endif
w = q * coef
!------------------------
! Deallocations
!------------------------
return
end function w
