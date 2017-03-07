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
! Program unit: BoundaryMassForceMatrix2D                                
! Description: Generation of the generalised boundary mass force matrix RN, on the base of the cosine matrix 
!              T and the parameter Fi.
!              (Di Monaco et al., 2011, EACFM)                        
!-------------------------------------------------------------------------------

subroutine BoundaryMassForceMatrix2D(T,RN,FiS,FiN) 
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
double precision :: FiS, FiN 
double precision,dimension(1:SPACEDIM,1:SPACEDIM) :: T,RN
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
RN(1,1) = FiS * T(1,1) * T(1,1) + FiN * T(3,1) * T(3,1)
RN(1,3) = (FiS - FiN) * T(1,1) * T(3,1) 
RN(3,1) = RN(1,3)         
RN(3,3) = FiS * T(3,1) * T(3,1) + FiN * T(1,1) * T(1,1)
!------------------------
! Deallocations
!------------------------
return
end subroutine BoundaryMassForceMatrix2D

