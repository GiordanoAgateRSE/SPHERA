!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)

! SPHERA authors and email contact are provided on SPHERA documentation.

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
! Program unit: BoundaryMassForceMatrix2D                                
! Description: Generation of the generalized boundary mass force matrix RN, on 
!              the base of the cosine matrix T and the parameter Fi. (Di Monaco 
!              et al., 2011, EACFM)
!-------------------------------------------------------------------------------
#ifdef SPACE_2D
subroutine BoundaryMassForceMatrix2D(FiS,FiN,TT,RN) 
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(in) :: FiS,FiN 
double precision,dimension(1:SPACEDIM,1:SPACEDIM),intent(in) :: TT
double precision,dimension(1:SPACEDIM,1:SPACEDIM),intent(inout) :: RN
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
RN(1,1) = FiS * TT(1,1) * TT(1,1) + FiN * TT(3,1) * TT(3,1)
RN(1,3) = (FiS - FiN) * TT(1,1) * TT(3,1) 
RN(3,1) = RN(1,3)         
RN(3,3) = FiS * TT(3,1) * TT(3,1) + FiN * TT(1,1) * TT(1,1)
!------------------------
! Deallocations
!------------------------
return
end subroutine BoundaryMassForceMatrix2D
#endif
