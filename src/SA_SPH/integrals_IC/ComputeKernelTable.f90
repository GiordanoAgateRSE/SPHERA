!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2020 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: ComputeKernelTable                                   
! Description: To compute and save in kerneltab(0:ktrows,0:ktcols) the 
!              following values (Di Monaco et al., 201, EACFM):
!                 kerneltab(0:ktrows,0) = rob = r_0b/h = q_0b
!                 kerneltab(0:ktrows,1) = Integral (W_norm * q**2 * dq) 
!                                         (from q_0b to 2)
!                 kerneltab(0:ktrows,2) = Integral (dW_norm/dq * q * dq) 
!                                         (from q_0b to 2)
!                 kerneltab(0:ktrows,3) = Integral (dW_norm/dq * q**2 * dq) 
!                                         (from q_0b to 2)
!                 kerneltab(0:ktrows,4) = Integral (dW/dq * q**3 * dq) 
!                                         (from q_0b to 2)
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine ComputeKernelTable
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nr
double precision :: rob
double precision,external :: IWro2dro,JdWsRn
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
ktdelta = two / INT_KERNELTABLE
!------------------------
! Statements
!------------------------
do nr=0,ktrows
  rob = ktdelta * nr
  kerneltab(nr,0) = rob
  kerneltab(nr,1) = IWro2dro(rob)
  kerneltab(nr,2) = JdWsRn(rob,1,1) * Unosusquareh
  kerneltab(nr,3) = JdWsRn(rob,2,1) * Unosuh
  kerneltab(nr,4) = JdWsRn(rob,3,1)
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ComputeKernelTable
#endif
