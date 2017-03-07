!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
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
! Program unit: BoundaryPressureGradientMatrix3D                                
! Description: To generate the pressure gradient matrix RRP, based on the cosine
!              matrix T and the parameter vector Psi.
!              (Di Monaco et al., 2011, EACFM)                        
!-------------------------------------------------------------------------------
subroutine BoundaryPressureGradientMatrix3D(T,RGP,Psi)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(INOUT),dimension(SPACEDIM)          :: Psi
double precision,intent(INOUT),dimension(SPACEDIM,SPACEDIM) :: T,RGP
Integer(4) :: i
double precision,dimension(SPACEDIM,SPACEDIM) :: Diag,PsiR,TTR
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
Diag = zero
!------------------------
! Statements
!------------------------
do i=1,SPACEDIM
   Diag(i,i) = Psi(i)
enddo
call MatrixTransposition(T,TTR,SPACEDIM,SPACEDIM)
call MatrixProduct(Diag,TTR,PsiR,SPACEDIM,SPACEDIM,SPACEDIM)
call MatrixProduct(T,PsiR,RGP,SPACEDIM,SPACEDIM,SPACEDIM)
!------------------------
! Deallocations
!------------------------
return
end subroutine BoundaryPressureGradientMatrix3D

