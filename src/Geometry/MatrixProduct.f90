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
! Program unit: MatrixProduct  
! Description: Returning in CC the product between matrices AA and BB.
!              nr: number of rows of AA and CC
!              nc: number of columns of BB and CC
!              nrc: number of columns of AA = number of rows of BB    
!-------------------------------------------------------------------------------
subroutine MatrixProduct(AA,BB,CC,nr,nrc,nc)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(IN) :: nr,nrc,nc
double precision,intent(IN),dimension(nr,nrc) :: AA
double precision,intent(IN),dimension(nrc,nc) :: BB
double precision,intent(INOUT),dimension(nr, nc) :: CC
integer(4) :: i,j,k
double precision :: sum
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
do i=1,nr
   do j=1,nc
      sum = zero
      do k=1,nrc
         sum = sum + AA(i,k) * BB(k,j)
      enddo
      CC(i,j) = sum
   enddo
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine MatrixProduct

