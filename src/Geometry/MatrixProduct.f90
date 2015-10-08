!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-) 
!      
!     
!   
!      
!  

! This file is part of SPHERA.
!  
!  
!  
!  
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
!  
!  
!  
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: MatrixProduct  
! Description: Returning in CC the product between matrices AA and BB.
!              nr: number of rows of AA and CC
!              nc: number of columns of BB and CC
!              nrc: number of columns of AA = number of rows of BB    
!----------------------------------------------------------------------------------------------------------------------------------

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
      end do
      CC(i,j) = sum
   end do
end do
!------------------------
! Deallocations
!------------------------
return
end subroutine MatrixProduct

