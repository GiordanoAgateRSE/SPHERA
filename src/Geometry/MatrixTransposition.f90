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
! Program unit: MatrixTransposition   
! Description: Returns in AAT(n,m) the transposed matrix of AA(m, n).     
!----------------------------------------------------------------------------------------------------------------------------------

subroutine MatrixTransposition(AA,AAT,m,n) 
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(IN) :: m,n
double precision,intent(IN),dimension(m,n) :: AA
double precision,intent(INOUT),dimension(n,m) :: AAT
integer(4) :: i,j
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
do i=1,n
    do j=1,m
        AAT(i,j) = AA(j,i)
    end do
end do
!------------------------
! Deallocations
!------------------------
return
end subroutine MatrixTransposition

