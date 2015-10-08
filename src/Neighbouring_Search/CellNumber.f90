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
! Program unit: CellNumber             
! Description: To return the ID of the cell of indices  (i,j,k).  
!----------------------------------------------------------------------------------------------------------------------------------

integer(4) function CellNumber(i,j,k)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(IN) :: i,j,k
integer(4) :: ni,nj,nk
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
ni = Grid%ncd(1)
nj = Grid%ncd(2)
nk = Grid%ncd(3)
!------------------------
! Statements
!------------------------
if ((i<1).or.(i>ni).or.(j<1).or.(j>nj).or.(k<1).or.(k>nk)) then
! .. the cell is outside the grid limits
   CellNumber = 0  
   else
! .. return the cell number
      CellNumber = ((k - 1) * nj + (j - 1)) * ni + i
end if
!------------------------
! Deallocations
!------------------------
end function CellNumber

