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
! Program unit: ParticleCellNumber              
! Description: To return the ID of the grid cell where particle np is located. If particle is outside of the grid, it returns -1 .      
!----------------------------------------------------------------------------------------------------------------------------------

integer(4) function ParticleCellNumber(coord)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(in) :: coord(3)
integer(4) :: i, j, k
double precision :: xp, yp, zp
integer(4),external :: CellNumber
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
xp = coord(1) - Grid%extr(1,1)
yp = coord(2) - Grid%extr(2,1)
zp = coord(3) - Grid%extr(3,1)
i = ceiling(xp / Grid%dcd(1))
k = ceiling(zp / Grid%dcd(3)) 
if (ncord==3) then
   j = ceiling(yp / Grid%dcd(2))
   else
      j = 1
end if
ParticleCellNumber = CellNumber(i, j, k)
!------------------------
! Deallocations
!------------------------
return
end Function ParticleCellNumber

