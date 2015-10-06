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
! Program unit: CalcVarp              
! Description: To calculate physical quantities at a monitoring point.     
!----------------------------------------------------------------------------------------------------------------------------------

subroutine CalcVarp 
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
integer(4) :: i,ii,jj,kk,pointcellnumber
double precision xp,yp,zp
type (TyCtlPoint) :: pglocal
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
if (Npointst<1) return
do i=1, Npointst
   pglocal%coord(:) = control_points(i)%coord(:)
   pglocal%pres = zero
   pglocal%dens = zero
   pglocal%vel(:) = zero
   pglocal%uni = zero
   xp = pglocal%coord(1) - Grid%extr(1,1)
   yp = pglocal%coord(2) - Grid%extr(2,1)
   zp = pglocal%coord(3) - Grid%extr(3,1)
   ii = ceiling(xp / Grid%dcd(1))
   jj = ceiling(yp / Grid%dcd(2))
   kk = ceiling(zp / Grid%dcd(3)) 
   pointcellnumber = CellNumber(ii, jj, kk)
   pglocal%cella = pointcellnumber
   call GetVarPart (pglocal)
   control_points(i)%pres = pglocal%pres
   control_points(i)%dens = pglocal%dens
   control_points(i)%vel(:) = pglocal%vel(:)
   control_points(i)%uni = pglocal%uni
   control_points(i)%cella  = pglocal%cella
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine CalcVarp

