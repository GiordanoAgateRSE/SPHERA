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
! Program unit: InterpolateBoundaryIntegrals2D                                       
! Description:  Interpolation in table "BoundIntegralTab(:,:)", defined in module "SA_SPH_module", the values in 
!               columns "Colmn(nc), nc=1, Ncols" corresponding to the input value "x" to be interpolated, in turn, in column 0. 
!               It returns:
!                  Func(nc), nc=1, Ncols : values interpolated in columns Col(nc), nc=1, Ncols
!               (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine InterpolateBoundaryIntegrals2D(x,Ncols,Colmn,Func)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
use SA_SPH_module
!------------------------
! Declarations
!------------------------
implicit none
double precision :: x 
integer(4) :: Ncols
integer(4),dimension(1:NUMCOLS_BIT) :: Colmn
double precision,dimension(1:NUMCOLS_BIT) :: Func
integer(4) :: nc,i,j
double precision :: xi,fi,fip1
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
if (x<BoundIntegralTab2D(1,0)) x = BoundIntegralTab2D(1,0)
i = Int((x-BoundIntegralTab2D(1,0))/DELTAX_BIT)+1
xi = BoundIntegralTab2D(i,0)
do nc=1,Ncols
   j = Colmn(nc)
   fi = BoundIntegralTab2D(i,j)
   fip1 = BoundIntegralTab2D(i+1,j)
   Func(nc) = fi+(fip1-fi)*(x-xi)/DELTAX_BIT
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine InterpolateBoundaryIntegrals2D

