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
! Program unit: WIntegr                                          
! Description: Computing the definite integral  
!2h
!S W(r,h)rdr
!ri
!               (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

double precision function WIntegr(ri,h) 
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
double precision,parameter :: a   =  0.5d0
double precision,parameter :: a2  = -0.375d0
double precision,parameter :: a3  =  0.15d0
double precision,parameter :: aa  = -0.125d0
double precision,parameter :: aa1 =  0.05d0
double precision,parameter :: b   =  0.35d0
double precision,intent(IN) :: ri,h
double precision :: WIntegr1,WIntegr2,ro,ro2,ro3,q1,q2,q4
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
! KERNELCONST2D = 0.454728408833987d0 =10/(7*pigreco), as defined in the modules
ro = ri / h
!------------------------
! Statements
!------------------------
if ((zero<=ro).and.(ro<one)) then
   ro2 = ro * ro
   ro3 = ro2 * ro
   WIntegr1 = (a + a2 * ro2 + a3 * ro3) * ro2
   WIntegr = KERNELCONST2D * (b - WIntegr1)
   elseif ((one<=ro).and.(ro<two)) then
      q1 = two - ro
      q2 = q1 * q1
      q4 = q2 * q2
      WIntegr2 = (aa + aa1 * q1) * q4
      WIntegr = -KERNELCONST2D * WIntegr2
      else
         WIntegr = zero
endif
!------------------------
! Deallocations
!------------------------
return
end function WIntegr

