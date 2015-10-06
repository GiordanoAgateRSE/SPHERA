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
! Program unit: IWro2dro                                        
! Description: Computes the definite integral
!2
!S W*(ro')ro2' dro'
!ro
!               (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

double precision function IWro2dro(ro)
!------------------------
! Modules
!------------------------ 
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
double precision,parameter :: a1 = 0.333333333333333d0 
double precision,parameter :: a2 = 0.3d0               
double precision,parameter :: a3 = 0.125d0             
double precision,parameter :: a4 = 0.25d0
double precision,parameter :: b1 = 0.0625d0             != 1 / 16
double precision,parameter :: b2 = 0.025d0              != 1 / 40
double precision,parameter :: b3 = 4.16666666666667d-03 != 1 / 240
double precision,intent(IN) :: ro
double precision :: ro2,ro3,duemro,duemro2,duemro4,IWro2dro1,IWro2dro2
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
if ((ro>=zero).and.(ro<one)) then
   ro2 = ro * ro
   ro3 = ro2 * ro
   IWro2dro1 = (a1 - a2 * ro2 + a3 * ro3) * ro3
   IWro2dro = KERNELCONST3D * (a4 - IWro2dro1)
   elseif ((ro>=one).and.(ro<two)) then
      ro2 = ro * ro
      duemro = two - ro
     duemro2 = duemro * duemro
     duemro4 = duemro2 * duemro2
     IWro2dro2 = -(b1 * ro2 + b2 * duemro * ro + b3 * duemro2) * duemro4
     IWro2dro = KERNELCONST3D * (-IWro2dro2)
     else
        IWro2dro = zero
endif
!------------------------
! Deallocations
!------------------------
return
end function IWro2dro

