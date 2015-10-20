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
! Program unit: Vector_Product    
! Description: To return in ww the cross product of vectors uu and vv.           
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Vector_Product(uu,VV,ww,SPACEDIM)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(IN) :: SPACEDIM
double precision,intent(IN),dimension(SPACEDIM) :: uu,VV
double precision,intent(INOUT),dimension(SPACEDIM) :: ww
integer(4) :: i,j,k
integer(4),dimension(3) :: iseg=(/2,3,1/)
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
do i=1,SPACEDIM
   j = iseg(i)
   k = iseg(j)
   ww(i) = uu(j) * VV(k) - uu(k) * VV(j)
end do
!------------------------
! Deallocations
!------------------------
return
end subroutine Vector_Product

