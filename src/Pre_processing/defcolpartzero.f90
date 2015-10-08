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
! Program unit: defcolpartzero                   
! Description:  On the particle colours for visualization purposes.               
!----------------------------------------------------------------------------------------------------------------------------------

subroutine defcolpartzero (ir,partz,pg)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(IN) :: ir
type (TyZone),intent(IN),dimension(NPartZone) :: partz
type (TyParticle),intent(INOUT) :: pg
integer(4) :: nbande, numbanda
double precision :: aldx
integer(4),dimension(5) :: iclnumb
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
iclnumb(1)=1
iclnumb(2)=2
iclnumb(3)=4
iclnumb(4)=5
iclnumb(5)=6
!------------------------
! Statements
!------------------------
if (partz(ir)%bend=="u") then        
! Uniform color 
   pg%icol = partz(ir)%icol
   elseif (partz(ir)%bend=="o")then   
! Colour based on external option 
      pg%icol = partz(ir)%icol
      elseif(partz(ir)%bend=="b") then
! Vertical strips 
         nbande = partz(ir)%icol
         aldx = (partz(ir)%coordMM(1,2) - partz(ir)%coordMM(1,1)) / nbande
         numbanda = int((pg%coord(1) - partz(ir)%coordMM(1,1)) / aldx) + 1
         numbanda = min(nbande,numbanda)
         numbanda = max(0,numbanda)
         pg%icol = iclnumb(numbanda)
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine defcolpartzero

