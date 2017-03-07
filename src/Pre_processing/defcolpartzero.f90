!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.8.0.
! SPHERA v.8.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with SPHERA. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: defcolpartzero                   
! Description: on the particle colours for visualization purposes.               
!-------------------------------------------------------------------------------
subroutine defcolpartzero(ir,partz,pg)
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

