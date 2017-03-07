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
! Program unit: AddBoundaryContribution_to_CE2D                                
! Description: To compute boundary terms for the 2D continuity equation 
!              (rodivV). Equation refers to particle npi. It performs implicit 
!              computation of gradPsuro. 
!              (Di Monaco et al., 2011, EACFM)                        
!-------------------------------------------------------------------------------
subroutine AddBoundaryContribution_to_CE2D(npi,IntNcbs,BCrodivV)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
use SA_SPH_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(IN) :: npi,IntNcbs
double precision,intent(INOUT) :: BCrodivV
integer(4) :: pd,icbs,iside,sidestr,ibdt,ibdp
double precision :: IntWds,roi,vin 
integer(4),dimension(1:PLANEDIM)    :: acix
double precision,dimension(1:PLANEDIM) :: IntLocXY,nnlocal,Dvel
type (TyBoundarySide) :: RifBoundarySide
character(4) :: strtype
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
! Active coordinate indexes
acix(1) = 1        
acix(2) = 3
BCrodivV = zero
if (IntNcbs<=0) return
roi = pg(npi)%dens
ibdt = BoundaryDataPointer(3,npi)
!------------------------
! Statements
!------------------------
do icbs=1,IntNcbs
   ibdp = ibdt + icbs - 1
   IntLocXY(1:PLANEDIM) = BoundaryDataTab(ibdp)%LocXYZ(1:PLANEDIM)
   iside = BoundaryDataTab(ibdp)%CloBoNum
   RifBoundarySide = BoundarySide(iside)
   sidestr = RifBoundarySide%stretch
   strtype = Tratto(sidestr)%tipo
   if (strtype=="fixe".OR.strtype=="tapi".OR.strtype=="velo".OR.               &
      strtype=="flow".OR.strtype=="sour") then 
      IntWdS = BoundaryDataTab(ibdp)%BoundaryIntegral(1)
      vin = zero
      do pd=1,PLANEDIM
         nnlocal(pd) = RifBoundarySide%T(acix(pd), acix(2))
      enddo
      select case (strtype)
         case ("fixe")
            do pd=1,PLANEDIM
               Dvel(pd) = pg(npi)%var(acix(pd))
            enddo
            vin = Dvel(1) * nnlocal(1) + Dvel(2) * nnlocal(2)
         case ("tapi")
            do pd=1,PLANEDIM
               Dvel(pd) = pg(npi)%var(acix(pd)) -                              &
                  RifBoundarySide%velocity(acix(pd))
            enddo
            vin = Dvel(1) * nnlocal(1) + Dvel(2) * nnlocal(2)
         case ("velo", "flow", "sour")
            if ((Domain%time_stage==1).or.(Domain%time_split==1)) then 
               pg(npi)%koddens = 2
               pg(npi)%densass = roi
               BCrodivV = zero
            endif
            return
      endselect
! Boundary contribution to the continuity equation 
      BCrodivV = BCrodivV + two * vin * roi * IntWdS
   endif 
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine AddBoundaryContribution_to_CE2D

