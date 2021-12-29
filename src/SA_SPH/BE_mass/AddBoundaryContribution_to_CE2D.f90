!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2021 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.9.0.0
! SPHERA v.9.0.0 is free software: you can redistribute it and/or modify
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
!              (rodivV) of a generic "npi" SPH particle. It performs implicit 
!              computation of "gradPsuro". In case of a neighbouring inlet 
!              section, the particle density and pressure do not change.
!              (Di Monaco et al., 2011, EACFM)                       
!-------------------------------------------------------------------------------
#ifdef SPACE_2D
subroutine AddBoundaryContribution_to_CE2D(npi,IntNcbs,BCrodivV)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use Memory_I_O_interface_module
use SA_SPH_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: npi,IntNcbs
double precision,intent(inout) :: BCrodivV
integer(4) :: pd,icbs,iside,sidestr,ibdt,ibdp
double precision :: IntWds,roi,vin 
integer(4),dimension(1:PLANEDIM) :: acix
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
BCrodivV = zero
if (IntNcbs<=0) return
! Active coordinate indices
acix(1) = 1        
acix(2) = 3
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
   if (strtype=="fixe".or.strtype=="tapi".or.strtype=="velo".or.               &
      strtype=="flow".or.strtype=="sour") then 
      IntWdS = BoundaryDataTab(ibdp)%BoundaryIntegral(1)
      vin = zero
      do pd=1,PLANEDIM
         nnlocal(pd) = RifBoundarySide%T(acix(pd),acix(2))
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
#endif
