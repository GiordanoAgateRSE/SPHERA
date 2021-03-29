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
! Program unit: SASPH_continuity
! Description: SASPH boundary term for the continuity equation
!-------------------------------------------------------------------------------
subroutine SASPH_continuity(npi                                                &
#ifdef SPACE_3D
                               ,Ncbf_Max                                       &
#endif
                            )
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Dynamic_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: npi
#ifdef SPACE_3D
integer(4),intent(inout) :: Ncbf_Max
#endif
#ifdef SPACE_3D
integer(4) :: Ncbf
double precision :: BCtorodivV
#elif defined SPACE_2D
integer(4) :: Ncbs,IntNcbs
double precision :: BCrodivV
#endif
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
! Seaching for the neighbouring faces/sides of the particle "npi"
#ifdef SPACE_3D
   BCtorodivV = zero
   Ncbf = BoundaryDataPointer(1,npi)
! Detecting the faces with actual contributions
   if (Ncbf>0) then
!$omp critical (omp_Ncbf_Max_2)
      Ncbf_Max = max(Ncbf_Max,Ncbf)
!$omp end critical (omp_Ncbf_Max_2)
      call AddBoundaryContribution_to_CE3D(npi,Ncbf,BCtorodivV)
   endif
#elif defined SPACE_2D
      BCrodivV = zero
      Ncbs = BoundaryDataPointer(1,npi)
      IntNcbs = BoundaryDataPointer(2,npi)
! Detecting the sides with actual contributions
      if ((Ncbs>0).and.(IntNcbs>0)) then
         call AddBoundaryContribution_to_CE2D(npi,IntNcbs,BCrodivV)
      endif
#endif
if (pg(npi)%koddens==0) then
#ifdef SPACE_3D
      pg(npi)%dden = pg(npi)%dden - BCtorodivV
#elif defined SPACE_2D
         pg(npi)%dden = pg(npi)%dden - BCrodivV
#endif
   else
      pg(npi)%dden = zero
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine SASPH_continuity
