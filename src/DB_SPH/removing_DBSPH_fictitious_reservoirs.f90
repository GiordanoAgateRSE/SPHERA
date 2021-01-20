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
! Program unit: removing_DBSPH_fictitious_reservoirs
! Description: Removing the fictitious reservoirs used for DB-SPH 
!              initialization (Amicarelli et al., 2013, CAF)    
!-------------------------------------------------------------------------------
subroutine removing_DBSPH_fictitious_reservoirs
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: npi
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
!$omp parallel do default(none)                                                &
!$omp shared(nag,pg,OpCount,Partz)                                             &
!$omp private(npi)
do npi=1,nag
! Fictitious air reservoirs
   if (Partz(pg(npi)%izona)%DBSPH_fictitious_reservoir_flag.eqv.(.true.)) then   
      OpCount(pg(npi)%imed) = OpCount(pg(npi)%imed) + 1    
      pg(npi)%cella = -1
   endif
#ifdef SPACE_3D
! Fictitious fluid reservoir top
   if (Partz(pg(npi)%izona)%IC_source_type==2) then
      if (pg(npi)%coord(3)>Partz(pg(npi)%izona)%H_res) then
         OpCount(pg(npi)%imed) = OpCount(pg(npi)%imed) + 1
         pg(npi)%cella = -1
      endif
   endif
#endif
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine removing_DBSPH_fictitious_reservoirs
