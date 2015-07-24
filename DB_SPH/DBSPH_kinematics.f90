!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-; SPHERA has been authored for RSE SpA by 
!    Andrea Amicarelli, Antonio Di Monaco, Sauro Manenti, Elia Bon, Daria Gatti, Giordano Agate, Stefano Falappi, 
!    Barbara Flamini, Roberto Guandalini, David Zuccalà).
! Main numerical developments of SPHERA: 
!    Amicarelli et al. (2015,CAF), Amicarelli et al. (2013,IJNME), Manenti et al. (2012,JHE), Di Monaco et al. (2011,EACFM). 
! Email contact: andrea.amicarelli@rse-web.it

! This file is part of SPHERA.
! SPHERA is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with SPHERA. If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: DBSPH_kinematics
! Description: Imposing input kinematics for the DB-SPH elements (linear interpolation of input data).             
!----------------------------------------------------------------------------------------------------------------------------------

subroutine DBSPH_kinematics
!------------------------
! Modules
!------------------------ 
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: j
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
! Loop over the kinematics records to provide a linear interpolation of the 
! imposed velocity to the first surface element.
if (DBSPH%n_w>0) then
   do j=1,DBSPH%n_kinematics_records
      if (DBSPH%kinematics(j,1)>=tempo) then
         if (DBSPH%kinematics(j,1)==tempo) then
            pg_w(1)%vel(:) = DBSPH%kinematics(j,2:4)
            else
               pg_w(1)%vel(:) = DBSPH%kinematics(j-1,2:4) +                    &
                                (DBSPH%kinematics(j,2:4) -                     &
                                DBSPH%kinematics(j-1,2:4)) /                   &
                                (DBSPH%kinematics(j,1) -                       &
                                DBSPH%kinematics(j-1,1))                       &
                                * (tempo - DBSPH%kinematics(j-1,1))
         endif
         exit
      endif
   enddo
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine DBSPH_kinematics

