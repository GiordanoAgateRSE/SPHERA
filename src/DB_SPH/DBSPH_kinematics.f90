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

