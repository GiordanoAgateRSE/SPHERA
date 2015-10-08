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
integer(4) :: j,i_a
double precision :: vel_aux(3),omega_aux(3),rel_pos(3)
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
            vel_aux(:) = DBSPH%kinematics(j,2:4)
            omega_aux(:) = DBSPH%kinematics(j,5:7)
            else
               vel_aux(:) = DBSPH%kinematics(j-1,2:4) +                        &
                            (DBSPH%kinematics(j,2:4) -                         &
                            DBSPH%kinematics(j-1,2:4)) /                       &
                            (DBSPH%kinematics(j,1) -                           &
                            DBSPH%kinematics(j-1,1)) * (tempo -                &
                            DBSPH%kinematics(j-1,1)) 
               omega_aux(:) = DBSPH%kinematics(j-1,5:7) +                      &
                              (DBSPH%kinematics(j,5:7) -                       &
                              DBSPH%kinematics(j-1,5:7)) /                     &
                              (DBSPH%kinematics(j,1) -                         &
                              DBSPH%kinematics(j-1,1)) * (tempo -              &
                              DBSPH%kinematics(j-1,1))                  
         endif
!$omp parallel do default(none)                                                &
!$omp shared(DBSPH,pg_w)                                                       &
!$omp private(i_a,rel_pos,vel_aux,omega_aux)
         do i_a=1,DBSPH%n_w
            rel_pos = pg_w(i_a)%coord(:) - DBSPH%rotation_centre(:)          
            call Vector_Product(omega_aux,rel_pos,pg_w(i_a)%vel(:),3)
            pg_w(i_a)%vel(:) = pg_w(i_a)%vel(:) + vel_aux(:)
         enddo
!$omp end parallel do
         exit
      endif
   enddo
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine DBSPH_kinematics

