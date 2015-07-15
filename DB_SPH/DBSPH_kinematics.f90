!AA601 the whole subroutine
!cfile BC_wall_elements.f90
!************************************************************************************
!                          S P H E R A 5.0.3 
!
!                Smoothed Particle Hydrodinamics Code
!
!************************************************************************************
!
! File name     : DBSPH_kinematics
!
! Creation      : Amicarelli A., 26Jan16
!
!************************************************************************************
! Module purpose : Imposed kinematics for the DBSPH frontier system (linear interpolation of input data)
!                  (method derived from time_integration_body_dynamics)
!
! Calling routine: Euler
!
! Called routines: /
!
!************************************************************************************

subroutine DBSPH_kinematics

! Modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE

! Declarations
implicit none
integer(4) :: j

! Interface blocks

! Allocations

! Initializations

! Statements
! Loop over kinematics records to provide linear interpolation of imposed velocity to the first surface element
if (DBSPH%n_w>0) then
do j=1,DBSPH%n_kinematics_records
   if (DBSPH%kinematics(j,1)>=tempo) then
      if (DBSPH%kinematics(j,1)==tempo) then
         pg_w(1)%vel(:) = DBSPH%kinematics(j,2:4)
         else
            pg_w(1)%vel(:) = DBSPH%kinematics(j-1,2:4) + (DBSPH%kinematics(j,2:4)-DBSPH%kinematics(j-1,2:4))/ &
                             (DBSPH%kinematics(j,1)-DBSPH%kinematics(j-1,1)) * (tempo-DBSPH%kinematics(j-1,1))
      endif
      exit
   endif
enddo
endif

! Deallocations

return
end subroutine DBSPH_kinematics
!---split

