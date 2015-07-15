!AA601 the whole subroutine
!cfile BC_wall_elements.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : DBSPH_inlet_outlet
!
! Creation      : Amicarelli, 26Jan15
!
!************************************************************************************
! Module purpose : Impose boundary conditions at the inlet and outlet sections (DBSPH) 
!
! Calling routines: inter_SmoothVelo_2D,inter_SmoothVelo_3D,Euler
!
! Called routines: /
!
!************************************************************************************
subroutine DBSPH_inlet_outlet(npi)

! Modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE

! Declarations 
implicit none
integer(4),intent(in)           :: npi

!Initializations

!Statements
! By hypothesis a fluid particle cannot be close to two or more inlet/outlet sections (otherwise it randomly takes an inlet velocity value)
! Inlet boundary conditions
if (pg(npi)%DBSPH_inlet_ID>0) pg(npi)%vel(:) = pg_w(pg(npi)%DBSPH_inlet_ID)%vel(:)
! Outlet boundary conditions
if (pg(npi)%DBSPH_outlet_ID>0) then
   pg(npi)%dens = Med(1)%den0
   pg(npi)%pres = pg_w(pg(npi)%DBSPH_outlet_ID)%pres
endif

!Deallocations

return
end subroutine DBSPH_inlet_outlet
!---split

