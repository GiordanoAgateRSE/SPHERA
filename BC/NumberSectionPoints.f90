!cfile NumberSectionPoints.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
integer(4) function NumberSectionPoints ( values, opt )
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. implicit declarations
  implicit none
!
!.. dummy arguments
double precision,dimension(3,2) :: values
character(1)                    :: opt
!
!.. local Scalars ..
integer(4) :: n
!
!.. local Arrays ..
integer(4),dimension(3) :: Nmesh
!
!.. executable statements
!
!creazione mesh di lato dd
 Nmesh = 1
 do n = 1,SPACEDIM  !Ncord
    if ( n == 1 .AND. opt == "x" ) cycle
    if ( n == 2 .AND. opt == "y" ) cycle
    if ( n == 3 .AND. opt == "z" ) cycle
    Nmesh(n)   = nint ( (values(n,2)-values(n,1)) / Domain%dd )
 end do

 NumberSectionPoints = Product(Nmesh)

return
end function NumberSectionPoints
!---split

