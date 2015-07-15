!cfile CreateSectionPoints.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine  CreateSectionPoints ( vp, values, opt, seccor )
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. implicit declarations
  implicit none
!
!.. dummy arguments
integer(4)   :: seccor
character(1) :: opt
double precision, dimension(3)   :: vp
double precision, dimension(3,2) :: values
!
!.. local scalars
integer(4)   :: i,j,k,n, npse
integer(4),      dimension(3)   :: Nmesh
double precision,dimension(3)   :: Cc,CcStart

character(80), external :: lcase
!
!.. executable statements
!
!  creazione mesh di lato dd
 Nmesh = 1
 npse = Control_Sections(seccor)%icont(1)
!
 CcStart(:) = vp(:) - Domain%dd
 do n = 1,SPACEDIM
   if ( n == 1 .AND. lcase(opt) == "x" ) cycle
   if ( n == 2 .AND. lcase(opt) == "y" ) cycle
   if ( n == 3 .AND. lcase(opt) == "z" ) cycle
   Nmesh(n)   = nint ( (values(n,2)-values(n,1)) / Domain%dd )
   CcStart(n) = values(n,1) - Domain%dd * half
 end do

 Cc(1) = CcStart(1)

 do i = 1, Nmesh(1)
   Cc(1) = Cc(1) + Domain%dd
   Cc(2) = CcStart(2)
   do j = 1, Nmesh(2)
     Cc(2) = Cc(2) + Domain%dd
     Cc(3) = CcStart(3)
     do k = 1, Nmesh(3)
       Cc(3) = Cc(3) + Domain%dd
       npse = npse + 1
       Control_Points(npse)%coord(:) = Cc(:)
     end do
   end do
 end do
 
return
end subroutine  CreateSectionPoints
!---split

