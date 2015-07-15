!cfile Vector_Product.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : Vector_Product
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!
!************************************************************************************
! Module purpose : Module to Compute and return in ww(1 to SPACEDIM) the components
!                  of the vector product of vectors uu(1 to SPACEDIM) and
!                  vv(1 to SPACEDIM)
!
!AA501b modified
!AA601 sub
! Calling routine: DefineLocalSystemVersors,RHS_body_dynamics,area_triangle
!
! Called routines: 
!
!************************************************************************************
!
subroutine Vector_Product ( uu, VV, ww, SPACEDIM )
!Computes and returns in ww(1 to SPACEDIM) the components of the vector product
!of vectors uu(1 to SPACEDIM) and vv(1 to SPACEDIM)
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
integer(4),      intent(IN)                        :: SPACEDIM
double precision,intent(IN),   dimension(SPACEDIM) :: uu
double precision,intent(IN),   dimension(SPACEDIM) :: VV
double precision,intent(INOUT),dimension(SPACEDIM) :: ww
!
!.. Local Scalars ..
integer(4) :: i, j, k
!
!.. Local Arrays ..
integer(4), dimension(3) :: iseg = (/ 2,3,1 /)
!
!.. Executable Statements ..
!
 do i = 1, SPACEDIM
   j = iseg(i)
   k = iseg(j)
   ww(i) = uu(j) * VV(k) - uu(k) * VV(j)
 end do
!
return
end subroutine Vector_Product
!---split

