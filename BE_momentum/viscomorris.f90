!cfile viscomorris.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : viscomorris
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
! Module purpose : Module to 
! 
! Calling routine: inter_EqMoto
!
! Called routines: 
!
!************************************************************************************
!
subroutine viscomorris (npi,npj,npartint,dervel,rvw)
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
integer(4),      intent(IN)  :: npi
integer(4),      intent(IN)  :: npj
integer(4),      intent(IN)  :: npartint
double precision,intent(IN)  :: dervel(3)
double precision,intent(OUT) :: rvw(3)
!
!.. Local Scalars ..
double precision :: amassj,rhotilde,anuitilde,factivis
!
!.. Executable Statements ..
!
 if ( pg(npj)%vel_type /= "std" ) then          !non part fix o altro
   amassj    = pg(npi)%mass
   rhotilde  = pg(npi)%dens
   anuitilde = two * pg(npi)%visc         !CONTROLLARE
 else
!   anuitilde = two * (pg(npi)%visc + pg(npj)%visc)
!   rhotilde  = (pg(npj)%dens + pg(npi)%dens)
   amassj = pg(npj)%mass
   rhotilde  = (pg(npi)%visc * pg(npi)%dens + pg(npj)%visc * pg(npj)%dens + 0.001d0)
   anuitilde = 4.0d0 * (pg(npi)%visc * pg(npj)%visc)      ! ATTENZIONE!! viscosita' cinematica
 end if
!
!.. formula utilizzata nel nov 2005 da rapporto Bon/Gatti/Zuccala/DiMonaco/Gallati
! eta     = 0.01d0 * hh
! gradd   = gradmod / (rij+eta)
! rvw(:)  =-gradd * dervel(:)
!.. formula utilizzata nel nov 2007 come da letteratura
! eta     = 0.01d0 * hh * hh
 factivis = amassj * anuitilde / rhotilde
 
 rvw(1:3) = factivis * &
            (-dervel(1:3) * PartKernel(2,npartint) * (rag(1,npartint) * rag(1,npartint) +  &
             rag(2,npartint) * rag(2,npartint) + rag(3,npartint) * rag(3,npartint)) )
!
return
end subroutine viscomorris
!---split

