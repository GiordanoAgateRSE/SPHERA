!cfile diffumorris.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : diffumorris
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
! Calling routine: inter_EqCont_2d, inter_EqCont_3d
!
! Called routines: 
!
!************************************************************************************
!
subroutine diffumorris (npi,npj,npartint,dervol,factdiff,rvw)
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
 integer(4)       :: npi,npj,npartint
 double precision :: dervol,factdiff,rvw
!
!.. Local Scalars ..
! double precision        :: gradd
! double precision        :: dvar(3),cuei,eta
 double precision        :: anuitilde,rhotilde,amassj,coei,coej
! double precision        :: coeffveli,coeffvelj,dist,distj,disti
! double precision        :: dvi(3),dvj(3)
!
!.. Executable Statements ..
!
! coeffvelj = (Dsqrt(pg(npj)%dv(1)*pg(npj)%dv(1)+pg(npj)%dv(3)*pg(npj)%dv(3)))
! coeffveli = (Dsqrt(pg(npi)%dv(1)*pg(npi)%dv(1)+pg(npi)%dv(3)*pg(npi)%dv(3)))
! coej=two*pg(npj)%coefdif * coeffvelj
! coei=two*pg(npi)%coefdif * coeffveli

 pg(npj)%coefdif = med(pg(npj)%imed)%codif
 pg(npi)%coefdif = med(pg(npi)%imed)%codif
 
 coej = pg(npj)%coefdif
 coei = pg(npi)%coefdif

!! Schema armonico
 amassj    = pg(npj)%mass
 anuitilde = 4.0d0 * (coei * coej) / (coei + coej + 0.00001d0)
 rhotilde  = pg(npj)%dens

 if (pg(npi)%visc == med(2)%numx .or. pg(npj)%visc == med(2)%numx) then
   anuitilde = zero
 end if
 
 if ( pg(npj)%vel_type /= "std" ) then          !non part fix o altro
    rhotilde  = pg(npi)%dens
    amassj    = pg(npi)%mass
    anuitilde = zero
 end if

!.. formula utilizzata nel nov 2005 da rapporto Bon/Gatti/Zuccala/DiMonaco/Gallati
! factdiff = amassj * anuitilde / rhotilde
!! eta      = 0.01d0 * hh
! gradd    = gradmod / (rij+eta)
! rvw      =-gradd * dervol

! disti = pg(npi)%coord(1) * pg(npi)%coord(1) + pg(npi)%coord(3) * pg(npi)%coord(3)
! distj = pg(npj)%coord(1) * pg(npj)%coord(1) + pg(npj)%coord(3) * pg(npj)%coord(3)
! dist  = distj - disti
! cuei  = -two * gradd * dist * amassj / pg(npj)%dens

!.. formula utilizzata nel nov 2007 come da letteratura
 factdiff = amassj * anuitilde / rhotilde
! eta     = 0.01d0 * hh * hh
 rvw     = -dervol * PartKernel(2,npartint) * &
            (rag(1,npartint)*rag(1,npartint) + rag(2,npartint)*rag(2,npartint) + rag(3,npartint)*rag(3,npartint)) 

return
end subroutine diffumorris
!---split

