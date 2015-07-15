!cfile viscomon.f90

!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : viscomon
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!AA504
! 03  Amicarelli        08Apr14        (v5.04) Modifications for Monaghan's viscosisty (active also for separating elements) and 
!                                      neglection of the molecular viscosity depending on the velocity divergence (momentum equation)
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
subroutine viscomon (npi,npj,npartint,dervel,rvwalfa,rvwbeta)

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
double precision,intent(OUT) :: rvwalfa(3), rvwbeta(3)
!
!.. Local Scalars ..
double precision :: celtilde,rhotilde,amassj,vrij,TermMon
!
!.. Executable Statements ..
!
 vrij = -dervel(1) * rag(1,npartint) - dervel(2) * rag(2,npartint) - dervel(3) * rag(3,npartint)
!
!AA504 removed part: Monaghan's term is always active even if particle detach themselves 
   if ( pg(npj)%vel_type /= "std" ) then          !non part fix o altro
     amassj   = pg(npi)%mass
     rhotilde = two * pg(npi)%dens
     celtilde = Med(pg(npi)%imed)%celerita + Med(pg(npi)%imed)%celerita
   else
     amassj   = pg(npj)%mass
     rhotilde = pg(npi)%dens + pg(npj)%dens
!........................................... 2011 mar 15
     if (esplosione) then  !introdotto per stabilita' flushing
       celtilde = pg(npi)%Csound + pg(npj)%Csound
     else
       celtilde = Med(pg(npi)%imed)%celerita + Med(pg(npj)%imed)%celerita
     end if
!........................................... 2011 mar 15
   end if
!
!.. formula utilizzata nel nov 2007 come da letteratura
   TermMon = Med(pg(npi)%imed)%alfaMon * celtilde * Domain%h / rhotilde

!AA504 Molecular viscosity term depending on velocity divergence (compressible flows) can be neglected (and may causes several problems when activated)
   rvwalfa(1:3) = amassj * TermMon * &
                  vrij * rag(1:3,npartint) * PartKernel(2,npartint)

!.. rvwbeta e' da verificare
!   rvwbeta(1:3)= - amassj * two * Med(pg(npi)%imed)%betaMon * vrij*vrij * squareh / rhotilde & 
!                 * (PartKernel(2,npartint)/PartKernel(1,npartint)*(PartKernel(2,npartint)/PartKernel(1,npartint)))
   rvwbeta(1:3) = zero

!
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!if (it_corrente == 2 .and. npi > 16900 .and. npi < 17200) then
!  write (99,*) ' it = ',it_corrente,' n. particella = ',npi,npj
!  write (99,'(a,4d23.15)') ' viscomon   vrij,amassj,rhotilde,celtilde ',vrij,amassj,rhotilde,celtilde
!  write (99,'(a,3d23.15)') ' viscomon   rvwalfa ',rvwalfa
!end if
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!
return
end subroutine viscomon
!---split

