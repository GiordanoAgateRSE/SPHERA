!cfile viscomorris_wall_elements.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : viscomorris_wall_elements
!
! Last updating : Amicarelli/Agate  30 Nov 2011
!
!************************************************************************************
! Module purpose : Wall element contributions to Morris' viscosity term 
! 
! Calling routine: inter_EqMoto
!
! Called routines: /
!
!************************************************************************************
!
subroutine viscomorris_wall_elements(npi,npj,npartint,dervel,rvw)
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
! Declarations
! Implicit
implicit none
! Formal Arguments
integer(4),      intent(IN)  :: npi,npj,npartint
double precision,intent(IN)  :: dervel(3)
double precision,intent(OUT) :: rvw(3)
! Local Scalars
double precision :: rhoWw,rhotilde,anuitilde,factivis,dis2
!
! Statements
 rhoWw = pg_w(npj)%dens * pg_w(npj)%weight * kernel_fw(1,npartint) 
 rhotilde  = (pg(npi)%visc * pg(npi)%dens + pg(npi)%visc * pg_w(npj)%dens + 0.001d0)
 anuitilde = 4.0d0 * (pg(npi)%visc * pg(npi)%visc)      ! ATTENZIONE!! viscosita' cinematica
 factivis = rhoWw * anuitilde / rhotilde
 dis2 = (rag_fw(1,npartint)*rag_fw(1,npartint)+rag_fw(2,npartint)*rag_fw(2,npartint)+rag_fw(3,npartint)*rag_fw(3,npartint))
 rvw(:) = factivis * dervel(:) * (rag_fw(:,npartint)*pg_w(npj)%normal(:)) / dis2
!
return
end subroutine viscomorris_wall_elements
!---split

!cfile viscomon_wall_elements.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : viscomon_wall_elements
!
! Last updating : Amicarelli/Agate  30 Nov 2011
!
!************************************************************************************
! Module purpose : Wall element contributions of Monaghan's viscosity term 
! 
! Calling routine: inter_EqMoto
!
! Called routines: /
!
!************************************************************************************
!
subroutine viscomon_wall_elements(npi,npj,npartint,dervel,rvwalfa,rvwbeta)
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
integer(4),      intent(IN)  :: npi,npj,npartint
double precision,intent(IN)  :: dervel(3)
double precision,intent(OUT) :: rvwalfa(3), rvwbeta(3)
!
!.. Local Scalars ..
double precision :: rhoWw,rhotilde,celtilde,vrij,TermMon,dis2
!
!.. Executable Statements ..
!
 vrij = -dervel(1) * rag_fw(1,npartint) - dervel(2) * rag_fw(2,npartint) - dervel(3) * rag_fw(3,npartint)
 dis2 = (rag_fw(1,npartint)*rag_fw(1,npartint)+rag_fw(2,npartint)*rag_fw(2,npartint)+rag_fw(3,npartint)*rag_fw(3,npartint))
 if (vrij > zero) then
   rvwalfa = zero
   rvwbeta = zero
   else
      rhotilde = pg(npi)%dens + pg_w(npj)%dens
      celtilde = 2. * Med(pg(npi)%imed)%celerita 
      rhoWw = pg_w(npj)%dens * pg_w(npj)%weight * kernel_fw(1,npartint) 
      TermMon = Med(pg(npi)%imed)%alfaMon * celtilde * Domain%h / rhotilde
      rvwalfa(1:3) = rhoWw * TermMon * vrij * pg_w(npj)%normal(:) / dis2
      rvwbeta(1:3) = zero
 end if
!
return
end subroutine viscomon_wall_elements
!---split

