!cfile InterFix.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : InterFix
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
! Module purpose : Module to accumulate the contributions of the particles that are
!                  in the sphere of influence of the particle considered
!
! Calling routine: NormFix
!
! Called routines: 
!
!************************************************************************************
!
subroutine InterFix (npi,appo,unity)
!* implementa il meccanismo di ricerca delle particelle che agiscono
!* su quella i-esima.
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
integer(4), parameter :: local_d = 500  ! num max part entro 2h
!
!.. Formal Arguments ..
integer(4),      intent(IN)    :: npi
double precision,intent(INOUT) :: unity
double precision,intent(INOUT),dimension(3) :: appo
!
!.. Local Scalars ..
integer(4) :: npj,contj,npartint   !!!,nfix
double precision :: rhoj,amassj,pesoj    !rhoi,
!
!.. Local Arrays ..
double precision,dimension(3) :: pesogradj
!
!.. Executable Statements ..
!
!* azzeramento quantita generali
 unity   = zero
 appo(:) = zero
!!! nfix  = zero
!
!*_______________________________________________________________
!*prima passata per trovare celle interagenti e memorizzazione
!
 do contj = 1, nPartIntorno(npi)
!
   npartint = (npi-1)* NMAXPARTJ + contj
   npj = PartIntorno(npartint)
!
   if ( pg(npj)%vel_type == "std" ) cycle      !non part fix o altro 
!!!   nfix = nfix + 1  
!
!   rhoi   = pg(npi)%dens
   rhoj   = pg(npj)%dens
   amassj = pg(npj)%mass
!
!* calcolo unita'
   pesoj = amassj * Partkernel(4,npartint) / rhoj
   pesogradj(1:3) = amassj * rag(1:3,npartint) * PartKernel(1,npartint) / rhoj
!
   unity = unity + pesoj  
   appo(:) = appo(:) + pesogradj(:)  
!
 end do
!
 appo(:) = -appo(:)
!
return
end subroutine InterFix 
!---split

