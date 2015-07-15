!cfile defcolpartzero.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : defcolpartzero
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
! Module purpose : Module for particle color definition as from input file
!
! Calling routine: GenerateSourceParticles_2D, GenerateSourceParticles_3D
!                  SetParticles
!
! Called routines: 
!
!************************************************************************************
!
subroutine defcolpartzero (ir,partz,pg)
!* definisce colore particella nello stato zero
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
integer(4),       intent(IN)                         :: ir
type (TyZone),    intent(IN),   dimension(NPartZone) :: partz
type (TyParticle),intent(INOUT)                      :: pg
!
!.. Local Scalars ..
integer(4)       :: nbande, numbanda
double precision :: aldx
!
!.. Local Arrays ..
integer(4), dimension(5) :: iclnumb
!
!.. Executable Statements ..
!
 iclnumb(1)=1
 iclnumb(2)=2
 iclnumb(3)=4
 iclnumb(4)=5
 iclnumb(5)=6

 if ( partz(ir)%bend == "u")then        ! color uniform
    pg%icol = partz(ir)%icol

 else if ( partz(ir)%bend == "o")then   ! color optional (solo su opzione esterna)
    pg%icol = partz(ir)%icol

 else if(partz(ir)%bend == "b")then     ! bande verticali
    nbande = partz(ir)%icol
    aldx   =( partz(ir)%coordMM(1,2)-partz(ir)%coordMM(1,1))/nbande
    numbanda = int((pg%coord(1)-partz(ir)%coordMM(1,1))/aldx)+1
    numbanda = min(nbande,numbanda)
    numbanda = max(0,numbanda)
    pg%icol = iclnumb(numbanda)
 end if

return
end subroutine defcolpartzero
!---split

