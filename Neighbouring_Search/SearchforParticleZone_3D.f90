!cfile SearchforParticleZone_3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : SearchforParticleZone_3D
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
! Module purpose : Module to search the highest pointer of the fluid
!
! Calling routine: PreSourceParticles_3D, Gest_Input
!
! Called routines: LocalNormalCoordinates
!
!************************************************************************************
!
  subroutine SearchforParticleZone_3D (partizone)   

!Restituisce in 'partizone' l'indice di zona più alto cui corrisponde un volume fluido (di particelle)
!Nel caso non esista una tale volume (zona) fluido assegna partzone = sourzone (zona dellasorgente)
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
  integer(4),intent(INOUT)   :: partizone
!
!.. Local Scalars ..
  integer(4)  :: iz, sourzone, mate
  character(4) :: tipo
!
!.. Executable Statements ..
!
  partizone = 0
  sourzone = 0
!
!.. 
  do iz = NPartZone, 1, -1
    tipo = Partz(iz)%tipo
    if (tipo /= "sour") then
      mate = Partz(iz)%Medium
      if (mate > 0) then
        partizone = iz
        exit
      end if
    else
      sourzone = iz
    end if
  end do
  if (partizone == 0) partizone = sourzone

  return
  end subroutine SearchforParticleZone_3D
!---split

