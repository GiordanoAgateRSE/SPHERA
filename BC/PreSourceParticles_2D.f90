!cfile PreSourceParticles_2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : PreSourceParticles_2D
!
! Last updating : November 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Agate             2011           varie
! 04  Amicarelli        23/11/2011     multiple inlet
!
!************************************************************************************
! Module purpose : Module to generate new source particles to simulate inlet fluid flow
!                  !!! ONLY in 2D and with one inlet !!!
!
! Calling routine: Loop_Irre_2D
!
! Called routines: 
!
!************************************************************************************
!
subroutine PreSourceParticles_2D
!Genera e immette nel campo nuove particelle di sorgente per simulare una portata fluida entrante
!Attenzione: questa routine funziona solo in presenza di un solo lato di contorno funzionante
!            da sorgente e solo nel caso bidimensionale!!!!
!
!.. assign modules
use GLOBAL_MODULE
use FILES_ENTITIES
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Scalars ..
!
!AA405 sub
integer(4)       :: nt, nA, isi, sd, ip, i_source
!
!double precision :: Time
double precision :: deltapart, linedist, sidelen, eps
!
!.. Local Arrays ..
double precision,dimension(1:SPACEDIM) :: A
double precision,dimension(1:SPACEDIM) :: ss
!
integer(4), external :: ParticleCellNumber
!
!.. Executable Statements ..
!
!Time = tempo
!
! Ricerca del numero del lato sorgente
  SourceSide = 0
  SpCount = 0
!
!AA405 
  i_source=0
!
  do isi = 1, NumBSides
    if ( BoundarySide(isi)%tipo == "sour" ) then
      SourceSide = isi
!
!AA405
      i_source=i_source+1
!
!AA405 rm start
!      exit
!    end if
!  end do
!  if ( SourceSide > 0 ) then
!AA405 rm end
!
    nt = BoundarySide(SourceSide)%stretch
    irz = Tratto(nt)%zone
    mat = partz(irz)%Medium 
    nA = BoundarySide(SourceSide)%Vertex(1)
    do sd = 1, SPACEDIM
      A(sd)  = Vertice(sd, nA)
      ss(sd) = BoundarySide(SourceSide)%T(sd, 1)
      nn(sd) = BoundarySide(SourceSide)%T(sd, 3)
    end do
    deltapart = Domain%dd
    sidelen = BoundarySide(SourceSide)%length
!
!AA405 sub
    NumPartperLine(i_source) = Int(sidelen / deltapart + 0.01d0)
!
    eps = -half
    yfila = eps * deltapart 
    linedist = -half * deltapart
!
!AA405 sub
    do ip = 1, NumPartperLine(i_source)
!
      linedist = linedist + deltapart
      do sd = 1, SPACEDIM
!AA405 sub
        PartLine(i_source, ip, sd) = A(sd) + linedist * ss(sd)
!
      end do
    end do
!    ParticleVolume = deltapart * deltapart
    ParticleVolume = Domain%PVolume
!
!AA405 sub
    RowPeriod = ParticleVolume * NumPartperLine(i_source) / Tratto(nt)%FlowRate 
!
!!    RowVelocity = partz(irz)%vel(1) * nn(1) + partz(irz)%vel(3) * nn(3)
!    RowVelocity = Tratto(nt)%NormVelocity
!AA405 sub
    RowVelocity(i_source) = Domain%dd / RowPeriod
!
!AA405 sub start
    Tratto(nt)%NormVelocity = RowVelocity(i_source)
    partz(irz)%vel(1) = RowVelocity(i_source) * nn(1)
    partz(irz)%vel(2) = RowVelocity(i_source) * nn(2)
    partz(irz)%vel(3) = RowVelocity(i_source) * nn(3)
!    if (RowVelocity == zero) exit  
!AA405 sub end
!
!    if (tempo < Tratto(nt)%trampa) RowVelocity = RowVelocity * tempo / tratto(nt)%trampa
!    RowPeriod = deltapart/ Abs(RowVelocity) 
    pinttimeratio = -1
!
!AA405 rm
!  end if
!
!AA405 start
  end if ! BoundarySide(isi)%tipo == "sour"
  end do ! isi
!AA405 end
!
return
end subroutine PreSourceParticles_2D
!---split

