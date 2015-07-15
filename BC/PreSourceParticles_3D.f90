!cfile PreSourceParticles_3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : PreSourceParticles_3D
!
! Last updating : November 23, 2011
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
!                  !!! ONLY in 3D and with one inlet face with four nodes !!!
!
! Calling routine: Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
!
subroutine PreSourceParticles_3D 
!Genera e immette nel campo nuove particelle di sorgente per simulare una portata fluida entrante
!Attenzione: questa routine funziona solo in presenza di una sola faccia di contorno fuzionante
!            da sorgente e solo nel caso che questa sia un parallelogramma (4 nodi)
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
!.. Local Parameters ..
!integer(4), parameter :: MAXPARTLINE = 2000
!
!.. Local Scalars ..
!
!AA405 sub
integer(4) :: nt, isi, sd, ip, i_source    !nA, 
!
integer(4) :: i, j, NumPartR, NumPartS, nodes
!double precision :: Time
double precision :: deltapart, eps
double precision :: deltaR, deltaS, LenR, LenS, distR, distS, csi, etalocal
!
!.. Executable Statements ..
!
!Time = tempo
!
! Ricerca del numero del lato sorgente
  SourceFace = 0
  SpCount = 0
!
!AA405 
  i_source=0
!
  do isi = 1, NumFacce
    nt = BoundaryFace(isi)%stretch
    if ( Tratto(nt)%tipo == "sour" ) then
      SourceFace = isi
!
!AA405
      i_source=i_source+1
!
!AA405 rm start
!      exit
!    end if
!  end do
!AA405 rm end
!
  if ( SourceFace > 0 ) then
  
!Cerca l'indice di zona 'partzone' da assegnare alle particelle generate
!Tale indice sara' il più alto tra quelli che identificano volumi di particelle presenti nel campo.
!Alle particelle generate dalla sorgente sara' assegnato (di ufficio) un indice di
!materiale uguale a quello ('partzone') della zona trovata
!Se una tale zona non è stata assegnata (casi con sola sorgente) allora l'indice di zona
!assegnato alle particelle generate sara' quello assegnato in input ('izone').

!    call SearchforParticleZone_3D(partzone)
  
    nt = BoundaryFace(SourceFace)%stretch
    mat = Tratto(nt)%Medium
    izone = Tratto(nt)%zone
    nodes = BoundaryFace(SourceFace)%nodes   !Nota: Inserire un controllo con mess box nel caso di nodes =/ 4 
    deltapart = Domain%dd
!    RowVelocity = zero
!    nA = BoundaryFace(SourceFace)%Node(nodes)%name
!
!AA601 comm
! LenR and LenS are the length scales of the inlet section: they are computed as the distance between the first and the last inlet vertices 
! and the third and the last inlet vertices, respectively
! Particles are aligned with Plast-P1 and Plast-P3, where P1 the first boundary vertex, ..., Plast the last boundary vertex. 
! In case of a triangular inlet, we have particles aligned with one direction: P3-P1.
! In case of a quadrilateral inlet, we have particles distributed along two directions: P4-P1 and P4-P3.
    LenR = zero
    LenS = zero
    do sd = 1, SPACEDIM
!???   A(sd) = BoundaryFace(SourceFace)%Node(nodes)%GX(sd)
!???   rr(sd) = BoundaryFace(SourceFace)%T(sd, 1)
!???   ss(sd) = BoundaryFace(SourceFace)%T(sd, 2)
!      nn(sd) = BoundaryFace(SourceFace)%T(sd, 3)
      LenR = LenR + (BoundaryFace(SourceFace)%Node(1)%GX(sd) - BoundaryFace(SourceFace)%Node(nodes)%GX(sd))**2
      LenS = LenS + (BoundaryFace(SourceFace)%Node(3)%GX(sd) - BoundaryFace(SourceFace)%Node(nodes)%GX(sd))**2
!      RowVelocity = RowVelocity + partz(izone)%vel(sd) * nn(sd) 
!      partz(izone)%vel(sd) = RowVelocity * nn(sd) 
    end do
    
!    if (RowVelocity == zero) return
    LenR = Dsqrt(LenR)
    LenS = Dsqrt(LenS)
    NumPartR = Int(LenR / deltapart + 0.01d0)
    NumPartS = Int(LenS / deltapart + 0.01d0)
    deltaR = LenR / NumPartR 
    deltaS = LenS / NumPartS 
    eps = -half
    zfila = eps * deltapart
    distR = -half * deltaR
    ip = 0
    do i = 1, NumPartR
      distR = distR + deltaR
      csi = distR / LenR
      distS = -half * deltaS
      do j = 1, NumPartS
        distS = distS + deltaS
        etalocal = distS / LenS
        ip = ip + 1
        do sd = 1, SPACEDIM
          P(sd) = BoundaryFace(SourceFace)%Node(4)%GX(sd) * (one - csi) + &
                  BoundaryFace(SourceFace)%Node(1)%GX(sd) * csi
          Q(sd) = BoundaryFace(SourceFace)%Node(3)%GX(sd) * (one - csi) + &
                  BoundaryFace(SourceFace)%Node(2)%GX(sd) * csi
!
!AA405sub
          PartLine(i_source,ip,sd) = P(sd) * (one - etalocal) + Q(sd) * etalocal
!
         end do
      end do
    end do
!
!AA405 sub
    NumPartFace(i_source) = ip
!
!    ParticleVolume = deltaR * deltaS * deltapart
!    RowPeriod = deltapart/ Abs(RowVelocity)
    ParticleVolume = Domain%PVolume
!
!AA405 sub
    RowPeriod = ParticleVolume * NumPartFace(i_source) / Tratto(nt)%FlowRate 
!
!    RowVelocity = Tratto(nt)%NormVelocity 
!
!AA405 sub start
    RowVelocity(i_source) = Domain%dd / RowPeriod
    Tratto(nt)%NormVelocity = RowVelocity(i_source)
!    partz(irz)%vel(1) = RowVelocity * nn(1)
!    partz(irz)%vel(2) = RowVelocity * nn(2)
!    partz(irz)%vel(3) = RowVelocity * nn(3)
!    if (RowVelocity == zero) return 
!AA405 sub end 

    pinttimeratio = -1
  end if
!
!AA405 start
  end if ! BoundarySide(isi)%tipo == "sour"
  end do ! isi
!AA405 end
!
return
end subroutine PreSourceParticles_3D
!---split

