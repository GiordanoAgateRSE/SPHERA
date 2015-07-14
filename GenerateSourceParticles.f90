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

!cfile GenerateSourceParticles_2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : GenerateSourceParticles_2D
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
! 04  Amicarelli/Agate  05/07/2011     Time stage parameters
! 05  Amicarelli        23/11/2011     multiple inlet
!
!************************************************************************************
! Module purpose : Module to generate new source particles to simulate inlet fluid flow
!                  !!! ONLY in 2D and with one inlet !!!
!
! Calling routine: Loop_Irre_2D
!
!AA601 sub
! Called routines: defcolpartzero, wavy_inlet
!
!************************************************************************************
!
subroutine GenerateSourceParticles_2D
!Genera e immette nel campo nuove particelle di sorgente per simulare una portata fluida entrante
!Attenzione: questa routine funziona solo in presenza di un solo lato di contorno funzionante
!            da sorgente e solo nel caso bidimensionale!!!!
!
!.. assign modules
use GLOBAL_MODULE
use FILES_ENTITIES
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Scalars ..
!
!AA405 sub
integer(4)       :: nt, sd, ip, inttimeratio, isi, i_source
!
double precision :: Time
!
!AA406 sub
double precision :: SourceTime, TimeFrac, DisplFrac,rnd
!
character(len=lencard) :: nomsub = "GenerateSourceParticles_2D"
!
integer(4), external :: ParticleCellNumber
!
!.. Executable Statements ..
!
Time = tempo
!
if ( SourceSide == 0 ) return
!if (FlowRate == zero) return
!
!AA405 rm
!if (RowVelocity == zero) return 
!
inttimeratio = Int(Time / RowPeriod)
if ( inttimeratio > pinttimeratio ) then
!
!AA406
  itime_jet = itime_jet + 1
!
  SourceTime = inttimeratio * RowPeriod 
  TimeFrac = Time - SourceTime 
!
!AA405 rm
!  DisplFrac = RowVelocity * TimeFrac 
!
  continue
!
!AA405 start
  i_source=0
  do isi = 1, NumBSides
  if ( BoundarySide(isi)%tipo == "sour" ) then
  SourceSide = isi
  i_source=i_source+1
  DisplFrac = RowVelocity(i_source) * TimeFrac 
!AA405 end
!
!AA406 sub
  do ip = 1,NumPartperLine(i_source)
!
! Genera una nuova fila di particelle
    nag = nag + 1
    SpCount(mat) = SpCount(mat) + 1
!
    if (nag > PARTICLEBUFFER) then
      call diagnostic (arg1=6,arg2=1,arg3=nomsub)
! Inserire segnalazione di superamento limite 
! oppure ridimensionamento automatico arrays  pg(), nPartintorno() e connesse 
    end if    
!
!.. initializes the attributes of the new particle
!
    Pg(nag) = PgZero
!
    if (Domain%RKscheme > 1) ts0_pg(nag) = ts_pgZero
!
    nt = BoundarySide(SourceSide)%stretch
    do sd = 1, SPACEDIM
      nn(sd) = BoundarySide(SourceSide)%T(sd, 3)
!      pg(nag)%coord(sd) = PartLine(ip, sd) + (yfila + DisplFrac) * nn(sd)
!
!AA405 sub
      pg(nag)%coord(sd) = PartLine(i_source, ip, sd) - (yfila + DisplFrac) * nn(sd)
!      pg(nag)%vel(sd) = partz(irz)%vel(sd) 
      pg(nag)%vel(sd) = Tratto(nt)%NormVelocity * nn(sd) 
      if (tempo < Tratto(nt)%trampa) pg(nag)%vel(sd) = pg(nag)%vel(sd) * tempo / tratto(nt)%trampa
      pg(nag)%var(sd) = pg(nag)%vel(sd) 
    end do
!AA601 sub
    if (Domain%tipo == "bsph") call wavy_inlet(isi)
!AA406
    if (Domain%tipo == "bsph") then
       pg(nag)%rhoSPH_new = zero
       pg(nag)%Gamma = 1.
       pg(nag)%uni = zero
       pg(nag)%sigma = zero
       pg(nag)%dShep = zero 
       pg(nag)%FS = 0 
!
    endif
!
    pg(nag)%izona= irz
    pg(nag)%mass = ParticleVolume * Med(Mat)%den0
!    pg(nag)%acc  = zero
!    pg(nag)%zer  = zero
    pg(nag)%imed = mat           ! indice mezzo 
    pg(nag)%visc = Med(mat)%visc
    pg(nag)%mu   = Med(mat)%visc * Med(Mat)%den0
!
    if (index(Med(mat)%tipo,"liquid") > 0 .or. index(Med(mat)%tipo,"smagorin") > 0) then
      pg(nag)%state  = "flu"
      pg(nag)%VolFra = VFmn
    else if (index(Med(mat)%tipo,"granular") > 0 .or. index(Med(mat)%tipo,"general") > 0) then
!!      pg(nag)%visc = Med(mat)%numx
!!      pg(nag)%mu   = Med(mat)%mumx
      pg(nag)%state  = "sol"
      pg(nag)%VolFra = VFmx
    else if (index(Med(mat)%tipo,"gas") > 0) then
      pg(nag)%state  = "flu"
      pg(nag)%VolFra = VFmn
!    else
!      pg(nag)%state = "   "
    end if
!
!    pg(nag)%secinv    = zero
    pg(nag)%vel_type  = partz(irz)%move        ! indice moto
    if ( partz(irz)%move /= "std" ) pg(nag)%visc = zero
    pg(nag)%slip  = partz(irz)%slip        ! condizione al contorno
!
!AA406 sub
    pg(nag)%cella = ParticleCellNumber(pg(nag)%coord)
!
!    pg(nag)%mno   = zero
!    pg(nag)%velmorr(:)= zero
!    pg(nag)%dudx      = zero
!    pg(nag)%dudy      = zero
!    pg(nag)%dvdx      = zero
!    pg(nag)%dvdy      = zero
! colore part
    call defcolpartzero ( irz,partz,pg(nag) )
!
! inizializzazione press e dens in relazione a quanto assegnato in input
! data p si ricava ro
!
    if (partz(irz)%pressure == "pa") then       ! pressione assegnata 
      pg(nag)%pres = partz(irz)%valp  
    else if (partz(irz)%pressure == "qp") then   ! quota piezometrica assegnata 
      pg(nag)%pres = Med(mat)%den0 * Domain%grav(3) * (pg(nag)%coord(3) - partz(irz)%valp)
    end if  
    pg(nag)%dens = Med(mat)%den0 * (one + pg(nag)%pres/Med(mat)%eps)       
!    pg(nag)%dden = zero                 
  end do
!
!AA405 start
  end if ! BoundarySide(isi)%tipo == "sour"
  end do ! isi
!AA405 end
!
  pinttimeratio = inttimeratio
!
end if
!
return
end subroutine GenerateSourceParticles_2D
!---split

!cfile GenerateSourceParticles_3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : GenerateSourceParticles_3D
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
! 04  Amicarelli/Agate  05/07/2011     Time stage parameters
! 05  Amicarelli        23/11/2011     multiple inlet
!
!************************************************************************************
! Module purpose : Module to generate new source particles to simulate inlet fluid flow
!                  !!! ONLY in 3D and with one inlet face with four nodes !!!
!
! Calling routine: Loop_Irre_3D
!
! Called routines: defcolpartzero
!
!************************************************************************************
!
subroutine GenerateSourceParticles_3D 
!Genera e immette nel campo nuove particelle di sorgente per simulare una portata fluida entrante
!Attenzione: questa routine funziona solo in presenza di una sola faccia di contorno fuzionante
!            da sorgente e solo nel caso che questa sia un parallelogramma (4 nodi)
!
!.. assign modules
use GLOBAL_MODULE
use FILES_ENTITIES
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE
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
integer(4) :: nt, sd, ip, inttimeratio, isi, i_source
!
double precision :: Time
!
!AA406 sub
double precision :: SourceTime, TimeFrac, DisplFrac,rnd1 !,rnd2
!
character(len=lencard) :: nomsub = "GenerateSourceParticles_3D"
!
integer(4), external :: ParticleCellNumber
!
!.. Executable Statements ..
!
Time = tempo

if ( SourceFace == 0 ) return
!
!AA405 rm
!if (RowVelocity == zero) return
!
inttimeratio = Int(Time / RowPeriod)
if ( inttimeratio > pinttimeratio ) then
!
!AA406
  itime_jet = itime_jet + 1
!
  SourceTime = inttimeratio * RowPeriod
  TimeFrac = Time - SourceTime
!
!AA405 rm
!  DisplFrac = RowVelocity * TimeFrac
!
!AA405 start
  i_source=0
  do isi = 1, NumFacce
  SourceFace = isi
  nt = BoundaryFace(SourceFace)%stretch
  if (tratto(nt)%tipo == "sour" ) then
  i_source=i_source+1
  DisplFrac = RowVelocity(i_source) * TimeFrac 
!AA405 end
!
!AA405 sub
  do ip = 1, NumPartFace(i_source)
!
!Genera una nuova fila di particelle
    nag = nag + 1
    SpCount(mat) = SpCount(mat) + 1
!    partz(partzone)%limit(2) = nag
  
    if (nag > PARTICLEBUFFER) then             !Attenzione!
      call diagnostic (arg1=6,arg2=2,arg3=nomsub)
!  Attenzione!
!Inserire segnalazione di superamento limite
!oppure ridimensionamento automatico arrays  pg(), nPartintorno()  e connesse
    end if
!
    pg(nag) = PgZero
!
    if (Domain%RKscheme > 1) ts0_pg(nag) = ts_pgZero
!
!AA405 rm
!    nt = BoundaryFace(SourceFace)%stretch
!
    do sd = 1, SPACEDIM
      nn(sd) = BoundaryFace(SourceFace)%T(sd, 3)
!
!AA405 sub
      pg(nag)%coord(sd) = PartLine(i_source, ip, sd) + (zfila + DisplFrac) * nn(sd)
!
!      pg(nag)%coord(sd) = PartLine(ip, sd) - (zfila + DisplFrac) * nn(sd)
!      pg(nag)%vel(sd) = partz(izone)%vel(sd)
      pg(nag)%vel(sd) = Tratto(nt)%NormVelocity * nn(sd)
      pg(nag)%var(sd) = pg(nag)%vel(sd)
    end do
!
!AA601 sub
if (Domain%tipo == "bsph") call wavy_inlet(isi)
!AA406 start
    if (Domain%tipo == "bsph") then
       pg(nag)%rhoSPH_new = zero
       pg(nag)%Gamma = 1.
       pg(nag)%uni = zero
       pg(nag)%sigma = zero
       pg(nag)%dShep = zero 
       pg(nag)%FS = 0 
    endif
!AA406 end
!
    pg(nag)%izona= izone        !Particle zone
    pg(nag)%mass = ParticleVolume * Med(Mat)%den0
!    pg(nag)%acc  = zero
!    pg(nag)%zer  = zero
    pg(nag)%imed = mat           ! indice mezzo 
    pg(nag)%visc = Med(mat)%visc
    pg(nag)%mu   = Med(mat)%visc * Med(Mat)%den0
!
    if (index(Med(mat)%tipo,"liquid") > 0 .or. index(Med(mat)%tipo,"smagorin") > 0) then
      pg(nag)%state  = "flu"
      pg(nag)%VolFra = VFmn
    else if (index(Med(mat)%tipo,"granular") > 0 .or. index(Med(mat)%tipo,"general") > 0) then
!!      pg(nag)%visc = Med(mat)%numx
!!      pg(nag)%mu   = Med(mat)%mumx
      pg(nag)%state  = "sol"
      pg(nag)%VolFra = VFmx
    else if (index(Med(mat)%tipo,"gas") > 0) then
      pg(nag)%state  = "flu"
      pg(nag)%VolFra = VFmn
!    else
!      pg(nag)%state = "   "
    end if
!
!    pg(nag)%secinv    = zero
    pg(nag)%vel_type  = partz(izone)%move        ! indice moto
    if (partz(izone)%move /= "std" ) pg(nag)%visc = zero
    pg(nag)%slip      = partz(izone)%slip        ! condizione al contorno
!
!AA406 sub
    pg(nag)%cella = ParticleCellNumber(pg(nag)%coord)
!
!    pg(nag)%mno       = zero
!    pg(nag)%velmorr(:)= zero
!    pg(nag)%dudx      = zero
!    pg(nag)%dudy      = zero
!    pg(nag)%dvdx      = zero
!    pg(nag)%dvdy      = zero
!colore part
    call defcolpartzero ( izone,partz,pg(nag) )
!
!inizializzazione press e dens in relazione a quanto assegnato in input
!data p si ricava ro
!
    if (partz(izone)%pressure == "pa") then       ! pressione assegnata
      pg(nag)%pres = partz(izone)%valp
    else if (partz(izone)%pressure == "qp") then   ! quota piezometrica assegnata
      pg(nag)%pres = Med(mat)%den0 * Domain%grav(3) * (pg(nag)%coord(3) - partz(izone)%valp)
    end if
    pg(nag)%dens = Med(mat)%den0 * (one + pg(nag)%pres/Med(mat)%eps)
!    pg(nag)%dden = zero
  end do
!
!
!AA405 start
  end if ! BoundarySide(isi)%tipo == "sour"
  end do ! isi
!AA405 end
!
  pinttimeratio = inttimeratio
!
end if
!
return
end subroutine GenerateSourceParticles_3D
!---split
