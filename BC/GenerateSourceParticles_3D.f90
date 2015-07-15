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

