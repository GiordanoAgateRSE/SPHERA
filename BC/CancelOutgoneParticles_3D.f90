!cfile CancelOutgoneParticles_3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : CancelOutgoneParticles_3D
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
! Module purpose : Module to count and delete the outgoing particles on boundaries
!                  of type "leve", "flow", "velo", "crit", "open"
!
! Calling routine: Loop_Irre_3D
!
! Called routines: LocalNormalCoordinates
!
!************************************************************************************
!
  subroutine CancelOutgoneParticles_3D    

!Conta e cancella le particelle uscite dalle facce di contorno di tipo "leve", "flow", "velo", "crit", "open"
!La cancellazione avviene con due diverse modalita'
!a) Se la particella appartiene alla zona di particelle 'maxzone' con indice piu alto,
!   l'unica in cui sono consentite sia una riduzione (per uscita dal contorno) che un aumento
!   (per generazione dalla sorgente)  del numero di particelle attive, allora la cancellazione consiste
!   nel sostituire, nell'array delle particelle "pg()", alla particella "npi" uscita l'ultima particella "nag"
!   e nel decrememtare il numero di particelle attive nag <== nag-1;
!   viene contestualmente modificato l'indirizzo Partz(maxzone)%limit(2) (indice dell'ultima particella
!   della zona) mentre l'indirizzo Partz(maxzone)%limit(1) (indice della prima particella della zona) rimane
!   immutato.
!b) Se la particella appartiene ad una delle zone di particelle Partz(izone) con indice izone < maxzone, nelle
!   quali è consentita solo una diminuzione del numero di particelle (per uscita dal contorno), allora
!   viene semplicemente posto pg(npi)%cella = 0 (particella fuori dalla griglia di riferimento). 
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Local Scalars ..
  integer(4)        :: nfc, iof, npi, sdi, sdj, fkod, nodes, i     !, izona, nt
  logical           :: esci
  double precision  :: deltax21,deltay21,deltaz21,deltax31,deltay31,deltaz31,deltax,deltay,deltaz
  double precision  :: DetPnew, DetPold
!
!.. Local Arrays
  double precision, dimension(1:SPACEDIM) :: LocXY, LocXYZnew,locXYZold,csi
!
!.. Executable Statements ..
!
!.. loops on the boundary opened faces
!
  do iof = 1, NumOpenFaces
!
    nfc = OpenFace(iof)
    nodes = BoundaryFace(nfc)%nodes
    deltax21 = BoundaryFace(nfc)%Node(2)%GX(1) - BoundaryFace(nfc)%Node(1)%GX(1)
    deltax31 = BoundaryFace(nfc)%Node(3)%GX(1) - BoundaryFace(nfc)%Node(1)%GX(1)
    deltay21 = BoundaryFace(nfc)%Node(2)%GX(2) - BoundaryFace(nfc)%Node(1)%GX(2)
    deltay31 = BoundaryFace(nfc)%Node(3)%GX(2) - BoundaryFace(nfc)%Node(1)%GX(2)
    deltaz21 = BoundaryFace(nfc)%Node(2)%GX(3) - BoundaryFace(nfc)%Node(1)%GX(3)
    deltaz31 = BoundaryFace(nfc)%Node(3)%GX(3) - BoundaryFace(nfc)%Node(1)%GX(3)
!
!$omp parallel do default(none) &
!$omp private(npi,i,deltax,deltay,deltaz,DetPnew,DetPold,LocXY,LocXYZnew,locXYZold,sdi,sdj,esci,fkod,csi) &
!$omp shared(nag,pg,BoundaryFace,nfc,nodes,deltax21,deltax31,deltay21,deltay31,deltaz21,deltaz31,OpCount)
!
!.. loops on the particles in the domain
!
    do npi = 1,nag
!
      if (pg(npi)%cella == 0) cycle
!
!.. consider the three nodes that identify the plane of the current face and evaluates the current and old
!.. volume of the tetrahedrons having the three points as base and the current and old particle positions 
!.. as vertices
!
      deltax = pg(npi)%coord(1) - BoundaryFace(nfc)%Node(1)%GX(1)
      deltay = pg(npi)%coord(2) - BoundaryFace(nfc)%Node(1)%GX(2)
      deltaz = pg(npi)%coord(3) - BoundaryFace(nfc)%Node(1)%GX(3)
      DetPnew = deltax21 * deltay31 * deltaz + deltax31 * deltay * deltaz21 + deltax * deltay21 * deltaz31 - &
                deltax21 * deltay * deltaz31 - deltax31 * deltay21 * deltaz - deltax * deltay31 * deltaz21
!
!.. if the current determinant is less than zero (volume oriented), the particle might be out of the boundary face,
!.. since the reference normal is oriented inside the domain
!
      if (DetPnew <= zero) then
!
!.. verifies the signes of determinants: if different, the particle passed through the plane of the boundary face
!
        deltax = pg(npi)%CoordOld(1) - BoundaryFace(nfc)%Node(1)%GX(1)
        deltay = pg(npi)%CoordOld(2) - BoundaryFace(nfc)%Node(1)%GX(2)
        deltaz = pg(npi)%CoordOld(3) - BoundaryFace(nfc)%Node(1)%GX(3)
        DetPold = deltax21 * deltay31 * deltaz + deltax31 * deltay * deltaz21 + deltax * deltay21 * deltaz31 - &
                  deltax21 * deltay * deltaz31 - deltax31 * deltay21 * deltaz - deltax * deltay31 * deltaz21
        if (sign(one,DetPnew) == sign(one,DetPold)) cycle
!
!.. verifies if the particle path crosses the boundary face area. As first step, evaluates all the coordinates
!.. in the local system of the face plane [x'.y',z'] = Tgl*[x-x3,y-y3,z-z3]
!        
        do sdi = 1, spacedim
!
          LocXYZnew(sdi) = zero
          LocXYZold(sdi) = zero
!
          do sdj = 1, spacedim
!
            LocXYZnew(sdi) = LocXYZnew(sdi) + BoundaryFace(nfc)%T(sdj, sdi) *  &
                             (pg(npi)%coord(sdj) - BoundaryFace(nfc)%Node(nodes)%GX(sdj))
            LocXYZold(sdi) = LocXYZold(sdi) + BoundaryFace(nfc)%T(sdj, sdi) *  &
                             (pg(npi)%CoordOld(sdj) - BoundaryFace(nfc)%Node(nodes)%GX(sdj))
!
          end do
        end do
!
!.. evaluates the coordinates of the intersection point in the local coordinate system; by definition, the z' value
!.. is equal to the constant of the plane equation z' = c
!
        deltaz =  -LocXYZold(3) / (LocXYZnew(3) - LocXYZold(3))
        locXY(1) =  deltaz * (LocXYZnew(1) - LocXYZold(1))
        locXY(2) =  deltaz * (LocXYZnew(2) - LocXYZold(2))
!
!.. transform the local coordinates on the face in the XYZ system into the local coordinates csi1,csi2,csi3 in the local system
!
        Call LocalNormalCoordinates (LocXY, csi, nfc)
!
!.. checks if the projected point falls inside the face (triangular=3 nodes, rectangular=4 nodes)
!
        esci = .true.
        fkod = 6 - BoundaryFace(nfc)%nodes 
        do i = 1,fkod
          if (csi(i) < zero .or. csi(i) > one) then
            esci = .false.
            exit
          end if
        end do
!
!.. the particle projection falls outside the face and therefore must be deleted 
!
        if (esci) then  
!
          OpCount(pg(npi)%imed) = OpCount(pg(npi)%imed) + 1    
          pg(npi)%cella = -1
!
        end if
      end if
    end do
!
!$omp end parallel do
!
  end do
!
!.. the particle array is compacted; if the npi-th particle is zero the last particle is moved into the npi-th array location
!
!!!!!!  npi = 0
!!!!!!  do while (nag > 0 .AND. npi < nag)
!!!!!!    npi = npi + 1
!!!!!!    if (pg(npi)%cella == 0 ) then
!!!!!!!.. the counter of the removed particles is increased by 1
!!!!!!      OpCount(pg(npi)%imed) = OpCount(pg(npi)%imed) + 1    
!!!!!!      pg(npi) = pg(nag)
!!!!!!      pg(nag)%cella = 0
!!!!!!      pg(nag) = PgZero
!!!!!!      nag = nag - 1
!!!!!!      npi = npi - 1
!!!!!!    end if
!!!!!!  end do
!
  return
  end subroutine CancelOutgoneParticles_3D
!---split

