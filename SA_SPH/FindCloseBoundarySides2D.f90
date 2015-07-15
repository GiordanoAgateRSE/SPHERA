!cfile FindCloseBoundarySides2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : FindCloseBoundarySides2D
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
! Module purpose : Module for Finds the "close" boundary sides
!
! Calling routine: ComputeBoundaryDataTab
!
! Called routines: diagnostic
!
!************************************************************************************
!
  subroutine FindCloseBoundarySides2D (npi, Ncbs, Cloboside, LocXY)
!
!.. Finds the "close" boundary sides, i.e. those sited at a distance from particle npi <= 2h (where h is the smoothing length)
!.. Returns:
!.. Ncbs                     = Number of close boundary sides (= 0, 1, 2)
!.. Cloboside(1:Ncbs)        = List of close boundary sides
!.. LocXY(1:PLANEDIM,1:Ncbs) = Local coordinates of particle npi with respect each boundary side (vertex V1)
!
!.. assign modules
  use FILES_ENTITIES
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
  use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
  integer(4),                                                 intent(in)  :: npi
  integer(4),                                                 intent(out) :: Ncbs
  integer(4),      dimension(1:MAXCLOSEBOUNDSIDES),           intent(out) :: Cloboside
  double precision,dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES),intent(out) :: LocXY
!
!.. Local Scalars ..
  integer(4)       :: icbs, icb, isous, isi, v1, v2, sd, iside1, iside2  !mate, 
  double precision :: xp, yp, sidel, xpmin, xpmax, ypmin, ypmax, xpq
  double precision :: Lmxpq, sidelen, ypmn, ypmx
  character(len=lencard) :: nomsub = "FindCloseBoundarySides2D"
!
!.. local Arrays ..
  integer(4),      dimension(1:PLANEDIM) :: acix
  double precision,dimension(1:PLANEDIM) :: Plocal, P1, P1P, sss, nnn
!
!.. Executable Statements ..
!
!.. initializations ..
!
  acix(1) = 1  
  acix(2) = 3
  Cloboside = 0
  LocXY     = zero
  Ncbs  = 0
  Plocal(:) = pg(npi)%coord(acix(:))
!  mate  = pg(npi)%imed
  ypmin = - doubleh
  ypmax = doubleh
  pg(npi)%CloseBcOut = 0
!
!.. loops on all the boundary sides of the domain
!
  side_loop: do isi = 1, NumBSides
!
!.. the perimeter and pool boundary types are skipped
!
    if (BoundarySide(isi)%tipo /= "peri" .AND. BoundarySide(isi)%tipo /= "pool" ) then
!
!.. loads the side coordinates P1, the local versors sss and nnn and
!.. evaluates the distances between the current particle coordinates Plocal and the 
!.. reference vertex V1
!    
      v1 = BoundarySide(isi)%Vertex(1)
      v2 = BoundarySide(isi)%Vertex(2)
!
      do sd = 1, PLANEDIM
!
        P1(sd) = Vertice(acix(sd),v1)
        P1P(sd) = Plocal(sd) - P1(sd)
        sss(sd) = BoundarySide(isi)%T(acix(sd), acix(1))
        nnn(sd) = BoundarySide(isi)%T(acix(sd), acix(2))
!
      end do
!
!.. evaluates the particle coordinates xp and yp in the side local reference system, where
!.. x in the direction of the side segment, y is the normal direction and V1 the origin
!
      sidel = BoundarySide(isi)%length
!
      xp = P1P(1) * sss(1) + P1P(2) * sss(2)
      yp = P1P(1) * nnn(1) + P1P(2) * nnn(2)
!
!.. set the interaction area in the neighborough of the side segment having a distance equal to
!.. +/- 2*h in the y local direction and a distance of 2h from the vertices V1 and V2 in the x local direction
! 
      xpmin = - doubleh
      xpmax = sidel + doubleh
!
!.. checks if the particle has local coordinates falling inside the interaction area
!        
      if ((xpmin < xp .AND. xp < xpmax) .AND. (ypmin < yp .AND. yp < ypmax)) then
!
!.. the particle falls in the extreme cicle around V1
!            
        if (xp < zero) then
!
          xpq = xp * xp
          ypmx = Dsqrt(doublesquareh - xpq)
          ypmn = -ypmx 
!
!.. the particle falls in the segment natural length
!
        else if (xp <= sidel) then
!
          ypmx = ypmax    !CONTROLLARE
          ypmn = ypmin 
!
!.. the particle falls in the extreme cicle around V2
!
        else if (xp < xpmax) then
!
          Lmxpq = (sidel - xp) * (sidel - xp)
          ypmx = Dsqrt(doublesquareh - Lmxpq)
          ypmn = -ypmx 
!
        end if
!
!.. the boundary must be considered for the current particle as close boundary
!            
        if (ypmn < yp .AND. yp < ypmx) then
!
!.. the number of close boundaries is increased
!
          Ncbs = Ncbs + 1
!
!.. checks against the maximum number allowed for the closest boundaries
!
          if (Ncbs > MAXCLOSEBOUNDSIDES) then
!.. The particle npi has more than two boundary sides: reduce dd and restart
            write (nout,'(1x,a,i12,a)')   ' ERROR! The particle ',npi,' has more than two boundary sides: reduce dd and restart.'
            write (nout,'(1x,a,3f15.10)') '        Coordinate: ',pg(npi)%coord(:)
            write (nscr,'(1x,a,i12,a)')   ' ERROR! The particle ',npi,' has more than two boundary sides: reduce dd and restart.'
            write (nscr,'(1x,a,3f15.10)') '        Coordinate: ',pg(npi)%coord(:)
            call diagnostic (arg1=8,arg2=8,arg3=nomsub)
          end if
!
!.. stores the boundary index and the local coordinates on the boundary segment
!
          Cloboside(Ncbs) = isi
          LocXY(1, Ncbs) = xp
          LocXY(2, Ncbs) = yp
!
        end if
      end if
    end if
  end do side_loop
!
!.. searches for a nearest source type side
!
  isous = 0
  do icbs = 1, Ncbs
    isi = Cloboside(icbs)
!
!.. incremento numero particelle vicine al contorno e calcolo massima quota
!$omp critical (numpa)
    BoundarySide(isi)%CloseParticles = BoundarySide(isi)%CloseParticles + 1
    if (BoundarySide(isi)%CloseParticles_maxQuota < pg(npi)%coord(3)) BoundarySide(isi)%CloseParticles_maxQuota = pg(npi)%coord(3)
!$omp end critical (numpa)
!
    if (BoundarySide(isi)%tipo == "sour") then
      XP = LocXY(1, icbs)
      sidel = BoundarySide(isi)%length
      if (XP > zero .AND. XP < sidel) then
        isous = icbs    !contrassegna la sorgente per eliminare gli altri eventuali lati vicini
        exit
      else                !cancella il lato sorgente perché ininfluente
        if (icbs < Ncbs) then
          do icb = icbs+1, Ncbs
            Cloboside(icb -1) = Cloboside(icb)
            LocXY(1, icb -1) = LocXY(1, icb)
            LocXY(2, icb -1) = LocXY(2, icb)
          end do
          Ncbs = Ncbs -1
          exit
        else if (icbs == Ncbs) then
          Ncbs = Ncbs -1
          exit
        end if
      end if
    end if
  end do
!
!.. a source has been found: the other nearest sides are deleted
!
  if (isous > 0) then   
    Ncbs = 1
    Cloboside(Ncbs) = Cloboside(isous)
    LocXY(1, Ncbs) = LocXY(1, isous)
    LocXY(2, Ncbs) = LocXY(2, isous)
  end if
!
!.. checks if more than two sides are close the particle
!
  if (Ncbs > LIMCLOSEBOUNDSIDES) then
    call diagnostic (arg1=8,arg2=9,arg3=nomsub)
  end if
!
!.. there are two close boundaries:
!
!Case of two close boundary sides
!Change of the origin of local abcsissa in one of the two adjacent boundary sides
!in order that the abscissa origin coincide with the common vertex in each side.
!
  if (Ncbs == 2) then
!
    iside1 = Cloboside(1)
    iside2 = Cloboside(2)
!
    if (BoundarySide(iside1)%previous_side == iside2) then
      sidelen = BoundarySide(iside2)%length
      LocXY(1, 2) = sidelen - LocXY(1, 2)
    else if (BoundarySide(iside2)%previous_side == iside1) then
      sidelen = BoundarySide(iside1)%length
      LocXY(1, 1) = sidelen - LocXY(1, 1)
    else
!     Sides not adjacent: this case should never occur if the minimum length
!      of boundary sides is >4h (diameter of influence circle)
    end if
  end if
!
  return
  end subroutine FindCloseBoundarySides2D
!---split

