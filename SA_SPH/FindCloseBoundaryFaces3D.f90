!cfile FindCloseBoundaryFaces3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name    : FindCloseBoundaryFaces3D
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
! Module purpose : Module Loop
!
! Calling routine: ComputeBoundaryDataTab
!
! Called routines: diagnostic
!                  CellIndices
!                  CellNumber
!                  IsPointInternal
!                  ParticleCellNumber
!                  LocalNormalCoordinates
!
!************************************************************************************
!
  subroutine FindCloseBoundaryFaces3D ( npi, Ncbf, Clobface, LocX, Nfzn )

!.. Finds the "close" boundary faces, i.e. those sited at a distance from the particle npi 
!.. less than or equal to 2h (where h is the smoothing length)
!
!.. Returns:
!.. Ncbf                   = Number of close boundary faces
!.. Clobface(1 to Ncbf)    = List of close boundary faces
!.. LocX(1:SPACEDIM, Ncbf) = Local coordinates of particle npi with respect each boundary side
!
!.. The algorithm looks for boundary faces intersected by the cell boxes of the reference frame
!.. sited all around particle npi, and cancels the repeated ones
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
  integer(4),      intent(IN)                                                :: npi
  integer(4),      intent(INOUT)                                             :: Ncbf, Nfzn
  integer(4),      intent(INOUT),dimension(1:Domain%MAXCLOSEBOUNDFACES)             :: Clobface
  double precision,intent(INOUT),dimension(1:SPACEDIM, 1:Domain%MAXCLOSEBOUNDFACES) :: LocX
!
!.. Local Scalars ..
  integer(4)       :: nc, ic, jc, kc, i, j, k, sdi, sdj, nodes, irestocell, fkod      !ni, nj, nk,  
  integer(4)       :: flpini, flp, flpfin, nfpercell, intbf, icbf, nbface, stretch
  double precision :: pin, pinmin, pinmax
  logical          :: Thereis
  character(len=lencard) :: nomsub = "FindCloseBoundaryFaces3D"
!
!.. Local Arrays ..
  double precision,dimension(1:SPACEDIM) :: PXLoc, csi
!
!.. External routines ..
  integer(4), external :: CellNumber, ParticleCellNumber, CellIndices
  logical,    external :: IsPointInternal
!
!.. Executable Statements ..
!
  Clobface = 0
  Ncbf = 0
  Nfzn = 0
  LocX = zero
  pg(npi)%CloseBcOut = 0
!
!.. find the cell number for the current particle
!
  nc = ParticleCellNumber(pg(npi)%coord)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  if (Domain%RKscheme /= 2) then
!    if (nc /= pg(npi)%cella) then
!      write (nout,*) "Iteration=",it_corrente,"nc =",nc,"different from pg%cella= ",pg(npi)%cella
!      call diagnostic (arg1=8,arg2=5,arg3=nomsub)
!    end if
!  end if
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  if (nc <= 0) return
!
!.. find the indices of the cell in the domain grid
!
!  ni = Grid%ncd(1)        ! inutile
!  nj = Grid%ncd(2)        ! inutile
!  nk = Grid%ncd(3)        ! inutile
  irestocell = CellIndices (nc, ic, jc, kc)
!
!.. loops on the cells surrounding the current one
!
  do i = ic - 1, ic + 1
!!!!!!!!!    if (1 <= i .and. i <= ni) then        ! inutile
!
      Do j = jc - 1, jc + 1
!!!!!!!!!        if (1 <= j .and. j <= nj) then        ! inutile
!
          do k = kc - 1, kc + 1
!!!!!!!!!            if (1 <= k .and. k <= nk) then        ! inutile
!
!.. find the cell number of the surrounding cells
!
              nc = CellNumber(i, j, k)
              if (nc == 0) cycle
!
!.. load the number of boundary faces cutting the cell and the initial and final pointers
!
              nfpercell = GCBFPointers(nc, 1)
              if (nfpercell > 0) then
                flpini = GCBFPointers(nc, 2)
                flpfin = flpini + nfpercell - 1
!
!.. loops on the cutting faces
!
                do flp = flpini, flpfin
!
!.. load the face index
!
                  intbf = GCBFVector(flp)
                  thereis = .false.
                  pinmin = zero
                  pinmax = doubleh
!                  pinmin = -doubleh  !prova 2 dic 2008
                  stretch = BoundaryFace(intbf)%stretch
!
                  if (Tratto(stretch)%tipo == "sour") pinmin = -doubleh   !prova 2 dic 2008
!
!.. checks if the face intbf is already included in the Clobface() array
!
                  if (Ncbf > 0) thereis = any(Clobface(1:Ncbf) == intbf)
!
!.. the face is not yet considered
!
                  if (.Not. Thereis) then         
!
!.. evaluate the normal coordinate of the npi particle with respect the face intbf
!
                    nodes = BoundaryFace(intbf)%nodes
                    do sdi = 1, SPACEDIM
                      PXLoc(sdi) = zero
                      do sdj = 1, SPACEDIM
                        PXLoc(sdi) = PXLoc(sdi) + BoundaryFace(intbf)%T(sdj, sdi) * &
                                     (pg(npi)%coord(sdj) - BoundaryFace(intbf)%Node(nodes)%GX(sdj))
                      end do
                    end do
                    pin = PXLoc(3)
                    call LocalNormalCoordinates (PXLoc, csi, intbf)
                    fkod = BoundaryFace(intbf)%nodes - 2
!
!.. the distance between the particle and the face is greater than zero and less then 2h (influence diameter)
!.. the face is considered and it is added to the Clobface() array
!
!!!!primaprova                    if (pin > pinmin .and. pin < pinmax) then !!!nuovo test
                    if (pin >= pinmin .and. pin < pinmax) then !!!nuovo test
!
                      if (Tratto(stretch)%tipo == "sour" .or. Tratto(stretch)%tipo == "velo" .or. &
                           Tratto(stretch)%tipo == "flow") then
                        if (IsPointInternal(fkod, csi))  then  !La proiezione normale della particella npi sul piano
                                                               ! della faccia "iface" è interna alla faccia "iface"
                          Ncbf = ncbf + 1 
                          if ( ncbf <= Domain%MAXCLOSEBOUNDFACES ) then
                            Clobface(Ncbf) = intbf
                            LocX(3, Ncbf) = pin
                            pg(npi)%CloseBcOut = 1
                          else
                            call diagnostic (arg1=8,arg2=6,arg3=nomsub)
                          end if
                        end if
                      else
                        Ncbf = ncbf + 1 
                        if ( ncbf <= Domain%MAXCLOSEBOUNDFACES ) then
                          Clobface(Ncbf) = intbf
                          LocX(3, Ncbf) = pin
                        else
                          call diagnostic (arg1=8,arg2=7,arg3=nomsub)
                        end if
                      end if
!
!!!!primaprova                    else if (pin <= pinmin) then   !AG 05dic2008  aggiungere test se la proiezione della particella sta all'interno della faccia
                    else if (pin < pinmin) then   !AG 05dic2008  aggiungere test se la proiezione della particella sta all'interno della faccia
                      if (IsPointInternal(fkod, csi)) then    !La proiezione normale della particella npi sul piano
                                                              ! della faccia "iface" è interna alla faccia "iface"
                         Nfzn = Nfzn + 1
                      end if
!????
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!                      write (778,*) " Nfzn >",Nfzn," it>",it_corrente
!!                      write (778,*) " doubleh>",doubleh," nfpercell>",nfpercell
!!                      write (778,*) " dcd  >",grid%dcd(1:3)
!!                      write (778,*) " particella>",npi," intbf>",intbf," cella>",nc
!!                      write (778,*) " pin  >",pin,"faccia>",BoundaryFace(intbf)%stretch
!!                      write (778,*) " coord>",pg(npi)%coord(1:3)
!!                      write (778,*) " T    >",BoundaryFace(intbf)%T(1:3, sdi)
!!                      write (778,*) " gx   >",BoundaryFace(intbf)%Node(nodes)%GX(1:3)
!!                      write (778,*) "==============================================================="
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    end if
                  end if
                end do
!
              end if
!!!!!!!!!            end if        ! inutile
          end do
!
!!!!!!!!!        end if        ! inutile
      end do
!
!!!!!!!!!    end if        ! inutile
  end do
!
!.. evaluates the tangent coordinates (r,s) of the particle npi in the local coordinates of the closest faces
!
  if (Ncbf > 0) then
!
!.. loops on the found faces
!
    do icbf = 1, Ncbf
      nbface = Clobface(icbf)
      nodes = BoundaryFace(nbface)%nodes
!
!.. incremento numero particelle vicine al contorno e calcolo massima quota
!$omp critical (numpa)
      BoundaryFace(nbface)%CloseParticles = BoundaryFace(nbface)%CloseParticles + 1
      if (BoundaryFace(nbface)%CloseParticles_maxQuota < pg(npi)%coord(3)) &
                              BoundaryFace(nbface)%CloseParticles_maxQuota = pg(npi)%coord(3)
!$omp end critical (numpa)
!
      do sdi = 1, PLANEDIM
        LocX(sdi, icbf) = zero
!
        do sdj = 1, SPACEDIM
          LocX(sdi, icbf) = LocX(sdi, icbf) + BoundaryFace(nbface)%T(sdj, sdi) * &
                            (pg(npi)%coord(sdj) - BoundaryFace(nbface)%Node(nodes)%GX(sdj))
        end do
      end do
    end do
!
  end if
!
  return
  end subroutine FindCloseBoundaryFaces3D
!---split

