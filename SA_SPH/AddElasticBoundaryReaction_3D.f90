!cfile AddElasticBoundaryReaction_3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : AddElasticBoundaryReaction_3D
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
! Module purpose : Module to compute the boundary integral IntWdS
!
! Calling routine: AddBoundaryContribution_to_CE3D
!                  AddBoundaryContributions_to_ME3D
!
! Called routines: 
!
!************************************************************************************
!
subroutine AddElasticBoundaryReaction_3D (npi, Ncbf, BoundReaction)

! Adds supplementari normal boundary reaction to reinforce insufficient
! pressure gradient in case of few neighbouring particles and presence of
! normal component of mass force (gravity)
! The normal reaction is computed with the formula R=(c0^2/d) ln(zi/d) [for zi<d],
! stemming from the compressible reaction of the fluid, where:
! c0^2 = E/ro0 is the square celerity of the fluid;
! zi is the distance of the particle Pi from the boundary face
! d is a reference distance from which the reaction is added

!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
integer(4),      intent(IN)    :: npi, Ncbf
double precision,intent(INOUT),dimension(1:SPACEDIM) :: BoundReaction
!
!.. Local Parameters ..
double precision,parameter :: zmincoeff = 0.25d0
double precision,parameter :: reafactor = 1.0d0
!
!.. Local Scalars ..
integer(4) :: sd, icbf, iface, nt, ibdt, ibdp, mate, fkod
integer(4) :: ne, NCloseEdgeF
double precision :: zi, zimin, celer02, vin, normreact
double precision :: scaprod, edgelen2, tau, edgedist2
!
!.. Local Arrays ..
double precision,dimension(1:SPACEDIM)   :: PXLoc, csi
double precision,dimension(1:SPACEDIM)   :: XQ, QP, QPcosdir
logical, dimension(1:Domain%MAXCLOSEBOUNDFACES) :: ReaFace
!
! External functions and subrotuines
logical, external    :: IsPointInternal
!
!.. Executable Statements ..
!
  BoundReaction = zero
  mate = pg(npi)%imed
  zimin = zmincoeff * Domain%dd
  celer02 = Med(mate)%eps / Med(mate)%den0
!
  ibdt = BoundaryDataPointer(3,npi)
  do icbf = 1, Ncbf
!
    ibdp = ibdt + icbf - 1
    iface = BoundaryDataTab(ibdp)%CloBoNum
    nt = BoundaryFace(iface)%stretch
!
    ReaFace(icbf) = .false.
!
    if (Tratto(nt)%tipo == "fixe" .or. Tratto(nt)%tipo == "tapi") then
!
      PXLoc(:) = BoundaryDataTab(ibdp)%LocXYZ(:)
      zi = PXLoc(3)
!
      if (zi < zimin) then
!
        call LocalNormalCoordinates (PXLoc, csi, iface)
        fkod = BoundaryFace(iface)%nodes - 2
!
        if (IsPointInternal(fkod, csi)) then    !La proiezione normale della particella npi sul piano
                                                ! della faccia "iface" è interna alla faccia "iface"
          vin = zero
          do sd = 1,SPACEDIM
            vin = vin + pg(npi)%var(sd) * BoundaryFace(iface)%T(sd, 3)
          end do
          if (vin < zero) then
!!!            normreact = -reafactor * celer02 * DLog(zi / zimin) / zimin
            normreact = -reafactor * celer02 * DLog((Domain%h + zi - zimin) / Domain%h) / Domain%h
            BoundReaction(:) = BoundReaction(:) + normreact * BoundaryFace(iface)%T(:, 3)
          end if
          ReaFace(icbf) = .true.
        else
          ReaFace(icbf) = .false.
        end if
      end if
    end if
  end do
!
!.. Reazione da eventuali spigoli (edges) vicini
  do ne = 1, NumBEdges
    !.. Verifica se la particella è vicina ad almeno una delle facce comuni allo spigolo ne
    NCloseEdgeF = 0
    
    ibdt = BoundaryDataPointer(3,npi)
    do icbf = 1, Ncbf
!
      ibdp = ibdt + icbf - 1
      iface = BoundaryDataTab(ibdp)%CloBoNum
      if (iface == BoundaryConvexEdge(ne)%face(1) .and. .Not. ReaFace(icbf)) then
        NCloseEdgeF = NCloseEdgeF + 1
      else if (iface == BoundaryConvexEdge(ne)%face(2) .and. .Not. ReaFace(icbf)) then
        NCloseEdgeF = NCloseEdgeF + 1
      end if
    end do
    if (NCloseEdgeF /= 2) return
    
    !..Distanza zi della particella npi dallo spigolo ne
    scaprod = zero
    do sd = 1,SPACEDIM
      scaprod = scaprod + (pg(npi)%Coord(sd) - BoundaryConvexEdge(ne)%node(1)%GX(sd)) * BoundaryConvexEdge(ne)%component(sd)
    end do
    edgelen2 = BoundaryConvexEdge(ne)%length * BoundaryConvexEdge(ne)%length
    tau = scaprod / edgelen2
    if (tau >= zero .and. tau <= one) then
      edgedist2 = zero
      XQ(:) = BoundaryConvexEdge(ne)%node(1)%GX(:) + BoundaryConvexEdge(ne)%component(:) * tau
      QP(:) = pg(npi)%Coord(:) - XQ(:)
      do sd = 1,SPACEDIM
        edgedist2 = edgedist2 + QP(sd) * QP(sd)
      end do
      zi = Dsqrt(edgedist2)
      if (zi < zimin) then
        vin = zero
        QPcosdir(:) = QP(:) / zi
        do sd = 1,SPACEDIM
          vin = vin + pg(npi)%var(sd) * QPcosdir(sd)
        end do
        if (vin < zero) then
          normreact = -reafactor * celer02 * DLog(zi / zimin) / zimin
          BoundReaction(:) = BoundReaction(:) +  normreact * QPcosdir(:)
        end if
      end if
    end if
  end do
!
!prova!
!  BoundReaction = floor(BoundReaction * azzeramento) / azzeramento
!  where (dabs(BoundReaction) < arrotondamento) BoundReaction = zero
!prova!
!
return 
end subroutine AddElasticBoundaryReaction_3D
!---split

