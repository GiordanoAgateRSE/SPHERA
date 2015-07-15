!cfile ComputeBoundaryVolumeIntegrals_P0.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ComputeBoundaryVolumeIntegrals_P0
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
! Module purpose : Module to compute boundary contributions to rodivV
!
! Calling routine: ComputeBoundaryDataTab
!
! Called routines: InterpolateTable
!                  IsPointInternal
!                  LocalNormalCoordinates
!
!************************************************************************************
!
subroutine ComputeBoundaryVolumeIntegrals_P0 (icbf, Clobface, LocX, IntWdV, IntdWrm1dV, IntGWZrm1dV, IntGWdV, IntGWrRdV)
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
  integer(4),      intent(IN)  :: icbf
  integer(4),      intent(IN), dimension(1:Domain%MAXCLOSEBOUNDFACES) :: Clobface
  double precision,intent(OUT) :: IntWdV, IntdWrm1dV, IntGWZrm1dV
  double precision,intent(IN), dimension(1:SPACEDIM,1:Domain%MAXCLOSEBOUNDFACES) :: LocX
  double precision,intent(OUT),dimension(1:SPACEDIM)                      :: IntGWdV
  double precision,intent(OUT),dimension(1:SPACEDIM,1:SPACEDIM)           :: IntGWrRdV
!
!.. Local Parameters ..
integer(4),parameter :: nicols = 4
double precision,parameter :: eps = 0.05d0
!
!Computes the boundary volume integral IntWdV
!
!.. Local Scalars ..
integer(4) :: SD, sdj, ipk, fkod
integer(4) :: iface
double precision :: deltaPiPx, PiPxdist, rob, delta_alpha, tsegnato, zpmin, RZ, RZ2, robpre
double precision :: JW3ro2dA, JdW3ro1dA, JdW3ro2dA, JdW3ro3dA
!
!.. Local Arrays ..
integer(4),dimension(1:nicols) :: icol
double precision,dimension(1:SPACEDIM) :: LocPi, LocPj, LocPiPj, PXLoc, csi, RLoccos
double precision,dimension(1:nicols) :: ivalue
!
! External functions and subrotuines
logical, external :: IsPointInternal
!
!.. Executable Statements ..
!
  icol(1) = 1
  icol(2) = 2
  icol(3) = 3
  icol(4) = 4
!
  IntWdV = zero
  IntdWrm1dV = zero
  IntGWZrm1dV = zero
  IntGWdV(:) = zero
  IntGWrRdV(:,:) = zero
!
  JW3ro2dA  = zero
  JdW3ro1dA = zero
  JdW3ro2dA = zero
  JdW3ro3dA = zero
!
  iface = Clobface(icbf)
  zpmin = eps * Domain%h
!
  LocPi(:) = LocX(:, icbf)                     !Coordinate locali della particella reale Pi
!
  if (LocPi(3) < zero)  return
!
  if (LocPi(3) < zpmin)  LocPi(3) = zpmin
!
  robpre = 2.1
  do ipk = 1,BITrows
    LocPiPj(:) = BoundIntegralTab(ipk, 1:3)      !Componenti locali del segmento orientato PiPj
    LocPj(:) = LocPi(:) + LocPiPj(:)             !Coordinate locali del punto Pj
    if (LocPj(3) <= 0) then   !Il punto Pj is si trova nel semispazio della
                              ! faccia "iface", opposto a quello della particella reale "Pi"
      !Coordinate locali del punto PX intersezione del segmento PiPj con la faccia "iface"
      tsegnato = LocPi(3) / (LocPi(3) - LocPj(3))
      PXLoc(:) = LocPi(:) + (LocPj(:) - LocPi(:)) * tsegnato
      PXLoc(3) = zero
!
      call LocalNormalCoordinates (PXLoc, csi, iface)
!
      fkod = BoundaryFace(iface)%nodes - 2
!
      if (IsPointInternal(fkod, csi)) then    !PX è interno alla faccia e quindi il punto
                                                !Pj da' contributo all'integrale di contorno

!

        !*******  Distance between Pi and Px  ***********************************
        PiPxdist = zero
        do SD = 1,SPACEDIM
          deltaPiPx = PXLoc(SD) - LocPi(SD)
          PiPxdist = PiPxdist + deltaPiPx * deltaPiPx
        end do
        PiPxdist = Dsqrt(PiPxdist)
!
        !*******  Normalised distance and related functions   *******************
        delta_alpha = BoundIntegralTab(ipk, 4)
        rob = PiPxdist / Domain%h
!
        if (Abs(rob - robpre) > 0.001) then
!
          ivalue = zero
          call InterpolateTable (rob, nicols, icol, ivalue)
!
          JW3ro2dA  = ivalue(1) * delta_alpha
          JdW3ro1dA = ivalue(2) * delta_alpha
          JdW3ro2dA = ivalue(3) * delta_alpha
          JdW3ro3dA = ivalue(4) * delta_alpha
          robpre = rob
!
        end if
!
        RLoccos(:) = BoundIntegralTab(ipk, 5:7)
        RZ = RLoccos(3)
        RZ2 = RZ * RZ
!
        !*******  Boundary volume integrals  ***********************************
        IntWdV = IntWdV + JW3ro2dA
        IntdWrm1dV = IntdWrm1dV + JdW3ro1dA
        IntGWZrm1dV = IntGWZrm1dV + RZ2 * JdW3ro1dA
        do SD = 1,SPACEDIM
          IntGWdV(SD) = IntGWdV(SD) + RLoccos(SD) * JdW3ro2dA
          do sdj = SD,SPACEDIM
            IntGWrRdV(SD, sdj) = IntGWrRdV(SD, sdj) + RLoccos(SD) * RLoccos(sdj) * JdW3ro3dA
          end do
        end do
      end if
    end if
  end do
!
! Completamento della matrice simmetrica IntGWrRdV(3, 3)
  IntGWrRdV(2,1) = IntGWrRdV(1,2)
  IntGWrRdV(3,1) = IntGWrRdV(1,3)
  IntGWrRdV(3,2) = IntGWrRdV(2,3)
!
  if (LocPi(3) == zpmin) then
    PXLoc(:) = LocPi(:)
    PXLoc(3) = zero
!
    call LocalNormalCoordinates (PXLoc, csi, iface)
!
    fkod = BoundaryFace(iface)%nodes - 2
!
    if (IsPointInternal(fkod, csi)) then    !La proiezione della particella è interna alla faccia
      IntWdV = half
    Else                                    !La proiezione della particella è esterna alla faccia
      IntWdV = zero
    end if
  end if
!
return
end subroutine ComputeBoundaryVolumeIntegrals_P0
!---split

