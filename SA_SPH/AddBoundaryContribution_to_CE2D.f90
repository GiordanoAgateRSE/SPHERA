!cfile AddBoundaryContribution_to_CE2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : AddBoundaryContribution_to_CE2D
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
! Calling routine: Loop_Irre_2D
!
! Called routines: /
!
!************************************************************************************
!
  subroutine AddBoundaryContribution_to_CE2D (npi, IntNcbs, BCrodivV)  !, Ncbs
!
!Computes boundary contributions to rodivV
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
  use Diagnostic_MODULE
  use BoundIntegralTab_Module
!
!.. Implicit Declarations ..
  implicit none
!
!  double precision,parameter :: eps = 0.005d0 !AdM 15-10-08
!
!.. Formal Arguments ..
  integer(4),      intent(IN)      :: npi
!  integer(4),      intent(IN)      :: Ncbs
  integer(4),      intent(IN)      :: IntNcbs
  double precision,intent(INOUT)   :: BCrodivV
!
!.. Local Scalars ..
  integer(4)       :: pd, icbs, iside, sidestr, ibdt, ibdp
  double precision :: IntWds, roi, vin   !,IntWdV,  xpi, ypi,xpmin, xpmax, interlen,  etalocal !AdM 15-10-08
  character(4)     :: strtype
  type (TyBoundarySide) :: RifBoundarySide
!
!.. Local Arrays ..
  double precision,dimension(1:PLANEDIM)    :: IntLocXY
  integer(4),      dimension(1:PLANEDIM)    :: acix
  double precision,dimension(1:PLANEDIM)    :: nnlocal
  double precision,dimension(1:PLANEDIM)    :: Dvel
!
!.. Executable Statements ..
!
  acix(1) = 1        !active coordinate indexes
  acix(2) = 3
!
  BCrodivV = zero
!
  if (IntNcbs <= 0) return
!
  roi = pg(npi)%dens
  ibdt = BoundaryDataPointer(3,npi)
!
  do icbs = 1, IntNcbs
!
    ibdp = ibdt + icbs - 1
    IntLocXY(1:PLANEDIM) = BoundaryDataTab(ibdp)%LocXYZ(1:PLANEDIM)
    iside = BoundaryDataTab(ibdp)%CloBoNum
!
    RifBoundarySide = BoundarySide(iside)
    sidestr = RifBoundarySide%stretch
    strtype = Tratto(sidestr)%tipo
!
    if (strtype == "fixe" .OR. strtype == "tapi" .OR. strtype == "velo" .OR. strtype == "flow" .OR. strtype == "sour") then 
!
      IntWdS = BoundaryDataTab(ibdp)%BoundaryIntegral(1)
!      IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(3) !AdM 15-10-08
      vin = zero
!
      do pd = 1, PLANEDIM
        nnlocal(pd) = RifBoundarySide%T(acix(pd), acix(2))
      end do
!
      select case (strtype)
!
        case ("fixe")
          do pd = 1, PLANEDIM
            Dvel(pd) = pg(npi)%var(acix(pd))
          end do
          vin = Dvel(1) * nnlocal(1) + Dvel(2) * nnlocal(2)
!          etalocal = eps * Domain%h !AdM 15-10-08
!          SIntWds2 = IntWdV / (ypi + etalocal) !AdM 15-10-08
!
        case ("tapi")
          do pd = 1, PLANEDIM
            Dvel(pd) = pg(npi)%var(acix(pd))-RifBoundarySide%velocity(acix(pd))
          end do
          vin = Dvel(1) * nnlocal(1) + Dvel(2) * nnlocal(2)
!          etalocal = eps * Domain%h !AdM 15-10-08
!          SIntWds2 = IntWdV / (ypi + etalocal) !AdM 15-10-08
!
        case ("velo", "flow", "sour")
!AA404 sub
          if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
            pg(npi)%koddens = 2
            pg(npi)%densass = roi
            BCrodivV = zero
          end if
          return
!
      end select
!
        !****  Boundary contribution to rodivV  ****
!
!      BCrodivV = BCrodivV + vin * roi * ( IntWds + SIntWds2 ) !AdM 15-10-08
      BCrodivV = BCrodivV + two * vin * roi * IntWdS
!
    end if 
!
  end do
!
  return
  end subroutine AddBoundaryContribution_to_CE2D
!---split

