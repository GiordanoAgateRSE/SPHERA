!cfile AddBoundaryContribution_to_CE3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : AddBoundaryContribution_to_CE3D
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
! Module purpose : Module to compute boundary contributions to rodivV, in 3D geometry,
!                  relative to particle npi
!                  Additional term due to a linear distribution of normal velocity
!                  has been introduced
!
! Calling routine: Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
!
  subroutine AddBoundaryContribution_to_CE3D (npi, Ncbf, BCtorodivV)
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
  integer(4),      intent(IN)    :: npi
  integer(4),      intent(IN)    :: Ncbf
  double precision,intent(INOUT) :: BCtorodivV
!
!.. Local Scalars ..
  integer(4)       :: sd, sdj, icbf, iface, ibdt, ibdp
  integer(4)       :: stretch        
  double precision :: roi, scaprod
  character(4)     :: boundtype            
!
!.. Local Arrays ..
  double precision,dimension(1:SPACEDIM) :: vb
  double precision,dimension(1:SPACEDIM) :: vi
  double precision,dimension(1:SPACEDIM) :: dvij
  double precision,dimension(1:SPACEDIM) :: LocPi, LocDvij
!
!.. Executable Statements ..
!
  roi = pg(npi)%dens
  vi(:) = pg(npi)%var(:)
  BCtorodivV = zero
!
  if (Ncbf <= 0) return
!
  ibdt = BoundaryDataPointer(3,npi)
!
  do icbf = 1, Ncbf
!
    ibdp = ibdt + icbf - 1
    LocPi(1:SPACEDIM) = BoundaryDataTab(ibdp)%LocXYZ(1:SPACEDIM)
    iface = BoundaryDataTab(ibdp)%CloBoNum
!
    stretch = BoundaryFace(iface)%stretch
    boundtype = Tratto(stretch)%tipo
!
    if (boundtype == "fixe" .OR. boundtype == "tapi") then
!
      if (LocPi(3) > zero) then          !La faccia "iface" interagisce con la particella Pi
!
        vb(:) = BoundaryFace(iface)%velocity(:)
        dvij(:) = two * (vi(:) - vb(:))
!
!!!        scaprod = zero
!!!        do sd = 1, SPACEDIM
!!!          Locdvij(sd) = zero  !componenti locali del vettore 2*(vi-vb)
!!!          do sdj = 1, SPACEDIM
!!!            Locdvij(sd) = Locdvij(sd) + dvij(sdj) * BoundaryFace(iface)%T(sdj,sd)
!!!          end do
!!!          scaprod = scaprod + Locdvij(sd) * BoundaryDataTab(ibdp)%BoundaryIntegral(3+sd)
!!!        end do
        scaprod = zero
        sd = 3
        Locdvij(sd) = zero  !componenti locali del vettore 2*(vi-vb)
        do sdj = 1, SPACEDIM
          Locdvij(sd) = Locdvij(sd) + dvij(sdj) * BoundaryFace(iface)%T(sdj,sd)
        end do
        scaprod = scaprod + Locdvij(sd) * BoundaryDataTab(ibdp)%BoundaryIntegral(3+sd)        
!
        !*******  Boundary contribution to divV  ***********************************
        BCtorodivV = BCtorodivV + roi * scaprod
!
      end if
!
    else if (boundtype == "velo" .or. boundtype == "flow" .or. boundtype == "sour") then
!
!AA404 sub
          if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
        pg(npi)%koddens = 2
        pg(npi)%densass = roi
        BCtorodivV = zero
      end if
      return
!
    end if
!
  end do
!
  return
  end subroutine AddBoundaryContribution_to_CE3D
!---split
