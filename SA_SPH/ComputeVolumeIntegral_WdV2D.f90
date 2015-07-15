!cfile ComputeVolumeIntegral_WdV2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ComputeVolumeIntegral_WdV2D
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
! Module purpose : Module for Computing the integral of WdV extented to the volume
!                  delimited by the influence circle (radius=2h) of the particle i,
!                  whose local coordinates are xpi=LocXY(1, icbs) and
!                  ypi=LocXY(2, icbs), and the adjacent boundary side icbs
 
!
! Calling routine: ComputeBoundaryDataTab
!
! Called routines: 
!
!************************************************************************************
!
subroutine ComputeVolumeIntegral_WdV2D (icbs, Ncbslocal, Cloboside, LocXY, BoundarySide, &
                                        xpmin, xpmax, interlen, VIntWdV)

!Computes the integral of WdV extented to the volume delimited by the influence circle (radius=2h)
!of the particle i, !whose local coordinates are xpi=LocXY(1, icbs) and ypi=LocXY(2, icbs),
!and the adjacent boundary side icbs
!
!.. assign modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
integer(4),      parameter :: Nalfadiv = 10
!double precision,parameter :: eps=0.01d0 
double precision,parameter :: eps=0.05d0
!
!.. Formal Arguments ..
integer(4),      intent(IN)    :: icbs
integer(4),      intent(IN)    :: Ncbslocal
integer(4),      intent(IN),dimension(1:MAXCLOSEBOUNDSIDES)            :: Cloboside
double precision,intent(IN),dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES) :: LocXY
type (TyBoundarySide),intent(IN),dimension(1:NumBSides)  :: BoundarySide
double precision,intent(IN)    :: xpmin
double precision,intent(IN)    :: xpmax
double precision,intent(IN)    :: interlen
double precision,intent(INOUT) :: VIntWdV
!
!.. Local Scalars ..
integer(4)       :: jcbs, ndiv, ipt
double precision :: xpi, ypi, yplimite, ypj, &             !ypiq, 
                    angle, dalfarif, Intalfa, tanalfa, ris, dalfa, &
                    alfaA, alfaB, alfa, csiPA, etaPA, csiPB, etaPB
character(len=lencard) :: nomsub = "ComputeVolumeIntegral_WdV2D"
!
! External functions and subrotuines
double precision,    external    :: WIntegr
!
!.. Executable Statements ..
!
  dalfarif = PIGRECO / Nalfadiv
!
  VIntWdV = zero
  if (interlen <= zero) return
  yplimite = eps * Domain%h
  xpi = LocXY(1, icbs)
  ypi = LocXY(2, icbs)
!
  if (ypi >= yplimite) then
!    ypiq = ypi * ypi
    csiPA = xpmin - xpi
    etaPA = ypi
    csiPB = xpmax - xpi
    etaPB = ypi
    alfaA = Atan2(csiPA, etaPA)
    alfaB = Atan2(csiPB, etaPB)
    Intalfa = alfaB - alfaA
    ndiv = Int(intalfa / dalfarif + half)
    if ( ndiv < 2 ) ndiv = 2
    dalfa = intalfa / ndiv
    alfa = alfaA - half * dalfa
    do ipt = 1, ndiv
      alfa = alfa + dalfa
      tanalfa = Tan(alfa)
      ris = ypi * Dsqrt(one + tanalfa * tanalfa)
      VIntWdV = VIntWdV + WIntegr(ris, Domain%h) * dalfa
    end do
!
  else if (ypi < yplimite .and. Ncbslocal == 1) then        !Il lato vicino e' solo uno cioè icbs
    if (xpmin >= xpi .and. xpi <= xpmax) then      
      VIntWdV = half
    else    
      VIntWdV = zero
    end if
!         
  else if (ypi < yplimite .and. Ncbslocal == 2) then
!                                           !I lati vicini alla particella "i" possono essere due
    jcbs = Ncbslocal + 1 - icbs             !indice del secondo lato vicino
    ypj = LocXY(2, jcbs)                    !distanza della particella npi dal secondo lato
    if (ypj <= yplimite) then               !la particella è molto vicina al vertice comune ai due lati adiacenti
      angle = zero
      if (BoundarySide(Cloboside(1))%previous_side == Cloboside(2)) then
        angle = BoundarySide(Cloboside(1))%angle
      else if (BoundarySide(Cloboside(2))%previous_side == Cloboside(1)) then
        angle = BoundarySide(Cloboside(2))%angle
      else
        write (nout,'(a,2i10)') 'ERROR!! Sides not consecutive',Cloboside(1),Cloboside(2)
        call diagnostic (arg1=8,arg2=4,arg3=nomsub)
      end if
!
      VIntWdV = half * (angle + PIGRECO) / (two * PIGRECO)
!
    else    !la particella è molto vicina solo al lato icbs
      if (xpmin >= xpi .and. xpi <= xpmax) then      
        VIntWdV = half
      else    
        VIntWdV = zero
      end if         
    end if
  end if
!
return
end subroutine ComputeVolumeIntegral_WdV2D
!---split

