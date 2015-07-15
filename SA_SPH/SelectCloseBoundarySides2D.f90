!cfile SelectCloseBoundarySides2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : SelectCloseBoundarySides2D
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
! Module purpose : Module for Selection from close boundary sides those that really
!                  give contribution to the equations of particle 'npi'.
!
! Calling routine: ComputeBoundaryDataTab
!
! Called routines: 
!
!************************************************************************************
!
  subroutine SelectCloseBoundarySides2D (npi, Ncbs, Cloboside, LocXY, IntNcbs, Intboside, IntLocXY)
!
!Returns:
!IntNcbs                     = Number of close boundary sides which give contribution (= 0, 1, 2)
!Intboside(1:IntNcbs)        = List of close boundary sides which give contribution
!IntLocXY(1:PLANEDIM,1:Ncbs) = Local coordinates of particle np with respect each boundary side
!                                which give contribution                                
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Local Parameters ..
  double precision, parameter :: eps = 0.501d0
!
!.. Formal Arguments ..
  integer(4),                                                  intent(in)  :: npi
  integer(4),                                                  intent(in)  :: Ncbs
  integer(4),                                                  intent(out) :: IntNcbs
  integer(4),       dimension(1:MAXCLOSEBOUNDSIDES),           intent(in)  :: Cloboside
  integer(4),       dimension(1:MAXCLOSEBOUNDSIDES),           intent(out) :: Intboside
  double precision, dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES),intent(in)  :: LocXY
  double precision, dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES),intent(out) :: IntLocXY
!
!.. Local Scalars ..
  integer(4)       :: icbs, isi, nt
  double precision :: sidelen, yxpmin
!
!.. Executable Statements ..
!
!.. initializations
!
  IntNcbs       = 0
  Intboside(:)  = 0
  IntLocXY(:,:) = zero
!
!.. loops on the close boundaries previously found
!
  do icbs = 1, Ncbs
!
    isi = Cloboside(icbs)
    sidelen = BoundarySide(isi)%length
    nt = BoundarySide(isi)%stretch
!
!.. the boundary is of type source:
!.. the side is considered only if the particle is in front of the side itself and
!.. inside the domain or very closest the side (< 0.5*Domain%dd)
!
    if (Tratto(nt)%tipo == "sour") then
!
      yxpmin = -eps * Domain%dd
      if (LocXY(2, icbs) > yxpmin .AND. LocXY(1, icbs) > zero .And. LocXY(1, icbs) < sidelen) then
        IntNcbs = IntNcbs + 1
        Intboside(IntNcbs) = Cloboside(icbs)
        IntLocXY(:, IntNcbs) = LocXY(:, icbs)
      end if
!
!.. the boundary is type outlet velocity
!
    else if (Tratto(nt)%tipo == "velo" .or. Tratto(nt)%tipo == "flow") then
!
      yxpmin = zero
      if (LocXY(2, icbs) > yxpmin .AND. LocXY(1, icbs) > zero .AND. LocXY(1, icbs) < sidelen) then
        IntNcbs = IntNcbs + 1
        Intboside(IntNcbs) = Cloboside(icbs)
        IntLocXY(:, IntNcbs) = LocXY(:, icbs)
!.. particella corrente vicino a contorno di uscita (vedi erosione crit_erosion)
        pg(npi)%CloseBcOut = 1
!
      end if
!
!.. the boundary is type outlet open
!
    else if (Tratto(nt)%tipo == "open") then
!
      yxpmin = zero
      if (LocXY(2, icbs) > yxpmin .AND. LocXY(1, icbs) > zero .AND. LocXY(1, icbs) < sidelen) then
        IntNcbs = IntNcbs + 1
        Intboside(IntNcbs) = Cloboside(icbs)
        IntLocXY(:, IntNcbs) = LocXY(:, icbs)
!.. particella corrente vicino a contorno di uscita (vedi erosione crit_erosion)
        pg(npi)%CloseBcOut = 1
!
      end if
!
!.. the boundary is not a source velo flow or open
!
    else
!
      yxpmin = zero          
!
!.. the current particle is inside the domain and in front of an opened (??) boundary condition
!
      if (LocXY(2, icbs) > yxpmin ) then
!
!.. the closest boundary is accounted for the calculation of rhs contribution 
!
        IntNcbs = IntNcbs + 1
        Intboside(IntNcbs) = Cloboside(icbs)
        IntLocXY(:, IntNcbs) = LocXY(:, icbs)
!
      end if
!
    end if
!
  end do
!        
  return
  end subroutine SelectCloseBoundarySides2D
!---split

