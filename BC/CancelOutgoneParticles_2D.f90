!cfile CancelOutgoneParticles_2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : CancelOutgoneParticles_2D
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
! Calling routine: Loop_Irre_2D
!
! Called routines: 
!
!************************************************************************************
!
  subroutine CancelOutgoneParticles_2D
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
  integer(4)        :: ios, v1, pd, npi, isi, v2
!  double precision  :: yp,xcross,ycross,m,mv
  double precision  :: detV1PnewV2, detV1PoldV2, detPoldV1Pnew, detPoldV2Pnew
!
!.. local arrays
  integer(4),       dimension(1:PLANEDIM)       :: acix
  double precision, dimension(1:PLANEDIM)       :: OP1,op2
!  double precision, dimension(1:PLANEDIM)       :: Onn
  double precision, dimension(1:PLANEDIM)       :: Plocalnew,Plocalold
!  double precision, dimension(1:PLANEDIM)       :: OP1P
!
!.. Executable Statements ..
!
!.. 2D scheme: uses the X and Z axes
!
  acix(1)=1  
  acix(2)=3
!
!.. all the opened boundary sides are considered and the relative position of all the particles 
!.. is checked
!
  do ios = 1, NumOpenSides
!
    isi = OpenSide(ios)
    v1 = BoundarySide(isi)%Vertex(1)
    v2 = BoundarySide(isi)%Vertex(2)
!
!.. load the reference vertex coordinates and local side versor (tangent and normal)
! 
    do pd = 1, PLANEDIM
!
      OP1(pd) = Vertice(acix(pd),v1)
      OP2(pd) = Vertice(acix(pd),v2)
!      Onn(pd) = BoundarySide(isi)%T(acix(pd),3)
!
    end do
!
!$omp parallel do default(none) &
!$omp private(npi,ios,Plocalnew,Plocalold,detV1PnewV2,detV1PoldV2,detPoldV1Pnew,detPoldV2Pnew) &
!$omp shared(nag,pg,acix,OP1,OP2,OpCount)
!
!.. loops on all the particles
!
    do npi = 1,nag
!
      if (pg(npi)%cella == 0) cycle
!
!.. load the particle coordinates
!
      Plocalnew(1:planedim) = pg(npi)%coord(acix(1:planedim))
      Plocalold(1:planedim) = pg(npi)%CoordOld(acix(1:planedim))
!
!.. evaluates the distance components in X and Z directions between the vertex V1 and the current particle
!
!      OP1P(:) = Plocalnew(:) - OP1(:)
!
!.. evaluates the normal distance component with respect the opened boundary side
!
!      yp = OP1P(1) * Onn(1) + OP1P(2) * Onn(2)
!
!.. evaluates the determinant of the outside triplet V1-Pnew-V2
!
      detV1PnewV2 = (Plocalnew(2)-op1(2)) * (op2(1)-op1(1)) - (plocalnew(1)-op1(1)) * (op2(2)-op1(2))
!
!.. if the normal component is less than zero, the particle might be out of the boundary side, since
!.. the reference normal is oriented inside the domain
!
      if (detV1PnewV2 <= zero) then
!      if (yp <= zero) then
!
!.. so the crossing point between the last path of the particle and the boundary segment is tested
!
        detV1PoldV2 = (Plocalold(2)-op1(2)) * (op2(1)-op1(1)) - (plocalold(1)-op1(1)) * (op2(2)-op1(2))
        if (sign(one,detV1PnewV2) == sign(one,detV1PoldV2)) cycle
        detPoldV1Pnew = (op1(2)-plocalold(2)) * (plocalnew(1)-plocalold(1)) - (op1(1)-plocalold(1)) * (plocalnew(2)-plocalold(2))
        detPoldV2Pnew = (op2(2)-plocalold(2)) * (plocalnew(1)-plocalold(1)) - (op2(1)-plocalold(1)) * (plocalnew(2)-plocalold(2))
        if (sign(one,detPoldV1Pnew) == sign(one,detPoldV2Pnew)) cycle
!
        OpCount(pg(npi)%imed) = OpCount(pg(npi)%imed) + 1 
        pg(npi)%cella = -1
!
!        m  = (Plocalnew(2) - Plocalold(2))/(Plocalnew(1) - Plocalold(1))
!        mv = (op2(2) - op1(2))/(op2(1) - op1(1))
!        xcross = (op1(2) - plocalold(2) + m * plocalold(1) - mv * op1(1))/ (m - mv)
!        ycross = plocalold(2) + m * (xcross - plocalold(1))                         
      end if
    end do
!
!$omp end parallel do
!
  end do
!
!.. the particle array is compacted; if the npi-th particle is zero the last particle is moved into the npi-th array location
!
!!!!!  npi = 0
!!!!!  do while (nag > 0 .AND. npi < nag)
!!!!!    npi = npi + 1
!!!!!    if (pg(npi)%cella == 0) then
!!!!!!.. the counter of the removed particles is increased by 1
!!!!!      OpCount(pg(npi)%imed) = OpCount(pg(npi)%imed) + 1 
!!!!!      pg(npi) = pg(nag)
!!!!!      pg(nag)%cella = 0
!!!!!      pg(nag) = PgZero
!!!!!      nag = nag - 1
!!!!!      npi = npi - 1
!!!!!    end if
!!!!!  end do
!
  return
  end subroutine CancelOutgoneParticles_2D
!---split

