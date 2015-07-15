!cfile DefineLocalSystemVersors.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : DefineLocalSystemVersors
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
! Module purpose : Module to define versor components (cosine directors) of the
!                  local reference system
!
! Calling routine: DefineBoundaryFaceGeometry3D
!
! Called routines: Vector_Product
!
!************************************************************************************
!
!AA504 sub
subroutine DefineLocalSystemVersors (Nf)
!Defines versor components (cosine directors) of the local reference system
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!AA504
use ALLOC_Module
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
!AA504 rm part
integer(4) :: Nf
!
!.. Local Scalars ..
integer(4) :: i, j, n, nnodes, sidek, refnode, nod
double precision :: U1len, W12len, LocX
!
!Dynamic Arrays
double precision, dimension(SPACEDIM) ::  RR, ss, nnlocal, U1, U2, W12
!
!integer(4), dimension(2, 3) :: ((nindex(i,j),j=1,3),i=1,2) = (/ 1,2,3, 1,3,4 /)
integer(4), dimension(2, 3) :: nindex ! = (/ 1,1, 2,3, 3,4 /) modifica per compatibilita xlf90
!
!
!.. Executable Statements ..
!
! modifica per compatibilita xlf90
    nindex(1, 1) = 1
    nindex(1, 2) = 2
    nindex(1, 3) = 3
    nindex(2, 1) = 1
    nindex(2, 2) = 3
    nindex(2, 3) = 4
! fine modifica
!
 nnodes = 4
 if ( BoundaryFace(Nf)%Node(4)%name <= 0 ) nnodes = 3
 BoundaryFace(Nf)%nodes = nnodes
 sidek = nnodes - 2                   !sidek=1 (triangle),  sidek=2 (parallelogram)
!
 do n = 1, nnodes                     !copies face node global coordinates
    nod = BoundaryFace(nf)%Node(n)%name
    do i = 1, SPACEDIM
       BoundaryFace(nf)%Node(n)%GX(i) = Vertice(i,nod)
    end do
!write(97,*) nod, BoundaryFace(nf)%Node(n)%GX(:)
 end do

!defines main sides components
 U1(1:SPACEDIM) = BoundaryFace(Nf)%Node(nindex(sidek, 1))%GX(1:SPACEDIM) - BoundaryFace(Nf)%Node(nindex(sidek, 3))%GX(1:SPACEDIM)
 U2(1:SPACEDIM) = BoundaryFace(Nf)%Node(nindex(sidek, 2))%GX(1:SPACEDIM) - BoundaryFace(Nf)%Node(nindex(sidek, 3))%GX(1:SPACEDIM)

 U1len = Dsqrt( U1(1)*U1(1) + U1(2)*U1(2) + U1(3)*U1(3) ) !length of side U1

 RR(1:SPACEDIM) = U1(1:SPACEDIM) / U1len      !computes versor components (director cosines) of side W12

 call Vector_Product ( U1, U2, W12, SPACEDIM ) !computes vector product W12=W12xU2

 W12len = Dsqrt( W12(1)*W12(1) + W12(2)*W12(2) + W12(3)*W12(3) ) !length of vector W12

 BoundaryFace(Nf)%Area = W12len / float(3 - sidek)        !define the area of the face nf (denominator=2 do triangle, =1 do parallelogram)

 nnlocal(1:SPACEDIM) = W12(1:SPACEDIM) / W12len      !computes versor components (director cosines) of the normal to face nf

 call Vector_Product ( nnlocal, RR, ss, SPACEDIM )   !computes vector product ss=rrxnn; ss is the versor of the third local axis

!defines cosine director matrix of face nf
 BoundaryFace(Nf)%T(1:SPACEDIM, 1) = RR(1:SPACEDIM)
 BoundaryFace(Nf)%T(1:SPACEDIM, 2) = ss(1:SPACEDIM)
 BoundaryFace(Nf)%T(1:SPACEDIM, 3) = nnlocal(1:SPACEDIM)

 refnode = nnodes        !defines local coordinates of face nodes
!
 do n = 1, nnodes - 1
    do i = 1, PLANEDIM
       LocX = zero
       do j = 1, SPACEDIM
          LocX = LocX + BoundaryFace(Nf)%T(j, i) * (BoundaryFace(Nf)%Node(n)%GX(j) - BoundaryFace(Nf)%Node(refnode)%GX(j))
       end do
       BoundaryFace(nf)%Node(n)%LX(i) = LocX
    end do
 end do
!
 do i = 1, PLANEDIM
    BoundaryFace(nf)%Node(refnode)%LX(i) = zero
 end do
!
return
end subroutine DefineLocalSystemVersors
!---split

