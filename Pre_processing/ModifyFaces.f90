!cfile ModifyFaces.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ModifyFaces
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
! Module purpose : Module 
!
! Calling routine: Gest_Input
!
! Called routines: 
!
!************************************************************************************
!
subroutine ModifyFaces ( NumberEntities )
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
integer(4),intent(IN),dimension(20) :: NumberEntities
!
!.. Local Scalars ..
integer(4)       :: n,i,new
double precision :: d13, d24
!
!.. Executable Statements ..
!
!genera i nuovi triangoli, divisione lungo diagonale minima
 new = NumberEntities(11)
 do n = 1, NumberEntities(11)
    if ( BoundaryFace(n)%Node(4)%name == 0 ) cycle
    new = new + 1
    d13 = zero
    d24 = zero
    do i =1, SPACEDIM
       d13 = d13 + (Vertice(i,BoundaryFace(n)%Node(1)%name)-Vertice(i,BoundaryFace(n)%Node(3)%name)) * &
                   (Vertice(i,BoundaryFace(n)%Node(1)%name)-Vertice(i,BoundaryFace(n)%Node(3)%name))
       d24 = d24 + (Vertice(i,BoundaryFace(n)%Node(2)%name)-Vertice(i,BoundaryFace(n)%Node(4)%name)) * &
                   (Vertice(i,BoundaryFace(n)%Node(2)%name)-Vertice(i,BoundaryFace(n)%Node(4)%name))
    end do
    if ( d13 < d24 ) then
       BoundaryFace(new) = BoundaryFace(n)
       BoundaryFace(n)%Node(4)%name   =-BoundaryFace(n)%Node(4)%name
       BoundaryFace(new)%Node(2)%name = BoundaryFace(new)%Node(3)%name
       BoundaryFace(new)%Node(3)%name = BoundaryFace(new)%Node(4)%name
       BoundaryFace(new)%Node(4)%name =-99999999
    else
       BoundaryFace(new) = BoundaryFace(n)
       i                 = BoundaryFace(n)%Node(1)%name
       BoundaryFace(n)%Node(1)%name   = BoundaryFace(n)%Node(2)%name
       BoundaryFace(n)%Node(2)%name   = BoundaryFace(n)%Node(3)%name
       BoundaryFace(n)%Node(3)%name   = BoundaryFace(n)%Node(4)%name
       BoundaryFace(n)%Node(4)%name   =-i
       BoundaryFace(new)%Node(3)%name = BoundaryFace(new)%Node(4)%name
       BoundaryFace(new)%Node(4)%name =-99999999
    end if
 end do

return
end subroutine ModifyFaces
!---split

