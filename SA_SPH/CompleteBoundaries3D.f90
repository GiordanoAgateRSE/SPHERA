!cfile CompleteBoundaries3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : CompleteBoundaries3D
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
! Module purpose : Module to Compute and store complemetary data starting
!                  from assigned data
!
! Calling routine: Gest_Input
!
! Called routines: 
!
!************************************************************************************
!
subroutine CompleteBoundaries3D
!Computes and stores complemetary data starting from assigned data
!
!.. assign modules
use GLOBAL_MODULE
use AdM_User_Type
use ALLOC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Scalars ..
integer(4) :: Kf, Nt, Nf
!
!.. Executable Statements ..
!
!Boundary face list
!
 Kf = 0
 BFaceList(1:NumFacce) = 0
!
  do Nt = 1, NumTratti
!
    Tratto(Nt)%inivertex   = 0
    Tratto(Nt)%numvertices = 0
!
    do Nf = 1, NumFacce
!
        if ( BoundaryFace(Nf)%stretch == Nt ) then
!
            Kf = Kf + 1
            if ( Tratto(Nt)%iniface == 0 ) then !Stores initial face index
               Tratto(Nt)%iniface = Kf
!write(97,"(a,i5,a,i5)") "tratto",nt,"  iniface=",tratto(nt)%iniface
            end if
            BFaceList(Kf) = Nf
!write(97,"(a,i5,a,i5)")  "BFaceList=(",kf,") =",Nf
            Tratto(Nt)%numvertices = Tratto(Nt)%numvertices + 1
        end if
!
    end do
!write(97,"(a,i5,a,i5)") "tratto",nt,"  numvertices=",tratto(nt)%numvertices
 end do
!
!write(97,*) "Tratto%NumVertices"
!write(97,"(10i8)") Tratto(1:NumTratti)%NumVertices
return
end subroutine CompleteBoundaries3D
!---split

