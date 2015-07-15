!cfile IsParticleInternal3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : IsParticleInternal3D
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
! Module purpose : Module check if a particle is internal to the domain 3D
!
! Calling routine: SetParticles
!
! Called routines: LocalNormalCoordinates
!                  IsPointInternal
!
!************************************************************************************
!
Logical Function IsParticleInternal3D ( mib, PX, IsopraS )
!
!Checks if point Px() is internal to the perimeter mib;
!in the affirmative returns 'true'; otherwise returns 'false'
!The perimeter can be both convex and concave !!
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
integer(4),parameter       :: intxy = 3
double precision,parameter :: eps = 0.001d0
!
!.. Formal Arguments ..
integer(4),      intent(IN)   :: mib
double precision,intent(IN),dimension(SPACEDIM) :: PX
integer(4),      intent(IN)   :: IsopraS
!
!.. Local Scalars ..
integer(4)       :: kf, nf, i, j, sd, nnodes, norig
integer(4)       :: Nints, IntSotto, IntSopra, fkod
double precision :: tpar
double precision,dimension(SPACEDIM) :: P1, Pint, LPint
double precision,dimension(3)        :: csi

!Dynamic Array
double precision,dimension(Tratto(mib)%numvertices) :: XYInts
!
!.. External Routines ..
logical, external :: IsPointInternal
!
!.. Executable Statements ..
!
 Nints    = 0
 IntSotto = 0
 IntSopra = 0
 IntSopra = IsopraS
 IsParticleInternal3D = .FALSE.

 do kf = Tratto(mib)%iniface, Tratto(mib)%iniface + Tratto(mib)%numvertices - 1

    nf     = BFaceList(kf)
    nnodes = 4
    if ( BoundaryFace(nf)%Node(4)%name <= 0 ) nnodes = 3
    norig  = nnodes                                      !nodo origine del sistema locale

    do sd = 1, SPACEDIM
       P1(sd) = Vertice(sd,BoundaryFace(nf)%Node(norig)%name)
    end do
    tpar = zero
    do sd = 1, SPACEDIM
       tpar = tpar + BoundaryFace(nf)%T(sd, 3) * (P1(sd) - PX(sd))
    end do

    if ( Abs(BoundaryFace(nf)%T(3, 3)) > eps ) then

        tpar = tpar / BoundaryFace(nf)%T(3, 3)

        do sd = 1, SPACEDIM                          !Pint()= coordinate globali del punto di intersezione
           Pint(sd) = PX(sd)
        end do
        Pint(3) = Pint(3) + tpar

        LPint = zero
        do sd = 1, PLANEDIM                          !LPint()= coordinate locali del punto di intersezione
           LPint(sd) = zero
           do j = 1, SPACEDIM
              LPint(sd) = LPint(sd) + BoundaryFace(nf)%T(j, sd) * (Pint(j) - P1(j))
           end do
        end do

        call LocalNormalCoordinates ( LPint, csi, nf )

        fkod = nnodes - 2
        if ( IsPointInternal ( fkod, csi ) ) then    !Il punto di intersezione è interno alla faccia;
                                                     !cioè la faccia interseca la verticale per Px.
                                                     !Memorizza coordinata z di intersezione
            Nints = Nints + 1
            XYInts(Nints) = Pint(3)
        end if
    end if

 end do

 if ( Nints > 0 ) then
    do i = 1, Nints
       if ( XYInts(i) <= PX(intxy) ) then
          IntSotto = IntSotto + 1
       Else
          IntSopra = IntSopra + 1
       end if
    end do
    if ( Mod(IntSotto,2) == 1 .AND. Mod(IntSopra,2) == 1 ) then
        IsParticleInternal3D = .TRUE.
    end if
 end if

return
End Function IsParticleInternal3D
!---split

