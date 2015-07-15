!AA601 the whole subroutine
!cfile BC_wall_elements.f90
!************************************************************************************
!                          S P H E R A 5.0.3 
!
!                Smoothed Particle Hydrodinamics Code
!
!************************************************************************************
!
! Subroutine : wavy_inlet
!
! Creation   : Amicarelli A., 26Jan15; Main features (see Purpose)
!
!************************************************************************************
! Purpose : It provides a wavy flow at the inlet section: each layer is staggered 
!           by 0.5dx with respect to the previous and the following ones, which are instead aligned.
!           This is a numerical feature to reduce the SPH truncation error (necessary for DBSPH jets).
!
! Calling routine: GenerateSourceParticles_2D,GenerateSourceParticles_3D
! Called subroutines : / 
!
!************************************************************************************
subroutine wavy_inlet(i_inlet)

! Modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE

! Declarations
implicit none
!AA601!!! sub
integer(4),intent(in) :: i_inlet
integer(4) :: npart1,npart2,i
!AA601!!! sub
double precision :: Length1,Length2,rnd
double precision,dimension(3) :: cos_dir_1

!Initializations
Length1 = 0.d0
Length2 = 0.d0
npart1 = 0
npart2 = 0

! Statements
if (ncord==2) then
   npart1 = NumPartperLine(i_inlet)
   cos_dir_1(:) = BoundarySide(i_inlet)%T(:,1)
   else
! Length1 and Length 2 are the length scales of the inlet section: they are computed as the distance between the first and the last inlet vertices 
! and the third and the last inlet vertices, respectively.
! Particles are aligned with Plast-P1 and Plast-P3, where P1 is the first boundary vertex, ... and Plast the last boundary vertex. 
! In case of a triangular inlet, we have particles aligned with one direction: P3-P1.
! In case of a quadrilateral inlet, we have particles distributed along two directions: P4-P1 and P4-P3.
      do i=1,3
         Length1 = Length1 + (BoundaryFace(i_inlet)%Node(1)%GX(i) - BoundaryFace(i_inlet)%Node(BoundaryFace(i_inlet)%nodes)%GX(i))**2
         Length2 = Length2 + (BoundaryFace(i_inlet)%Node(3)%GX(i) - BoundaryFace(i_inlet)%Node(BoundaryFace(i_inlet)%nodes)%GX(i))**2
      end do
      Length1 = Dsqrt(Length1)
      Length2 = Dsqrt(Length2)
      npart1 = Int(Length1/Domain%dd+0.01d0)
      npart2 = Int(Length2/Domain%dd+0.01d0)
      npart1 = npart1*npart2 
      cos_dir_1(:) = BoundaryFace(i_inlet)%T(:,1)
   endif
! (ID particle = nag) indicates the last generated particle (of the on-going inlet section)
   select case (mod(itime_jet,4))
      case (1,3)  
         pg(nag)%coord(:) = pg(nag)%coord(:) - 0.25d0 * Domain%dd * cos_dir_1(:)
      case (2,0) 
         pg(nag)%coord(:) = pg(nag)%coord(:) + 0.25d0 * Domain%dd * cos_dir_1(:)
   end select
call random_number(rnd)
pg(nag)%coord(:) = pg(nag)%coord(:) + (two * rnd - one) * 0.1d0 * Domain%dd * cos_dir_1(:)

return
end subroutine wavy_inlet
!---split

