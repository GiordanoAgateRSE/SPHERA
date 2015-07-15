!cfile IsParticleInternal2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : IsParticleInternal2D
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
! Module purpose : Module check if a particle is internal to the domain 2D
!
! Calling routine: SetParticles
!
! Called routines: diagnostic
!
!************************************************************************************
!
  logical function IsParticleInternal2D (Nt, PX)
!
!.. assign modules
  use FILES_ENTITIES
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
  use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
  integer(4),intent(IN)                           :: Nt
  double precision,intent(IN),dimension(SPACEDIM) :: PX
!
!.. Local Scalars ..
  integer(4)       :: inizio,fine,iv, ni,n,n2
  double precision :: xa,za, xba,zba,  xi, zs
  character(len=lencard)  :: nomsub = "IsParticleInternal2D"
!
!.. Executable Statements ..
!
!.. initializations
!
  IsParticleInternal2D = .FALSE.
  ni = 0
  inizio = Tratto(Nt)%inivertex
  fine   = Tratto(Nt)%inivertex + Tratto(Nt)%numvertices - 2
!
!.. loops on the line vertices
!
  do iv = inizio, fine
!
!. checks for the storage limits
!
    if (iv > NumBVertices) then
      call diagnostic (arg1=10,arg2=2,arg3=nomsub)
    end if
!
!.. set the pointer to the current vertex
!
    n = BoundaryVertex(iv)
!
    if (n > NumVertici) then
      call diagnostic (arg1=10,arg2=3,arg3=nomsub)
    end if
!
!.. set the coordinates of the current vertex
!
    xa = Vertice(1,n)
    za = Vertice(3,n)
!
!. set the next vertex and its coordinates
!
    n2 = BoundaryVertex(iv+1)    
    xba = Vertice(1,n2)
    zba = Vertice(3,n2)
!
!.. the current segment is not vertical
!
    if ( abs(xa-xba) >= xyz_tolerance ) then
!
!.. the current segment is not horizontal
!
      if ( abs(za-zba) >= xyz_tolerance ) then
!
!.. evaluates the x coordinate of the particle projection on the segment along X
!
        xi = xa + (xba - xa) * (px(3) - za)/(zba - za)
!
      else
!
!.. the segment is horizontal: the X value of the mean point is assumed
!
          xi = half * (xa+xba)
!
      end if
    else
!
!.. the segment is vertical: the X value of the vertices is assumed
!
       xi = xa
!
    end if
!
!.. order the vertices in order to have the first one with the lower Z value
!
    if ( za > zba ) then
      zs  = za
      za  = zba
      zba = zs
    end if
!
!.. the Z value of the particle is inside the segment Z values
!
    if ( (PX(3)- za) > xyz_tolerance .AND. (PX(3) - zba) < xyz_tolerance ) then
!
      if ( xi > PX(1) ) then
        ni = ni + 1
      end if
    end if
!

  end do
!
  if ( MOD(ni,2) == 1 ) then
    IsParticleInternal2D = .TRUE.
  end if
!
  return
  end function IsParticleInternal2D
!---split

