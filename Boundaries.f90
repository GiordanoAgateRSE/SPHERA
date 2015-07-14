!cfile BoundaryMassForceMatrix2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : BoundaryMassForceMatrix2D
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
! Module purpose : Module for generation of generalized boundary mass force matrix
!
! Calling routine: DefineBoundarySideGeometry2D
!
! Called routines: 
!
!************************************************************************************
!
subroutine BoundaryMassForceMatrix2D (T, RN, FiS, FiN) 

!Generates the generalised boundary mass force matrix RN(1 To SPACEDIM, 1 To SPACEDIM)
!on the base of the cosine matrix T(1 To SPACEDIM, 1 To SPACEDIM) and the parameter Fi
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
double precision                                   :: FiS, FiN 
double precision,dimension(1:SPACEDIM, 1:SPACEDIM) :: T
double precision,dimension(1:SPACEDIM, 1:SPACEDIM) :: RN
!
!.. Executable Statements ..
!
RN(1, 1) = FiS * T(1, 1) * T(1, 1) + FiN * T(3, 1) * T(3, 1)
RN(1, 3) = (FiS - FiN) * T(1, 1) * T(3, 1) 
RN(3, 1) = RN(1, 3)         
RN(3, 3) = FiS * T(3, 1) * T(3, 1) + FiN * T(1, 1) * T(1, 1)
!
return
!
end subroutine BoundaryMassForceMatrix2D
!---split

!cfile BoundaryMassForceMatrix3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : BoundaryMassForceMatrix3D
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
! Module purpose : Module to generate the boundary mass force matrix
!                  RN(1 To SPACEDIM, 1 To SPACEDIM) on the base of the cosine
!                  matrix T(1 To SPACEDIM, 1 To SPACEDIM) and the parameters Fi()
!
! Calling routine: DefineBoundaryFaceGeometry3D
!
! Called routines: MatrixTransposition
!                  MatrixProduct
!
!************************************************************************************
!
subroutine BoundaryMassForceMatrix3D ( T, RMF, Fi )
!Generates the boundary mass force matrix RN(1 To SPACEDIM, 1 To SPACEDIM)
!on the base of the cosine matrix T(1 To SPACEDIM, 1 To SPACEDIM) and the parameters Fi()
!
!.. assign modules
use GLOBAL_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
double precision,intent(INOUT),dimension(SPACEDIM,SPACEDIM) :: T
double precision,intent(INOUT),dimension(SPACEDIM,SPACEDIM) :: RMF
double precision,intent(INOUT),dimension(SPACEDIM)          :: Fi
!
!.. Local Scalars ..
integer(4) :: i
!
!.. Local Arrays ..
double precision, dimension(SPACEDIM,SPACEDIM) :: Diag, FiR, TTR
!
!.. Executable Statements ..
!
 Diag = zero
 do i = 1, SPACEDIM
   Diag(i, i) = Fi(i)
 end do
!
 call MatrixTransposition ( T, TTR, SPACEDIM, SPACEDIM )
!
 call MatrixProduct ( Diag, TTR, FiR, SPACEDIM, SPACEDIM, SPACEDIM )
!
 call MatrixProduct ( T, FiR, RMF, SPACEDIM, SPACEDIM, SPACEDIM )
!
return
!
end subroutine BoundaryMassForceMatrix3D
!---split

!cfile BoundaryPressureGradientMatrix3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : BoundaryPressureGradientMatrix3D
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
! Module purpose : Module to generate the pressure gradient matrix
!                  RRP(1 To SPACEDIM, 1 To SPACEDIM) on the base of the cosine
!                  matrix T(1 To SPACEDIM, 1 To SPACEDIM) and the parameters Psi()
!
! Calling routine: DefineBoundaryFaceGeometry3D
!
! Called routines: MatrixTransposition
!                  MatrixProduct
!
!************************************************************************************
!
subroutine BoundaryPressureGradientMatrix3D ( T, RGP, Psi )
!Generates the pressure gradient matrix RRP(1 To SPACEDIM, 1 To SPACEDIM) on the base of
!the cosine matrix T(1 To SPACEDIM, 1 To SPACEDIM) and the parameters Psi()
!
!.. assign modules
use GLOBAL_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
double precision,intent(INOUT),dimension(SPACEDIM,SPACEDIM) :: T
double precision,intent(INOUT),dimension(SPACEDIM,SPACEDIM) :: RGP
double precision,intent(INOUT),dimension(SPACEDIM)          :: Psi
!
!.. Local Scalars ..
Integer(4) :: i
!
!.. Local Arrays ..
double precision,dimension(SPACEDIM,SPACEDIM) :: Diag, PsiR, TTR
!
!.. Executable Statements ..
!
 Diag = zero
 do i = 1, SPACEDIM
   Diag(i, i) = Psi(i)
 end do

 call MatrixTransposition ( T, TTR, SPACEDIM, SPACEDIM )

 call MatrixProduct ( Diag, TTR, PsiR, SPACEDIM, SPACEDIM, SPACEDIM )

 call MatrixProduct ( T, PsiR, RGP, SPACEDIM, SPACEDIM, SPACEDIM )
!
return
!
end subroutine BoundaryPressureGradientMatrix3D
!---split

!cfile BoundaryReflectionMatrix2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : BoundaryReflectionMatrix2D
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
! Module purpose : Module for generation of generalized reflection matrix
!
! Calling routine: DefineBoundarySideGeometry2D
!
! Called routines: 
!
!************************************************************************************
!
subroutine BoundaryReflectionMatrix2D (T, R, PsiS, PsiN)

!Generates the generalised reflection matrix R(1 To SPACEDIM, 1 To SPACEDIM) on the base of
!the cosine matrix T(1 To SPACEDIM, 1 To SPACEDIM) and the parameters PsiS, PsiN
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
double precision   :: PsiS, PsiN 
double precision,dimension(1:SPACEDIM, 1:SPACEDIM)  :: T
double precision,dimension(1:SPACEDIM, 1:SPACEDIM)  :: R
!
!.. Executable Statements ..
!
R(1, 1) = PsiS * T(1, 1) * T(1, 1) + PsiN * T(3, 1) * T(3, 1) 
R(1, 3) = (PsiS - PsiN) * T(1, 1) * T(3, 1)
R(3, 1) = R(1, 3)
R(3, 3) = PsiS * T(3, 1) * T(3, 1) + PsiN * T(1, 1) * T(1, 1)
!
return
!
end subroutine BoundaryReflectionMatrix2D
!---split

!cfile BoundaryVolumeIntegrals2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : BoundaryVolumeIntegrals2D
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
! Module purpose : Module to compute the boundary volume integrals IntWdV only
!
! Calling routine: ComputeBoundaryDataTab
!
! Called routines: 
!
!************************************************************************************
!
subroutine BoundaryVolumeIntegrals2D (icbs, LocXY, xpmin, xpmax, interlen, IntWdV, IntDpWdV)
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
integer(4),      parameter :: Nalfadiv = 10
double precision,parameter :: eps = 0.05d0
!
!.. Formal Arguments ..
integer(4),      intent(IN)    :: icbs
double precision,intent(IN),   dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES) :: LocXY
double precision,intent(IN)    :: xpmin
double precision,intent(IN)    :: xpmax
double precision,intent(IN)    :: interlen
double precision,intent(INOUT) :: IntWdV
double precision,intent(INOUT),dimension(1:2) :: IntDpWdV
!
!.. Local Scalars ..
!integer(4)       :: NumBSides, iside
integer(4)       :: k, ndiv
double precision :: xpi, ypi, yplimite, &             !ypiq, 
                    dalfarif, Intalfa, dalfa, &
                    alfaA, alfaB, alfa_k, csiPA, csiPB, sinalfa, cosalfa, rb, rob, mult
!
!.. Local Arrays ..
!
! External functions and subrotuines
double precision,external :: WIntegr, J2Wro2
!
!.. Executable Statements ..
!
  dalfarif = PIGRECO / Nalfadiv
!
  IntWdV = zero
  IntDpWdV = zero
!
  if (interlen <= zero) return
!
  yplimite = eps * Domain%h
  xpi = LocXY(1, icbs)
  ypi = LocXY(2, icbs)
!
  if (ypi < yplimite) ypi = yplimite
!
!  ypiq = ypi * ypi
  csiPA = xpmin - xpi
  csiPB = xpmax - xpi
  alfaA = Atan2(csiPA, ypi)
  alfaB = Atan2(csiPB, ypi)
  Intalfa = alfaB - alfaA
  ndiv = Int(intalfa / dalfarif + half)
  if ( ndiv < 1 ) ndiv = 1
  dalfa = intalfa / ndiv
  alfa_k = alfaA - half * dalfa
!
  do k = 1, ndiv
    alfa_k = alfa_k + dalfa
    sinalfa = sin(alfa_k)
    cosalfa = cos(alfa_k)
!    if (cosalfa == zero) then 
!      rb = zero
!    else
      rb = ypi / cosalfa
!    end if
    rob = rb / Domain%h
    mult = Domain%h * J2Wro2(rob) * dalfa
    IntDpWdV(1) = IntDpWdV(1) + sinalfa * mult
    IntDpWdV(2) = IntDpWdV(2) - cosalfa * mult
    IntWdV = IntWdV + WIntegr(rb, Domain%h) * dalfa
  end do
!
! controllo se la particella è molto vicina allo spigolo interno 
!
  if (ypi == yplimite) then
!    if ( ) then ! la proiezione della particella è interna alla faccia
      IntWdV = half
!    else
!      IntWdV = zero
  end if
!
return
end subroutine BoundaryVolumeIntegrals2D
!---split

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

!cfile ComputeBoundaryIntegralTab.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ComputeBoundaryIntegralTab
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
! Module purpose : Module to compute local coordinates x, y, z of a grid of points,
!                  regularly distributed on the semisphere
!                  z<0 (radius = 2h), whose centre is the origin O of local axis.
!
! Calling routine: Gest_Trans
!
! Called routines: 
!
!************************************************************************************
!
subroutine ComputeBoundaryIntegralTab

!Computes local coordinates x, y, z of a grid of points, regularly distributed
!on the semisphere z<0 (radius = 2h), whose centre is the origin O of local axis.
!The semisphere will be superposed to the influence sphere of the
!generic particle near a plane boundary face, and oriented in such a way that
!the axis x, y, z coincide with the face local axes r, s, n.
!In the first three columns of the array BoundaryIntegralTab() the coordinates
!x, y, z, of each point is stored; in the forth column then relative d_alpha
!(portion of solid angle relative to the point, necessary fo integrations)
!is stored
!BITcols = 4
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Scalars ..
integer(4) :: i, j, points
double precision ::  delta_fi, delta_teta, fi, teta, xp, yp, zp, delta_alpha, &
                     delta_fimz, duesendelta_fimz
!
!.. Executable Statements ..
!
  delta_fi = FI_INTERVAL / FI_STEPS
  delta_fimz = half * delta_fi
  duesendelta_fimz = two * Sin(delta_fimz)
  delta_teta = TETA_INTERVAL / TETA_STEPS
!
  points = 0
  fi = -half * delta_fi
  Do i = 1, FI_STEPS
    fi = fi + delta_fi
    teta = -half * delta_teta
    zp = -doubleh * Cos(fi)
    Do j = 1, TETA_STEPS
      teta = teta + delta_teta
      xp = doubleh * Sin(fi) * Cos(teta)
      yp = doubleh * Sin(fi) * Sin(teta)
      delta_alpha = delta_teta * duesendelta_fimz * Sin(fi)
      points = points + 1
      BoundIntegralTab(points, 1) = xp
      BoundIntegralTab(points, 2) = yp
      BoundIntegralTab(points, 3) = zp
      BoundIntegralTab(points, 4) = delta_alpha
!.. coseni direttori del raggio OP di componenti (xp,yp,zp)
      BoundIntegralTab(points, 5) = xp / doubleh
      BoundIntegralTab(points, 6) = yp / doubleh
      BoundIntegralTab(points, 7) = zp / doubleh
    end do
  end do
!
return
end subroutine ComputeBoundaryIntegralTab
!---split

!cfile ComputeBoundaryDataTab.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ComputeBoundaryDataTab
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! 00  Agate/Guandalini  13/11/08       creation
!
!************************************************************************************
! Module purpose : Module to calculate array to store close boundary and integrals
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: diagnostic
!                  BoundaryVolumeIntegrals2D
!                  ComputeBoundaryVolumeIntegrals_P0
!                  ComputeSurfaceIntegral_WdS2D
!                  ComputeVolumeIntegral_WdV2D
!                  FindBoundaryIntersection2D
!                  FindCloseBoundarySides2D
!                  FindCloseBoundaryFaces3D
!                  InterpolateBoundaryIntegrals2D
!                  SelectCloseBoundarySides2D
!
!************************************************************************************
!
subroutine ComputeBoundaryDataTab 
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
!
!.. Local Scalars ..
  integer(4)       :: npi, Ncbs, IntNcbs, icbs, Ncbf, icbf, ibdt, Nfzn, Ncols
  double precision :: IntWdS, IntWdV, IntWdV1, IntdWrm1dV, IntGWZrm1dV, &
                      IntWd1s0, IntWd3s0, IntWd1s2
  double precision :: deltai, ypi, xpmin, xpmax, interlen     !xpi, 
  character(len=lencard) :: nomsub = "ComputeBoundaryDataTab"
  
!
!.. Local Arrays ..
  integer(4),      dimension(1:NUMCOLS_BIT)                   :: Colmn
  integer(4),      dimension(1:MAXCLOSEBOUNDSIDES)            :: Cloboside, Intboside
  double precision,dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES) :: LocXY, IntLocXY
  double precision,dimension(1:PLANEDIM)                      :: IntDpWdV
  integer(4),      dimension(1:Domain%MAXCLOSEBOUNDFACES)     :: Cloboface
  double precision,dimension(1:SPACEDIM,1:Domain%MAXCLOSEBOUNDFACES) :: LocX
  double precision,dimension(1:SPACEDIM)                      :: IntGWdV
  double precision,dimension(1:SPACEDIM,1:SPACEDIM)           :: IntGWrRdV
  double precision,dimension(1:NUMCOLS_BIT)                   :: Func
!  type (TyBoundaryData),dimension(:),allocatable :: buffer1
!
!.. Executable Statements ..
!
  BoundaryDataPointer = 0
!
  if (ncord == 2) then
!
!..  azzeramento contatore numero particelle vicine al contorno
    BoundarySide(:)%CloseParticles = 0
    BoundarySide(:)%CloseParticles_maxQuota = const_m_9999
!
!$omp parallel do default(none) &
!$omp private(npi,Ncbs,Cloboside,LocXY,IntNcbs,Intboside,IntLocXY,ibdt,icbs) &
!$omp private(xpmin,xpmax,interlen,Ncols,Colmn,deltai,Func,ypi) &
!$omp private(IntWdS,IntWdV,IntDpWdV,IntWdV1,IntWd1s0,IntWd3s0,IntWd1s2) &
!$omp shared(nag,pg,Domain,BoundaryDataTab,BoundaryDataPointer,MaxNcbs,nomsub,nout,nscr,BoundarySide,squareh)

    do npi = 1,nag
!
      if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" ) cycle
!
!.. searches for the boundary sides nearest the npi-th current particle

      call FindCloseBoundarySides2D (npi,Ncbs, Cloboside, LocXY)

!
!.. some nearest boundary has been detected
      if (Ncbs > 0) then
!
!.. selects the boundary sides that effectively contribute to the motion field terms
        call SelectCloseBoundarySides2D (npi, Ncbs, Cloboside, LocXY, IntNcbs, Intboside, IntLocXY)
!
        if (IntNcbs > 0) then
!
          BoundaryDataPointer(1,npi) = Ncbs
          BoundaryDataPointer(2,npi) = IntNcbs
          ibdt = MAXCLOSEBOUNDSIDES * (npi-1)
          BoundaryDataPointer(3,npi) = ibdt+1
!
          do icbs = 1,IntNcbs
!
            ibdt = ibdt + 1
!.. controllo dimensioni vettori
            if (ibdt > MaxNcbs) then
              call diagnostic (arg1=8,arg2=1,arg3=nomsub)
!
!              MaxNcbs = ibdt * 1.2
!              write(nout,'(a,i15)') " New Max num particles*BoundaryCloseSides the new value is: MaxNcbs = ",MaxNcbs
!              allocate (buffer1(1:MaxNcbs), stat = ier)
!              if (ier /= 0) then
!                write (nout,'(1x,a,i2)') "  REDIMENSION Array BUFFER1 not allocated. Error code: ",ier
!                call diagnostic (arg1=4,arg3=nomsub)
!              else
!                write (nout,'(1x,a)') "  REDIMENSION Array BUFFER1 successfully allocated "
!              end if 
!              buffer1(1:ibdt-1) = BoundaryDataTab(1:ibdt-1)
!              deallocate (BoundaryDataTab)
!              allocate (BoundaryDataTab(1:MaxNcbs), stat = ier)
!              if (ier /= 0) then
!                write (nout,'(1x,a,i2)') "  REDIMENSION Array BoundaryDataTab not allocated. Error code: ",ier
!                call diagnostic (arg1=4,arg3=nomsub)
!              else
!                write (nout,'(1x,a)') "  REDIMENSION Array BoundaryDataTab successfully allocated "
!              end if 
!              BoundaryDataTab(1:ibdt-1) = buffer1(1:ibdt-1)
!              deallocate (buffer1)
            end if
!
            BoundaryDataTab(ibdt)%CloBoNum  = Intboside(icbs)
            BoundaryDataTab(ibdt)%LocXYZ(1:PLANEDIM) = IntLocXY(1:PLANEDIM,icbs)
            BoundaryDataTab(ibdt)%LocXYZ(3) = zero
            Func = zero
!
            if (IntNcbs == 2) then  !Calcolo numerico degli integrali IntWds e IntWdV
!.. Calcolo integrali
!.. 2D computation of the boundary integrals (semy-analitic approach, volumetric method)
!.. to find the intersection between the boundary and the kernel support
              call FindBoundaryIntersection2D (icbs, Intboside, IntLocXY, BoundarySide, xpmin, xpmax, interlen)
!.. computation of the 1D integrals
              call ComputeSurfaceIntegral_WdS2D (icbs, IntLocXY, xpmin, interlen, IntWdS)  !Intboside, xpmax, 
!.. computation of the 2D integrals
              call BoundaryVolumeIntegrals2D (icbs, IntLocXY, xpmin, xpmax, interlen, IntWdV, IntDpWdV)
              call ComputeVolumeIntegral_WdV2D (icbs, IntNcbs, Intboside, IntLocXY, BoundarySide, xpmin, xpmax, interlen, IntWdV1)
!.. Interpolation of some integrals (from tables)
!.. Interpolazione da tabella degli integrali IntWd1s0, IntWd3s0, IntWd1s2 
              Ncols=3
              Colmn(1)=3
              Colmn(2)=4
              Colmn(3)=5
              ypi = IntLocXY(2, icbs)
              deltai = ypi / Domain%h
              call InterpolateBoundaryIntegrals2D (deltai, Ncols, Colmn, Func)
              IntWd1s0 = Func(1) / squareh
              IntWd3s0 = Func(2) / squareh
              IntWd1s2 = Func(3) / squareh
!
              BoundaryDataTab(ibdt)%BoundaryIntegral(1) = IntWdS
              BoundaryDataTab(ibdt)%BoundaryIntegral(2) = IntWdV
              BoundaryDataTab(ibdt)%BoundaryIntegral(3) = IntWdV1
              BoundaryDataTab(ibdt)%BoundaryIntegral(4:5) = IntDpWdV(1:2)
              BoundaryDataTab(ibdt)%BoundaryIntegral(6) = IntWd1s0
              BoundaryDataTab(ibdt)%BoundaryIntegral(7) = IntWd3s0
              BoundaryDataTab(ibdt)%BoundaryIntegral(8) = IntWd1s2
! 
            else if (IntNcbs == 1) then
 ! Interpolazione da tabella degli integrali IntWdS,IntWdV, IntWd1s0, IntWd3s0, IntWd1s2 
              Ncols=5
              Colmn(1) = 1
              Colmn(2) = 2
              Colmn(3) = 3
              Colmn(4) = 4
              Colmn(5) = 5
!              xpi = IntLocXY(1, icbs)
              ypi = IntLocXY(2, icbs)
              deltai = ypi / Domain%h
              call InterpolateBoundaryIntegrals2D (deltai, Ncols, Colmn, Func)
              IntWdS = Func(1) / Domain%h
              IntWdV = Func(2)
              IntWd1s0 = Func(3) / squareh
              IntWd3s0 = Func(4) / squareh
              IntWd1s2 = Func(5) / squareh
              BoundaryDataTab(ibdt)%BoundaryIntegral(1) = IntWds
              BoundaryDataTab(ibdt)%BoundaryIntegral(2) = IntWdV
              BoundaryDataTab(ibdt)%BoundaryIntegral(3) = IntWdV
              BoundaryDataTab(ibdt)%BoundaryIntegral(4:5) = zero
              BoundaryDataTab(ibdt)%BoundaryIntegral(6) = IntWd1s0
              BoundaryDataTab(ibdt)%BoundaryIntegral(7) = IntWd3s0
              BoundaryDataTab(ibdt)%BoundaryIntegral(8) = IntWd1s2
            end if
          end do
        end if
      end if
    end do
!
!$omp end parallel do
!
!----------------------
!.. 3D
!----------------------
  else
!
!..  azzeramento contatore numero particelle vicine al contorno
    BoundaryFace(:)%CloseParticles = 0
    BoundaryFace(:)%CloseParticles_maxQuota = const_m_9999
!
    
!AA504sub
!$omp parallel do default(none) &
!$omp private(npi,Ncbf,Cloboface,LocX,Nfzn,icbf,ibdt) &
!$omp private(IntWdV,IntdWrm1dV,IntGWZrm1dV,IntGWdV,IntGWrRdV) &
!$omp shared(nag,pg,BoundaryDataTab,BoundaryDataPointer,EpCount,MaxNcbf,nout,nscr,nomsub,Domain)

    loop_particle:  do npi = 1,nag
!
      if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" ) cycle loop_particle
!
!.. searches for the boundary faces nearest the npi-th current particle

      call FindCloseBoundaryFaces3D (npi, Ncbf, Cloboface, LocX, Nfzn)

!
      if (Ncbf == 0) then
!
        BoundaryDataPointer(:,npi) = 0
!
      else
!
!.. some nearest boundary has been detected
!
!.. conteggio particelle sfuggite dal contorno ma ancora all'interno della griglia (EpCount)
!..  coordinate normali alle faces tutte negative
        if (Nfzn == Ncbf) then
          EpCount(pg(npi)%imed) = EpCount(pg(npi)%imed) + 1
!.. eliminazione della particella
!          write (699,'(a,f10.4,a,i10,a,3f15.6)') 'time ',tempo,' escaped particle n.',npi,' coordinates ', pg(npi)%coord
!          pg(nag)%cella = -2
        end if
!
        BoundaryDataPointer(1,npi) = Ncbf
        BoundaryDataPointer(2,npi) = 0
        ibdt = Domain%MAXCLOSEBOUNDFACES * (npi-1)
        BoundaryDataPointer(3,npi) = ibdt + 1
!
!.. controllo dimensioni vettori
        if (ibdt > MaxNcbf) then
          call diagnostic (arg1=8,arg2=2,arg3=nomsub)
        end if
!
        do icbf = 1,Ncbf
!
          call ComputeBoundaryVolumeIntegrals_P0 (icbf, Cloboface, LocX, IntWdV, IntdWrm1dV, IntGWZrm1dV, IntGWdV, IntGWrRdV)
!
          ibdt = ibdt + 1
!
          BoundaryDataTab(ibdt)%CloBoNum  = CloboFace(icbf)
          BoundaryDataTab(ibdt)%LocXYZ(1:SPACEDIM)    = LocX(1:SPACEDIM,icbf)
          BoundaryDataTab(ibdt)%BoundaryIntegral(1)   = zero
          BoundaryDataTab(ibdt)%BoundaryIntegral(2)   = IntWdV
          BoundaryDataTab(ibdt)%BoundaryIntegral(3)   = IntdWrm1dV
          BoundaryDataTab(ibdt)%BoundaryIntegral(4:6) = IntGWdV(1:3)
          BoundaryDataTab(ibdt)%BoundaryIntegral(7)   = IntGWZrm1dV
          BoundaryDataTab(ibdt)%BoundaryIntegral(8)   = zero
          BoundaryDataTab(ibdt)%IntGiWrRdV(:,:)       = IntGWrRdV(:,:)
!
        end do
!
      end if
!
    end do loop_particle
!
!AA504 sub
!$omp end parallel do
!
  end if
!
return
end subroutine ComputeBoundaryDataTab
!---split

!cfile ComputeBoundaryVolumeIntegrals_P0.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ComputeBoundaryVolumeIntegrals_P0
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
! Module purpose : Module to compute boundary contributions to rodivV
!
! Calling routine: ComputeBoundaryDataTab
!
! Called routines: InterpolateTable
!                  IsPointInternal
!                  LocalNormalCoordinates
!
!************************************************************************************
!
subroutine ComputeBoundaryVolumeIntegrals_P0 (icbf, Clobface, LocX, IntWdV, IntdWrm1dV, IntGWZrm1dV, IntGWdV, IntGWrRdV)
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
  integer(4),      intent(IN)  :: icbf
  integer(4),      intent(IN), dimension(1:Domain%MAXCLOSEBOUNDFACES) :: Clobface
  double precision,intent(OUT) :: IntWdV, IntdWrm1dV, IntGWZrm1dV
  double precision,intent(IN), dimension(1:SPACEDIM,1:Domain%MAXCLOSEBOUNDFACES) :: LocX
  double precision,intent(OUT),dimension(1:SPACEDIM)                      :: IntGWdV
  double precision,intent(OUT),dimension(1:SPACEDIM,1:SPACEDIM)           :: IntGWrRdV
!
!.. Local Parameters ..
integer(4),parameter :: nicols = 4
double precision,parameter :: eps = 0.05d0
!
!Computes the boundary volume integral IntWdV
!
!.. Local Scalars ..
integer(4) :: SD, sdj, ipk, fkod
integer(4) :: iface
double precision :: deltaPiPx, PiPxdist, rob, delta_alpha, tsegnato, zpmin, RZ, RZ2, robpre
double precision :: JW3ro2dA, JdW3ro1dA, JdW3ro2dA, JdW3ro3dA
!
!.. Local Arrays ..
integer(4),dimension(1:nicols) :: icol
double precision,dimension(1:SPACEDIM) :: LocPi, LocPj, LocPiPj, PXLoc, csi, RLoccos
double precision,dimension(1:nicols) :: ivalue
!
! External functions and subrotuines
logical, external :: IsPointInternal
!
!.. Executable Statements ..
!
  icol(1) = 1
  icol(2) = 2
  icol(3) = 3
  icol(4) = 4
!
  IntWdV = zero
  IntdWrm1dV = zero
  IntGWZrm1dV = zero
  IntGWdV(:) = zero
  IntGWrRdV(:,:) = zero
!
  JW3ro2dA  = zero
  JdW3ro1dA = zero
  JdW3ro2dA = zero
  JdW3ro3dA = zero
!
  iface = Clobface(icbf)
  zpmin = eps * Domain%h
!
  LocPi(:) = LocX(:, icbf)                     !Coordinate locali della particella reale Pi
!
  if (LocPi(3) < zero)  return
!
  if (LocPi(3) < zpmin)  LocPi(3) = zpmin
!
  robpre = 2.1
  do ipk = 1,BITrows
    LocPiPj(:) = BoundIntegralTab(ipk, 1:3)      !Componenti locali del segmento orientato PiPj
    LocPj(:) = LocPi(:) + LocPiPj(:)             !Coordinate locali del punto Pj
    if (LocPj(3) <= 0) then   !Il punto Pj is si trova nel semispazio della
                              ! faccia "iface", opposto a quello della particella reale "Pi"
      !Coordinate locali del punto PX intersezione del segmento PiPj con la faccia "iface"
      tsegnato = LocPi(3) / (LocPi(3) - LocPj(3))
      PXLoc(:) = LocPi(:) + (LocPj(:) - LocPi(:)) * tsegnato
      PXLoc(3) = zero
!
      call LocalNormalCoordinates (PXLoc, csi, iface)
!
      fkod = BoundaryFace(iface)%nodes - 2
!
      if (IsPointInternal(fkod, csi)) then    !PX è interno alla faccia e quindi il punto
                                                !Pj da' contributo all'integrale di contorno
!
        !*******  Distance between Pi and Px  ***********************************
        PiPxdist = zero
        do SD = 1,SPACEDIM
          deltaPiPx = PXLoc(SD) - LocPi(SD)
          PiPxdist = PiPxdist + deltaPiPx * deltaPiPx
        end do
        PiPxdist = Dsqrt(PiPxdist)
!
        !*******  Normalised distance and related functions   *******************
        delta_alpha = BoundIntegralTab(ipk, 4)
        rob = PiPxdist / Domain%h
!
        if (Abs(rob - robpre) > 0.001) then
!
          ivalue = zero
          call InterpolateTable (rob, nicols, icol, ivalue)
!
          JW3ro2dA  = ivalue(1) * delta_alpha
          JdW3ro1dA = ivalue(2) * delta_alpha
          JdW3ro2dA = ivalue(3) * delta_alpha
          JdW3ro3dA = ivalue(4) * delta_alpha
          robpre = rob
!
        end if
!
        RLoccos(:) = BoundIntegralTab(ipk, 5:7)
        RZ = RLoccos(3)
        RZ2 = RZ * RZ
!
        !*******  Boundary volume integrals  ***********************************
        IntWdV = IntWdV + JW3ro2dA
        IntdWrm1dV = IntdWrm1dV + JdW3ro1dA
        IntGWZrm1dV = IntGWZrm1dV + RZ2 * JdW3ro1dA
        do SD = 1,SPACEDIM
          IntGWdV(SD) = IntGWdV(SD) + RLoccos(SD) * JdW3ro2dA
          do sdj = SD,SPACEDIM
            IntGWrRdV(SD, sdj) = IntGWrRdV(SD, sdj) + RLoccos(SD) * RLoccos(sdj) * JdW3ro3dA
          end do
        end do
      end if
    end if
  end do
!
! Completamento della matrice simmetrica IntGWrRdV(3, 3)
  IntGWrRdV(2,1) = IntGWrRdV(1,2)
  IntGWrRdV(3,1) = IntGWrRdV(1,3)
  IntGWrRdV(3,2) = IntGWrRdV(2,3)
!
  if (LocPi(3) == zpmin) then
    PXLoc(:) = LocPi(:)
    PXLoc(3) = zero
!
    call LocalNormalCoordinates (PXLoc, csi, iface)
!
    fkod = BoundaryFace(iface)%nodes - 2
!
    if (IsPointInternal(fkod, csi)) then    !La proiezione della particella è interna alla faccia
      IntWdV = half
    Else                                    !La proiezione della particella è esterna alla faccia
      IntWdV = zero
    end if
  end if
!
return
end subroutine ComputeBoundaryVolumeIntegrals_P0
!---split

!cfile ComputeKernelTable.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ComputeKernelTable
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
! Module purpose : Module to compute boundary contributions to rodivV
!
! Calling routine: Gest_Trans
!
! Called routines: diagnostic
!
!************************************************************************************
!
subroutine ComputeKernelTable
!
!Precomputes and stores in kerneltab(0:ktrows, 0:ktcols) the following values
!kerneltab(0:ktrows, 0) = rob = rb/h
!kerneltab(0:ktrows, 1) = Int W* ro2 dro         (from rob to 2)
!kerneltab(0:ktrows, 2) = Int dW*/dro ro dro     (from rob to 2)
!kerneltab(0:ktrows, 3) = Int dW*/dro ro^2 dro   (from rob to 2)
!kerneltab(0:ktrows, 4) = Int dW*/dro ro^3 dro   (from rob to 2)
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Scalars ..
integer(4)        :: nr
double precision  :: rob
character(len=lencard) :: nomsub = "ComputeKernelTable"
!
!.. External routines ..
double precision, external :: IWro2dro
double precision, external :: JdWsRn
!
!.. Executable Statements ..
!
ktdelta = two / INT_KERNELTABLE
!
if (ncord == 3) then
  do nr = 0, ktrows
    rob = ktdelta * nr
    kerneltab(nr,0) = rob
    kerneltab(nr,1) = IWro2dro(rob)
    kerneltab(nr,2) = JdWsRn(rob,3,1,1) * Unosusquareh
    kerneltab(nr,3) = JdWsRn(rob,3,2,1) * Unosuh
    kerneltab(nr,4) = JdWsRn(rob,3,3,1)
  end do
!
else if (ncord == 2) then
  call diagnostic (arg1=8,arg2=3,arg3=nomsub)
!  do nr = 0,ktrows
!    rob = ktdelta * nr
!    kerneltab(nr,0) = rob
!    kerneltab(nr,1) = IWro2dro(rob)
!    kerneltab(nr,2) = JdWsRn(rob,2,1,1) * Unosusquareh
!    kerneltab(nr,3) = JdWsRn(rob,2,2,1) * Unosuh
!    kerneltab(nr,4) = JdWsRn(rob,2,3,1)
!  end do
end if
!
return
end subroutine ComputeKernelTable
!---split

!cfile ComputeSurfaceIntegral_WdS2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ComputeSurfaceIntegral_WdS2D
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
! Module purpose : Module for Computing the surface integral of kernel W along the
!                  segments intercepted by the influence circle (radius=2h) of the
!                  particle i, whose local coordinates are xpi=LocXY(1, icbs) and
!                  ypi=LocXY(2, icbs), on the adjacent boundary side icbs.
 
!
! Calling routine: ComputeBoundaryDataTab
!
! Called routines: 
!
!************************************************************************************
!
subroutine ComputeSurfaceIntegral_WdS2D (icbs, LocXY, xpmin, interlen, SIntWds) !Cloboside, xpmax, 

!Computes the surface integral of kernel W along the segments intercepted
!by the influence circle (radius=2h) of the particle i, whose local coordinates
!are xpi=LocXY(1, icbs) and ypi=LocXY(2, icbs), on the adjacent boundary side icbs.
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
integer(4), parameter :: ndivrif = 2
!
!.. Formal Arguments ..
integer(4),      intent(IN)    :: icbs
!integer(4),      intent(IN),dimension(1:MAXCLOSEBOUNDSIDES)            :: Cloboside
double precision,intent(IN),dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES) :: LocXY
double precision,intent(IN)    :: xpmin
!double precision,intent(IN)    :: xpmax
double precision,intent(IN)    :: interlen
double precision,intent(INOUT) :: SIntWds
!
!.. Local Scalars ..
integer(4)       :: ndiv, ipt    !iside, 
double precision :: xpi, ypi, ypiq, dsrif, deltas, xpip, dxpip, ris, risq, Wds
!
!.. Local Arrays ..
!integer(4),dimension(1:PLANEDIM)  :: acix
!
! External functions and subrotuines
double precision,    external    :: w
!
!.. Executable Statements ..
!
!  acix(1) = 1        !active coordinate indexes
!  acix(2) = 3

  SIntWds = zero
!  iside = Cloboside(icbs)
  dsrif = Domain%h / ndivrif

  if (interlen > zero) then
    xpi = LocXY(1, icbs)
    ypi = LocXY(2, icbs)
    if (ypi < zero) then      !if the particle is out of the boundary operates as if the particle
        ypi = zero            !where on the boundary
    end if    
    ypiq = ypi * ypi
    ndiv = Int(interlen / dsrif + half)
    if (ndiv < ndivrif) then
        ndiv = ndivrif
    end if
    deltas = interlen / ndiv
    xpip = xpmin - half * deltas
    do ipt = 1, ndiv
      xpip = xpip + deltas
      dxpip = xpip - xpi
      risq = dxpip * dxpip + ypiq
      ris = Dsqrt(risq)
      Wds = w(ris, Domain%h, Domain%coefke) * deltas
      SIntWds = SIntWds + Wds
    end do
  end if
!
return
end subroutine ComputeSurfaceIntegral_WdS2D
!---split

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

!cfile DefineBoundaryFaceGeometry3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : DefineBoundaryFaceGeometry3D
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
! Module purpose : Module to define boundary face from geometry 3D
!
! Calling routine: Gest_Input
!
! Called routines: BoundaryMassForceMatrix3D
!                  BoundaryPressureGradientMatrix3D
!                  DefineLocalSystemVersors
!
!************************************************************************************
!
subroutine DefineBoundaryFaceGeometry3D
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
integer(4) :: Kf, Nf, Nt
!
!.. Local Arrays ..
!double precision, dimension(SPACEDIM)          :: Psi
double precision, dimension(SPACEDIM)          :: Fi
double precision, dimension(SPACEDIM,SPACEDIM) :: TT, RGP, RMF

! 20081028  Data Psi /1.0d0, 1.0d0, 0.0d0/
! 20081028  Data Fi /0.0d0, 0.0d0, 1.0d0/
!Data Psi /0.0d0, 0.0d0, 0.0d0/
Data Fi /1.0d0, 1.0d0, 1.0d0/
!
!.. Executable Statements ..
!
  do Kf = 1, NumFacce, 1
!
    Nf = BFaceList(Kf)
!
    if ( Nf == 0 ) cycle
!
    Nt = BoundaryFace(Nf)%stretch

!AA504 sub
    call DefineLocalSystemVersors (Nf)
!
    !Psi(1:SPACEDIM) = Tratto(nt)%PsiCoeff(1:SPACEDIM)
    !Fi (1:SPACEDIM) = Tratto(nt)%FiCoeff(1:SPACEDIM)
    TT (1:SPACEDIM, 1:SPACEDIM) = BoundaryFace(nf)%T(1:SPACEDIM, 1:SPACEDIM)
!
    RGP = zero
    RMF = zero
! 20081028     call BoundaryPressureGradientMatrix3D ( TT, RGP, Psi )
    call BoundaryMassForceMatrix3D ( TT, RMF, Fi )
!
! 20081028    BoundaryFace(nf)%RPsi(1:SPACEDIM, 1:SPACEDIM) = RGP(1:SPACEDIM, 1:SPACEDIM)
    BoundaryFace(nf)%RFi (1:SPACEDIM, 1:SPACEDIM) = RMF(1:SPACEDIM, 1:SPACEDIM)
!
    if ( Tratto(nt)%tipo == "tapi" ) then
       BoundaryFace(nf)%velocity(1:SPACEDIM) = Tratto(nt)%velocity(1:SPACEDIM)
    else
       BoundaryFace(nf)%velocity(1:SPACEDIM) = zero
    end if
   !nv = nv + 1
!
  end do
!
return
end subroutine DefineBoundaryFaceGeometry3D
!---split

!cfile DefineBoundarySideGeometry2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : DefineBoundarySideGeometry2D
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
! Module purpose : Module for definition of boundary side 
!
! Calling routine: Gest_Input
!
! Called routines: BoundaryReflectionMatrix2D
!                  BoundaryMassForceMatrix2D
!                  DefineBoundarySideRelativeAngles2D
!
!************************************************************************************
!
  subroutine DefineBoundarySideGeometry2D
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
  integer(4)       :: NumBS, nt, ns, nst, nv, v1, v2
  double precision :: Dx, Dz, L, sx, sz, PsiS, PsiN, FiS, FiN  
  character(4)     :: tipo          
!
!.. Local Arrays ..
  double precision,dimension(1:SPACEDIM, 1:SPACEDIM) :: TT, RR, RRN
!
!.. Executable Statements ..
!
!PsiN = Coefficiente di proiezione tale che (dp/dn)b=PsiN*(dp/dn)i                
!PsiS = Coefficiente di proiezione tale che (dp/ds)b=PsiS*(dp/ds)i                   
!FiN = Coefficiente di proiezione tale che (dp/dn)b=FiN*roi*gn                 
!FiS = Coefficiente di proiezione tale che (dp/ds)b=FiS*roi*gs                  
!
  NumBS = 0
!
  do nt = 1, NumTratti
!
    nst = Tratto(nt)%numvertices - 1
    nv  = Tratto(nt)%inivertex
    Tratto(nt)%iniside = NumBS + 1
!
    do ns = 1, nst
      v1 = BoundaryVertex(nv)
      v2 = BoundaryVertex(nv + 1)
      NumBS = NumBS + 1
!
      Dx = Vertice(1,v2) - Vertice(1,v1)
      Dz = Vertice(3,v2) - Vertice(3,v1)
      L = Dsqrt(Dx * Dx + Dz * Dz)
      sx = Dx / L
      sz = Dz / L
!
      BoundarySide(NumBS)%stretch = nt
      BoundarySide(NumBS)%tipo = Tratto(nt)%tipo
      BoundarySide(NumBS)%Vertex(1) = v1
      BoundarySide(NumBS)%Vertex(2) = v2
      BoundarySide(NumBS)%length = L
!
      BoundarySide(NumBS)%T(1,1) = sx
      !BoundarySide(NumBS)%T(2,1) = zero
      BoundarySide(NumBS)%T(3,1) = sz
      !BoundarySide(NumBS)%T(1,2) = zero
      !BoundarySide(NumBS)%T(2,2) = one
      !BoundarySide(NumBS)%T(3,2) = zero
      BoundarySide(NumBS)%T(1,3) =-sz
      !BoundarySide(NumBS)%T(2,3) = zero
      BoundarySide(NumBS)%T(3,3) = sx
!
      TT=BoundarySide(NumBS)%T
!
      tipo = Tratto(Nt)%tipo   
      select case ( tipo )
!
        case ( "fixe", "tapi" )
!            PsiS = one
!            PsiN = one
!            FiS = zero
!            FiN = zero
!.. per far funzionare le modifiche 'AdM 15-10-08' in AddBoundaryContributions_to_ME2D e 3D
!..   vanno modificati PsiS PsiN FiS e FiN come segue. 
            PsiS = zero
            PsiN = zero
            FiS = one
            FiN = one 
!        case ( "leve", "velo", "crit", "open" )
!            PsiS = one  
!            PsiN = zero
!            FiS = zero
!            FiN = one
        case ( "leve" )
            PsiS = zero
            PsiN = zero  
            FiS = one
            FiN = one
        case ( "velo", "flow" )
            PsiS = zero
            PsiN = zero
            FiS = zero        
            FiN = zero
        case ( "crit" )
            PsiS = zero
            PsiN = zero
            FiS = one
            FiN = one
        case ( "open" )
            PsiS = zero
            PsiN = zero  
            FiS = zero 
            FiN = zero   
        case ( "sour" )
            PsiS = zero
            PsiN = zero
!            FiS = one
!            FiN = one
            FiS = zero
            FiN = zero
      end select                                                        
       
!     RR = zero
!     RRN = zero

      call BoundaryReflectionMatrix2D (TT, RR, PsiS, PsiN)
!      call BoundaryReflectionMatrix2D (TT, RR, Psi)  

      call BoundaryMassForceMatrix2D (TT, RRN, FiS, FiN) 
!      call BoundaryMassForceMatrix2D (TT, RRN, Fi)     

      BoundarySide(NumBS)%R  = RR
      BoundarySide(NumBS)%RN = RRN

      if (Tratto(nt)%tipo == "tapi") then
        BoundarySide(NumBS)%velocity = Tratto(nt)%velocity
      else
        BoundarySide(NumBS)%velocity(:) = zero
      end if
!
      nv = nv + 1
!
    end do
!
  end do
!
  call DefineBoundarySideRelativeAngles2D
!
  return
  end subroutine DefineBoundarySideGeometry2D
!---split

!cfile DefineBoundarySideRelativeAngles2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : DefineBoundarySideRelativeAngles2D
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
! Module purpose : Module for calculation of previous adjacent side and 
!                  relative angle for each boundary side
!
! Calling routine: DefineBoundarySideGeometry2D
!
! Called routines: 
!
!************************************************************************************
!
  subroutine DefineBoundarySideRelativeAngles2D
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
  integer(4)       :: nti, ntj, isi, jsi, ksi
  double precision :: sinangle, cosangle
!
!.. local arrays ..
  double precision, dimension(1:SPACEDIM, 1:SPACEDIM) ::    Tmx,  PTmx
!
!.. Executable Statements ..
!
!.. loops on all the domain boundary sides
!
  do isi = 1, NumBSides
!
    nti = BoundarySide(isi)%stretch
!
!.. skips the perimeter and pool types
!
    if (Tratto(nti)%tipo /= "peri" .AND. Tratto(nti)%tipo /= "pool" ) then
!
!.. loops on all the other sides
!
      ksi=0
      do jsi = 1, NumBSides
!
       ntj = BoundarySide(jsi)%stretch
       if (Tratto(ntj)%tipo == "peri" .or. Tratto(ntj)%tipo == "pool" ) cycle
!
!.. checks if the sides are adjacents; !ksi is the previous adjacent side
!
        if (BoundarySide(jsi)%Vertex(2) == BoundarySide(isi)%Vertex(1)) then
          ksi = jsi                    
          exit
        end if
!
      end do
!
      BoundarySide(isi)%angle = zero
!      BoundarySide(isi)%previous_side = 0
!
!.. an adjacent side has been found
!
      if (ksi > 0) then
!
!.. evaluates the angle between the two sides in radiants
!
        ntj = BoundarySide(ksi)%stretch
!
! errore corretto 25set08       if (Tratto(nti)%tipo /= "peri" .AND. Tratto(nti)%tipo /= "pool" ) then
        if (Tratto(ntj)%tipo /= "peri" .AND. Tratto(ntj)%tipo /= "pool" ) then
!            
          Tmx  = BoundarySide(isi)%T
! errore corretto 25set08         PTmx = BoundarySide(jsi)%T
          PTmx = BoundarySide(ksi)%T
!
          sinangle = PTmx(1, 1) * Tmx(3, 1) - PTmx(3, 1) * Tmx(1, 1)
          cosangle = PTmx(1, 1) * Tmx(1, 1) + PTmx(3, 1) * Tmx(3, 1)
          BoundarySide(isi)%angle = Atan2(sinangle, cosangle)
!
        end if
! errore corretto 06ott08
!        BoundarySide(isi)%previous_side = ksi
!
      end if
!
      BoundarySide(isi)%previous_side = ksi
!
    end if
!
  end do
!
  return
  end subroutine DefineBoundarySideRelativeAngles2D
!---split

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

!cfile EvaluateBER_TimeStep.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : EvaluateBER_TimeStep
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
! Module purpose : Module to compute boundary contributions to rodivV
!
! Calling routine: Gest_Trans
!
! Called routines: 
!
!************************************************************************************
!
  subroutine EvaluateBER_TimeStep
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
  integer(4) :: imed
  double precision :: DTlocal
!
!.. Executable Statements ..
!
  DTminBER = 1.0d30
  do imed = 1,NMedium
    DTlocal = Domain%h / Med(imed)%celerita
    if (DTlocal < DTminBER) DTminBER = DTlocal
  end do
  DTminBER = sqrttwo * DTminBER
!
  return
  end subroutine EvaluateBER_TimeStep
!---split

!cfile FindBoundaryConvexEdges3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : FindBoundaryConvexEdges3D
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!AA504
! 03  Amicarelli        08/04/2014     (v5.04) omp parallelization
!
!************************************************************************************
! Module purpose : Module to compute the boundary integral IntWdS
!
! Calling routine: Gest_Input
!
! Called routines: 
!
!************************************************************************************
!
subroutine FindBoundaryConvexEdges3D
!
!.. Cerca eventuali spigoli (edges) convessi del contorno e ne memorizza i dati geometrici
!.. nell'array BoundaryConvexEdge() as TyBoundaryConvexEdge
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Scalars ..
integer(4) :: kf, nf, sd, nodes, nodrif, nt, fk, i, j, nodi, nodj
double precision :: delta, length2, zitakk
integer(4) :: kf1, nf1, fk1, nt1, ii, jj, kk, nodii, nodjj, nodkk, NBFfin, NBFini, nodes1
logical    :: EdgeFound
character(len=lencard) :: nomsub = "FindBoundaryConvexEdges3D"
!
!.. Local Arrays ..
integer(4), dimension(2, 4) :: iseg
!
!.. Executable Statements ..
!
  iseg(1, 1) = 2
  iseg(1, 2) = 3
  iseg(1, 3) = 1
  iseg(1, 4) = 0
  iseg(2, 1) = 2
  iseg(2, 2) = 3
  iseg(2, 3) = 4
  iseg(2, 4) = 1
!
  NumBEdges = 0
  NBFfin = NumFacce - 1
!AA504 omp directives
!$omp parallel do default(none) &
!$omp shared(NBFfin,BFaceList,BoundaryFace,NumFacce,Tratto,Vertice,NumBEdges,Domain,nomsub,BoundaryConvexEdge,iseg) &
!$omp private(kf,nf,nt,nodes,fk,NBFini,kf1,nf1,nt1,EdgeFound,nodes1,fk1,i,j,nodi,nodj,ii,jj,nodii,nodjj,kk,nodkk,zitakk,length2,nodrif,sd,delta)
 do kf = 1, NBFfin
    nf = BFaceList(kf)
    nt = BoundaryFace(nf)%stretch
    if (Tratto(nt)%tipo == "fixe" .or. Tratto(nt)%tipo == "tapi") then
      nodes = BoundaryFace(nf)%nodes
      fk = nodes - 2
      NBFini = kf + 1
      do kf1 = NBFini, NumFacce
        nf1 = BFaceList(kf1)
        nt1 = BoundaryFace(nf1)%stretch
        EdgeFound = .False.
        if (Tratto(nt1)%tipo == "fixe" .or. Tratto(nt1)%tipo == "tapi") then
          nodes1 = BoundaryFace(nf1)%nodes
          fk1 = nodes1 - 2
          do i = 1, nodes
            j = iseg(fk, i)
            nodi = BoundaryFace(nf)%node(i)%name 
            nodj = BoundaryFace(nf)%node(j)%name 
            do ii = 1, nodes1
              jj = iseg(fk1, ii)
              nodii = BoundaryFace(nf1)%node(ii)%name 
              nodjj = BoundaryFace(nf1)%node(jj)%name 
!
              if ((nodi == nodii .and. nodj == nodjj) .or. (nodi == nodjj .and. nodj == nodii)) then
!.. Trovato lato comune alle facce nf e nf1
!.. Verifica se si tratta di uno spigolo convesso
                kk = iseg(fk1, jj)             !Nodo della faccia nf1 diverso da nodii e nodjj
                nodkk = BoundaryFace(nf1)%node(kk)%name 
                zitakk = zero              !Coordinata normale del nodo nodkk rispetto alla faccia nf
                nodrif = BoundaryFace(nf)%nodes
                do sd = 1, SPACEDIM
                  zitakk = zitakk + BoundaryFace(nf)%T(sd, 3) * (Vertice(sd,nodkk) - BoundaryFace(nf)%node(nodrif)%GX(sd))   
                end do
                if (zitakk < zero) then     !Lo spigolo di estremi (nodi,nodj) è convesso
!AA504 omp directives
!$omp critical (omp_FBCE3D)
                  NumBEdges = NumBEdges + 1
                  if (NumBEdges > Domain%MAXNUMCONVEXEDGES) call diagnostic (arg1=8,arg2=10,arg3=nomsub)
                  BoundaryConvexEdge(NumBEdges)%face(1) = nf
                  BoundaryConvexEdge(NumBEdges)%face(2) = nf1
                  BoundaryConvexEdge(NumBEdges)%node(1)%name = nodi
                  BoundaryConvexEdge(NumBEdges)%node(2)%name = nodj
                  length2 = zero
                  do sd = 1, SPACEDIM
                    BoundaryConvexEdge(NumBEdges)%node(1)%GX(sd) = Vertice(sd,nodi)
                    BoundaryConvexEdge(NumBEdges)%node(2)%GX(sd) = Vertice(sd,nodj)
                    delta = Vertice(sd,nodj) - Vertice(sd,nodi)
                    BoundaryConvexEdge(NumBEdges)%component(sd) = delta
                    length2 = length2 + delta * delta
                  end do
                  BoundaryConvexEdge(NumBEdges)%length = Dsqrt(length2)
!$omp end critical (omp_FBCE3D)
!AA504 omp directives                  
                  EdgeFound = .True.
                  Exit 
                end if
              end if
            end do
            if (EdgeFound) Exit
          end do
          if (EdgeFound) Exit
        end if
      end do
    end if
 end do
!AA504 omp directives
!$omp end parallel do  

return
end subroutine FindBoundaryConvexEdges3D
!---split

!cfile FindBoundaryIntersection2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : FindBoundaryIntersection2D
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
! Module purpose : Module for Find the intersection segment between the circle of
!                  influence of particle i, whose local coordinates are
!                  xpi=LocXY(1, icbs) and ypi=LocXY(2, icbs), and the straight
!                  boundary side iside=Cloboside(icbs), which lies on the local
!                  x-axis and extends from x=0 to bsidelen = BoundarySide(iside)%Length 
!
! Calling routine: ComputeBoundaryDataTab
!
! Called routines: 
!
!************************************************************************************
!
subroutine FindBoundaryIntersection2D (icbs, Cloboside, LocXY, BoundarySide, &
                                       xpmin, xpmax, interlen)

!Finds the intersection segment between the circle of influence of particle i,
!whose local coordinates are xpi=LocXY(1, icbs) and ypi=LocXY(2, icbs), and the
!straight boundary side iside=Cloboside(icbs), which lies on the local x-axis and extends
!from x=0 to bsidelen = BoundarySide(iside)%Length 
!Returns:
!xpmin =        Minimum abscissa of intersected segment
!xpmax =        Maximum abscissa of intersected segment
!interlen =        Length of intersected segment
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
integer(4)                                                       :: icbs
integer(4),           dimension(1:MAXCLOSEBOUNDSIDES)            :: Cloboside
double precision,     dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES) :: LocXY
type (TyBoundarySide),dimension(1:NumBSides)                     :: BoundarySide
double precision                                                 :: xpmin
double precision                                                 :: xpmax
double precision                                                 :: interlen
!
!.. Local Scalars ..
integer(4)       :: iside
double precision :: xpi, ypi, ypiq, bsidelen, halfchord, minlen, eps, XA, XB
!
!.. Local Arrays ..
!integer(4),dimension(1:PLANEDIM) :: acix
!
!.. Executable Statements ..
!
!  acix(1)=1        !active coordinate indexes
!  acix(2)=3
!
  eps = 0.01d0
  minlen = eps*Domain%h
  interlen = zero
!
  xpi = LocXY(1, icbs)
  ypi = LocXY(2, icbs)
!
  if (ypi < zero) then    !if the particle is out of the boundary operates as if the particle
    ypi = zero            !where on the boundary
  end if    
  ypiq = ypi * ypi
  halfchord = Dsqrt(doublesquareh - ypiq)
  XA = xpi - halfchord
  XB = xpi + halfchord
  iside = Cloboside(icbs)
  bsidelen = BoundarySide(iside)%Length
  xpmin = zero
  if (XA > xpmin) then
    xpmin = XA
  end if
  xpmax = bsidelen
  if (XB < xpmax) then
    xpmax = XB
  end if
  interlen = xpmax - xpmin
  if (interlen <= minlen) then
    interlen = zero               !Intersection is too small
  end if
!
return
end subroutine FindBoundaryIntersection2D
!---split

!cfile FindCloseBoundaryFaces3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name    : FindCloseBoundaryFaces3D
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
! Module purpose : Module Loop
!
! Calling routine: ComputeBoundaryDataTab
!
! Called routines: diagnostic
!                  CellIndices
!                  CellNumber
!                  IsPointInternal
!                  ParticleCellNumber
!                  LocalNormalCoordinates
!
!************************************************************************************
!
  subroutine FindCloseBoundaryFaces3D ( npi, Ncbf, Clobface, LocX, Nfzn )

!.. Finds the "close" boundary faces, i.e. those sited at a distance from the particle npi 
!.. less than or equal to 2h (where h is the smoothing length)
!
!.. Returns:
!.. Ncbf                   = Number of close boundary faces
!.. Clobface(1 to Ncbf)    = List of close boundary faces
!.. LocX(1:SPACEDIM, Ncbf) = Local coordinates of particle npi with respect each boundary side
!
!.. The algorithm looks for boundary faces intersected by the cell boxes of the reference frame
!.. sited all around particle npi, and cancels the repeated ones
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
  integer(4),      intent(IN)                                                :: npi
  integer(4),      intent(INOUT)                                             :: Ncbf, Nfzn
  integer(4),      intent(INOUT),dimension(1:Domain%MAXCLOSEBOUNDFACES)             :: Clobface
  double precision,intent(INOUT),dimension(1:SPACEDIM, 1:Domain%MAXCLOSEBOUNDFACES) :: LocX
!
!.. Local Scalars ..
  integer(4)       :: nc, ic, jc, kc, i, j, k, sdi, sdj, nodes, irestocell, fkod      !ni, nj, nk,  
  integer(4)       :: flpini, flp, flpfin, nfpercell, intbf, icbf, nbface, stretch
  double precision :: pin, pinmin, pinmax
  logical          :: Thereis
  character(len=lencard) :: nomsub = "FindCloseBoundaryFaces3D"
!
!.. Local Arrays ..
  double precision,dimension(1:SPACEDIM) :: PXLoc, csi
!
!.. External routines ..
  integer(4), external :: CellNumber, ParticleCellNumber, CellIndices
  logical,    external :: IsPointInternal
!
!.. Executable Statements ..
!
  Clobface = 0
  Ncbf = 0
  Nfzn = 0
  LocX = zero
  pg(npi)%CloseBcOut = 0
!
!.. find the cell number for the current particle
!
  nc = ParticleCellNumber(pg(npi)%coord)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  if (Domain%RKscheme /= 2) then
!    if (nc /= pg(npi)%cella) then
!      write (nout,*) "Iteration=",it_corrente,"nc =",nc,"different from pg%cella= ",pg(npi)%cella
!      call diagnostic (arg1=8,arg2=5,arg3=nomsub)
!    end if
!  end if
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  if (nc <= 0) return
!
!.. find the indices of the cell in the domain grid
!
!  ni = Grid%ncd(1)        ! inutile
!  nj = Grid%ncd(2)        ! inutile
!  nk = Grid%ncd(3)        ! inutile
  irestocell = CellIndices (nc, ic, jc, kc)
!
!.. loops on the cells surrounding the current one
!
  do i = ic - 1, ic + 1
!!!!!!!!!    if (1 <= i .and. i <= ni) then        ! inutile
!
      Do j = jc - 1, jc + 1
!!!!!!!!!        if (1 <= j .and. j <= nj) then        ! inutile
!
          do k = kc - 1, kc + 1
!!!!!!!!!            if (1 <= k .and. k <= nk) then        ! inutile
!
!.. find the cell number of the surrounding cells
!
              nc = CellNumber(i, j, k)
              if (nc == 0) cycle
!
!.. load the number of boundary faces cutting the cell and the initial and final pointers
!
              nfpercell = GCBFPointers(nc, 1)
              if (nfpercell > 0) then
                flpini = GCBFPointers(nc, 2)
                flpfin = flpini + nfpercell - 1
!
!.. loops on the cutting faces
!
                do flp = flpini, flpfin
!
!.. load the face index
!
                  intbf = GCBFVector(flp)
                  thereis = .false.
                  pinmin = zero
                  pinmax = doubleh
!                  pinmin = -doubleh  !prova 2 dic 2008
                  stretch = BoundaryFace(intbf)%stretch
!
                  if (Tratto(stretch)%tipo == "sour") pinmin = -doubleh   !prova 2 dic 2008
!
!.. checks if the face intbf is already included in the Clobface() array
!
                  if (Ncbf > 0) thereis = any(Clobface(1:Ncbf) == intbf)
!
!.. the face is not yet considered
!
                  if (.Not. Thereis) then         
!
!.. evaluate the normal coordinate of the npi particle with respect the face intbf
!
                    nodes = BoundaryFace(intbf)%nodes
                    do sdi = 1, SPACEDIM
                      PXLoc(sdi) = zero
                      do sdj = 1, SPACEDIM
                        PXLoc(sdi) = PXLoc(sdi) + BoundaryFace(intbf)%T(sdj, sdi) * &
                                     (pg(npi)%coord(sdj) - BoundaryFace(intbf)%Node(nodes)%GX(sdj))
                      end do
                    end do
                    pin = PXLoc(3)
                    call LocalNormalCoordinates (PXLoc, csi, intbf)
                    fkod = BoundaryFace(intbf)%nodes - 2
!
!.. the distance between the particle and the face is greater than zero and less then 2h (influence diameter)
!.. the face is considered and it is added to the Clobface() array
!
!!!!primaprova                    if (pin > pinmin .and. pin < pinmax) then !!!nuovo test
                    if (pin >= pinmin .and. pin < pinmax) then !!!nuovo test
!
                      if (Tratto(stretch)%tipo == "sour" .or. Tratto(stretch)%tipo == "velo" .or. &
                           Tratto(stretch)%tipo == "flow") then
                        if (IsPointInternal(fkod, csi))  then  !La proiezione normale della particella npi sul piano
                                                               ! della faccia "iface" è interna alla faccia "iface"
                          Ncbf = ncbf + 1 
                          if ( ncbf <= Domain%MAXCLOSEBOUNDFACES ) then
                            Clobface(Ncbf) = intbf
                            LocX(3, Ncbf) = pin
                            pg(npi)%CloseBcOut = 1
                          else
                            call diagnostic (arg1=8,arg2=6,arg3=nomsub)
                          end if
                        end if
                      else
                        Ncbf = ncbf + 1 
                        if ( ncbf <= Domain%MAXCLOSEBOUNDFACES ) then
                          Clobface(Ncbf) = intbf
                          LocX(3, Ncbf) = pin
                        else
                          call diagnostic (arg1=8,arg2=7,arg3=nomsub)
                        end if
                      end if
!
!!!!primaprova                    else if (pin <= pinmin) then   !AG 05dic2008  aggiungere test se la proiezione della particella sta all'interno della faccia
                    else if (pin < pinmin) then   !AG 05dic2008  aggiungere test se la proiezione della particella sta all'interno della faccia
                      if (IsPointInternal(fkod, csi)) then    !La proiezione normale della particella npi sul piano
                                                              ! della faccia "iface" è interna alla faccia "iface"
                         Nfzn = Nfzn + 1
                      end if
!????
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!                      write (778,*) " Nfzn >",Nfzn," it>",it_corrente
!!                      write (778,*) " doubleh>",doubleh," nfpercell>",nfpercell
!!                      write (778,*) " dcd  >",grid%dcd(1:3)
!!                      write (778,*) " particella>",npi," intbf>",intbf," cella>",nc
!!                      write (778,*) " pin  >",pin,"faccia>",BoundaryFace(intbf)%stretch
!!                      write (778,*) " coord>",pg(npi)%coord(1:3)
!!                      write (778,*) " T    >",BoundaryFace(intbf)%T(1:3, sdi)
!!                      write (778,*) " gx   >",BoundaryFace(intbf)%Node(nodes)%GX(1:3)
!!                      write (778,*) "==============================================================="
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    end if
                  end if
                end do
!
              end if
!!!!!!!!!            end if        ! inutile
          end do
!
!!!!!!!!!        end if        ! inutile
      end do
!
!!!!!!!!!    end if        ! inutile
  end do
!
!.. evaluates the tangent coordinates (r,s) of the particle npi in the local coordinates of the closest faces
!
  if (Ncbf > 0) then
!
!.. loops on the found faces
!
    do icbf = 1, Ncbf
      nbface = Clobface(icbf)
      nodes = BoundaryFace(nbface)%nodes
!
!.. incremento numero particelle vicine al contorno e calcolo massima quota
!$omp critical (numpa)
      BoundaryFace(nbface)%CloseParticles = BoundaryFace(nbface)%CloseParticles + 1
      if (BoundaryFace(nbface)%CloseParticles_maxQuota < pg(npi)%coord(3)) &
                              BoundaryFace(nbface)%CloseParticles_maxQuota = pg(npi)%coord(3)
!$omp end critical (numpa)
!
      do sdi = 1, PLANEDIM
        LocX(sdi, icbf) = zero
!
        do sdj = 1, SPACEDIM
          LocX(sdi, icbf) = LocX(sdi, icbf) + BoundaryFace(nbface)%T(sdj, sdi) * &
                            (pg(npi)%coord(sdj) - BoundaryFace(nbface)%Node(nodes)%GX(sdj))
        end do
      end do
    end do
!
  end if
!
  return
  end subroutine FindCloseBoundaryFaces3D
!---split

!cfile FindCloseBoundarySides2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : FindCloseBoundarySides2D
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
! Module purpose : Module for Finds the "close" boundary sides
!
! Calling routine: ComputeBoundaryDataTab
!
! Called routines: diagnostic
!
!************************************************************************************
!
  subroutine FindCloseBoundarySides2D (npi, Ncbs, Cloboside, LocXY)
!
!.. Finds the "close" boundary sides, i.e. those sited at a distance from particle npi <= 2h (where h is the smoothing length)
!.. Returns:
!.. Ncbs                     = Number of close boundary sides (= 0, 1, 2)
!.. Cloboside(1:Ncbs)        = List of close boundary sides
!.. LocXY(1:PLANEDIM,1:Ncbs) = Local coordinates of particle npi with respect each boundary side (vertex V1)
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
  integer(4),                                                 intent(in)  :: npi
  integer(4),                                                 intent(out) :: Ncbs
  integer(4),      dimension(1:MAXCLOSEBOUNDSIDES),           intent(out) :: Cloboside
  double precision,dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES),intent(out) :: LocXY
!
!.. Local Scalars ..
  integer(4)       :: icbs, icb, isous, isi, v1, v2, sd, iside1, iside2  !mate, 
  double precision :: xp, yp, sidel, xpmin, xpmax, ypmin, ypmax, xpq
  double precision :: Lmxpq, sidelen, ypmn, ypmx
  character(len=lencard) :: nomsub = "FindCloseBoundarySides2D"
!
!.. local Arrays ..
  integer(4),      dimension(1:PLANEDIM) :: acix
  double precision,dimension(1:PLANEDIM) :: Plocal, P1, P1P, sss, nnn
!
!.. Executable Statements ..
!
!.. initializations ..
!
  acix(1) = 1  
  acix(2) = 3
  Cloboside = 0
  LocXY     = zero
  Ncbs  = 0
  Plocal(:) = pg(npi)%coord(acix(:))
!  mate  = pg(npi)%imed
  ypmin = - doubleh
  ypmax = doubleh
  pg(npi)%CloseBcOut = 0
!
!.. loops on all the boundary sides of the domain
!
  side_loop: do isi = 1, NumBSides
!
!.. the perimeter and pool boundary types are skipped
!
    if (BoundarySide(isi)%tipo /= "peri" .AND. BoundarySide(isi)%tipo /= "pool" ) then
!
!.. loads the side coordinates P1, the local versors sss and nnn and
!.. evaluates the distances between the current particle coordinates Plocal and the 
!.. reference vertex V1
!    
      v1 = BoundarySide(isi)%Vertex(1)
      v2 = BoundarySide(isi)%Vertex(2)
!
      do sd = 1, PLANEDIM
!
        P1(sd) = Vertice(acix(sd),v1)
        P1P(sd) = Plocal(sd) - P1(sd)
        sss(sd) = BoundarySide(isi)%T(acix(sd), acix(1))
        nnn(sd) = BoundarySide(isi)%T(acix(sd), acix(2))
!
      end do
!
!.. evaluates the particle coordinates xp and yp in the side local reference system, where
!.. x in the direction of the side segment, y is the normal direction and V1 the origin
!
      sidel = BoundarySide(isi)%length
!
      xp = P1P(1) * sss(1) + P1P(2) * sss(2)
      yp = P1P(1) * nnn(1) + P1P(2) * nnn(2)
!
!.. set the interaction area in the neighborough of the side segment having a distance equal to
!.. +/- 2*h in the y local direction and a distance of 2h from the vertices V1 and V2 in the x local direction
! 
      xpmin = - doubleh
      xpmax = sidel + doubleh
!
!.. checks if the particle has local coordinates falling inside the interaction area
!        
      if ((xpmin < xp .AND. xp < xpmax) .AND. (ypmin < yp .AND. yp < ypmax)) then
!
!.. the particle falls in the extreme cicle around V1
!            
        if (xp < zero) then
!
          xpq = xp * xp
          ypmx = Dsqrt(doublesquareh - xpq)
          ypmn = -ypmx 
!
!.. the particle falls in the segment natural length
!
        else if (xp <= sidel) then
!
          ypmx = ypmax    !CONTROLLARE
          ypmn = ypmin 
!
!.. the particle falls in the extreme cicle around V2
!
        else if (xp < xpmax) then
!
          Lmxpq = (sidel - xp) * (sidel - xp)
          ypmx = Dsqrt(doublesquareh - Lmxpq)
          ypmn = -ypmx 
!
        end if
!
!.. the boundary must be considered for the current particle as close boundary
!            
        if (ypmn < yp .AND. yp < ypmx) then
!
!.. the number of close boundaries is increased
!
          Ncbs = Ncbs + 1
!
!.. checks against the maximum number allowed for the closest boundaries
!
          if (Ncbs > MAXCLOSEBOUNDSIDES) then
!.. The particle npi has more than two boundary sides: reduce dd and restart
            write (nout,'(1x,a,i12,a)')   ' ERROR! The particle ',npi,' has more than two boundary sides: reduce dd and restart.'
            write (nout,'(1x,a,3f15.10)') '        Coordinate: ',pg(npi)%coord(:)
            write (nscr,'(1x,a,i12,a)')   ' ERROR! The particle ',npi,' has more than two boundary sides: reduce dd and restart.'
            write (nscr,'(1x,a,3f15.10)') '        Coordinate: ',pg(npi)%coord(:)
            call diagnostic (arg1=8,arg2=8,arg3=nomsub)
          end if
!
!.. stores the boundary index and the local coordinates on the boundary segment
!
          Cloboside(Ncbs) = isi
          LocXY(1, Ncbs) = xp
          LocXY(2, Ncbs) = yp
!
        end if
      end if
    end if
  end do side_loop
!
!.. searches for a nearest source type side
!
  isous = 0
  do icbs = 1, Ncbs
    isi = Cloboside(icbs)
!
!.. incremento numero particelle vicine al contorno e calcolo massima quota
!$omp critical (numpa)
    BoundarySide(isi)%CloseParticles = BoundarySide(isi)%CloseParticles + 1
    if (BoundarySide(isi)%CloseParticles_maxQuota < pg(npi)%coord(3)) BoundarySide(isi)%CloseParticles_maxQuota = pg(npi)%coord(3)
!$omp end critical (numpa)
!
    if (BoundarySide(isi)%tipo == "sour") then
      XP = LocXY(1, icbs)
      sidel = BoundarySide(isi)%length
      if (XP > zero .AND. XP < sidel) then
        isous = icbs    !contrassegna la sorgente per eliminare gli altri eventuali lati vicini
        exit
      else                !cancella il lato sorgente perché ininfluente
        if (icbs < Ncbs) then
          do icb = icbs+1, Ncbs
            Cloboside(icb -1) = Cloboside(icb)
            LocXY(1, icb -1) = LocXY(1, icb)
            LocXY(2, icb -1) = LocXY(2, icb)
          end do
          Ncbs = Ncbs -1
          exit
        else if (icbs == Ncbs) then
          Ncbs = Ncbs -1
          exit
        end if
      end if
    end if
  end do
!
!.. a source has been found: the other nearest sides are deleted
!
  if (isous > 0) then   
    Ncbs = 1
    Cloboside(Ncbs) = Cloboside(isous)
    LocXY(1, Ncbs) = LocXY(1, isous)
    LocXY(2, Ncbs) = LocXY(2, isous)
  end if
!
!.. checks if more than two sides are close the particle
!
  if (Ncbs > LIMCLOSEBOUNDSIDES) then
    call diagnostic (arg1=8,arg2=9,arg3=nomsub)
  end if
!
!.. there are two close boundaries:
!
!Case of two close boundary sides
!Change of the origin of local abcsissa in one of the two adjacent boundary sides
!in order that the abscissa origin coincide with the common vertex in each side.
!
  if (Ncbs == 2) then
!
    iside1 = Cloboside(1)
    iside2 = Cloboside(2)
!
    if (BoundarySide(iside1)%previous_side == iside2) then
      sidelen = BoundarySide(iside2)%length
      LocXY(1, 2) = sidelen - LocXY(1, 2)
    else if (BoundarySide(iside2)%previous_side == iside1) then
      sidelen = BoundarySide(iside1)%length
      LocXY(1, 1) = sidelen - LocXY(1, 1)
    else
!     Sides not adjacent: this case should never occur if the minimum length
!      of boundary sides is >4h (diameter of influence circle)
    end if
  end if
!
  return
  end subroutine FindCloseBoundarySides2D
!---split

!cfile GridCellBoundaryFacesIntersections3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : GridCellBoundaryFacesIntersections3D
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!AA504
! 03  Amicarelli        08/04/2014     (v5.04) omp parallelization and optimization              
!
!************************************************************************************
! Module purpose : Module to find the boundary faces intercepted by each frame cell
!
! Calling routine: Gest_Trans
!
! Called routines: 
!
!************************************************************************************
!
  subroutine GridCellBoundaryFacesIntersections3D ( NumCellmax )

!Ricerca le facce di contorno intersecate da ciascuna cella della griglia di riferimento nc [1, NumCells]
!Nella riga generica nc del vettore CFBFPointers(1 to NumCells,1 to 2) mette
!nella prima colonna il numero delle facce intersecate,
!nella seconda colonna il puntatore alla posizione del vettore CFBFVector()
!dove inizia la lista degli indici delle facce intersecate.
!La ricerca si base su un principio di esclusione e si svolge in due fasi:
!nella prima fase per ogni cella esclude (come possibili intersecate) le facce i cui nodi giacciono
!tutti in uno dei semispazi, definiti dai piani delle facce della cella, che non includono la cella stessa;
!nella seconda fase per ciascuna faccia, non esclusa nella prima fase, verifica che tutti gli otto
!vertici della cella siano tutti contenuti in uno dei semispazi definiti dal piano della faccia,
!nel qual caso la faccia viene esclusa (come possibile intersecata).

!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
  use DIAGNOSTIC_MODULE
!AA504
  use files_entities

!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
  integer(4),intent(IN) :: NumCellmax
!
!.. Local Scalars ..
!AA504 sub
  integer(4)       :: nc, nf, kf, i, j, k, i0, j0, k0, flpointer, nodes, no, sd, ier, flpointer_cell, i_flpointer
  integer(4)       :: irestocell
  integer(4)       :: ii, jj, kk, nv
  double precision :: XYZnod, CellNodeZita, deltaXYZ
  logical          :: Found
  character(len=lencard)  :: nomsub = "GridCellBoundaryFacesIntersections3D"
!
!.. Local Arrays
  double precision,dimension(1:SPACEDIM, 1:2) :: CellXYZ
  double precision,dimension(1:8, 1:SPACEDIM) :: CellNodeXYZ
  integer(4),      dimension(-2:2)            :: Signcount
!AA504 
  integer(4),dimension(:),allocatable         :: GCBFVector_aux,GCBFVector_cell

!.. external declarations
  integer(4), external :: CellIndices

!AA504 start
!Allocating GCBFVector_aux 
  allocate(GCBFVector_aux(GCBFVecDim))    
!AA504 end   
 
!.. Executable Statements ..

!AA504 rm comments

  if (NumCellmax < Grid%nmax) call diagnostic (arg1=7,arg3=nomsub)            
!
  flpointer = 0
!
!.. loops on all the cells of the grid
!

!AA504 omp directives
!$omp parallel do default(none) &
!$omp shared(Grid,GCBFPointers,NumFacce,BFaceList,Tratto,BoundaryFace,flpointer,GCBFVector_aux,nomsub,Domain,nout) &
!$omp private(nc,irestocell,i,i0,j,j0,k,k0,CellXYZ,nv,ii,jj,kk,CellNodeXYZ,kf,nf,nodes,Found,sd) &
!$omp private(Signcount,no,XYZnod,CellNodeZita,deltaXYZ,i_flpointer,flpointer_cell,GCBFVector_cell)
  Do nc = 1, Grid%nmax
!
!AA504 sub      
    GCBFPointers(nc, 1) = zero
!
!.. detects the indices of the cell in the grid
!
    irestocell = CellIndices (nc, i, j, k)        
!
!.. evaluates the minimum and maximum coordinates of the cell in the XYZ system
!
    i0 = i - 1
    j0 = j - 1
    k0 = k - 1
    CellXYZ(1, 1) = Grid%extr(1,1) + Grid%dcd(1)*i0
    CellXYZ(1, 2) = Grid%extr(1,1) + Grid%dcd(1)*i
    CellXYZ(2, 1) = Grid%extr(2,1) + Grid%dcd(2)*j0
    CellXYZ(2, 2) = Grid%extr(2,1) + Grid%dcd(2)*j
    CellXYZ(3, 1) = Grid%extr(3,1) + Grid%dcd(3)*k0
    CellXYZ(3, 2) = Grid%extr(3,1) + Grid%dcd(3)*k
!
!.. evaluates the coordinate triplets of all the cell vertices
!
    nv = 0
    do ii = 1, 2
      do jj = 1, 2
        do kk = 1, 2
          nv = nv + 1
          CellNodeXYZ(nv, 1) = CellXYZ(1, ii)
          CellNodeXYZ(nv, 2) = CellXYZ(2, jj)
          CellNodeXYZ(nv, 3) = CellXYZ(3, kk)
        end do
      end do
    end do
    
!AA504 start
!Allocating GCBFVector_cell and initializing the private variable flpointer_cell
  allocate(GCBFVector_cell(Domain%MAXCLOSEBOUNDFACES))
  flpointer_cell = 0
!AA504 end 
    
!
!.. loops on all the boundary faces assigned to the domain 
!
    Do kf = 1, NumFacce                                    
!
      nf = BFaceList(kf)
!
!.. skip the perimeter and pool conditions
!
      if (Tratto(BoundaryFace(nf)%stretch)%tipo /= "peri" .AND. Tratto(BoundaryFace(nf)%stretch)%tipo /= "pool") then
!
        nodes = BoundaryFace(nf)%nodes
!
            !Prima fase di esclusione

        Found = .True.
!
!.. loops on the system coordinates
!
        Do sd = 1, SPACEDIM
!
          Signcount(-2:2) = 0
!
!.. loops on the nodes of the face (3 or 4)
!
          Do no = 1, nodes
!
            XYZnod = BoundaryFace(nf)%Node(no)%GX(sd)
            if (XYZnod < CellXYZ(sd, 1)) then
              Signcount(-2) = Signcount(-2) + 1
            else if (XYZnod == CellXYZ(sd, 1)) then
              Signcount(-1) = Signcount(-1) + 1
            else if (XYZnod < CellXYZ(sd, 2)) then
              Signcount(0) = Signcount(0) + 1
            else if (XYZnod == CellXYZ(sd, 2)) then
              Signcount(1) = Signcount(1) + 1
            else if (XYZnod > CellXYZ(sd, 2)) then
              Signcount(2) = Signcount(2) + 1
            end if
!
          end do
!
          if ((Signcount(-2) + Signcount(-1) == nodes) .And..Not. (Signcount(-1) == nodes)) then
            Found = .False.
            Exit
          else if ((Signcount(2) + Signcount(1) == nodes) .And..Not. (Signcount(1) == nodes)) then
            Found = .False.
            Exit
          end if
!
        end do
!
!.. Seconda fase di esclusione
!
        if (Found) then
          SignCount(-2:2) = 0
!
!.. loops on the ??
!
          do nv = 1, 8
!
!.. evaluates the normal component of the coordinates of each vertex of the cell nc-th referred to the face nf
!
            CellNodeZita = zero
!
            do sd = 1, SPACEDIM
               deltaXYZ = CellNodeXYZ(nv, sd) - BoundaryFace(nf)%node(nodes)%GX(sd)
               CellNodeZita = CellNodeZita + BoundaryFace(nf)%T(sd, 3) * deltaXYZ
            end do
!
!.. check for the normal components
!
            if (CellNodeZita < zero) then
              SignCount(-1) = SignCount(-1) + 1
            else if (CellNodeZita == zero) then
              SignCount(0) = SignCount(0) + 1
            else
              SignCount(1) = SignCount(1) + 1
            end if
!
          end do
!                
          if ((SignCount(-1) + SignCount(0) == 8) .And. .Not. (SignCount(0) == 4)) then
            Found = .False.
          else if ((SignCount(1) + SignCount(0) == 8) .And. .Not. (SignCount(0) == 4)) then
            Found = .False.
          end if
!
!.. final check. It found = .true. the face nf cut the cell nc-th and it is counted
!
          if (Found) then        

!AA504 sub start
             flpointer_cell = flpointer_cell + 1
             GCBFVector_cell(flpointer_cell) = nf
             GCBFPointers(nc, 1) = GCBFPointers(nc, 1) + 1
             if (flpointer_cell>Domain%MAXCLOSEBOUNDFACES) then
                 write (nout,'(1x,a)') " Too many faces crossing a given cell. Please increase the parameter MAXCLOSEBOUNDFACES. "
                 call diagnostic (arg1=4,arg3=nomsub)
             endif    
!AA504 sub end            

          end if
        end if
      end if
    end do

!AA504 sub start
!$omp critical
    do i_flpointer=1,flpointer_cell    
       flpointer = flpointer + 1
       GCBFVector_aux(flpointer) = GCBFVector_cell(i_flpointer)
       if (i_flpointer==1) GCBFPointers(nc,2) = flpointer
    end do
!$omp end critical
!Deallocate auxiliary array
    deallocate(GCBFVector_cell)
!AA504 end 

  end do
!$omp end parallel do 
!AA504 omp directives
!
!.. save the dimension of the GCBFVector array
!
  GCBFVecDim = flpointer

!AA504 start
!Allocating GCBFVector 
  allocate(GCBFVector(GCBFVecDim),stat=ier)    
  if (ier/=0) then
     write (nout,'(1x,a,i2)') "   Array GCBFVector not allocated. Error code: ",ier
     call diagnostic (arg1=4,arg3=nomsub)
     else
        write (nout,'(1x,a)') "   Array GCBFVector successfully allocated "
  end if
!Loop over the reference GCBFVector array
!$omp parallel do default(none) shared(GCBFVector,GCBFVector_aux,GCBFVecDim) private(i)
  do i=1,GCBFVecDim
!Copy auxiliary array in reference array
     GCBFVector(i) = GCBFVector_aux(i)
  end do
!$omp end parallel do  
!Deallocate auxiliary array
  deallocate(GCBFVector_aux)
!AA504 end 
          
  return
  end subroutine GridCellBoundaryFacesIntersections3D
!---split

!cfile InterpolateBoundaryIntegrals2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : InterpolateBoundaryIntegrals2D
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
! Module purpose : Module for Interpolation in table "BoundIntegralTab(:,:)", 
!                  defined in module "BoundIntegralTab_Module", the values in 
!                  columns "Colmn(nc), nc=1, Ncols" corresponding to the input
!                  value "x" to be interpolated, in turn, in column 0.
!
! Calling routine: ComputeBoundaryDataTab
!
! Called routines: 
!
!************************************************************************************
!
subroutine InterpolateBoundaryIntegrals2D (x, Ncols, Colmn, Func)

!Interpolates in table "BoundIntegralTab(:,:)", defined in module "BoundIntegralTab_Module",  
!the values in columns "Colmn(nc), nc=1, Ncols" corresponding to the input value "x" to be
!interpolated, in turn, in column 0.
!Returns:
!Func(nc), nc=1, Ncols    =    Values interpolated in columns Col(nc), nc=1, Ncols
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use BoundIntegralTab_Module
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
double precision                          :: x 
integer(4)                                :: Ncols
integer(4),      dimension(1:NUMCOLS_BIT) :: Colmn
double precision,dimension(1:NUMCOLS_BIT) :: Func
!
!.. Local Scalars ..
integer(4)       :: nc, i, j
double precision :: xi, fi, fip1   !xip1, 
!
!.. Executable Statements ..
!
  if (x < BoundIntegralTab2D(1,0)) x = BoundIntegralTab2D(1,0)
!
  i = Int((x-BoundIntegralTab2D(1,0))/DELTAX_BIT)+1
  xi = BoundIntegralTab2D(i,0)
!  xip1 = BoundIntegralTab2D(i+1,0)

  do nc = 1, Ncols
    j = Colmn(nc)
    fi = BoundIntegralTab2D(i,j)
    fip1 = BoundIntegralTab2D(i+1,j)
    Func(nc) = fi+(fip1-fi)*(x-xi)/DELTAX_BIT
  end do
!
return
end subroutine InterpolateBoundaryIntegrals2D
!---split

!cfile InterpolateTable.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine InterpolateTable (xval, nicols, icol, ivalue)
!
!Interpolates values in the array "Table()" with "nrows" rows and "ncols" columns
!Independent variable is in column 0 of Table()
!nicols = number of colums of dependent variables to be interpolated
!icol() = listof columns of dependent variables to be interpolated
!ivalue() = list of the "nicols" interpolated values
!
use Global_MODULE
!
implicit none
!
integer(4)       :: nicols, nr0, nc, ic, nr1
double precision :: xval, xval0, csi, deltaX
!
integer(4),      dimension(1:nicols)         :: icol
double precision,dimension(1:nicols)         :: ivalue
!
deltaX = ktdelta
xval0 = kerneltab(0, 0)
nr0 = Int((xval - xval0) / deltaX)
!
if (nr0 <= 0) then
  Do ic = 1, nicols
    nc = icol(ic)
    ivalue(ic) = kerneltab(0, nc)
  end do
else if (nr0 >= ktrows) then
  Do ic = 1, nicols
    nc = icol(ic)
    ivalue(ic) = kerneltab(ktrows, nc)
  end do
Else
  nr1 = nr0 + 1
  csi = (xval - kerneltab(nr0, 0)) / deltaX
  Do ic = 1, nicols
    nc = icol(ic)
    ivalue(ic) = kerneltab(nr0, nc) + csi * (kerneltab(nr1, nc) - kerneltab(nr0, nc))
  end do
end if
!
return
end subroutine InterpolateTable
!---split

!cfile IWro2dro.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
double precision Function IWro2dro(ro)

!Computes the definite integral
!2
!S W*(ro')ro2' dro'
!ro
!
!.. assign modules
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
double precision,parameter :: a1 = 0.333333333333333d0 != 1 / 3
double precision,parameter :: a2 = 0.3d0               != 3 / 10
double precision,parameter :: a3 = 0.125d0             != 1 / 8
double precision,parameter :: a4 = 0.25d0
double precision,parameter :: b1 = 0.0625d0             != 1 / 16
double precision,parameter :: b2 = 0.025d0              != 1 / 40
double precision,parameter :: b3 = 4.16666666666667d-03 != 1 / 240
!
!.. Formal Arguments ..
double precision,intent(IN) :: ro
!
!.. Local Scalars ..
double precision :: ro2,ro3,duemro,duemro2,duemro4,IWro2dro1,IWro2dro2
!
!.. Executable Statements ..
!
if (ro >= zero .And. ro < one) then
!    0.333333333333333 = 1 / 3
!    0.3               = 3 / 10
!    0.125             = 1 / 8
  ro2 = ro * ro
  ro3 = ro2 * ro
  IWro2dro1 = (a1 - a2 * ro2 + a3 * ro3) * ro3
  IWro2dro = KERNELCONST3D * (a4 - IWro2dro1)
else if (ro >= one .And. ro < two) then
!    0.0625               = 1 / 16
!    0.025                = 1 / 40
!    4.16666666666667d-03 = 1 / 240
  ro2 = ro * ro
  duemro = two - ro
  duemro2 = duemro * duemro
  duemro4 = duemro2 * duemro2
  IWro2dro2 = -(b1 * ro2 + b2 * duemro * ro + b3 * duemro2) * duemro4
  IWro2dro = KERNELCONST3D * (-IWro2dro2)
Else
  IWro2dro = zero
end if
!
return
End Function IWro2dro
!---split

!cfile J2Wro2.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : J2Wro2
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
! Module purpose : Module to compute the boundary volume integrals IntWdV only
!
! Calling routine: BoundaryVolumeIntegrals2D
!
! Called routines: 
!
!************************************************************************************
!
double precision function J2Wro2 (ro)
!
!.. assign modules
use GLOBAL_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
double precision,parameter :: a1 = 0.00833333333333333d0 != 1 / 120
double precision,parameter :: a2 = 0.26666666666666667d0 != 8 / 30
!
!.. Formal Arguments ..
double precision,intent(IN) :: ro
!
!.. Local Scalars ..
double precision :: ro2, ro3
!
!.. Executable Statements ..
!
  ro2 = ro * ro
  ro3 = ro2 * ro
!
  if (zero <= ro .and. ro < one) then
    J2Wro2 = KERNELCONST2D * (0.25d0 - (a1 * (40.0d0 - 36.0d0 * ro2 + 15.0d0 * ro3) * ro3))
  else if (one <= ro .and. ro < two) then
    J2Wro2 = KERNELCONST2D * (a2 - (a1 * (80.0d0 - 90.0d0 * ro + 36.0d0 * ro2 - 5.0d0 * ro3) * ro3))
  else
    J2Wro2 = zero
  end if
!
return
End function J2Wro2
!---split

!cfile JdWsRn.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : JdWsRn
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
! Module purpose : Module to compute boundary contributions to rodivV
!
! Calling routine: ComputeKernelTable
!
! Called routines: 
!
!************************************************************************************
!
double precision Function JdWsRn(ro, SD, n, kernel)
!
!.. assign modules
  use GLOBAL_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
  integer(4),      intent(IN)      :: SD, n, kernel
  double precision,intent(IN)      :: ro
!
!.. Parameters
  double precision,parameter :: KS2D = 0.454728408833987d0 != 10 / (7*pigreco)
  double precision,parameter :: KA2D = 0.09947192d0        !=5./(16*pigreco)
  double precision,parameter :: KS3D = 0.31830989d0        != 1 / pigreco
  double precision,parameter :: KA3D = 0.07460394d0        !=15./(16*pigreco)
!
!..  0.1875                     = 3/16
!..  0.5625                     = 9/16
!
!.. Local Scalars ..
  double precision :: ro2, ro3, ro4, ro5, duemro
!
!.. Executable Statements ..
!
  ro2 = ro * ro
  ro3 = ro2 * ro
  ro4 = ro2 * ro2
  ro5 = ro4 * ro
!
  JdWsRn = zero
  Select Case (kernel)
    Case (1)              !kernel spline cubica (Monagan)
    
        Select Case (n)
            Case (0)      !n = 0
                if (ro >= zero .And. ro < one) then
                  JdWsRn = -one + (1.5d0 - 0.75d0 * ro) * ro2
                else if (ro >= one .And. ro < two) then
                  duemro = two - ro
                  JdWsRn = -0.25d0 * duemro * duemro * duemro
                end if
            Case (1)      !n = 1
                if (ro >= zero .And. ro < one) then
                  JdWsRn = -0.75d0 + (one - 0.5625d0 * ro) * ro3
                else if (ro >= one .And. ro < two) then
                  JdWsRn = -one + (1.5d0 - one * ro + 0.1875d0 * ro2) * ro2
                end if
            Case (2)      !n = 2
                if (ro >= zero .And. ro < one) then
                  JdWsRn = -0.7d0 + (0.75d0 - 0.45d0 * ro) * ro4
                else if (ro >= one .And. ro < two) then
                  JdWsRn = -0.8d0 + (one - 0.75d0 * ro + 0.15d0 * ro2) * ro3
                end if
            Case (3)      !n = 3
                if (ro >= zero .And. ro < one) then
                  JdWsRn = -0.75d0 + (0.6d0 - 0.375d0 * ro) * ro5
                else if (ro >= one .And. ro < two) then
                  JdWsRn = -0.8d0 + (0.75d0 - 0.6d0 * ro + 0.125d0 * ro2) * ro4
                end if
            case default
                JdWsRn = zero
        End Select
        
        Select Case (SD)
            Case (2)      !geometria 2D
                JdWsRn = JdWsRn * KS2D
            Case (3)      !geometria 3D
                JdWsRn = JdWsRn * KS3D
            case default
                JdWsRn = zero
        End Select
        
    Case (2)              !kernel anticluster cubica (Gallati)
    
        Select Case (n)
            Case (1)      !n = 1
                if (ro < two) then
                  JdWsRn = -4.0d0 + (6.0d0 - 4.0d0 * ro + 0.75d0 * ro2) * ro2
                end if
            Case (2)      !n = 2
                if (ro < two) then
                  JdWsRn = -3.2d0 + (4.0d0 - 3.0d0 * ro + 0.6d0* ro2) * ro3
                end if
            Case (3)      !n = 3
                if (ro < two) then
                  JdWsRn = -3.2d0 + (3.0d0 - 2.4d0 * ro + 0.5d0 * ro2) * ro4
                end if
            case default
                JdWsRn = zero
        End Select
        
        Select Case (SD)
            Case (2)      !geometria 2D
                JdWsRn = JdWsRn * KA2D
            Case (3)      !geometria 3D
                JdWsRn = JdWsRn * KA3D
            case default
                JdWsRn = zero
        End Select

    case default
        
  End Select
!
return
End Function JdWsRn
!---split

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

!cfile w.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!* kernel standard
double precision function w (r,h,coef)
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
double precision,parameter :: a1 = 0.666666667d0
!
!.. Local Scalars ..
double precision :: r,h,s,q,dms,coef
!
!.. Executable Statements ..
!
 s = r / h

 if ( s <= 1.0d0 ) then
    q = a1 + s * s * (s * 0.5d0 - 1.0d0)
 else if (s >= 2.0d0) then
    q = 0.0d0
 else   ! if (s>1.0d0 .and. s<2.0d0) then
    dms = 2.0d0 - s
    q = dms * dms * dms / 6.0d0
 end if

! coef = 0.682093d0/(h*h)

 w = q * coef

return
end function w
!---split

!cfile WIntegr.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : WIntegr
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
! Module purpose : Module for Computing the definite integral
!
! Calling routine: ComputeVolumeIntegral_WdV2D
!
! Called routines: 
!
!************************************************************************************
!
double precision Function WIntegr (ri, h) 

!Computes the definite integral
!2h
!S W(r,h)rdr
!ri
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
double precision,parameter :: a   =  0.5d0
double precision,parameter :: a2  = -0.375d0
double precision,parameter :: a3  =  0.15d0
double precision,parameter :: aa  = -0.125d0
double precision,parameter :: aa1 =  0.05d0
double precision,parameter :: b   =  0.35d0
!
!.. Formal Arguments ..
double precision,intent(IN) :: ri
double precision,intent(IN) :: h
!
!.. Local Scalars ..
double precision :: WIntegr1, WIntegr2, ro, ro2, ro3, q1, q2, q4
!
!.. Executable Statements ..
!

!KERNELCONST2D = 0.454728408833987d0     =10/(7*pigreco) definita in AdM_USER_TYPE

  ro = ri / h

  if ( zero <= ro .AND. ro < one ) then
    ro2 = ro * ro
    ro3 = ro2 * ro
    WIntegr1 = (a + a2 * ro2 + a3 * ro3) * ro2
    WIntegr = KERNELCONST2D * (b - WIntegr1)
  else if ( one <= ro .AND. ro < two ) then
    q1 = two - ro
    q2 = q1 * q1
    q4 = q2 * q2
    WIntegr2 = (aa + aa1 * q1) * q4
    WIntegr = -KERNELCONST2D * WIntegr2
  Else
    WIntegr = zero
  end if

return
End Function WIntegr
!---split
