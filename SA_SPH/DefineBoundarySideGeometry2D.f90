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

