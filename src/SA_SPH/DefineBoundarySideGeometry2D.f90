!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-) 
!      
!     
!   
!      
!  

! This file is part of SPHERA.
!  
!  
!  
!  
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
!  
!  
!  
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: DefineBoundarySideGeometry2D                                     
! Description: Definition of the boundary sides. 
!              (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine DefineBoundarySideGeometry2D
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: NumBS,nt,ns,nst,nv,v1,v2
double precision :: Dx,Dz,L,sx,sz,PsiS,PsiN,FiS,FiN  
double precision,dimension(1:SPACEDIM,1:SPACEDIM) :: TT,RR,RRN
character(4) :: tipo
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
NumBS = 0
!------------------------
! Statements
!------------------------
do nt=1,NumTratti
   nst = Tratto(nt)%numvertices - 1
   nv = Tratto(nt)%inivertex
   Tratto(nt)%iniside = NumBS + 1
   do ns=1,nst
      v1 = BoundaryVertex(nv)
      v2 = BoundaryVertex(nv + 1)
      NumBS = NumBS + 1
      Dx = Vertice(1,v2) - Vertice(1,v1)
      Dz = Vertice(3,v2) - Vertice(3,v1)
      L = Dsqrt(Dx * Dx + Dz * Dz)
      sx = Dx / L
      sz = Dz / L
      BoundarySide(NumBS)%stretch = nt
      BoundarySide(NumBS)%tipo = Tratto(nt)%tipo
      BoundarySide(NumBS)%Vertex(1) = v1
      BoundarySide(NumBS)%Vertex(2) = v2
      BoundarySide(NumBS)%length = L
      BoundarySide(NumBS)%T(1,1) = sx
      BoundarySide(NumBS)%T(3,1) = sz
      BoundarySide(NumBS)%T(1,3) =-sz
      BoundarySide(NumBS)%T(3,3) = sx
      TT=BoundarySide(NumBS)%T
      tipo = Tratto(Nt)%tipo   
      select case (tipo)
         case ("fixe","tapi")
            PsiS = zero
            PsiN = zero
            FiS = one
            FiN = one 
         case ("leve")
            PsiS = zero
            PsiN = zero  
            FiS = one
            FiN = one
         case ("velo","flow")
            PsiS = zero
            PsiN = zero
            FiS = zero        
            FiN = zero
         case ("crit")
            PsiS = zero
            PsiN = zero
            FiS = one
            FiN = one
         case ("open")
            PsiS = zero
            PsiN = zero  
            FiS = zero 
            FiN = zero   
         case ("sour")
            PsiS = zero
            PsiN = zero
            FiS = zero
            FiN = zero
      end select                                                        
      call BoundaryReflectionMatrix2D(TT,RR,PsiS,PsiN)
      call BoundaryMassForceMatrix2D(TT,RRN,FiS,FiN) 
      BoundarySide(NumBS)%R = RR
      BoundarySide(NumBS)%RN = RRN
      if (Tratto(nt)%tipo=="tapi") then
         BoundarySide(NumBS)%velocity = Tratto(nt)%velocity
         else
            BoundarySide(NumBS)%velocity(:) = zero
      endif
      nv = nv + 1
   enddo
enddo
call DefineBoundarySideRelativeAngles2D
!------------------------
! Deallocations
!------------------------
return
end subroutine DefineBoundarySideGeometry2D

