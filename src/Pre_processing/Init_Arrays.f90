!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.8.0.
! SPHERA v.8.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with SPHERA. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: Init_Arrays                   
! Description:                 
!-------------------------------------------------------------------------------
subroutine Init_Arrays
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
integer(4) :: i,j,n,i1
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
do i=1,size(Med)
   Med(i)%tipo = "Empty   "
   Med(i)%index = 0
   Med(i)%NIterSol = 0
   Med(i)%den0 = zero
   Med(i)%eps = zero
   Med(i)%celerita = zero
   Med(i)%alfaMon = zero
   Med(i)%betaMon = zero
   Med(i)%visc = zero
   Med(i)%coes = zero
   Med(i)%numx = zero
   Med(i)%mumx = zero
   Med(i)%taucri = zero
   Med(i)%cuin = zero
   Med(i)%phi = zero
   Med(i)%cons = zero
   Med(i)%Cs = zero
   Med(i)%RoughCoef = zero
   Med(i)%d50 = zero
   Med(i)%SettlingCoef = zero
   Med(i)%Codif = zero
   Med(i)%Gamma = zero
   Med(i)%InitialIntEn = zero
enddo
do i=1,size(Partz)
   Partz(i)%label = "        "
   Partz(i)%tipo = "    "
   Partz(i)%shape = " "
   Partz(i)%bend = " "
   Partz(i)%pressure = "  "
   Partz(i)%move = "   "
   Partz(i)%slip = " "
   Partz(i)%ipool = 0
   Partz(i)%npoints = 0
   Partz(i)%icol = 0
   Partz(i)%Medium = 0
   Partz(i)%npointv = 0
   Partz(i)%Indix(1) = 0
   Partz(i)%Indix(2) = 0
   Partz(i)%limit(1) = 1
   Partz(i)%limit(2) = 0
   Partz(i)%pool = zero
   do j=1,3
      Partz(i)%coordMM(j,1) = zero
      Partz(i)%coordMM(j,2) = zero
      Partz(i)%vel(j) = zero
   enddo
   do j=0,3
      do n=1,MAXPOINTSVLAW
         Partz(i)%vlaw(j,n) = zero
      enddo
   enddo
   Partz(i)%trampa = zero
   Partz(i)%valp = zero
enddo
do i=1,size(Control_Points)
   Control_Points(i)%cella = 0
   Control_Points(i)%coord(1) = zero
   Control_Points(i)%coord(2) = zero
   Control_Points(i)%coord(3) = zero
   Control_Points(i)%vel(1) = zero
   Control_Points(i)%vel(2) = zero
   Control_Points(i)%vel(3) = zero
   Control_Points(i)%pres = zero
   Control_Points(i)%dens = zero
   Control_Points(i)%uni = zero
   Control_Points(i)%dist = zero
enddo
do i=1,size(Section_Points)
   Section_Points(i)%cella = 0
   Section_Points(i)%coord(1) = zero
   Section_Points(i)%coord(2) = zero
   Section_Points(i)%coord(3) = zero
   Section_Points(i)%vel(1) = zero
   Section_Points(i)%vel(2) = zero
   Section_Points(i)%vel(3) = zero
   Section_Points(i)%pres = zero
   Section_Points(i)%dens = zero
   Section_Points(i)%uni = zero
   Section_Points(i)%dist = zero
enddo
Control_Sections(0)%Label = "Domain  "
Control_Sections(1:NSections)%Label = "        "
Control_Sections(NSections+1)%Label = "One more"
do i1=1,size(Control_Sections)
   i = i1 - 1
   Control_Sections(i)%Tipo = "  "
   Control_Sections(i)%Icont(1) = 0
   Control_Sections(i)%Icont(2) = 0
   Control_Sections(i)%ColorCode = 1
   Control_Sections(i)%Constant(1) = zero
   Control_Sections(i)%Constant(2) = zero
   Control_Sections(i)%Constant(3) = zero
   Control_Sections(i)%XYZRange(1,1) = zero
   Control_Sections(i)%XYZRange(2,1) = zero
   Control_Sections(i)%XYZRange(3,1) = zero
   Control_Sections(i)%XYZRange(1,2) = zero
   Control_Sections(i)%XYZRange(2,2) = zero
   Control_Sections(i)%XYZRange(3,2) = zero
   do j=1,SPACEDIM
      do n=1,SPACEDIM
         Control_Sections(i)%TGLsection(j,n) = zero
      enddo
   enddo
   Control_Sections(i)%TGLsection(1,2) = one
   Control_Sections(i)%TGLsection(2,1) =-one
   Control_Sections(i)%TGLsection(3,3) = one
enddo
do i=1,size(Control_Lines)
   Control_Lines(i)%label = "Empty   "
   Control_Lines(i)%icont(1) = 0
   Control_Lines(i)%icont(2) = 0
enddo
do j=1,SPACEDIM
   Vertice(j,1:NumVertici) = zero
enddo
do i=1,size(BoundaryFace)
   do j=1,MAXFACENODES
      BoundaryFace(i)%Node(j)%name = 0
      do n=1,SPACEDIM
         BoundaryFace(i)%Node(j)%GX(n) = zero
         BoundaryFace(i)%Node(j)%LX(n) = zero
      enddo
   enddo
   BoundaryFace(i)%nodes = 0
   BoundaryFace(i)%stretch = 0
   BoundaryFace(i)%CloseParticles = 0
   BoundaryFace(i)%CloseParticles_maxQuota = const_m_9999
   BoundaryFace(i)%area = zero
   do j=1,SPACEDIM
      do n=1,SPACEDIM
         BoundaryFace(i)%T(j,n) = zero
         BoundaryFace(i)%RPsi(j,n) = zero
         BoundaryFace(i)%RFi(j,n) = zero
      enddo
      BoundaryFace(i)%velocity(j) = zero
   enddo
enddo
if (allocated(BFaceList)) then
   do i=1,size(BFaceList)
      BFaceList(i) = 0
   enddo
endif
do i=1,size(BoundaryVertex)
   BoundaryVertex(i) = 0
enddo
do i=1,size(Tratto)
   Tratto(i)%tipo = "    "
   Tratto(i)%ColorCode = 0
   Tratto(i)%numvertices = 0
   Tratto(i)%inivertex = 0
   Tratto(i)%iniside = 0
   Tratto(i)%iniface = 0
   Tratto(i)%medium = 0
   Tratto(i)%zone = 0
   Tratto(i)%NormVelocity = zero
   Tratto(i)%trampa = zero
   Tratto(i)%ShearCoeff = zero
   do j=1,SPACEDIM
      Tratto(i)%velocity(j) = zero
      Tratto(i)%PsiCoeff(j) = zero
      Tratto(i)%FiCoeff(j)  = zero
   enddo
enddo
do i=1,size(BoundarySide)
   BoundarySide(i)%tipo = "    "
   BoundarySide(i)%stretch = 0
   BoundarySide(i)%previous_side = 0
   BoundarySide(i)%vertex(1:SPACEDIM-1) = 0
   BoundarySide(i)%CloseParticles = 0
   BoundarySide(i)%length = zero
   BoundarySide(i)%CloseParticles_maxQuota = const_m_9999
   do n=1,SPACEDIM
      BoundarySide(i)%T(n,1:SPACEDIM) = zero
      BoundarySide(i)%R(n,1:SPACEDIM) = zero
      BoundarySide(i)%RN(n,1:SPACEDIM) = zero
   enddo
   BoundarySide(i)%angle = zero
   do j=1,SPACEDIM
      BoundarySide(i)%velocity(j) = zero
   enddo
enddo
! In case of restart, it does not zero "domain" and "grid"
if (Restart) return
do j=1,3
   Grid%ncd(j) = 0
   Grid%dcd(j) = zero
   Grid%extr(j,1) = zero
   Grid%extr(j,2) = zero
enddo
Grid%nmax = 0
Domain%tipo = "semi"
Domain%file = "                                                                                "
Domain%Psurf = " "
Domain%RandomPos = " "
Domain%iplot_fr = 0
Domain%imemo_fr = 0
Domain%irest_fr = 0
Domain%icpoi_fr = 0
Domain%ipllb_fr = 0
Domain%ipllb_md = 0
Domain%istart = 0
Domain%ioutopt = 0
Domain%itmax = 0
do j=1,3
   Domain%coord(j,1) = zero
   Domain%coord(j,2) = zero
   Domain%grav(j) = zero
enddo
Domain%tmax = zero
Domain%dx = zero
Domain%trunc = zero
Domain%coefke = zero
Domain%coefkacl = zero
Domain%CFL = zero
Domain%vsc_coeff = zero
Domain%prif = zero
Domain%plot_fr = zero
Domain%memo_fr = zero
Domain%rest_fr = zero
Domain%cpoi_fr = zero
Domain%pllb_fr = zero
Domain%TetaP = zero
Domain%TetaV = zero
Domain%prif = zero
Domain%h = zero
Domain%start = zero
Domain%NormFix = .false.
Domain%Slip = .false.
!------------------------
! Deallocations
!------------------------
return
end subroutine Init_Arrays

