!cfile Init_arrays.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : Init_arrays
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
! Module purpose : Initialize the storage arrays for calculation
!
! Calling routine: Gest_input
!
! Called routines: none
!
!************************************************************************************

subroutine Init_Arrays
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
integer(4) :: i, j, n, i1
!
!.. Executable Statements ..
!
!Initialize Arrays

 do i = 1,size(Med)
   Med(i)%tipo         = "Empty   "
   Med(i)%modelloerosione = "        "
   Med(i)%index        = 0
   Med(i)%NIterSol     = 0
   Med(i)%den0         = zero
   Med(i)%eps          = zero
   Med(i)%celerita     = zero
   Med(i)%alfaMon      = zero
   Med(i)%betaMon      = zero
   Med(i)%visc         = zero
   Med(i)%coes         = zero
   Med(i)%numx         = zero
   Med(i)%mumx         = zero
   Med(i)%taucri       = zero
   Med(i)%cuin         = zero
   Med(i)%phi          = zero
   Med(i)%cons         = zero
   Med(i)%Cs           = zero
   Med(i)%RoughCoef    = zero
   Med(i)%D50          = zero
   Med(i)%SettlingCoef = zero
   Med(i)%Codif        = zero
   Med(i)%Gamma        = zero
   Med(i)%InitialIntEn = zero
 end do


 do i = 1,size(Partz)
    Partz(i)%label    = "        "
    Partz(i)%tipo     = "    "
    Partz(i)%shape    = " "
    Partz(i)%bend     = " "
    Partz(i)%pressure = "  "
    Partz(i)%move     = "   "
    Partz(i)%slip     = " "

    Partz(i)%ipool    = 0
    Partz(i)%npoints  = 0
    Partz(i)%icol     = 0
    Partz(i)%Medium   = 0
    Partz(i)%npointv  = 0
    Partz(i)%Indix(1) = 0
    Partz(i)%Indix(2) = 0
    Partz(i)%limit(1) = 1
    Partz(i)%limit(2) = 0
    Partz(i)%pool           = zero
    do j = 1,3
      Partz(i)%coordMM(j,1) = zero
      Partz(i)%coordMM(j,2) = zero
      Partz(i)%vel(j)       = zero
    end do
    do j = 0,3
      do n = 1,MAXPOINTSVLAW
        Partz(i)%vlaw(j,n) = zero
      end do
    end do
    Partz(i)%trampa = zero
    Partz(i)%valp = zero
 end do


 do i = 1,size(Control_Points)
   Control_Points(i)%cella    = 0
   Control_Points(i)%coord(1) = zero
   Control_Points(i)%coord(2) = zero
   Control_Points(i)%coord(3) = zero
   Control_Points(i)%vel  (1) = zero
   Control_Points(i)%vel  (2) = zero
   Control_Points(i)%vel  (3) = zero
   Control_Points(i)%pres     = zero
   Control_Points(i)%dens     = zero
   Control_Points(i)%uni      = zero
   Control_Points(i)%dist     = zero
 end do


 do i = 1,size(Section_Points)
   Section_Points(i)%cella    = 0
   Section_Points(i)%coord(1) = zero
   Section_Points(i)%coord(2) = zero
   Section_Points(i)%coord(3) = zero
   Section_Points(i)%vel  (1) = zero
   Section_Points(i)%vel  (2) = zero
   Section_Points(i)%vel  (3) = zero
   Section_Points(i)%pres     = zero
   Section_Points(i)%dens     = zero
   Section_Points(i)%uni      = zero
   Section_Points(i)%dist     = zero
 end do


 Control_Sections(0            )%Label = "Domain  "
 Control_Sections(1:NSections  )%Label = "        "
 Control_Sections(  NSections+1)%Label = "One more"
 do i1 = 1,size(Control_Sections)
   i = i1 - 1
   Control_Sections(i)%Tipo          = "  "
   Control_Sections(i)%Icont(1)      = 0
   Control_Sections(i)%Icont(2)      = 0
   Control_Sections(i)%ColorCode     = 1
   Control_Sections(i)%Constant(1)   = zero
   Control_Sections(i)%Constant(2)   = zero
   Control_Sections(i)%Constant(3)   = zero
   Control_Sections(i)%XYZRange(1,1) = zero
   Control_Sections(i)%XYZRange(2,1) = zero
   Control_Sections(i)%XYZRange(3,1) = zero
   Control_Sections(i)%XYZRange(1,2) = zero
   Control_Sections(i)%XYZRange(2,2) = zero
   Control_Sections(i)%XYZRange(3,2) = zero
   do j = 1,SPACEDIM
     do n = 1,SPACEDIM
       Control_Sections(i)%TGLsection(j,n) = zero
     end do
   end do
   Control_Sections(i)%TGLsection(1,2) = one
   Control_Sections(i)%TGLsection(2,1) =-one
   Control_Sections(i)%TGLsection(3,3) = one
 end do


 do i = 1,size(Control_Lines)
   Control_Lines (i)%label    = "Empty   "
   Control_Lines (i)%icont(1) = 0
   Control_Lines (i)%icont(2) = 0
 end do


 do j = 1,SPACEDIM
   Vertice(j,1:NumVertici) = zero
 end do

 do i = 1,size(BoundaryFace)
    do j = 1,MAXFACENODES
       BoundaryFace(i)%Node(j)%name = 0
       do n = 1,SPACEDIM
          BoundaryFace(i)%Node(j)%GX(n) = zero
          BoundaryFace(i)%Node(j)%LX(n) = zero
       end do
    end do
    BoundaryFace(i)%nodes   = 0
    BoundaryFace(i)%stretch  = 0
    BoundaryFace(i)%CloseParticles = 0
    BoundaryFace(i)%CloseParticles_maxQuota = const_m_9999
    BoundaryFace(i)%area    = zero
    do j = 1,SPACEDIM
      do n = 1,SPACEDIM
        BoundaryFace(i)%T(j,n)    = zero
        BoundaryFace(i)%RPsi(j,n) = zero
        BoundaryFace(i)%RFi(j,n)  = zero
      end do
      BoundaryFace(i)%velocity(j) = zero
    end do
 end do


 if (allocated(BFaceList)) then
   do i = 1,size(BFaceList)
     BFaceList(i) = 0
   end do
 end if


 do i = 1,size(BoundaryVertex)
   BoundaryVertex(i) = 0
 end do


 do i = 1,size(Tratto)
   Tratto(i)%tipo         = "    "
   Tratto(i)%ColorCode    = 0
   Tratto(i)%numvertices  = 0
   Tratto(i)%inivertex    = 0
   Tratto(i)%iniside      = 0
   Tratto(i)%iniface      = 0
   Tratto(i)%medium       = 0
   Tratto(i)%zone         = 0
   Tratto(i)%NormVelocity = zero
   Tratto(i)%trampa       = zero
   Tratto(i)%ShearCoeff   = zero
   do j = 1,SPACEDIM
     Tratto(i)%velocity(j) = zero
     Tratto(i)%PsiCoeff(j) = zero
     Tratto(i)%FiCoeff(j)  = zero
   end do
 end do



 do i = 1,size(BoundarySide)
   BoundarySide(i)%tipo           = "    "
   BoundarySide(i)%stretch        = 0
   BoundarySide(i)%previous_side  = 0
   BoundarySide(i)%vertex(1)      = 0
   BoundarySide(i)%vertex(2)      = 0
   BoundarySide(i)%CloseParticles = 0
   BoundarySide(i)%length         = zero
   BoundarySide(i)%CloseParticles_maxQuota = const_m_9999
   do n = 1,SPACEDIM
     BoundarySide(i)%T(n,1:SPACEDIM)  = zero
     BoundarySide(i)%R(n,1:SPACEDIM)  = zero
     BoundarySide(i)%RN(n,1:SPACEDIM) = zero
   end do
   BoundarySide(i)%angle = zero
   do j = 1,SPACEDIM
     BoundarySide(i)%velocity(J) = zero
   end do
 end do

! caso di restart non azzero domian e grid
if (Restart) return


 do j = 1,3
   Grid%ncd(j)    = 0
   Grid%dcd(j)    = zero
   Grid%extr(j,1) = zero
   Grid%extr(j,2) = zero
 end do
 Grid%nmax = 0


 Domain%tipo      = "semi"
 Domain%file      = "                                                                                "
 Domain%Psurf     = " "
 Domain%RandomPos = " "
 Domain%iplot_fr  = 0
 Domain%imemo_fr  = 0
 Domain%irest_fr  = 0
 Domain%icpoi_fr  = 0
 Domain%ipllb_fr  = 0
 Domain%ipllb_md  = 0
 Domain%istart    = 0
 Domain%ioutopt   = 0
 Domain%itmax     = 0
 do j = 1,3
   Domain%coord(j,1) = zero
   Domain%coord(j,2) = zero
   Domain%grav(j)    = zero
 end do
 Domain%tmax     = zero
 Domain%dd       = zero
 Domain%trunc    = zero
 Domain%coefke   = zero
 Domain%coefkacl = zero
! Domain%cote     = zero
 Domain%CFL      = zero
 Domain%prif     = zero
 Domain%plot_fr  = zero
 Domain%memo_fr  = zero
 Domain%rest_fr  = zero
 Domain%cpoi_fr  = zero
 Domain%pllb_fr  = zero
! Domain%TetaX    = zero
 Domain%TetaP    = zero
 Domain%TetaV    = zero
 Domain%pre      = zero
 Domain%h        = zero
 Domain%start    = zero
 Domain%NormFix  = .false.
 Domain%Slip     = .false.

return
end subroutine Init_Arrays
!---split

