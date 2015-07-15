!cfile FindLine.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
  subroutine FindLine (Xmin, Xmax, Nt)
!
!..Finds extremes of the rectangular frame which contains the boundary mib
!
!.. assign modules ..
  use FILES_ENTITIES
  use GLOBAL_MODULE
  use AdM_User_Type
  use ALLOC_MODULE
  use DIAGNOSTIC_MODULE
!
!.. implicit declarations ..
  implicit none
!
!.. dummy arguments ..
  integer(4),      intent(IN)                                   :: Nt
  double precision,intent(INOUT), dimension(SPACEDIM,NumTratti) :: Xmin, Xmax
!
!.. local Scalars ..
  integer(4) :: i, iv, nf
  character(len=lencard)  :: nomsub = "FindLine"
!
!.. executable statements
!
!.. loops on the segments of the line
!
  do iv = Tratto(Nt)%inivertex, Tratto(Nt)%inivertex + Tratto(Nt)%numvertices - 1
!
!.. checks for the maximum boundary lines storage
!
    if (iv > NumBVertices) then
      call diagnostic (arg1=10,arg2=2,arg3=nomsub)
    end if
!
!.. evaluates the minimum and maximum coordinates with respect the current vertex nf in the line nt
!
    nf = BoundaryVertex(iv)
!
    do i = 1, SPACEDIM
      if ( Vertice(i,nf) < Xmin(i,Nt) ) Xmin(i,Nt) = Vertice(i,nf)
      if ( Vertice(i,nf) > Xmax(i,Nt) ) Xmax(i,Nt) = Vertice(i,nf)
    end do
  end do
!
  return
  end subroutine FindLine
!---split

