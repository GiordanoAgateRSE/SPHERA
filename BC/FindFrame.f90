!cfile FindFrame.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
  subroutine FindFrame (Xmin, Xmax, Nt)
!
!Finds extremes of the rectangular frame which contains the boundary mib
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_User_Type
  use ALLOC_MODULE
!
!.. implicit declarations
  implicit none
!
!.. dummy arguments
  integer(4),intent(IN) :: Nt
  double precision,intent(INOUT),dimension(SPACEDIM,NumFacce) :: Xmin, Xmax
!
!.. local scalars
  integer(4) :: i, n, iv, nf, nod
!
!.. executable statements
!
  do iv = Tratto(Nt)%iniface, Tratto(Nt)%iniface + Tratto(Nt)%numvertices - 1
    nf = BFaceList(iv)
    do n = 1, 4
      nod = BoundaryFace(nf)%Node(n)%name
      if ( nod <= 0 ) cycle
      do i = 1, Ncord
        if ( Vertice(i,nod) < Xmin(i,Nt) ) Xmin(i,Nt) = Vertice(i,nod)
        if ( Vertice(i,nod) > Xmax(i,Nt) ) Xmax(i,Nt) = Vertice(i,nod)
      end do
    end do
  end do
!
  return
  end subroutine FindFrame
!---split

