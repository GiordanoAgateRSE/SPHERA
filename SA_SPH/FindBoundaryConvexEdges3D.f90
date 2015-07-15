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

