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

