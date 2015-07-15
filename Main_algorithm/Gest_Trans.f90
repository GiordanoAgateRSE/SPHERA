!cfile gest_trans.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : gest_trans
!
! Last updating : May 08, 2012
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Amicarelli/Agate  30Nov11        BSPH approach
! 04  Amicarelli/Agate  13nov12        (AA501) bBody dynamics 
! 05  Amicarelli/Agate  18apr13        add maximum number of neighbours read in input
!AA504
! 06  Amicarelli        08Apr14        (v5.04) Modifications for maximum depth and specific flow rate
!
!************************************************************************************
! Module purpose : Module for loop management
!
! Calling routine: sphera
!
! Called routines: diagnostic
!                  ComputeBoundaryIntegralTab
!                  ComputeKernelTable
!                  EvaluateBER_TimeStep
!                  GridCellBoundaryFacesIntersections3D
!                  loop_irre_2D
!                  loop_irre_3D
!                  s_ctime
!                  start_and_stop
!
!************************************************************************************
!
subroutine Gest_Trans 
!
!.. assign modules
!AA503
use AdM_USER_TYPE
use FILES_ENTITIES
use GLOBAL_MODULE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
integer(4),parameter :: ner0 = 0
!
!.. Local Scalars ..
integer(4)             :: npi,NumCellmax,i,ier,k,kk,k1,k2,nlinee,nvalori
character(len=lencard) :: nomsub = "GEST_TRANS"
character(len=lencard) :: filename,stringa,prefix,filevtk
!AA504 sub
character(len=200)     :: cargo
!
!.. Local Arrays ..
double precision,dimension(:),allocatable :: verticecolore
!
!.. External functions and subroutines ..
!
character(80), external :: lcase
logical,       external :: check_files2
integer(4),    external :: stepdata
!
!.. Executable Statements ..
!
 write ( nout,"(/)" )
 write (nout,*) "Initial number of particles      NAG= ",nag
 write ( nout,"(/)" )

!! scrittura su file di ouput delle particelle di campo
 if ( Domain%ioutopt < 0 ) then
   write (nout,*) 
   write (nout,*) "======== PARTICLES COORDINATES =========="
   do npi = 1,nag
     write (nout,"(i10,4f14.5)") npi, pg(npi)%coord, pg(npi)%tstop  !, pg(npi)%vel
   end do

 end if

!
!.. initialization arrays for vtkconverter
 nblocchi = 0
 blocchi = 0
 block = -1
 Time_Block = zero

!.. allocation arrays to store close boundaries and integrals for the loop
!!!! MaxNcbs = int(COEFNMAXPARTJ * MAXCLOSEBOUNDSIDES * Nag)
!!!! MaxNcbf = int(COEFNMAXPARTJ * MAXCLOSEBOUNDFACES * Nag)
!! MaxNcbs = int(COEFNMAXPARTJ * MAXCLOSEBOUNDSIDES * PARTICLEBUFFER)
!! MaxNcbf = int(COEFNMAXPARTJ * MAXCLOSEBOUNDFACES * PARTICLEBUFFER)
 MaxNcbs = int(MAXCLOSEBOUNDSIDES * PARTICLEBUFFER)
!AA504sub 
 MaxNcbf = int(Domain%MAXCLOSEBOUNDFACES * PARTICLEBUFFER)

!AA406
 if ((Domain%tipo == "semi") .or. (Domain%tipo == "bsph"))  then
!
 allocate (BoundaryDataPointer(1:3,1:PARTICLEBUFFER), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays BoundaryDataPointer not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays BoundaryDataPointer successfully allocated "
 end if
!
 if (Ncord == 2) then
   write(nout,'(a,i15)') "     Max num of close boundary sides: MaxNcbs = ",MaxNcbs
   allocate (BoundaryDataTab(1:MaxNcbs), stat = ier)
   if (ier /= 0) then
     write (nout,'(1x,a,i2)') "   Arrays BoundaryDataTab not allocated. Error code: ",ier
     call diagnostic (arg1=4,arg3=nomsub)
   else
     write (nout,'(1x,a)') "   Arrays BoundaryDataTab successfully allocated "
   end if
 else
   write(nout,'(a,i15)') "     Max num of close boundary faces: MaxNcbf = ",MaxNcbf
   allocate (BoundaryDataTab(1:MaxNcbf), stat = ier)
   if (ier /= 0) then
     write (nout,'(1x,a,i2)') "   Arrays BoundaryDataTab not allocated. Error code: ",ier
     call diagnostic (arg1=4,arg3=nomsub)
   else
     write (nout,'(1x,a)') "   Arrays BoundaryDataTab successfully allocated "
   end if
 end if
!
!AA406
 endif
!
!.. allocation arrays to store values for the loop
 NMAXPARTJ = Domain%COEFNMAXPARTJ * (Domain%h * four / Domain%dd) ** Ncord
 write(nout,'(a,i15)') "     Max num particles surrounding the current particle: NMAXPARTJ = ",NMAXPARTJ
!
 allocate (Array_Flu(1:PARTICLEBUFFER), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays Array_Flu not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays Array_Flu successfully allocated "
 end if
!
! allocate (Array_Sol(1:PARTICLEBUFFER), stat = ier)
! if (ier /= 0) then
!   write (nout,'(1x,a,i2)') "   Arrays Array_Sol not allocated. Error code: ",ier
!   call diagnostic (arg1=4,arg3=nomsub)
! else
!   write (nout,'(1x,a)') "   Arrays Array_Sol successfully allocated "
! end if
!
 allocate (nPartIntorno(1:PARTICLEBUFFER), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays NPARTINTORNO not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays NPARTINTORNO successfully allocated "
 end if
!
 allocate (PartIntorno(1:NMAXPARTJ*PARTICLEBUFFER), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays PARTINTORNO not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays PARTINTORNO successfully allocated "
 end if
 allocate (PartKernel(1:4,1:NMAXPARTJ*PARTICLEBUFFER), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays PARTKERNEL not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays PARTKERNEL successfully allocated "
 end if
 allocate (rag(1:3,1:NMAXPARTJ*PARTICLEBUFFER), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays RAG not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays RAG successfully allocated "
 end if
!
!AA406 start
 if (Domain%tipo == "bsph") then
 allocate (nPartIntorno_fw(1:PARTICLEBUFFER), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays NPARTINTORNO_fw not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays NPARTINTORNO_fw successfully allocated "
 end if
!
 allocate (PartIntorno_fw(1:NMAXPARTJ*PARTICLEBUFFER), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays PARTINTORNO_fw not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays PARTINTORNO_fw successfully allocated "
 end if
!
!AA406test
!AA501 sub
 if (Domain%tipo == "bsph") allocate (kernel_fw(2,1:NMAXPARTJ*PARTICLEBUFFER), stat = ier)
!
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays kernel_fw not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays kernel_fw successfully allocated "
 end if
 allocate (rag_fw(1:3,1:NMAXPARTJ*PARTICLEBUFFER), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays RAG_fw not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays RAG_fw successfully allocated "
 end if
 endif
!AA406 end
!

!AA501b start
 if (n_bodies > 0) then
 allocate (nPartIntorno_bp_f(n_body_part), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays nPartIntorno_bp_f not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays nPartIntorno_bp_f successfully allocated "
 end if
 allocate (PartIntorno_bp_f(NMAXPARTJ*n_body_part), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays PartIntorno_bp_f not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays PartIntorno_bp_f successfully allocated "
 end if
 allocate (KerDer_bp_f_cub_spl(NMAXPARTJ*n_body_part), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays KerDer_bp_f_cub_spl not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays KerDer_bp_f_cub_spl successfully allocated "
 end if
 allocate (KerDer_bp_f_Gal(NMAXPARTJ*n_body_part), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays KerDer_bp_f_Gal not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays KerDer_bp_f_Gal successfully allocated "
 end if 
 allocate (rag_bp_f(3,NMAXPARTJ*n_body_part), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays rag_bp_f not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays rag_bp_f successfully allocated "
 end if
 allocate (nPartIntorno_bp_bp(n_surf_body_part), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays nPartIntorno_bp_bp not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays nPartIntorno_bp_bp successfully allocated "
 end if
 allocate (PartIntorno_bp_bp(n_surf_body_part*NMAXPARTJ), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays PartIntorno_bp_bp not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays PartIntorno_bp_bp successfully allocated "
 end if
 allocate (rag_bp_bp(3,n_surf_body_part*NMAXPARTJ), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays rag_bp_bp not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays rag_bp_bp successfully allocated "
 end if
 if (ncord==2) allocate (impact_vel(n_surf_body_part,(n_bodies+NumBSides)), stat = ier)
 if (ncord==3) allocate (impact_vel(n_surf_body_part,(n_bodies+NumFacce)), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Array impact_vel not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Array impact_vel successfully allocated "
 end if
 impact_vel = 0.
 endif
!AA501b end

 write (nout,'(1x,a)') "..."
 write (nout,'(a,i15)') " Max number of particles  : PARTICLEBUFFER = ",PARTICLEBUFFER
 write (nout,*) " Size # of elements in array pg                  : ",size(pg)
 write (nout,*) " Size # of elements in array BoundaryDataTab     : ",size(BoundaryDataTab)
 write (nout,*) " Size # of elements in array BoundaryDataPointer : ",size(BoundaryDataPointer)
 write (nout,*) " Size # of elements in array Array_Flu           : ",size(Array_Flu)
! write (nout,*) " Size # of elements in array Array_Sol           : ",size(Array_Sol)
 write (nout,*) " Size # of elements in array nPartIntorno        : ",size(nPartIntorno)
 write (nout,*) " Size # of elements in array PartIntorno         : ",size(PartIntorno)
 write (nout,*) " Size # of elements in array PartKernel          : ",size(PartKernel)
 write (nout,*) " Size # of elements in array rag                 : ",size(rag)
!
!AA406 start
 if ((Domain%tipo == "bsph").and.(DBSPH%n_w > 0)) then
 write (nout,*) " Size # of elements in array pg_w                : ",size(pg_w)
 write (nout,*) " Size # of elements in array nPartIntorno_fw     : ",size(nPartIntorno_fw)
 write (nout,*) " Size # of elements in array PartIntorno_fw      : ",size(PartIntorno_fw)
 write (nout,*) " Size # of elements in array kernel_fw           : ",size(kernel_fw)
 write (nout,*) " Size # of elements in array rag_fw              : ",size(rag_fw)
 endif
!AA406 end

!AA501b start
 if (n_bodies > 0) then
 write (nout,*) " Size # of elements in array nPartIntorno_bp_f   : ",size(nPartIntorno_bp_f)
 write (nout,*) " Size # of elements in array PartIntorno_bp_f    : ",size(PartIntorno_bp_f)
 write (nout,*) " Size # of elements in array KerDer_bp_f_cub_spl : ",size(KerDer_bp_f_cub_spl)
 write (nout,*) " Size # of elements in array KerDer_bp_f_Gal     : ",size(KerDer_bp_f_Gal) 
 write (nout,*) " Size # of elements in array rag_bp_f            : ",size(rag_bp_f)
 write (nout,*) " Size # of elements in array surf_body_part      : ",size(surf_body_part) 
 write (nout,*) " Size # of elements in array nPartIntorno_bp_bp  : ",size(nPartIntorno_bp_bp)
 write (nout,*) " Size # of elements in array PartIntorno_bp_bp   : ",size(PartIntorno_bp_bp)
 write (nout,*) " Size # of elements in array rag_bp_bp           : ",size(rag_bp_bp)
 endif
!AA501b end

!
 write (nout,'(1x,a)') "..."
 write (nout,*) " Size in bytes of array pg                       : ",sizeof(pg)
 write (nout,*) " Size in bytes of array BoundaryDataTab          : ",sizeof(BoundaryDataTab)
 write (nout,*) " Size in bytes of array BoundaryDataPointer      : ",sizeof(BoundaryDataPointer)
 write (nout,*) " Size in bytes of array Array_Flu                : ",sizeof(Array_Flu)
! write (nout,*) " Size in bytes of array Array_Sol                : ",sizeof(Array_Sol)
 write (nout,*) " Size in bytes of array nPartIntorno             : ",sizeof(nPartIntorno)
 write (nout,*) " Size in bytes of array PartIntorno              : ",sizeof(PartIntorno)
 write (nout,*) " Size in bytes of array PartKernel               : ",sizeof(PartKernel)
 write (nout,*) " Size in bytes of array rag                      : ",sizeof(rag)
!
!AA406 start
 if ((Domain%tipo == "bsph").and.(DBSPH%n_w > 0)) then
 write (nout,*) " Size in bytes of array pg_w                     : ",sizeof(pg_w)
 write (nout,*) " Size in bytes of array nPartIntorno_fw          : ",sizeof(nPartIntorno_fw)
 write (nout,*) " Size in bytes of array PartIntorno_fw           : ",sizeof(PartIntorno_fw)
 write (nout,*) " Size in bytes of array kernel_fw                : ",sizeof(kernel_fw)
 write (nout,*) " Size in bytes of array rag_fw                   : ",sizeof(rag_fw)
 endif
!AA406 end

!AA501b start
 if (n_bodies > 0.) then
 write (nout,*) " Size in bytes of array nPartIntorno_bp_f        : ",sizeof(nPartIntorno_bp_f)
 write (nout,*) " Size in bytes of array PartIntorno_bp_f         : ",sizeof(PartIntorno_bp_f)
 write (nout,*) " Size in bytes of array KerDer_bp_f_cub_spl      : ",sizeof(KerDer_bp_f_cub_spl)
 write (nout,*) " Size in bytes of array KerDer_bp_f_Gal          : ",sizeof(KerDer_bp_f_Gal)
 write (nout,*) " Size in bytes of array rag_bp_f                 : ",sizeof(rag_bp_f)
 write (nout,*) " Size in bytes of array surf_body_part           : ",sizeof(surf_body_part) 
 write (nout,*) " Size in bytes of array nPartIntorno_bp_bp       : ",sizeof(nPartIntorno_bp_bp)
 write (nout,*) " Size in bytes of array PartIntorno_bp_bp        : ",sizeof(PartIntorno_bp_bp)
 write (nout,*) " Size in bytes of array rag_bp_bp                : ",sizeof(rag_bp_bp)
 endif
!AA501b end

!
 write (nout,'(1x,a)') "..."
!
 write (nout,'(1x,a)') " "
 write (nout,'(1x,a)') "   end allocation step. "
 write (nout,'(1x,a)') " "
!
!scrittura su file di output del numero cella e delle  particelle che gravano su di essa
 if ( Domain%ioutopt < 0 ) then
   write (nout,*) 
   write (nout,*) "Number of cells           NCELLS= ",grid%nmax
   write (nout,*) 
   write (nout,*) "======== CELLS AND RELATED PARTICLES =========="
   do i = 1,grid%nmax   
     if ( Icont(i+1) > Icont(i) ) then  !20051230
       write (nout,"(3(a,i5),a)") " cell", i," from",Icont(i)," to",Icont(i+1),"  particles:"
       write (nout,"(5i8)") NPartOrd(Icont(i):Icont(i+1)-1)  !20051230
     end if
   end do
 end if
!
 write (nout,*) 
 write (nout,*) 
 call s_ctime( nout )
 write (nout,*) 
 write (nscr,*) "Transient loop begins..."
 write (nout,*) "Transient loop begins..."
 write (nout,*)
!
!inizializza file post-processing
 if ( Domain%imemo_fr > 0 .OR. Domain%memo_fr > zero ) then
   open ( nres, file=nomefile(2) &
        , status="unknown"      &
        , access="sequential"   &
        , form  ="unformatted"  )
 else
   nres = -nres
 end if
!
!inizializza file restart in output
! if ( Domain%irest_fr > 0 .OR. Domain%rest_fr > zero ) then
!   open ( nsav, file=nomefile(3) &
!        , status="unknown"      &
!        , access="sequential"   &
!        , form  ="unformatted"  )
! else
!   nsav = -nsav
! end if
!
!inizializza file pelo libero in output
 if ( Domain%ipllb_fr > 0 .OR. Domain%pllb_fr > zero ) then
   open ( nplb, file=nomefile(4) &
        , status="unknown"      &
        , access="sequential"   &
        , form  ="formatted"  )
   write (nplb,"(a)") "tempo         free_surface_quota"
 else
   nplb = -nplb
 end if
!
!inizializza file fronte in output
 if ( Domain%imemo_fr > 0 .OR. Domain%memo_fr > zero ) then
   open ( nfro, file=nomefile(5) &
        , status="unknown"      &
        , access="sequential"   &
        , form  ="formatted"  )
   write (nfro,"(a)") "tempo         x fronte      (y fronte)    z fronte"
 else
   nfro = -nfro
 end if
!
! write (nres) menouno,const_m_9999,nag,ncord
! write (nres) domain, grid, &
!              NPartZone, NMedium, NPointst, NLines, NSections, &
!              NumVertici, NumFacce, NumTratti, NumBVertices, NumBSides, &
!              Vertice(1:SPACEDIM,1:NumVertici), &
!              BoundaryFace(1:NumFacce), &
!              Tratto (1:NumTratti),  &
!              BoundaryVertex(1:NumBVertices), & 
!              BoundarySide  (1:NumBSides), &
!              Partz(1:NPartZone), &
!              Med(1:nmedium), &
!              Control_Sections(0:NSections+1), &
!              Control_Points(1:NPointst), &
!              Control_Lines(1:NLines)

!
!.. create file for vtk to design the boundaries of the domain
!
 if (vtkconv) then
!
   prefix = nomecaso
   filevtk = "VTKConverter_"//prefix(1:len_trim(prefix))//"_domain.vtk"
   open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',status='unknown')
   write(unitvtk,'(a)') '# vtk DataFile Version 2.0'
   write(unitvtk,'(a,a)') 'Domain limits for the case:',prefix(1:len_trim(prefix))
   write(unitvtk,'(a)') 'ASCII'
   write(unitvtk,'(a)') 'DATASET POLYDATA'
   write(unitvtk,'(a,i8,a)') 'POINTS ',numvertici,' float'

   do i = 1,numvertici,4
     k1 = i
     k2 = k1 + 3
     if (k2 > numvertici) k2 = numvertici
     write (stringa,'(4(3(e12.5,1x)))') (vertice(1,k),vertice(2,k),vertice(3,k),k=k1,k2)
     stringa = adjustl(trim(stringa))
     write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
   end do
!
!.. 2D case
!
   if (ncord == 2) then
!
     nlinee = 0
     nvalori = 0
     do i = 1,numbsides
       if (boundaryside(i)%TIPO == 'peri') cycle
       nlinee = nlinee + 1
       nvalori = nvalori + 1
!!??       nvalori = nvalori + (boundaryside(i)%VERTEX(2) - boundaryside(i)%VERTEX(1) + 1)
       nvalori = nvalori + 2
     end do
     write(unitvtk,'(a,2i8)') 'LINES ', nlinee, nvalori
!
     allocate (verticecolore(numvertici))
     verticecolore = zero
     do i = 1,numbsides
       if (boundaryside(i)%TIPO == 'peri') cycle
       stringa = ' '
       k1 = boundaryside(i)%VERTEX(1) - 1
       k2 = boundaryside(i)%VERTEX(2) - 1
       write (stringa,'(a,2i10)') '2',k1,k2
       stringa = adjustl(trim(stringa))
       write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
       if (boundaryside(i)%TIPO == 'velo' .or. boundaryside(i)%TIPO == 'flow' .or. boundaryside(i)%TIPO == 'open') then
         verticecolore(k1+1) = 2.0
         verticecolore(k2+1) = 2.0
       else if (boundaryside(i)%TIPO == 'sour') then
         verticecolore(k1+1) = 1.0
         verticecolore(k2+1) = 1.0
       end if
     end do
!
     write(unitvtk,'(a,i8)') 'POINT_DATA ', numvertici
     write(unitvtk,'(a,a,a)') 'SCALARS ', prefix(1:len_trim(prefix)),' float 1'
     write(unitvtk,'(a)') 'LOOKUP_TABLE mytable'
!
     do i = 1,numvertici,8
       k1 = i
       k2 = k1 + 7
       if (k2 > numvertici) k2 = numvertici
       write (stringa,'(8(f10.3,1x))') (verticecolore(k),k=k1,k2)
       stringa = adjustl(trim(stringa))
       write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
     end do
     deallocate (verticecolore)
!
     write(unitvtk,'(a)') 'LOOKUP_TABLE mytable 3'
     write(unitvtk,'(a)') '0.0 0.0 0.0 1.0'
     write(unitvtk,'(a)') '0.0 0.0 1.0 1.0'
     write(unitvtk,'(a)') '1.0 0.0 0.0 1.0'
!
!.. 3D case
!
   else
!
     nlinee = 0
     nvalori = 0
     do i = 1,numfacce
       if (tratto(BoundaryFace(i)%stretch)%tipo == 'peri') cycle
       nlinee = nlinee + 1
       nvalori = nvalori + 1 + BoundaryFace(i)%nodes + 1
     end do
     write(unitvtk,'(a,2i8)') 'LINES ', nlinee, nvalori
!
     allocate (verticecolore(numvertici))
     verticecolore = 0.0
     do i = 1,numfacce
       if (tratto(BoundaryFace(i)%stretch)%tipo == 'peri') cycle
       stringa = ' '
       do k = 1,BoundaryFace(i)%nodes,8
         k1 = k 
         k2 = k1 + 7
         if (k2 > BoundaryFace(i)%nodes) k2 = BoundaryFace(i)%nodes
         write (stringa,'(10i10)') BoundaryFace(i)%nodes+1,(BoundaryFace(i)%node(kk)%name-1,kk=k1,k2),BoundaryFace(i)%node(1)%name-1
         stringa = adjustl(trim(stringa))
         write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
       end do
       if (tratto(BoundaryFace(i)%stretch)%tipo == 'velo' .or. tratto(BoundaryFace(i)%stretch)%tipo == 'flow' .or. &
           tratto(BoundaryFace(i)%stretch)%tipo == 'open') then
         do kk = k1,k2
           verticecolore(BoundaryFace(i)%node(kk)%name) = 2.0
           verticecolore(BoundaryFace(i)%node(kk)%name) = 2.0
         end do
       else if (tratto(BoundaryFace(i)%stretch)%tipo == 'sour') then
         do kk = k1,k2
           verticecolore(BoundaryFace(i)%node(kk)%name) = 1.0
           verticecolore(BoundaryFace(i)%node(kk)%name) = 1.0
         end do
       end if
     end do
!
     write(unitvtk,'(a,i8)') 'POINT_DATA ', numvertici
     write(unitvtk,'(a,a,a)') 'SCALARS ', prefix(1:len_trim(prefix)),' float 1'
     write(unitvtk,'(a)') 'LOOKUP_TABLE mytable'
!
     do i = 1,numvertici,8

       k1 = i
       k2 = k1 + 7
       if (k2 > numvertici) k2 = numvertici
       write (stringa,'(8(f10.3,1x))') (verticecolore(k),k=k1,k2)
       stringa = adjustl(trim(stringa))
       write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
     end do
     deallocate (verticecolore)
!
     write(unitvtk,'(a)') 'LOOKUP_TABLE mytable 3'
     write(unitvtk,'(a)') '0.0 0.0 0.0 1.0'
     write(unitvtk,'(a)') '0.0 0.0 1.0 1.0'
     write(unitvtk,'(a)') '1.0 0.0 0.0 1.0'
!
   end if
!
   flush(unitvtk)
   close (unitvtk)
 end if
!
!..
!
!AA406 sub
   if ( (Domain%tipo == "semi") .or. (Domain%tipo == "bsph") ) then
!
!.. calcolo DTmin reazione elastica boundary
   call EvaluateBER_TimeStep
!
   if (ncord==2) then
!       
     call start_and_stop(2,5)
     call loop_irre_2D 
     call start_and_stop(3,5)
!
   else if (ncord==3) then
!
     NumCellmax = Grid%nmax
     allocate ( GCBFPointers(NumCellmax,2), stat = ier )
     if (ier /= 0) then
       write (nout,'(1x,a,i2)') "   Array GCBFPointers not allocated. Error code: ",ier
       call diagnostic (arg1=4,arg3=nomsub)
     else
       write (nout,'(1x,a)') "   Array GCBFPointers successfully allocated "
     end if
!
!AA406
     if (Domain%tipo == "semi") then
!AA504 rm spart         
!.. To select the grid cells intercepting a boundary faces (Boundaries.f90)
     call GridCellBoundaryFacesIntersections3D (NumCellmax)
!
!.. To compute the local coordinates, solid angle and solid normal, relative to a boundary element (semi-analytic approach; Boundaries.f90)
     call ComputeBoundaryIntegralTab
!
!.. Computation of the boundary contributions for the continuity equation (Boundaries.f90)
     call ComputeKernelTable
!
!AA406 
     endif
     
!AA504 start
!Allocation and initialization of the Z_fluid_max array
!loop over the zones
     do i=1,NPartZone
        if (Partz(i)%IC_source_type == 2) then
            allocate(Z_fluid_max(Grid%ncd(1)*Grid%ncd(2)))
            Z_fluid_max = -999.
            allocate(q_max(Grid%ncd(1)*Grid%ncd(2)))
            q_max = 0.d0
            exit
        endif
     end do     
!AA504 end
     call start_and_stop(2,5)
!.. main loop
     call loop_irre_3D
     call start_and_stop(3,5)

!AA504 start 
!Writing the h_max array and deallocation of the Z_fluid_max array
     if (allocated(Z_fluid_max)) then
        call write_h_max        
        if (allocated(Z_fluid_max)) deallocate(Z_fluid_max)
     endif
!AA504 end
     
   end if 
!
 else

   call diagnostic (arg1=10,arg2=5,arg3=nomsub)
 end if
!
!!!!! deallocate ( NPartOrd,Icont, stat = ier )
!!!!! if (ier /= 0) then
!!!!!   write (nout,'(1x,a,i2)') "   Arrays NPartOrd,Icont  not deallocated. Error code: ",ier
!!!!!   call diagnostic (arg1=4,arg3=nomsub)
!!!!! else
!!!!!   write (nout,'(1x,a)') "   Arrays NPartOrd,Icont successfully deallocated "
!!!!! end if
!
!.. create file for vtk
!
  if (vtkconv) then
    prefix = nomecaso
    filevtk = "VTKConverter_"//prefix(1:len_trim(prefix))//".pvd"
    open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',status='unknown')
!
    write(unitvtk,'(a)') '<?xml version="1.0"?>'
    write(unitvtk,'(a)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
    write(unitvtk,'(a)') '  <Collection>'
    do i = 1,nblocchi
      stringa = ' '
      write(cargo,'(i10)') blocchi(i)
      cargo = adjustl(trim(cargo))
      filename = "VTKConverter_"//prefix(1:len_trim(prefix))//"_block_"//cargo(1:len_trim(cargo))//".vtu"
      stringa = '  <DataSet timestep="'
!AA504 sub     
      write(cargo,'(f15.6)')  Time_Block(i)
      cargo = adjustl(trim(cargo))
      stringa = stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" group="" part="'
      write(cargo,'(i10)') i
      cargo = adjustl(trim(cargo))
      stringa = stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" file="'//filename(1:len_trim(filename))//'"/>'
      write(unitvtk,'(a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '  </Collection>'
    write(unitvtk,'(a)') '</VTKFile>'
!
    flush(unitvtk)
    close (unitvtk)
  end if
!
!AA406 start
!.. create file for vtk (wall elements)
!
 if ((Domain%tipo == "bsph").and.(DBSPH%n_w > 0)) then
  if (vtkconv) then
    filevtk = "VTKConverter_wall_"//prefix(1:len_trim(prefix))//".pvd"
    open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',status='unknown')
    write(unitvtk,'(a)') '<?xml version="1.0"?>'
    write(unitvtk,'(a)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
    write(unitvtk,'(a)') '  <Collection>'
    do i = 1,nblocchi
      stringa = ' '
      write(cargo,'(i10)') blocchi(i)
      cargo = adjustl(trim(cargo))
      filename = "VTKConverter_"//prefix(1:len_trim(prefix))//"_block_wall_"//cargo(1:len_trim(cargo))//".vtu"
      stringa = '  <DataSet timestep="'
!AA504 sub       
      write(cargo,'(f15.6)')  Time_Block(i)
      cargo = adjustl(trim(cargo))
      stringa = stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" group="" part="'
      write(cargo,'(i10)') i
      cargo = adjustl(trim(cargo))
      stringa = stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" file="'//filename(1:len_trim(filename))//'"/>'
      write(unitvtk,'(a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '  </Collection>'
    write(unitvtk,'(a)') '</VTKFile>'
    flush(unitvtk)
    close (unitvtk)
  end if
  endif
!AA406 end

!AA501b start
 if (n_bodies > 0) then
 
! creation of the .pvd file for body particles
  if (vtkconv) then
    filevtk = "VTKConverter_body-part_"//prefix(1:len_trim(prefix))//".pvd"
    open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',status='unknown')
    write(unitvtk,'(a)') '<?xml version="1.0"?>'
    write(unitvtk,'(a)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
    write(unitvtk,'(a)') '  <Collection>'
    do i = 1,nblocchi
      stringa = ' '
      write(cargo,'(i10)') blocchi(i)
      cargo = adjustl(trim(cargo))
      filename = "VTKConverter_"//prefix(1:len_trim(prefix))//"_block_body-part_"//cargo(1:len_trim(cargo))//".vtu"
      stringa = '  <DataSet timestep="'
!AA504 sub 
      write(cargo,'(f15.6)')  Time_Block(i)
      cargo = adjustl(trim(cargo))
      stringa = stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" group="" part="'
      write(cargo,'(i10)') i
      cargo = adjustl(trim(cargo))
      stringa = stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" file="'//filename(1:len_trim(filename))//'"/>'
      write(unitvtk,'(a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '  </Collection>'
    write(unitvtk,'(a)') '</VTKFile>'
    flush(unitvtk)
    close (unitvtk)
  end if

! creation of the .pvd file for bodies
    if (vtkconv) then
    filevtk = "VTKConverter_body_"//prefix(1:len_trim(prefix))//".pvd"
    open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',status='unknown')
    write(unitvtk,'(a)') '<?xml version="1.0"?>'
    write(unitvtk,'(a)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
    write(unitvtk,'(a)') '  <Collection>'
    do i = 1,nblocchi
      stringa = ' '
      write(cargo,'(i10)') blocchi(i)
      cargo = adjustl(trim(cargo))
      filename = "VTKConverter_"//prefix(1:len_trim(prefix))//"_block_body_"//cargo(1:len_trim(cargo))//".vtu"
      stringa = '  <DataSet timestep="'
!AA504 sub       
      write(cargo,'(f15.6)')  Time_Block(i)
      cargo = adjustl(trim(cargo))
      stringa = stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" group="" part="'
      write(cargo,'(i10)') i
      cargo = adjustl(trim(cargo))
      stringa = stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" file="'//filename(1:len_trim(filename))//'"/>'
      write(unitvtk,'(a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '  </Collection>'
    write(unitvtk,'(a)') '</VTKFile>'
    flush(unitvtk)
    close (unitvtk)
  end if
  
  endif
!AA501b end

!
return
end subroutine Gest_Trans
!---split

