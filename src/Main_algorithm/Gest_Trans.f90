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
! Program unit: Gest_Trans         
! Description: Introductory procedure for the main algorithm.                    
!-------------------------------------------------------------------------------
subroutine Gest_Trans 
!------------------------
! Modules
!------------------------ 
use Hybrid_allocation_module
use I_O_file_module
use Static_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),parameter :: ner0 = 0
integer(4) :: npi,NumCellmax,i,ier,k,kk,k1,k2,nlinee,nvalori,alloc_stat
integer(4) :: aux_integer
character(len=lencard) :: nomsub = "GEST_TRANS"
character(len=lencard) :: filename,stringa,prefix,filevtk
character(len=200)     :: cargo
double precision,dimension(:),allocatable :: verticecolore
logical,external :: check_files2
integer(4),external :: stepdata
character(100), external :: lcase
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
write(nout,"(/)")
write(nout,*) "Initial number of particles      NAG= ",nag
write( nout,"(/)" )
if (Domain%ioutopt<0) then
   write(nout,*)
   write(nout,*) "======== PARTICLES COORDINATES =========="
   do npi=1,nag
      write(nout,"(i10,4f14.5)") npi,pg(npi)%coord,pg(npi)%tstop
   enddo
endif
! Array initialization for vtkconverter
nblocchi = 0
blocchi = 0
block = - 1
Time_Block = zero
MaxNcbs = int(MAXCLOSEBOUNDSIDES * PARTICLEBUFFER)
MaxNcbf = int(Domain%MAXCLOSEBOUNDFACES * PARTICLEBUFFER)
if ((Domain%tipo=="semi").or.(Domain%tipo=="bsph"))  then
   allocate(BoundaryDataPointer(1:3,1:PARTICLEBUFFER),stat=ier)
   if (ier/=0) then
      write(nout,'(1x,a,i2)')                                                  &
         "   Array BoundaryDataPointer not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(nout,'(1x,a)')                                                  &
            "   Array BoundaryDataPointer successfully allocated "
   endif
   if (Ncord==2) then
      write(nout,'(a,i15)') "     Max num of close boundary sides: MaxNcbs = ",&
         MaxNcbs
      allocate(BoundaryDataTab(1:MaxNcbs),stat=ier)
      if (ier/=0) then
         write(nout,'(1x,a,i2)')                                               &
            "   Array BoundaryDataTab not allocated. Error code: ",ier
          call diagnostic(arg1=4,arg3=nomsub)
         else
            write(nout,'(1x,a)')                                               &
               "   Array BoundaryDataTab successfully allocated "
      endif
      else
         write(nout,'(a,i15)')                                                 &
            "     Max num of close boundary faces: MaxNcbf = ",MaxNcbf
         allocate(BoundaryDataTab(1:MaxNcbf),stat=ier)
         if (ier/=0) then
            write(nout,'(1x,a,i2)')                                            &
               "   Array BoundaryDataTab not allocated. Error code: ",ier
            call diagnostic(arg1=4,arg3=nomsub)
            else
               write(nout,'(1x,a)')                                            &
                  "   Array BoundaryDataTab successfully allocated"
         endif
   endif
endif
NMAXPARTJ = Domain%COEFNMAXPARTJ * (Domain%h * four / Domain%dx) ** Ncord
write(nout,'(2a,i15)') "     Maximum number of neighbouring particles: ",      &
   "NMAXPARTJ = ",NMAXPARTJ
allocate(Array_Flu(1:PARTICLEBUFFER),stat=ier)
if (ier/=0) then
   write(nout,'(1x,a,i2)') "   Array Array_Flu not allocated. Error code: ",   &
      ier
   call diagnostic(arg1=4,arg3=nomsub)
   else
      write(nout,'(1x,a)') "   Array Array_Flu successfully allocated "
endif
allocate(nPartIntorno(1:PARTICLEBUFFER),stat=ier)
if (ier/=0) then
   write(nout,'(1x,a,i2)')"   Array NPARTINTORNO not allocated. Error code: "  &
   ,ier
   call diagnostic(arg1=4,arg3=nomsub)
   else
      write(nout,'(1x,a)') "   Array NPARTINTORNO successfully allocated "
endif
allocate(PartIntorno(1:NMAXPARTJ*PARTICLEBUFFER),stat=ier)
if (ier/=0) then
   write(nout,'(1x,a,i2)') "   Array PARTINTORNO not allocated. Error code: "  &
      ,ier
   call diagnostic(arg1=4,arg3=nomsub)
   else
      write(nout,'(1x,a)') "   Array PARTINTORNO successfully allocated "
endif
allocate(PartKernel(1:4,1:NMAXPARTJ*PARTICLEBUFFER),stat=ier)
if (ier/=0) then
   write(nout,'(1x,a,i2)') "   Array PARTKERNEL not allocated. Error code: ",  &
      ier
   call diagnostic(arg1=4,arg3=nomsub)
   else
      write(nout,'(1x,a)') "   Array PARTKERNEL successfully allocated "
endif
allocate(rag(1:3,1:NMAXPARTJ*PARTICLEBUFFER),stat=ier)
if (ier/=0) then
   write(nout,'(1x,a,i2)') "   Array RAG not allocated. Error code: ",ier
   call diagnostic(arg1=4,arg3=nomsub)
   else
      write(nout,'(1x,a)') "   Array RAG successfully allocated "
endif
if (Domain%tipo=="bsph") then
   allocate(nPartIntorno_fw(1:PARTICLEBUFFER),stat=ier)
   if (ier/=0) then
      write(nout,'(1x,a,i2)')                                                  &
         "   Array NPARTINTORNO_fw not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(nout,'(1x,a)')                                                  &
            "   Array NPARTINTORNO_fw successfully allocated "
   endif
   allocate(PartIntorno_fw(1:NMAXPARTJ*PARTICLEBUFFER),stat=ier)
   if (ier/=0) then
      write(nout,'(1x,a,i2)')                                                  &
         "   Array PARTINTORNO_fw not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(nout,'(1x,a)') "   Array PARTINTORNO_fw successfully allocated"
   endif
   allocate(grad_vel_VSL_fw(1:3,1:NMAXPARTJ*PARTICLEBUFFER),stat=ier)
   if (ier/=0) then
      write(nout,'(1x,a,i2)')                                                  &
         "   Array grad_vel_VSL_fw not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(nout,'(1x,a)') "   Array grad_vel_VSL_fw successfully allocated "
   endif
   allocate(kernel_fw(2,1:NMAXPARTJ*PARTICLEBUFFER),stat=ier)
   if (ier/=0) then
      write(nout,'(1x,a,i2)')                                                  &
         "   Array kernel_fw not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(nout,'(1x,a)') "   Array kernel_fw successfully allocated "
   endif
   allocate(rag_fw(1:3,1:NMAXPARTJ*PARTICLEBUFFER),stat=ier)
   if (ier/=0) then
      write(nout,'(1x,a,i2)') "   Array RAG_fw not allocated. Error code: ",   &
         ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(nout,'(1x,a)') "   Array RAG_fw successfully allocated "
   endif
endif
if (n_bodies>0) then
   allocate(nPartIntorno_bp_f(n_body_part),stat=ier)
   if (ier/=0) then
      write(nout,'(1x,a,i2)')                                                  &
         "   Array nPartIntorno_bp_f not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(nout,'(1x,a)')                                                  &
            "   Array nPartIntorno_bp_f successfully allocated "
   endif
   allocate(PartIntorno_bp_f(NMAXPARTJ*n_body_part),stat=ier)
   if (ier/=0) then
      write(nout,'(1x,a,i2)')                                                  &
         "   Array PartIntorno_bp_f not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(nout,'(1x,a)')                                                  &
            "   Array PartIntorno_bp_f successfully allocated "
   endif
   allocate(KerDer_bp_f_cub_spl(NMAXPARTJ*n_body_part),stat=ier)
   if (ier/=0) then
      write(nout,'(1x,a,i2)')                                                  &
         "   Array KerDer_bp_f_cub_spl not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(nout,'(1x,a)')                                                  &
            "   Array KerDer_bp_f_cub_spl successfully allocated "
   endif
   allocate(KerDer_bp_f_Gal(NMAXPARTJ*n_body_part),stat=ier)
   if (ier/=0) then
      write(nout,'(1x,a,i2)')                                                  &
         "   Array KerDer_bp_f_Gal not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(nout,'(1x,a)')                                                  &
            "   Array KerDer_bp_f_Gal successfully allocated "
   endif 
   allocate(rag_bp_f(3,NMAXPARTJ*n_body_part),stat=ier)
   if (ier/=0) then
      write(nout,'(1x,a,i2)') "   Array rag_bp_f not allocated. Error code: "  &
         ,ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(nout,'(1x,a)') "   Array rag_bp_f successfully allocated "
   endif
   allocate(nPartIntorno_bp_bp(n_surf_body_part),stat=ier)
   if (ier/=0) then
      write(nout,'(1x,a,i2)')                                                  &
         "   Array nPartIntorno_bp_bp not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
      write(nout,'(1x,a)')                                                     &
         "   Array nPartIntorno_bp_bp successfully allocated "
   endif
   allocate(PartIntorno_bp_bp(n_surf_body_part*NMAXPARTJ),stat=ier)
   if (ier/=0) then
      write(nout,'(1x,a,i2)')                                                  &
         "   Array PartIntorno_bp_bp not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(nout,'(1x,a)')                                                  &
            "   Array PartIntorno_bp_bp successfully allocated "
   endif
   allocate(rag_bp_bp(3,n_surf_body_part*NMAXPARTJ),stat=ier)
   if (ier/=0) then
      write(nout,'(1x,a,i2)')                                                  &
         "   Array rag_bp_bp not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(nout,'(1x,a)') "   Array rag_bp_bp successfully allocated "
   endif
   if (ncord==2)                                                               &
      allocate(impact_vel(n_surf_body_part,(n_bodies+NumBSides)),stat=ier)
   if (ncord==3)                                                               &
      allocate(impact_vel(n_surf_body_part,(n_bodies+NumFacce)),stat=ier)
   if (ier/=0) then
      write(nout,'(1x,a,i2)')                                                  &
         "   Array impact_vel not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(nout,'(1x,a)') "   Array impact_vel successfully allocated "
   endif
   impact_vel = 0.d0
endif
write(nout,'(1x,a)') "..."
write(nout,'(a,i15)') " Max number of particles  : PARTICLEBUFFER = ",         &
   PARTICLEBUFFER
write(nout,*) " Size # of elements in array pg                  : ",size(pg)
write(nout,*) " Size # of elements in array BoundaryDataTab     : ",           &
   size(BoundaryDataTab)
write(nout,*) " Size # of elements in array BoundaryDataPointer : ",           &
   size(BoundaryDataPointer)
write(nout,*) " Size # of elements in array Array_Flu           : ",           &
   size(Array_Flu)
write(nout,*) " Size # of elements in array nPartIntorno        : ",           &
   size(nPartIntorno)
write(nout,*) " Size # of elements in array PartIntorno         : ",           &
   size(PartIntorno)
write(nout,*) " Size # of elements in array PartKernel          : ",           &
   size(PartKernel)
write(nout,*) " Size # of elements in array rag                 : ",           &
   size(rag)
if ((Domain%tipo=="bsph").and.(DBSPH%n_w>0)) then
   write(nout,*) " Size # of elements in array pg_w                : ",        &
      size(pg_w)
   write(nout,*) " Size # of elements in array nPartIntorno_fw     : ",        &
      size(nPartIntorno_fw)
   write(nout,*) " Size # of elements in array PartIntorno_fw      : ",        &
      size(PartIntorno_fw)
   write(nout,*) " Size # of elements in array kernel_fw           : ",        &
      size(kernel_fw)
   write(nout,*) " Size # of elements in array rag_fw              : ",        &
      size(rag_fw)
   write(nout,*) " Size # of elements in array grad_vel_VSL_fw     : ",        &
      size(grad_vel_VSL_fw)
endif
if (n_bodies>0) then
   write(nout,*) " Size # of elements in array nPartIntorno_bp_f   : ",        &
      size(nPartIntorno_bp_f)
   write(nout,*) " Size # of elements in array PartIntorno_bp_f    : ",        &
      size(PartIntorno_bp_f)
   write(nout,*) " Size # of elements in array KerDer_bp_f_cub_spl : ",        &
      size(KerDer_bp_f_cub_spl)
   write(nout,*) " Size # of elements in array KerDer_bp_f_Gal     : ",        &
      size(KerDer_bp_f_Gal) 
   write(nout,*) " Size # of elements in array rag_bp_f            : ",        &
      size(rag_bp_f)
   write(nout,*) " Size # of elements in array surf_body_part      : ",        &
      size(surf_body_part) 
   write(nout,*) " Size # of elements in array nPartIntorno_bp_bp  : ",        &
      size(nPartIntorno_bp_bp)
   write(nout,*) " Size # of elements in array PartIntorno_bp_bp   : ",        &
      size(PartIntorno_bp_bp)
   write(nout,*) " Size # of elements in array rag_bp_bp           : ",        &
      size(rag_bp_bp)
endif
write(nout,'(1x,a)') "..."
write(nout,*) " Size in bytes of array pg                       : ",sizeof(pg)
write(nout,*) " Size in bytes of array BoundaryDataTab          : ",           &
   sizeof(BoundaryDataTab)
write(nout,*) " Size in bytes of array BoundaryDataPointer      : ",           &
   sizeof(BoundaryDataPointer)
write(nout,*) " Size in bytes of array Array_Flu                : ",           &
   sizeof(Array_Flu)
write(nout,*) " Size in bytes of array nPartIntorno             : ",           &
   sizeof(nPartIntorno)
write(nout,*) " Size in bytes of array PartIntorno              : ",           &
   sizeof(PartIntorno)
write(nout,*) " Size in bytes of array PartKernel               : ",           &
   sizeof(PartKernel)
write(nout,*) " Size in bytes of array rag                      : ",sizeof(rag)
if ((Domain%tipo=="bsph").and.(DBSPH%n_w>0)) then
   write(nout,*) " Size in bytes of array pg_w                     : ",        &
      sizeof(pg_w)
   write(nout,*) " Size in bytes of array nPartIntorno_fw          : ",        &
      sizeof(nPartIntorno_fw)
   write(nout,*) " Size in bytes of array PartIntorno_fw           : ",        &
      sizeof(PartIntorno_fw)
   write(nout,*) " Size in bytes of array kernel_fw                : ",        &
      sizeof(kernel_fw)
   write(nout,*) " Size in bytes of array rag_fw                   : ",        &
      sizeof(rag_fw)
   write(nout,*) " Size in bytes of array grad_vel_VSL_fw          : ",        &
      sizeof(grad_vel_VSL_fw)      
endif
if (n_bodies>0) then
   write(nout,*) " Size in bytes of array nPartIntorno_bp_f        : ",        &
      sizeof(nPartIntorno_bp_f)
   write(nout,*) " Size in bytes of array PartIntorno_bp_f         : ",        &
      sizeof(PartIntorno_bp_f)
   write(nout,*) " Size in bytes of array KerDer_bp_f_cub_spl      : ",        &
      sizeof(KerDer_bp_f_cub_spl)
   write(nout,*) " Size in bytes of array KerDer_bp_f_Gal          : ",        &
      sizeof(KerDer_bp_f_Gal)
   write(nout,*) " Size in bytes of array rag_bp_f                 : ",        &
      sizeof(rag_bp_f)
   write(nout,*) " Size in bytes of array surf_body_part           : ",        &
      sizeof(surf_body_part) 
   write(nout,*) " Size in bytes of array nPartIntorno_bp_bp       : ",        &
      sizeof(nPartIntorno_bp_bp)
   write(nout,*) " Size in bytes of array PartIntorno_bp_bp        : ",        &
      sizeof(PartIntorno_bp_bp)
   write(nout,*) " Size in bytes of array rag_bp_bp                : ",        &
      sizeof(rag_bp_bp)
endif
write(nout,'(1x,a)') "..."
write(nout,'(1x,a)') " "
write(nout,'(1x,a)') "   end allocation step. "
write(nout,'(1x,a)') " "
! Writing the cell number and particles within the cell 
if (Domain%ioutopt<0) then
   write(nout,*) 
   write(nout,*) "Number of cells           NCELLS= ",grid%nmax
   write(nout,*) 
   write(nout,*) "======== CELLS AND RELATED PARTICLES =========="
   do i=1,grid%nmax   
      if (Icont(i+1)>Icont(i)) then  
         write(nout,"(3(a,i5),a)") " cell", i," from",Icont(i)," to",          &
            Icont(i+1),"  particles:"
         write(nout,"(5i8)") NPartOrd(Icont(i):Icont(i+1)-1) 
      endif
   enddo
endif
write(nout,*) 
write(nout,*) 
call s_ctime(nout)
write(nout,*) 
write(nout,*) "Transient loop begins..."
write(nout,*)
! To initialize the post-processing file 
if ((Domain%imemo_fr>0).OR.(Domain%memo_fr>zero)) then
   open (nres,file=nomefile(2),status="unknown",access="sequential"            &
      ,form="unformatted")
   else
      nres = -nres
endif
if ((Domain%ipllb_fr>0).OR.(Domain%pllb_fr>zero)) then
   open (nplb,file=nomefile(4),status="unknown",access="sequential",           &
      form="formatted")
   write(nplb,"(a)") "time          free_surface_quota"
   else
      nplb = - nplb
endif
if ((Domain%imemo_fr>0).OR.(Domain%memo_fr>zero)) then
   open (nfro,file=nomefile(5),status="unknown",access="sequential"            &
      ,form="formatted")
   write(nfro,"(a)") "time          x fronte      (y fronte)    z fronte"
   else
      nfro = - nfro
endif
! To create vtk file to design the boundaries of the domain
if (vtkconv) then
   prefix = nomecaso
   filevtk = "VTKConverter_"//prefix(1:len_trim(prefix))//"_domain.vtk"
   open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',       &
      status='unknown')
   write(unitvtk,'(a)') '# vtk DataFile Version 2.0'
   write(unitvtk,'(a,a)') 'Domain limits for the case:',                       &
      prefix(1:len_trim(prefix))
   write(unitvtk,'(a)') 'ASCII'
   write(unitvtk,'(a)') 'DATASET POLYDATA'
   write(unitvtk,'(a,i8,a)') 'POINTS ',numvertici,' float'
   do i=1,numvertici,4
      k1 = i
      k2 = k1 + 3
      if (k2>numvertici) k2 = numvertici
      write(stringa,'(4(3(e12.5,1x)))') (vertice(1,k),vertice(2,k),            &
         vertice(3,k),k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
   enddo
! 2D case
   if (ncord==2) then
      nlinee = 0
      nvalori = 0
      do i=1,numbsides
         if (boundaryside(i)%TIPO=='peri') cycle
         nlinee = nlinee + 1
         nvalori = nvalori + 1
         nvalori = nvalori + 2
      enddo
      write(unitvtk,'(a,2i8)') 'LINES ', nlinee, nvalori
      allocate(verticecolore(numvertici))
      verticecolore = zero
      do i=1,numbsides
         if (boundaryside(i)%TIPO=='peri') cycle
         stringa = ' '
         k1 = boundaryside(i)%VERTEX(1) - 1
         k2 = boundaryside(i)%VERTEX(2) - 1
         write(stringa,'(a,2i10)') '2',k1,k2
         stringa = adjustl(trim(stringa))
         write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
         if ((boundaryside(i)%TIPO=='velo').or.(boundaryside(i)%TIPO=='flow')  &
            .or.(boundaryside(i)%TIPO=='open')) then
            verticecolore(k1+1) = 2.0
            verticecolore(k2+1) = 2.0
            else if (boundaryside(i)%TIPO=='sour') then
               verticecolore(k1+1) = 1.0
               verticecolore(k2+1) = 1.0
         endif
      enddo
      write(unitvtk,'(a,i8)') 'POINT_DATA ', numvertici
      write(unitvtk,'(a,a,a)') 'SCALARS ', prefix(1:len_trim(prefix)),' float 1'
      write(unitvtk,'(a)') 'LOOKUP_TABLE mytable'
      do i=1,numvertici,8
         k1 = i
         k2 = k1 + 7
         if (k2>numvertici) k2 = numvertici
         write(stringa,'(8(f10.3,1x))') (verticecolore(k),k=k1,k2)
         stringa = adjustl(trim(stringa))
         write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
      enddo
      deallocate(verticecolore)
      write(unitvtk,'(a)') 'LOOKUP_TABLE mytable 3'
      write(unitvtk,'(a)') '0.0 0.0 0.0 1.0'
      write(unitvtk,'(a)') '0.0 0.0 1.0 1.0'
      write(unitvtk,'(a)') '1.0 0.0 0.0 1.0'
! 3D case
      else
         nlinee = 0
         nvalori = 0
         do i=1,numfacce
            if (tratto(BoundaryFace(i)%stretch)%tipo=='peri') cycle
            nlinee = nlinee + 1
            nvalori = nvalori + 1 + BoundaryFace(i)%nodes + 1
         enddo
         write(unitvtk,'(a,2i8)') 'LINES ',nlinee,nvalori
         allocate(verticecolore(numvertici))
         verticecolore = 0.0
         do i=1,numfacce
            if (tratto(BoundaryFace(i)%stretch)%tipo=='peri') cycle
            stringa = ' '
            do k=1,BoundaryFace(i)%nodes,8
               k1 = k 
               k2 = k1 + 7
               if (k2>BoundaryFace(i)%nodes) k2 = BoundaryFace(i)%nodes
               write(stringa,'(10i10)') BoundaryFace(i)%nodes+1,               &
                 (BoundaryFace(i)%node(kk)%name-1,kk=k1,k2),                   &
                 BoundaryFace(i)%node(1)%name-1
               stringa = adjustl(trim(stringa))
               write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
            enddo
            if ((tratto(BoundaryFace(i)%stretch)%tipo=='velo').or.             &
               (tratto(BoundaryFace(i)%stretch)%tipo=='flow').or.              &
               (tratto(BoundaryFace(i)%stretch)%tipo=='open')) then
               do kk=k1,k2
                  verticecolore(BoundaryFace(i)%node(kk)%name) = 2.0
                  verticecolore(BoundaryFace(i)%node(kk)%name) = 2.0
               enddo
               else if (tratto(BoundaryFace(i)%stretch)%tipo=='sour') then
                  do kk=k1,k2
                     verticecolore(BoundaryFace(i)%node(kk)%name) = 1.0
                     verticecolore(BoundaryFace(i)%node(kk)%name) = 1.0
                  enddo
            endif
         enddo
         write(unitvtk,'(a,i8)') 'POINT_DATA ',numvertici
         write(unitvtk,'(a,a,a)') 'SCALARS ',prefix(1:len_trim(prefix)),       &
            ' float 1'
         write(unitvtk,'(a)') 'LOOKUP_TABLE mytable'
         do i=1,numvertici,8
            k1 = i
            k2 = k1 + 7
            if (k2>numvertici) k2 = numvertici
            write(stringa,'(8(f10.3,1x))') (verticecolore(k),k=k1,k2)
            stringa = adjustl(trim(stringa))
            write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
         enddo
         deallocate(verticecolore)
         write(unitvtk,'(a)') 'LOOKUP_TABLE mytable 3'
         write(unitvtk,'(a)') '0.0 0.0 0.0 1.0'
         write(unitvtk,'(a)') '0.0 0.0 1.0 1.0'
         write(unitvtk,'(a)') '1.0 0.0 0.0 1.0'
   endif
   flush(unitvtk)
   close (unitvtk)
endif
if ((Domain%tipo=="semi").or.(Domain%tipo=="bsph")) then
! To estimate "DTmin" for boundary elastic reaction 
   call EvaluateBER_TimeStep
! In case of bed-load transport
   if (Granular_flows_options%ID_erosion_criterion>0) then        
! Allocation of the 2D array of the minimum saturation flag (bed-load transport)
      if (.not.allocated(                                                      &
         Granular_flows_options%minimum_saturation_flag))then
         allocate(Granular_flows_options%minimum_saturation_flag(              &
            Grid%ncd(1),Grid%ncd(2)),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(nout,*) 'Allocation of ',                                    &
               'Granular_flows_options%minimum_saturation_flag ',              &
               'failed; the program stops here. '
            call diagnostic(arg1=4,arg2=1,arg3=nomsub)
            stop 
            else
               write (nout,*) 'Allocation of ',                                &
                  'Granular_flows_options%minimum_saturation_flag is',         &
                  ' successfully completed.'
         endif
      endif 
! Allocation of the 2D array of the maximum saturation flag (bed-load transport)
      if (.not.allocated(                                                      &
         Granular_flows_options%maximum_saturation_flag)) then
         allocate(Granular_flows_options%maximum_saturation_flag(              &
            Grid%ncd(1),Grid%ncd(2)),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(nout,*) 'Allocation of ',                                    &
               'Granular_flows_options%maximum_saturation_flag ',              &
               'failed; the program stops here. '
            call diagnostic(arg1=4,arg2=1,arg3=nomsub)
            stop 
            else
               write (nout,*) 'Allocation of ',                                &
                  'Granular_flows_options%maximum_saturation_flag is',         &
                  ' successfully completed.'
         endif
      endif
! Allocation of the 2D array of the saturation conditions
      if (.not.allocated(Granular_flows_options%saturation_conditions)&
         ) then
         allocate(Granular_flows_options%saturation_conditions(                &
            Grid%ncd(1),Grid%ncd(2)),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(nout,*) 'Allocation of ',                                    &
               'Granular_flows_options%saturation_conditions failed; '         &
               ,'the program stops here. '
            call diagnostic(arg1=4,arg2=1,arg3=nomsub)
            stop 
            else
               write (nout,*) 'Allocation of ',                                &
                  'Granular_flows_options%saturation_conditions is ',          &
                  'successfully completed.'
         endif
      endif
   endif
   if (ncord==2) then
      call start_and_stop(2,5)
      call loop_irre_2D 
      call start_and_stop(3,5)
      elseif (ncord==3) then
         NumCellmax = Grid%nmax
         if ((restart.eqv..false.).or.(Domain%tipo=="bsph")) then
            allocate(GCBFPointers(NumCellmax,2),stat=ier)
            if (ier/=0) then
               write(nout,'(1x,a,i2)')                                         &
                  "   Array GCBFPointers not allocated. Error code: ",ier
               call diagnostic(arg1=4,arg3=nomsub)
               else
                  write(nout,'(1x,a)')                                         &
                     "   Array GCBFPointers successfully allocated"
            endif
         endif
         if (Domain%tipo=="semi") then
! To select the grid cells intercepting a boundary faces
            if (restart.eqv..false.)                                           &
               call GridCellBoundaryFacesIntersections3D(NumCellmax)
! To compute the local coordinates, solid angle and solid normal, relative 
! to a boundary element (SA-SPH)
            call ComputeBoundaryIntegralTab
! Computation of the boundary contributions for the continuity equation (SA-SPH)
            call ComputeKernelTable
         endif
! Allocation and initialization of the Z_fluid_max array
! Loop over the zones
         do i=1,NPartZone
            if (Partz(i)%IC_source_type==2) then
               if (.not.allocated(Z_fluid_max)) then
                  allocate(Z_fluid_max(Grid%ncd(1)*Grid%ncd(2)),STAT=alloc_stat)
                  if (alloc_stat/=0) then
                     write(nout,*)                                             &
                     'Allocation of Z_fluid_max in Gest_Trans failed;',        &
                     ' the program terminates here.'
                     stop ! Stop the main program
                     else
                        write(nout,*)                                          &
                           'Allocation of Z_fluid_max in Gest_Trans ',         &
                           'successfully completed.'
                  endif
                  Z_fluid_max(:) = -999.d0
               endif
               aux_integer = Partz(i)%ID_last_vertex -                         &
                             Partz(i)%ID_first_vertex + 1  
               if (.not.allocated(q_max)) then
                  allocate(q_max(aux_integer),STAT=alloc_stat)
                  if (alloc_stat/=0) then
                     write(nout,*)                                             &
                     'Allocation of q_max in Gest_Trans failed;',              &
                     ' the program terminates here.'
                     stop ! Stop the main program
                     else
                        write(nout,*)                                          &
                           'Allocation of q_max in Gest_Trans successfully ',  &
                           'completed.'
                  endif
                  q_max(:) = 0.d0
               endif
               exit
            endif
         enddo  
         call start_and_stop(2,5)
! Main loop
         call loop_irre_3D
         call start_and_stop(3,5)
! Writing the h_max array and deallocation of the Z_fluid_max array
         if (allocated(Z_fluid_max)) then
            call write_h_max      
            if (allocated(Z_fluid_max)) deallocate(Z_fluid_max)
         endif
   endif 
   else
      call diagnostic(arg1=10,arg2=5,arg3=nomsub)
endif
! To create vtk file
if (vtkconv) then
   prefix = nomecaso
   filevtk = "VTKConverter_"//prefix(1:len_trim(prefix))//".pvd"
   open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',       &
      status='unknown')
   write(unitvtk,'(a)') '<?xml version="1.0"?>'
   write(unitvtk,'(2a)') '<VTKFile type="Collection" version="0.1" ',          &
      'byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
   write(unitvtk,'(a)') '  <Collection>'
   do i=1,nblocchi
      stringa = ' '
      write(cargo,'(i10)') blocchi(i)
      cargo = adjustl(trim(cargo))
      filename =                                                               &
"VTKConverter_"//prefix(1:len_trim(prefix))//"_block_"//cargo(1:len_trim(cargo))//".vtu"
      stringa = '  <DataSet timestep="'
      write(cargo,'(f15.6)')  Time_Block(i)
      cargo = adjustl(trim(cargo))
      stringa =                                                                &
stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" group="" part="'
      write(cargo,'(i10)') i
      cargo = adjustl(trim(cargo))
      stringa =                                                                &
stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" file="'//filename(1:len_trim(filename))//'"/>'
      write(unitvtk,'(a)') stringa(1:len_trim(stringa))
   enddo
   write(unitvtk,'(a)') '  </Collection>'
   write(unitvtk,'(a)') '</VTKFile>'
   flush(unitvtk)
   close (unitvtk)
endif
! To create vtk file for wall elements (DB-SPH)
if ((Domain%tipo=="bsph").and.(DBSPH%n_w>0)) then
   if (vtkconv) then
      filevtk = "VTKConverter_wall_"//prefix(1:len_trim(prefix))//".pvd"
      open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',    &
         status='unknown')
      write(unitvtk,'(a)') '<?xml version="1.0"?>'
      write(unitvtk,'(2a)') '<VTKFile type="Collection" version="0.1" ',       &
         'byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
      write(unitvtk,'(a)') '  <Collection>'
      do i=1,nblocchi
         stringa = ' '
         write(cargo,'(i10)') blocchi(i)
         cargo = adjustl(trim(cargo))
         filename =                                                            &
"VTKConverter_"//prefix(1:len_trim(prefix))//"_block_wall_"//cargo(1:len_trim(cargo))//".vtu"
         stringa = '  <DataSet timestep="'
         write(cargo,'(f15.6)')  Time_Block(i)
         cargo = adjustl(trim(cargo))
         stringa =                                                             &
stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" group="" part="'
         write(cargo,'(i10)') i
         cargo = adjustl(trim(cargo))
         stringa =                                                             &
stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" file="'//filename(1:len_trim(filename))//'"/>'
         write(unitvtk,'(a)') stringa(1:len_trim(stringa))
      enddo
      write(unitvtk,'(a)') '  </Collection>'
      write(unitvtk,'(a)') '</VTKFile>'
      flush(unitvtk)
      close (unitvtk)
   endif
endif
if (n_bodies>0) then
! Creation of the .pvd file for body particles (Body Transport)
   if (vtkconv) then
      filevtk = "VTKConverter_body-part_"//prefix(1:len_trim(prefix))//".pvd"
      open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',    &
         status='unknown')
      write(unitvtk,'(a)') '<?xml version="1.0"?>'
      write(unitvtk,'(2a)') '<VTKFile type="Collection" version="0.1" ',       &
         'byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
      write(unitvtk,'(a)') '  <Collection>'
      do i=1,nblocchi
         stringa = ' '
         write(cargo,'(i10)') blocchi(i)
         cargo = adjustl(trim(cargo))
         filename =                                                            &
"VTKConverter_"//prefix(1:len_trim(prefix))//"_block_body-part_"//cargo(1:len_trim(cargo))//".vtu"
         stringa = '  <DataSet timestep="'
         write(cargo,'(f15.6)')  Time_Block(i)
         cargo = adjustl(trim(cargo))
         stringa =                                                             &
stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" group="" part="'
         write(cargo,'(i10)') i
         cargo = adjustl(trim(cargo))
         stringa =                                                             &
stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" file="'//filename(1:len_trim(filename))//'"/>'
         write(unitvtk,'(a)') stringa(1:len_trim(stringa))
      enddo
      write(unitvtk,'(a)') '  </Collection>'
      write(unitvtk,'(a)') '</VTKFile>'
      flush(unitvtk)
      close (unitvtk)
   endif
! Creation of the .pvd file for bodies (Body Transport)
   if (vtkconv) then
      filevtk = "VTKConverter_body_"//prefix(1:len_trim(prefix))//".pvd"
      open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',    &
         status='unknown')
      write(unitvtk,'(a)') '<?xml version="1.0"?>'
      write(unitvtk,'(2a)') '<VTKFile type="Collection" version="0.1" ',       &
         'byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
      write(unitvtk,'(a)') '  <Collection>'
      do i=1,nblocchi
         stringa = ' '
         write(cargo,'(i10)') blocchi(i)
         cargo = adjustl(trim(cargo))
         filename =                                                            &
"VTKConverter_"//prefix(1:len_trim(prefix))//"_block_body_"//cargo(1:len_trim(cargo))//".vtu"
         stringa = '  <DataSet timestep="'
         write(cargo,'(f15.6)')  Time_Block(i)
         cargo = adjustl(trim(cargo))
         stringa =                                                             &
stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" group="" part="'
         write(cargo,'(i10)') i
         cargo = adjustl(trim(cargo))
         stringa =                                                             &
stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" file="'//filename(1:len_trim(filename))//'"/>'
         write(unitvtk,'(a)') stringa(1:len_trim(stringa))
      enddo
      write(unitvtk,'(a)') '  </Collection>'
      write(unitvtk,'(a)') '</VTKFile>'
      flush(unitvtk)
      close (unitvtk)
   endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine Gest_Trans

