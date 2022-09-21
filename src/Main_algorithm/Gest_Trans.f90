!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.9.0.0
! SPHERA v.9.0.0 is free software: you can redistribute it and/or modify
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
! Description: Introductory procedure for the main algorithm. Writing of the  
!              ".vtk" geometry file.
!-------------------------------------------------------------------------------
subroutine Gest_Trans 
!------------------------
! Modules
!------------------------
use Hybrid_allocation_module
use I_O_file_module
use Static_allocation_module
use Dynamic_allocation_module
use Memory_I_O_interface_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),parameter :: ner0 = 0
integer(4) :: npi,i,ier,k,k1,k2,nlinee,nvalori,alloc_stat
#ifdef SPACE_3D
integer(4) :: n_vertices_main_wall,NumCellmax,kk
character(100) :: array_name
#endif
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
interface
#ifdef SPACE_3D
   subroutine main_wall_info(n_vertices_main_wall,ID_main_wall)
      integer(4),intent(out) :: n_vertices_main_wall
      integer(4),intent(out),optional :: ID_main_wall
   end subroutine main_wall_info
#endif
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
write(ulog,"(/)")
write(ulog,*) "Initial number of particles      NAG= ",nag
write( ulog,"(/)" )
if (Domain%ioutopt<0) then
   write(ulog,*)
   write(ulog,*) "======== PARTICLES COORDINATES =========="
   do npi=1,nag
      write(ulog,"(i10,4f14.5)") npi,pg(npi)%coord,pg(npi)%tstop
   enddo
endif
! Array initialization for vtkconverter
nblocchi = 0
blocchi = 0
block = - 1
Time_Block = zero
#ifdef SPACE_3D
MaxNcbf = int(input_any_t%MAXCLOSEBOUNDFACES * PARTICLEBUFFER)
#elif defined SPACE_2D
MaxNcbs = int(MAXCLOSEBOUNDSIDES * PARTICLEBUFFER)
#endif
if ((Domain%tipo=="semi").or.(Domain%tipo=="bsph"))  then
   allocate(BoundaryDataPointer(1:3,1:PARTICLEBUFFER),stat=ier)
   if (ier/=0) then
      write(ulog,'(1x,a,i2)')                                                  &
         "   Array BoundaryDataPointer not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)')                                                  &
            "   Array BoundaryDataPointer successfully allocated "
   endif
#ifdef SPACE_2D
      write(ulog,'(a,i15)') "     Max num of close boundary sides: MaxNcbs = ",&
         MaxNcbs
      allocate(BoundaryDataTab(1:MaxNcbs),stat=ier)
      if (ier/=0) then
         write(ulog,'(1x,a,i2)')                                               &
            "   Array BoundaryDataTab not allocated. Error code: ",ier
          call diagnostic(arg1=4,arg3=nomsub)
         else
            write(ulog,'(1x,a)')                                               &
               "   Array BoundaryDataTab successfully allocated "
      endif
#elif defined SPACE_3D
         write(ulog,'(a,i15)')                                                 &
            "     Max num of close boundary faces: MaxNcbf = ",MaxNcbf
         allocate(BoundaryDataTab(1:MaxNcbf),stat=ier)
         if (ier/=0) then
            write(ulog,'(1x,a,i2)')                                            &
               "   Array BoundaryDataTab not allocated. Error code: ",ier
            call diagnostic(arg1=4,arg3=nomsub)
            else
               write(ulog,'(1x,a)')                                            &
                  "   Array BoundaryDataTab successfully allocated"
         endif
#endif
endif
allocate(Array_Flu(1:PARTICLEBUFFER),stat=ier)
if (ier/=0) then
   write(ulog,'(1x,a,i2)') "   Array Array_Flu not allocated. Error code: ",   &
      ier
   call diagnostic(arg1=4,arg3=nomsub)
   else
      write(ulog,'(1x,a)') "   Array Array_Flu successfully allocated "
endif
allocate(nPartIntorno(1:PARTICLEBUFFER),stat=ier)
if (ier/=0) then
   write(ulog,'(1x,a,i2)')"   Array NPARTINTORNO not allocated. Error code: "  &
   ,ier
   call diagnostic(arg1=4,arg3=nomsub)
   else
      write(ulog,'(1x,a)') "   Array NPARTINTORNO successfully allocated "
endif
allocate(PartIntorno(1:NMAXPARTJ*PARTICLEBUFFER),stat=ier)
if (ier/=0) then
   write(ulog,'(1x,a,i2)') "   Array PARTINTORNO not allocated. Error code: "  &
      ,ier
   call diagnostic(arg1=4,arg3=nomsub)
   else
      write(ulog,'(1x,a)') "   Array PARTINTORNO successfully allocated "
endif
allocate(PartKernel(1:4,1:NMAXPARTJ*PARTICLEBUFFER),stat=ier)
if (ier/=0) then
   write(ulog,'(1x,a,i2)') "   Array PARTKERNEL not allocated. Error code: ",  &
      ier
   call diagnostic(arg1=4,arg3=nomsub)
   else
      write(ulog,'(1x,a)') "   Array PARTKERNEL successfully allocated "
endif
allocate(rag(1:3,1:NMAXPARTJ*PARTICLEBUFFER),stat=ier)
if (ier/=0) then
   write(ulog,'(1x,a,i2)') "   Array RAG not allocated. Error code: ",ier
   call diagnostic(arg1=4,arg3=nomsub)
   else
      write(ulog,'(1x,a)') "   Array RAG successfully allocated "
endif
if (Domain%tipo=="bsph") then
   allocate(nPartIntorno_fw(1:PARTICLEBUFFER),stat=ier)
   if (ier/=0) then
      write(ulog,'(1x,a,i2)')                                                  &
         "   Array NPARTINTORNO_fw not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)')                                                  &
            "   Array NPARTINTORNO_fw successfully allocated "
   endif
   allocate(PartIntorno_fw(1:NMAXPARTJ*PARTICLEBUFFER),stat=ier)
   if (ier/=0) then
      write(ulog,'(1x,a,i2)')                                                  &
         "   Array PARTINTORNO_fw not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)') "   Array PARTINTORNO_fw successfully allocated"
   endif
   allocate(grad_vel_VSL_fw(1:3,1:NMAXPARTJ*PARTICLEBUFFER),stat=ier)
   if (ier/=0) then
      write(ulog,'(1x,a,i2)')                                                  &
         "   Array grad_vel_VSL_fw not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)') "   Array grad_vel_VSL_fw successfully allocated "
   endif
   allocate(kernel_fw(2,1:NMAXPARTJ*PARTICLEBUFFER),stat=ier)
   if (ier/=0) then
      write(ulog,'(1x,a,i2)')                                                  &
         "   Array kernel_fw not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)') "   Array kernel_fw successfully allocated "
   endif
   allocate(rag_fw(1:3,1:NMAXPARTJ*PARTICLEBUFFER),stat=ier)
   if (ier/=0) then
      write(ulog,'(1x,a,i2)') "   Array RAG_fw not allocated. Error code: ",   &
         ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)') "   Array RAG_fw successfully allocated "
   endif
endif
#ifdef SOLID_BODIES
   allocate(nPartIntorno_bp_f(n_body_part),stat=ier)
   if (ier/=0) then
      write(ulog,'(1x,a,i2)')                                                  &
         "   Array nPartIntorno_bp_f not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)')                                                  &
            "   Array nPartIntorno_bp_f successfully allocated "
   endif
   allocate(PartIntorno_bp_f(NMAXPARTJ*n_body_part),stat=ier)
   if (ier/=0) then
      write(ulog,'(1x,a,i2)')                                                  &
         "   Array PartIntorno_bp_f not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)')                                                  &
            "   Array PartIntorno_bp_f successfully allocated "
   endif
   allocate(proxy_normal_bp_f(NMAXPARTJ*n_body_part),stat=ier)
   if (ier/=0) then
      write(ulog,'(1x,a,i2)')                                                  &
         "   Array proxy_normal_bp_f not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)')                                                  &
            "   Array proxy_normal_bp_f successfully allocated "
   endif
   allocate(KerDer_bp_f_cub_spl(NMAXPARTJ*n_body_part),stat=ier)
   if (ier/=0) then
      write(ulog,'(1x,a,i2)')                                                  &
         "   Array KerDer_bp_f_cub_spl not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)')                                                  &
            "   Array KerDer_bp_f_cub_spl successfully allocated "
   endif
   allocate(KerDer_bp_f_Gal(NMAXPARTJ*n_body_part),stat=ier)
   if (ier/=0) then
      write(ulog,'(1x,a,i2)')                                                  &
         "   Array KerDer_bp_f_Gal not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)')                                                  &
            "   Array KerDer_bp_f_Gal successfully allocated "
   endif 
   allocate(rag_bp_f(3,NMAXPARTJ*n_body_part),stat=ier)
   if (ier/=0) then
      write(ulog,'(1x,a,i2)') "   Array rag_bp_f not allocated. Error code: "  &
         ,ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)') "   Array rag_bp_f successfully allocated "
   endif
   allocate(nPartIntorno_bp_bp(n_surf_body_part),stat=ier)
   if (ier/=0) then
      write(ulog,'(1x,a,i2)')                                                  &
         "   Array nPartIntorno_bp_bp not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
      write(ulog,'(1x,a)')                                                     &
         "   Array nPartIntorno_bp_bp successfully allocated "
   endif
   allocate(PartIntorno_bp_bp(n_surf_body_part*NMAXPARTJ),stat=ier)
   if (ier/=0) then
      write(ulog,'(1x,a,i2)')                                                  &
         "   Array PartIntorno_bp_bp not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)')                                                  &
            "   Array PartIntorno_bp_bp successfully allocated "
   endif
   allocate(rag_bp_bp(3,n_surf_body_part*NMAXPARTJ),stat=ier)
   if (ier/=0) then
      write(ulog,'(1x,a,i2)')                                                  &
         "   Array rag_bp_bp not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)') "   Array rag_bp_bp successfully allocated "
   endif
#ifdef SPACE_3D
      allocate(impact_vel(n_surf_body_part,(n_bodies+NumFacce)),stat=ier)
#elif defined SPACE_2D
         allocate(impact_vel(n_surf_body_part,(n_bodies+NumBSides)),stat=ier)
#endif
   if (ier/=0) then
      write(ulog,'(1x,a,i2)')                                                  &
         "   Array impact_vel not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)') "   Array impact_vel successfully allocated "
   endif
   impact_vel(:,:) = 0.d0
   if (FSI_free_slip_conditions) then
      array_name = "nPartIntorno_f_sbp"
      call allocate_de_int4_r1(.true.,nPartIntorno_f_sbp,PARTICLEBUFFER,       &
         array_name,ulog_flag=.true.)
      array_name = "PartIntorno_f_sbp"
      call allocate_de_int4_r1(.true.,PartIntorno_f_sbp,                       &
         extent_1=NMAXPARTJ*PARTICLEBUFFER,array_name=array_name,              &
         ulog_flag=.true.)
      array_name = "dis_f_sbp"
      call allocate_de_dp_r1(.true.,dis_f_sbp,                                 &
         extent_1=NMAXPARTJ*PARTICLEBUFFER,array_name=array_name,              &
         ulog_flag=.true.)
      array_name = "closest_f_sbp"
      call allocate_de_int4_r1(.true.,closest_f_sbp,                           &
         extent_1=NMAXPARTJ*PARTICLEBUFFER,array_name=array_name,              &
         ulog_flag=.true.)
   endif
#endif
write(ulog,'(1x,a)') "..."
write(ulog,'(a,i15)') " Max number of particles  : PARTICLEBUFFER = ",         &
   PARTICLEBUFFER
write(ulog,*) " Size # of elements in array pg                  : ",size(pg)
write(ulog,*) " Size # of elements in array BoundaryDataTab     : ",           &
   size(BoundaryDataTab)
write(ulog,*) " Size # of elements in array BoundaryDataPointer : ",           &
   size(BoundaryDataPointer)
write(ulog,*) " Size # of elements in array Array_Flu           : ",           &
   size(Array_Flu)
write(ulog,*) " Size # of elements in array nPartIntorno        : ",           &
   size(nPartIntorno)
write(ulog,*) " Size # of elements in array PartIntorno         : ",           &
   size(PartIntorno)
write(ulog,*) " Size # of elements in array PartKernel          : ",           &
   size(PartKernel)
write(ulog,*) " Size # of elements in array rag                 : ",           &
   size(rag)
if ((Domain%tipo=="bsph").and.(DBSPH%n_w>0)) then
   write(ulog,*) " Size # of elements in array pg_w                : ",        &
      size(pg_w)
   write(ulog,*) " Size # of elements in array nPartIntorno_fw     : ",        &
      size(nPartIntorno_fw)
   write(ulog,*) " Size # of elements in array PartIntorno_fw      : ",        &
      size(PartIntorno_fw)
   write(ulog,*) " Size # of elements in array kernel_fw           : ",        &
      size(kernel_fw)
   write(ulog,*) " Size # of elements in array rag_fw              : ",        &
      size(rag_fw)
   write(ulog,*) " Size # of elements in array grad_vel_VSL_fw     : ",        &
      size(grad_vel_VSL_fw)
endif
#ifdef SOLID_BODIES
   write(ulog,*) " Size # of elements in array nPartIntorno_bp_f   : ",        &
      size(nPartIntorno_bp_f)
   write(ulog,*) " Size # of elements in array PartIntorno_bp_f    : ",        &
      size(PartIntorno_bp_f)
   write(ulog,*) " Size # of elements in array proxy_normal_bp_f   : ",        &
      size(proxy_normal_bp_f)
   write(ulog,*) " Size # of elements in array KerDer_bp_f_cub_spl : ",        &
      size(KerDer_bp_f_cub_spl)
   write(ulog,*) " Size # of elements in array KerDer_bp_f_Gal     : ",        &
      size(KerDer_bp_f_Gal) 
   write(ulog,*) " Size # of elements in array rag_bp_f            : ",        &
      size(rag_bp_f)
   write(ulog,*) " Size # of elements in array surf_body_part      : ",        &
      size(surf_body_part) 
   write(ulog,*) " Size # of elements in array nPartIntorno_bp_bp  : ",        &
      size(nPartIntorno_bp_bp)
   write(ulog,*) " Size # of elements in array PartIntorno_bp_bp   : ",        &
      size(PartIntorno_bp_bp)
   write(ulog,*) " Size # of elements in array rag_bp_bp           : ",        &
      size(rag_bp_bp)
   if (FSI_free_slip_conditions) then
      write(ulog,*) " Size # of elements in array nPartIntorno_f_sbp  : ",     &
         size(nPartIntorno_f_sbp)
      write(ulog,*) " Size # of elements in array PartIntorno_f_sbp   : ",     &
         size(PartIntorno_f_sbp)
      write(ulog,*) " Size # of elements in array dis_f_sbp           : ",     &
         size(dis_f_sbp)
      write(ulog,*) " Size # of elements in array closest_f_sbp       : ",     &
         size(closest_f_sbp)
   endif
#endif
write(ulog,'(1x,a)') "..."
write(ulog,*) " Size in bytes of array pg                       : ",sizeof(pg)
write(ulog,*) " Size in bytes of array BoundaryDataTab          : ",           &
   sizeof(BoundaryDataTab)
write(ulog,*) " Size in bytes of array BoundaryDataPointer      : ",           &
   sizeof(BoundaryDataPointer)
write(ulog,*) " Size in bytes of array Array_Flu                : ",           &
   sizeof(Array_Flu)
write(ulog,*) " Size in bytes of array nPartIntorno             : ",           &
   sizeof(nPartIntorno)
write(ulog,*) " Size in bytes of array PartIntorno              : ",           &
   sizeof(PartIntorno)
write(ulog,*) " Size in bytes of array PartKernel               : ",           &
   sizeof(PartKernel)
write(ulog,*) " Size in bytes of array rag                      : ",sizeof(rag)
if ((Domain%tipo=="bsph").and.(DBSPH%n_w>0)) then
   write(ulog,*) " Size in bytes of array pg_w                     : ",        &
      sizeof(pg_w)
   write(ulog,*) " Size in bytes of array nPartIntorno_fw          : ",        &
      sizeof(nPartIntorno_fw)
   write(ulog,*) " Size in bytes of array PartIntorno_fw           : ",        &
      sizeof(PartIntorno_fw)
   write(ulog,*) " Size in bytes of array kernel_fw                : ",        &
      sizeof(kernel_fw)
   write(ulog,*) " Size in bytes of array rag_fw                   : ",        &
      sizeof(rag_fw)
   write(ulog,*) " Size in bytes of array grad_vel_VSL_fw          : ",        &
      sizeof(grad_vel_VSL_fw)      
endif
#ifdef SOLID_BODIES
   write(ulog,*) " Size in bytes of array nPartIntorno_bp_f        : ",        &
      sizeof(nPartIntorno_bp_f)
   write(ulog,*) " Size in bytes of array PartIntorno_bp_f         : ",        &
      sizeof(PartIntorno_bp_f)
   write(ulog,*) " Size in bytes of array proxy_normal_bp_f        : ",        &
      sizeof(proxy_normal_bp_f)
   write(ulog,*) " Size in bytes of array KerDer_bp_f_cub_spl      : ",        &
      sizeof(KerDer_bp_f_cub_spl)
   write(ulog,*) " Size in bytes of array KerDer_bp_f_Gal          : ",        &
      sizeof(KerDer_bp_f_Gal)
   write(ulog,*) " Size in bytes of array rag_bp_f                 : ",        &
      sizeof(rag_bp_f)
   write(ulog,*) " Size in bytes of array surf_body_part           : ",        &
      sizeof(surf_body_part) 
   write(ulog,*) " Size in bytes of array nPartIntorno_bp_bp       : ",        &
      sizeof(nPartIntorno_bp_bp)
   write(ulog,*) " Size in bytes of array PartIntorno_bp_bp        : ",        &
      sizeof(PartIntorno_bp_bp)
   write(ulog,*) " Size in bytes of array rag_bp_bp                : ",        &
      sizeof(rag_bp_bp)
   if (FSI_free_slip_conditions) then
      write(ulog,*) " Size in bytes of array nPartIntorno_f_sbp       : ",     &
         sizeof(nPartIntorno_f_sbp)
      write(ulog,*) " Size in bytes of array PartIntorno_f_sbp        : ",     &
         sizeof(PartIntorno_f_sbp)
      write(ulog,*) " Size in bytes of array dis_f_sbp                : ",     &
         sizeof(dis_f_sbp)
      write(ulog,*) " Size in bytes of array closest_f_sbp            : ",     &
         sizeof(closest_f_sbp)
   endif
#endif
write(ulog,'(1x,a)') "..."
write(ulog,'(1x,a)') " "
write(ulog,'(1x,a)') "   end allocation step. "
write(ulog,'(1x,a)') " "
! Writing the cell number and particles within the cell 
if (Domain%ioutopt<0) then
   write(ulog,*) 
   write(ulog,*) "Number of cells           NCELLS= ",grid%nmax
   write(ulog,*) 
   write(ulog,*) "======== CELLS AND RELATED PARTICLES =========="
   do i=1,grid%nmax   
      if (Icont(i+1)>Icont(i)) then  
         write(ulog,"(3(a,i5),a)") " cell", i," from",Icont(i)," to",          &
            Icont(i+1),"  particles:"
         write(ulog,"(5i8)") NPartOrd(Icont(i):Icont(i+1)-1) 
      endif
   enddo
endif
write(ulog,*) 
write(ulog,*) 
call s_ctime
write(ulog,*) 
write(ulog,*) "Transient loop begins..."
write(ulog,*)
! To initialize the post-processing file 
if ((Domain%imemo_fr>0).or.(input_any_t%memo_fr>zero)) then
   open(nres,file=nomefile(2),status="unknown",access="sequential"             &
      ,form="unformatted")
   else
      nres = -nres
endif
if ((Domain%ipllb_fr>0).or.(input_any_t%pllb_fr>zero)) then
   open(nplb,file=nomefile(4),status="unknown",access="sequential",            &
      form="formatted")
   write(nplb,"(a)") "time          free_surface_quota"
   open(uzlft,file=nomefile(6),status="unknown",access="sequential",           &
      form="formatted")
   write(uzlft,"(a)") "time          lower_fluid_top_height"
   else
      nplb = - nplb
      uzlft = - uzlft
endif
if ((Domain%imemo_fr>0).or.(input_any_t%memo_fr>zero)) then
   open(nfro,file=nomefile(5),status="unknown",access="sequential"             &
      ,form="formatted")
   write(nfro,"(a)") "time          x fronte      (y fronte)    z fronte"
   else
      nfro = - nfro
endif
! To create vtk file to design the boundaries of the domain
if (vtkconv) then
   prefix = nomecaso(1:len(prefix))
   filevtk = "VTKConverter_"//prefix(1:len_trim(prefix))//"_domain.vtk"
   open(unit=unitvtk,file=filevtk,form='formatted',access='sequential',        &
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
#ifdef SPACE_2D
      nlinee = 0
      nvalori = 0
      do i=1,NumBSides
         if (boundaryside(i)%TIPO=='peri') cycle
         nlinee = nlinee + 1
         nvalori = nvalori + 1
         nvalori = nvalori + 2
      enddo
      write(unitvtk,'(a,2i11)') 'LINES ', nlinee, nvalori
      allocate(verticecolore(numvertici))
      verticecolore = zero
      do i=1,NumBSides
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
#elif defined SPACE_3D
         nlinee = 0
         nvalori = 0
         do i=1,NumFacce
            if (tratto(BoundaryFace(i)%stretch)%tipo=='peri') cycle
            nlinee = nlinee + 1
            nvalori = nvalori + 1 + BoundaryFace(i)%nodes + 1
         enddo
         write(unitvtk,'(a,2i11)') 'LINES ',nlinee,nvalori
         allocate(verticecolore(numvertici))
         verticecolore = 0.0
         do i=1,NumFacce
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
#endif
   flush(unitvtk)
   close (unitvtk)
endif
if ((Domain%tipo=="semi").or.(Domain%tipo=="bsph")) then
! To estimate "DTmin" for boundary elastic reaction 
   call EvaluateBER_TimeStep
! In case of bed-load transport
   if (Granular_flows_options%KTGF_config>0) then        
! Allocation of the 2D array of the minimum saturation flag (bed-load transport)
      if (.not.allocated(                                                      &
         Granular_flows_options%minimum_saturation_flag))then
         allocate(Granular_flows_options%minimum_saturation_flag(              &
            Grid%ncd(1),Grid%ncd(2)),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(ulog,*) 'Allocation of ',                                    &
               'Granular_flows_options%minimum_saturation_flag ',              &
               'failed; the program stops here. '
            call diagnostic(arg1=4,arg2=1,arg3=nomsub)
            stop 
            else
               write(ulog,*) 'Allocation of ',                                &
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
            write(ulog,*) 'Allocation of ',                                    &
               'Granular_flows_options%maximum_saturation_flag ',              &
               'failed; the program stops here. '
            call diagnostic(arg1=4,arg2=1,arg3=nomsub)
            stop 
            else
               write(ulog,*) 'Allocation of ',                                &
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
            write(ulog,*) 'Allocation of ',                                    &
               'Granular_flows_options%saturation_conditions failed; '         &
               ,'the program stops here. '
            call diagnostic(arg1=4,arg2=1,arg3=nomsub)
            stop 
            else
               write(ulog,*) 'Allocation of ',                                 &
                  'Granular_flows_options%saturation_conditions is ',          &
                  'successfully completed.'
         endif
      endif
   endif
#ifdef SPACE_2D
      call start_and_stop(2,5)
      call time_step_loop 
      call start_and_stop(3,5)
#elif defined SPACE_3D
         NumCellmax = Grid%nmax
         if ((restart.eqv..false.).or.(Domain%tipo=="bsph")) then
            allocate(GCBFPointers(NumCellmax,2),stat=ier)
            if (ier/=0) then
               write(ulog,'(1x,a,i2)')                                         &
                  "   Array GCBFPointers not allocated. Error code: ",ier
               call diagnostic(arg1=4,arg3=nomsub)
               else
                  write(ulog,'(1x,a)')                                         &
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
! Allocation and initialization of the arrays: Z_fluid_max, Z_fluid_step, 
! q_max, U_max
         call main_wall_info(n_vertices_main_wall)
         if (n_vertices_main_wall>0) then
! In the absence of walls, there is no allocation
            if (.not.allocated(Z_fluid_max)) then
               allocate(Z_fluid_max(Grid%ncd(1)*Grid%ncd(2),2),STAT=alloc_stat)
               if (alloc_stat/=0) then
                  write(ulog,*)                                                &
                  'Allocation of Z_fluid_max in Gest_Trans failed;',           &
                  ' the program terminates here.'
                  stop
                  else
                     write(ulog,*)                                             &
                        'Allocation of Z_fluid_max in Gest_Trans ',            &
                        'successfully completed.'
               endif
               Z_fluid_max(:,:) = -999.d0
            endif
            if (.not.allocated(Z_fluid_step)) then
               allocate(Z_fluid_step(Grid%ncd(1)*Grid%ncd(2),2),STAT=alloc_stat)
               if (alloc_stat/=0) then
                  write(ulog,*)                                                &
                  'Allocation of "Z_fluid_step" in the subroutine ',           &
                  '"Gest_Trans" failed; the execution terminates here.'
                  stop
                  else
                     write(ulog,*)                                             &
                        'Allocation of "Z_fluid_step" in the ',                &
                        'subroutine "Gest_Trans" is successfully completed.'
               endif
               Z_fluid_step(:,:) = -999.d0
            endif
            if (.not.allocated(q_max)) then
               allocate(q_max(n_vertices_main_wall),STAT=alloc_stat)
               if (alloc_stat/=0) then
                  write(ulog,*)                                                &
                  'Allocation of q_max in Gest_Trans failed;',                 &
                  ' the program terminates here.'
                  stop
                  else
                     write(ulog,*)                                             &
                        'Allocation of q_max in Gest_Trans successfully ',     &
                        'completed.'
               endif
               q_max(:) = 0.d0
            endif
            if (.not.allocated(U_max)) then
               allocate(U_max(n_vertices_main_wall),STAT=alloc_stat)
               if (alloc_stat/=0) then
                  write(ulog,*)                                                &
                  'Allocation of U_max in Gest_Trans failed;',                 &
                  ' the program terminates here.'
                  stop
                  else
                     write(ulog,*)                                             &
                        'Allocation of U_max in Gest_Trans successfully ',     &
                        'completed.'
               endif
               U_max(:) = 0.d0
            endif
            array_name = "z_topog_max"
            call allocate_de_dp_r1(.true.,z_topog_max,Grid%ncd(1)*Grid%ncd(2), &
               array_name,ulog_flag=.true.)
         endif
         call start_and_stop(2,5)
! Main loop
         call time_step_loop
         call start_and_stop(3,5)
! Writing the h_max array and deallocation of the arrays Z_fluid_max and 
! Z_fluid_step
         if (allocated(Z_fluid_max)) then
            call write_h_max
            deallocate(Z_fluid_max,STAT=alloc_stat)
            if (alloc_stat/=0) then
               write(uerr,*) 'Subroutine "Gest_Trans". Deallocation of the ',  &
                  'array "Z_fluid_max" failed; the simulation terminates ',    &
                  ' here.'
               stop
               else
               write(ulog,*) 'Subroutine "Gest_Trans". Deallocation of the ',  &
                  'array "Z_fluid_max" is successfully completed.'
            endif
            deallocate(Z_fluid_step,STAT=alloc_stat)
            if (alloc_stat/=0) then
               write(uerr,*) 'Subroutine "Gest_Trans". Deallocation of the ',  &
                  'array "Z_fluid_step" failed; the simulation terminates ',   &
                  ' here.'
               stop
               else
                  write(ulog,*) 'Subroutine "Gest_Trans". Deallocation of ',   &
                     'the array "Z_fluid_step" is successfully completed.'
            endif
         endif
         array_name = "z_topog_max"
         call allocate_de_dp_r1(.false.,z_topog_max,array_name=array_name,     &
            ulog_flag=.true.)
#endif
   else
      call diagnostic(arg1=10,arg2=5,arg3=nomsub)
endif
! To create the ".pvd" file
if (vtkconv) then
   prefix = nomecaso(1:len(prefix))
   filevtk = "VTKConverter_"//prefix(1:len_trim(prefix))//".pvd"
   open(unit=unitvtk,file=filevtk,form='formatted',access='sequential',        &
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
      write(cargo,'(f15.9)') Time_Block(i)
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
      open(unit=unitvtk,file=filevtk,form='formatted',access='sequential',     &
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
         write(cargo,'(f15.9)')  Time_Block(i)
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
#ifdef SOLID_BODIES
! Creation of the .pvd file for body particles (Body Transport)
   if (vtkconv) then
      filevtk = "VTKConverter_body-part_"//prefix(1:len_trim(prefix))//".pvd"
      open(unit=unitvtk,file=filevtk,form='formatted',access='sequential',     &
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
         write(cargo,'(f15.9)')  Time_Block(i)
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
      open(unit=unitvtk,file=filevtk,form='formatted',access='sequential',     &
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
         write(cargo,'(f15.9)')  Time_Block(i)
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
#endif
!------------------------
! Deallocations
!------------------------
return
end subroutine Gest_Trans
