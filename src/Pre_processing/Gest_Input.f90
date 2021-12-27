!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2021 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: Gest_Input                   
! Description: Input check and management.                 
!-------------------------------------------------------------------------------
subroutine Gest_Input
!------------------------
! Modules
!------------------------
use I_O_file_module
use Static_allocation_module
use Dynamic_allocation_module
use Hybrid_allocation_module
use Memory_I_O_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4), parameter :: ner0 = 0
#ifdef SPACE_3D
! Flag for the transformation of a quadrilateral face in two triangular faces 
logical :: OnlyTriangle
integer(4) :: IC_loop,nt,nfc
#elif defined SPACE_2D
integer(4) :: isi
#endif
integer(4) :: npi,ier,i,n,nrecords,InputErr,alloc_stat
integer(4) :: machine_Julian_day,machine_hour,machine_minute,machine_second
integer(4),dimension(20) :: NumberEntities 
double precision :: eps_f
character(len=lencard) :: nomsub = "GEST_INPUT"
character(len=lencard) :: ainp,msg_err
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
write(ulog,'(1x,a)') ">> Input data management starts... "
! Array deallocation
call deallocation_sequence
! This line seems redundant, but it is useful to distinguish between the first 
! and the second execution of the reading program units. "dx" is the first 
! variable to be read from the main input file (beyond the strings).
Domain%dx = zero
NumberEntities = 0
Domain%istart = 0    
Domain%start = zero            
Domain%file = " " 
Domain%NormFix = .false.         
Domain%Slip = .false.
#ifdef SPACE_3D
   OnlyTriangle = .false.
#endif
restart = .false.
simulation_time = zero
dt = zero
it_start = 0
!------------------------
! Statements
!------------------------
if (exetype=="linux") then
   call system("date +%j%H%M%S>date_0.txt")
   open(unit_time_elapsed,file='date_0.txt',status="unknown",form="formatted")
   read(unit_time_elapsed,'(i3,i2,i2,i2)') machine_Julian_day,machine_hour,    &
      machine_minute,machine_second
   close(unit_time_elapsed)
   Domain%t0 = machine_Julian_day * 24 * 60 * 60 + machine_hour * 60 * 60 +    &
               machine_minute * 60 + machine_second
endif
! Allocations of temporary arrays
write(ulog,'(1x,a)') ">> Temporary storage allocation in routine "//trim(nomsub)
#ifdef SPACE_3D
allocate(Vertice(SPACEDIM,1),BoundaryFace(1),Tratto(1),BoundaryVertex(1),      &
   stat=ier)
#elif defined SPACE_2D
allocate(Vertice(SPACEDIM,1),Tratto(1),BoundaryVertex(1),stat=ier)
#endif
   if (ier/=0) then
      write(ulog,'(1x,a,i2)')                                                  &
"    Arrays VERTICE, possibly BoundaryFace, TRATTO, BOUNDARYVERTEX not allocated. Error code: "&
         ,ier
   call diagnostic(arg1=4,arg3=nomsub)
   else
      write(ulog,'(1x,a)')                                                     &
"    Arrays VERTICE, possibly BoundaryFace, TRATTO, BOUNDARYVERTEX successfully allocated "
endif
allocate(Partz(1),Med(1),Control_Points(1),Control_Lines(1),stat=ier)
   if (ier/=0) then
      write(ulog,'(1x,a,i2)')                                                  &
"    Arrays PARTZ, MED, Control_Points, Control_Lines not allocated. Error code: "&
         ,ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)')                                                  &
"    Arrays PARTZ, MED, Control_Points, Control_Lines successfully allocated "
endif
! Input data reading
! To read for the first time the input file to get the parameters for 
! array sizing
NumVertici = 1
#ifdef SPACE_3D
   NumFacce = 1
#endif
NumTratti = 1
NumBVertices = 1
NPartZone = 1
NMedium = 1
npoints = 1
NLines = 1
input_second_read = .false.
#ifdef SPACE_3D
call ReadInput(NumberEntities,OnlyTriangle,InputErr,ainp)
#elif defined SPACE_2D
call ReadInput(NumberEntities,InputErr,ainp)
#endif
input_second_read = .true.
! An error was detected in the input data. Execution fails.
msg_err = trim("dimensioning")
if (InputErr/=0) then
   InputErr = InputErr + 300
   call diagnostic(arg1=5,arg2=InputErr,arg3=msg_err)
endif
write(ulog,'(1x,a)')                                                           &
   ">Data are read from an ASCII input file in the routine ReadInput"
! Deallocations of temporary arrays
write(ulog,'(1x,a)') ">> Deallocation of temporary arrays "
#ifdef SPACE_3D
deallocate(Vertice,BoundaryFace,Tratto,BoundaryVertex,stat=ier)
#elif defined SPACE_2D
deallocate(Vertice,Tratto,BoundaryVertex,stat=ier)
#endif
if (ier/=0) then
   write(ulog,'(1x,a,i2)')                                                     &
"    Arrays VERTICE, possibly BoundaryFace, TRATTO, BOUNDARYVERTEX not deallocated. Error code: "&
      ,ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)')                                                  &
"    Arrays VERTICE, possibly BoundaryFace, TRATTO, BOUNDARYVERTEX successfully deallocated "
endif
deallocate(Partz,Med,Control_Points,Control_Lines,stat=ier)
if (ier/=0) then
   write(ulog,'(1x,a,i2)')                                                     &
"    Arrays PARTZ, MED, Control_Points, Control_Lines not deallocated. Error code: "&
      ,ier
   call diagnostic(arg1=4,arg3=nomsub)
   else
      write(ulog,'(1x,a)')                                                     &
"    Arrays PARTZ, MED, Control_Points, Control_Lines successfully deallocated "
endif
! A restart procedure has been invoked: restart positioning 
! (step number / step time)
if ((Domain%istart>0).or.(Domain%start>zero)) then
   restart = .true.
! To open the restart file from which restart data will be restored
   open(unit=nsav,file=trim(Domain%file),form="unformatted",status="old",      &
      iostat=ier)
   if (ier/=0) then
      ainp = Domain%file
      call diagnostic(arg1=5,arg2=201,arg3=trim(ainp))
      else
         write(ulog,'(1x,a)')                                                  &
">Data are read from the restart file "//trim(Domain%file)//" in the routine ReadRestartFile"
   endif
! To restore data from the restart file
! During the first reading of the restart file, only few parameters are read
   call ReadRestartFile(trim("heading"),ier,nrecords)
   msg_err = trim("heading")
   if (ier/=0) then
      ier = ier + 200
      call diagnostic(arg1=5,arg2=ier,arg3=msg_err)
   endif
   npoints = NumberEntities(4)
   NPointsl = NumberEntities(6)
   npointst = NumberEntities(4) + NumberEntities(6) + NumberEntities(13)
   NLines = NumberEntities(5)
   NPointse = NumberEntities(13)
   else
! No restart: standard initialization
      NMedium = NumberEntities(2)
      NPartZone = NumberEntities(3)
      npoints = NumberEntities(4)
      NPointsl = NumberEntities(6)
      npointst = NumberEntities(4) + NumberEntities(6) + NumberEntities(13)
      NLines = NumberEntities(5)
      NumVertici = NumberEntities(7)
      NumTratti = NumberEntities(8)
      NumBVertices = NumberEntities(9)
#ifdef SPACE_3D
         NumFacce = NumberEntities(11)
         if (OnlyTriangle) NumFacce = NumFacce + NumberEntities(18)
#elif defined SPACE_2D
            NumBSides = NumberEntities(10)
#endif
      NPointse = NumberEntities(13)
      if (NumberEntities(19)==1) Domain%Slip = .true.
      if (NumberEntities(20)==1) Domain%NormFix = .true.
endif
! Array allocations 
write(ulog,'(1x,a)') ">> Final storage allocation in routine "//trim(nomsub)
#ifdef SPACE_3D
allocate(Vertice(SPACEDIM,max(1,NumVertici)),BoundaryFace(max(1,NumFacce)),    &
   Tratto(max(1,NumTratti)),BoundaryVertex(max(1,NumBVertices)),stat=ier)
#elif defined SPACE_2D
allocate(Vertice(SPACEDIM,max(1,NumVertici)),Tratto(max(1,NumTratti)),         &
   BoundaryVertex(max(1,NumBVertices)),BoundarySide(max(1,NumBSides)),stat=ier)
#endif
if (ier/=0) then
   write(ulog,'(1x,a,i2)')                                                     &
"    Arrays VERTICE, BoundaryFace/BoundarySide, TRATTO, BOUNDARYVERTEX not allocated. Error code: "&
      ,ier
   call diagnostic(arg1=4,arg3=nomsub)
   else
      write(ulog,'(1x,a)')                                                     &
"    Arrays VERTICE, BoundaryFace/BoundarySide, TRATTO, BOUNDARYVERTEX successfully allocated "
endif
allocate(Partz(NPartZone),Med(NMedium),OpCount(NMedium),SpCount(NMedium),      &
   EpCount(NMedium),EpOrdGrid(NMedium),Control_Points(max(1,npointst)),        &
   Control_Lines(max(1,NLines)),stat=ier)
if (ier/=0) then
   write(ulog,'(1x,a,i2)')                                                     &
"    Arrays PARTZ, MED, Control_Points, Control_Lines, not allocated. Error code: "&
      ,ier
   call diagnostic(arg1=4,arg3=nomsub)
   else
      write(ulog,'(1x,a)')                                                     &
"    Arrays PARTZ, MED, Control_Points, Control_Lines, successfully allocated "
endif
rewind(ninp)
! Array initializations
call Init_Arrays
! No restart
if (.not.restart) then
#ifdef SPACE_3D
   call ReadInput(NumberEntities,OnlyTriangle,InputErr,ainp)
#elif defined SPACE_2D
   call ReadInput(NumberEntities,InputErr,ainp)
#endif
   msg_err = trim("readinput")
   if (InputErr/=0) then
      InputErr = InputErr + 300
      call diagnostic(arg1=5,arg2=InputErr,arg3=msg_err)
   endif
   close(ninp)
   nag = 0
#ifdef SPACE_2D
! 10 / (7 * pigreco) *(3./2.) /(h**2)
      Domain%coefke = 0.682093d0 / (Domain%h ** ncord)
! 5 / (16 * pigreco)/(h**2)
      Domain%coefkacl = 0.099472d0 / (Domain%h ** ncord)
#elif defined SPACE_3D
! 1 / pigreco/(h**3)
         Domain%coefke = 0.477464d0 / (Domain%h ** ncord) 
! 15 / (64 * pigreco)/(h**3)
         Domain%coefkacl = 0.074603d0 / (Domain%h ** ncord)
#endif
! Particle volume
   Domain%PVolume = Domain%dx ** ncord
! An irregular domain is considered 
   if ((Domain%tipo=="semi").or.(Domain%tipo=="bsph")) then
#ifdef SPACE_2D
         call DefineBoundarySideGeometry2D
#elif defined SPACE_3D
! To replace 4-sided geometries with 3-sided geometries
            if (OnlyTriangle) call ModifyFaces(NumberEntities)
            allocate(BFaceList(NumFacce),stat=ier) 
            if (ier/=0) then
               write(ulog,'(1x,a,i2)')                                         &
                  "    Array BFACELIST not allocated. Error code: ",ier
               call diagnostic(arg1=4,arg3=nomsub)
               else
                  write(ulog,'(1x,a)')                                         &
                     "    Array BFACELIST successfully allocated "
            endif
            call CompleteBoundaries3D
            call DefineBoundaryFaceGeometry3D
            allocate(BoundaryConvexEdge(1:input_any_t%MAXNUMCONVEXEDGES),      &
               stat=ier)
            if (ier/=0) then
               write(ulog,'(1x,a,i2)')                                         &
                  "   Array BoundaryConvexEdge not allocated. Error code: ",ier
               call diagnostic(arg1=4,arg3=nomsub)
               else
                  write(ulog,'(1x,a)')                                         &
                     "   Array BoundaryConvexEdge successfully allocated "
            endif
            call FindBoundaryConvexEdges3D
#endif
      nagpg = 0
      if (Granular_flows_options%KTGF_config>0) then
         do i=1,NMedium
            if (index(Med(i)%tipo,"granular")>0) then
               Med(i)%den0_s = Med(i)%den0
               if (Med(i)%saturated_medium_flag.eqv..true.) then
                  eps_f = 1.d0 - Med(i)%gran_vol_frac_max
                  else
                     eps_f = 0.d0
               endif
               Med(i)%den0 = eps_f *                                           &
                             Med(Granular_flows_options%ID_main_fluid)%den0 +  &
                             Med(i)%gran_vol_frac_max * Med(i)%den0_s 
               Med(i)%eps = Med(i)%eps * Med(i)%den0 / Med(i)%den0_s
            endif
         enddo
      endif
#ifdef SPACE_3D
      IC_loop = 1
      call GeneratePart(IC_loop)
#elif defined SPACE_2D
      call GeneratePart
#endif
      if (.not.(allocated(pg))) then      
! To assess the number of particles. Total number of particles is 
! allocated depending on the value "nag". 
         if (nag<100) then
! Initial domain empty (inlet section)
            PARTICLEBUFFER = int(INIPARTICLEBUFFER * Domain%COEFNMAXPARTI) + 1
            else
               PARTICLEBUFFER = int(nag * Domain%COEFNMAXPARTI) + 1
         endif
         if (((Domain%tipo=="semi").or.(Domain%tipo=="bsph"))) then
            allocate(pg(PARTICLEBUFFER),stat=ier)
            else
               call diagnostic(arg1=10,arg2=5,arg3=nomsub)
         endif   
         if (ier/=0) then
            write(ulog,'(1x,a,i2)') "    Array PG not allocated. Error code: ",&
               ier
            call diagnostic(arg1=4,arg3=nomsub)
            else
               write(ulog,'(1x,a)') "    Array PG successfully allocated "
               Pg(:) = PgZero
         endif
      endif
      if (Domain%RKscheme>1) then 
         if (Domain%tipo=="semi") then 
            allocate(ts0_pg(PARTICLEBUFFER),stat=ier)
            else
               call diagnostic(arg1=10,arg2=5,arg3=nomsub)
         endif   
         if (ier/=0) then
            write(ulog,'(1x,a,i2)')                                            &
               "    Array ts0_pg not allocated. Error code: ",ier
            call diagnostic(arg1=4,arg3=nomsub)
            else
               write(ulog,'(1x,a)') "    Array ts0_pg successfully allocated "
               ts0_pg(:) = ts_pgZero
         endif
      endif
! The background positioning grid is generated
      call CreaGrid
#ifdef SPACE_3D
! z0 (CLC class)
      if (CLC_flag.eqv..true.) then
         call CLC_pre_processing
         call z0_CLC
      endif
#endif
! Particles are created and initialized
#ifdef SPACE_3D
      IC_loop = 2
      call GeneratePart(IC_loop)
#elif defined SPACE_2D
      call GeneratePart
#endif
      else
         call diagnostic(arg1=10,arg2=5,arg3=nomsub)
   endif
   else
! A restart option is active
#ifdef SPACE_3D
      call ReadInput(NumberEntities,OnlyTriangle,InputErr,ainp)
#elif defined SPACE_2D
      call ReadInput(NumberEntities,InputErr,ainp)
#endif
      msg_err = trim("restart reading?")
      if (InputErr/=0) then
         InputErr = InputErr + 300
         call diagnostic(arg1=5,arg2=InputErr,arg3=msg_err)
      endif
      if (nag<100) then
! Initial domain is empty (inlet section)
         PARTICLEBUFFER = int(INIPARTICLEBUFFER * Domain%COEFNMAXPARTI) + 1
         else
            PARTICLEBUFFER = int(nag * Domain%COEFNMAXPARTI) + 1
      endif
      if ((Domain%tipo=="semi").or.(Domain%tipo=="bsph")) then
         allocate(pg(PARTICLEBUFFER),stat=ier)  
         else
            call diagnostic(arg1=10,arg2=5,arg3=nomsub)
      endif   
      if (ier/=0) then
         write(ulog,'(1x,a,i2)') "    Array PG not allocated. Error code: ",ier
         call diagnostic(arg1=4,arg3=nomsub)
         else
            write(ulog,'(1x,a)') "    Array PG successfully allocated "
            Pg(:) = PgZero
      endif
      if (Domain%tipo=="bsph") then
! DB-SPH pre-processing
         call Import_ply_surface_meshes
         call DBSPH_IC_surface_elements
         if (.not.allocated(NPartOrd_w)) then
            allocate(NPartOrd_w(DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet),       &
               Icont_w(grid%nmax+1),stat=ier)
            if (ier/=0) then
               write(ulog,*)                                                   &
                  'Error! Allocation of NPartOrd_w or Icont_w Gest_Input failed'
               call diagnostic(arg1=5,arg2=340)
               else
                  write(ulog,'(1x,a)')                                         &
                     "Arrays NPARTORD_w and ICONT_w successfully allocated."
                  NPartOrd_w(:) = 0
                  Icont_w(:) = 0
            endif
         endif
! Allocation of the array of the DB-SPH wall elements and semi-particles
         if (.not.allocated(pg_w)) then
            allocate(pg_w(DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet),STAT=        &
               alloc_stat)
            if (alloc_stat/=0) then
               write(ulog,*)                                                   &
                  'Allocation of pg_w in Gest_Input failed;',                  &
                  ' the program terminates here.'
               stop ! Stop the main program
               else
                  write(ulog,*)                                                &
                     "Allocation of pg_w in Gest_Input successfully completed."
            endif
         endif
      endif
! Allocation of the array of the bodies
      if (n_bodies>0) then
         if (.not.allocated(body_arr)) then
            allocate(body_arr(n_bodies),STAT=alloc_stat)
            if (alloc_stat/=0) then
               write(ulog,*)                                                   &
                  'Allocation of body_arr in Gest_Input failed;',              &
                  ' the program terminates here.'
               stop ! Stop the main program
               else
                  write(ulog,*)                                                &
                     "Allocation of body_arr in Gest_Input well completed. "
            endif
         endif
! Management of body dynamics input
         call Input_Body_Dynamics
! Allocation of the array of the body particles
         if (.not.allocated(bp_arr)) then
            allocate(bp_arr(n_body_part),STAT=alloc_stat)
            if (alloc_stat/=0) then
               write(ulog,*)                                                   &
                  'Allocation of bp_arr in Gest_Input failed;',                &
                  ' the program terminates here.'
               stop ! Stop the main program
               else
                  write(ulog,*)                                                &
                     "Allocation of bp_arr in Gest_Input successfully completed"
            endif
         endif
! Allocation of the array of the surface body particles
         if (.not.allocated(surf_body_part)) then
            allocate(surf_body_part(n_surf_body_part),STAT=alloc_stat)
            if (alloc_stat/=0) then
               write(ulog,*)                                                   &
                  'Allocation of surf_body_part in Gest_Input failed;',        &
                  ' the program terminates here.'
               stop ! Stop the main program
               else
                  write(ulog,*)                                                &
                     'Allocation of surf_body_part in Gest_Input',             &
                     ' successfully completed.'
            endif
         endif
      endif
      if (Domain%RKscheme>1) then
         if (Domain%tipo=="semi") then
           allocate(ts0_pg(PARTICLEBUFFER),stat=ier)  
           else
              call diagnostic(arg1=10,arg2=5,arg3=nomsub)
         endif
         if (ier/=0) then
            write(ulog,'(1x,a,i2)')                                            &
               "    Array ts0_pg not allocated. Error code: ",ier
            call diagnostic(arg1=4,arg3=nomsub)
            else
               write(ulog,'(1x,a)') "    Array ts0_pg successfully allocated "
               ts0_pg(:) = ts_pgZero
         endif
      endif
#ifdef SPACE_3D
      allocate(BFaceList(NumFacce),stat=ier)
      if (ier/=0) then
         write(ulog,'(1x,a,i2)')                                               &
            "    Array BFACELIST not allocated. Error code: ",ier
         call diagnostic(arg1=4,arg3=nomsub)
         else
            write(ulog,'(1x,a)') "    Array BFACELIST successfully allocated "
      endif
#endif
      call ReadRestartFile(trim("reading"),ier,nrecords)
      msg_err = trim("reading")
      if (ier/=0) then
         ier = ier + 200
         call diagnostic(arg1=5,arg2=ier,arg3=msg_err)
      endif
      close(nsav)
! To save current time for "result_converter"
      val_time = simulation_time
      close(ninp)
endif
! Writing on the log file 
if (Domain%ioutopt<0) then
   write(ulog,*) 
   write(ulog,*) "======== PARTICLES COORDINATES =========="
   do n=1,NPartZone
      write(ulog,*) 
      write(ulog,"(a,i5,3x,a)") "ZONE",n,Partz(n)%label
      do npi=Partz(n)%limit(1),Partz(n)%limit(2)
         write(ulog,"(i10,4f14.5)") npi,pg(npi)%coord,pg(npi)%tstop  
      enddo
   enddo
endif
! Management of body dynamics input
if ((n_bodies>0).and.(.not.restart)) call Input_Body_Dynamics
! Memory allocation for the particle ordering arrays
if ((Domain%tipo=="semi").or.(Domain%tipo=="bsph")) then
   allocate(NPartOrd(PARTICLEBUFFER),Icont(grid%nmax+1),stat=ier) 
   else
      call diagnostic(arg1=10,arg2=5,arg3=nomsub)
endif    
if (ier/=0) then
   write(ulog,'(1x,a,i2)')                                                     &
      "    Array NPARTORD,ICONT not allocated. Error code: ",ier
   call diagnostic(arg1=4,arg3=nomsub)
   else
      write(ulog,'(1x,a)') "    Array NPARTORD,ICONT successfully allocated "
      NPartOrd(:) = 0
      Icont(:) = 0
endif
if (n_bodies>0) then
   allocate(NPartOrd_bp(n_body_part),Icont_bp(grid%nmax+1),stat=ier) 
   if (ier/=0) then
      write(ulog,'(1x,a,i2)')                                                  &
         "    Arrays NPARTORD_bp and ICONT_bp not allocated. Error code: ",ier
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)')                                                  &
            "    Arrays NPARTORD_bp and ICONT_bp successfully allocated "
         NPartOrd_bp(:) = 0
         Icont_bp(:) = 0
   endif
endif
if ((Domain%tipo=="bsph").and.(.not.restart)) then
   call Import_ply_surface_meshes
   call DBSPH_IC_surface_elements
   if (.not.allocated(NPartOrd_w)) then
      allocate(NPartOrd_w(DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet),             &
         Icont_w(grid%nmax+1),stat=ier) 
      if (ier/=0) then
         write(ulog,*)                                                         &
            'Error! Allocation of NPartOrd_w or Icont_w Gest_Input failed.'           
         call diagnostic(arg1=5,arg2=340)
         else
            write(ulog,'(1x,a)')                                               &
               "Arrays NPARTORD_w and ICONT_w successfully allocated."
            NPartOrd_w(:) = 0
            Icont_w(:) = 0
      endif
   endif
endif
call OrdGrid1
if ((Domain%tipo=="bsph").and.(.not.restart)) then
   call DBSPH_find_close_faces 
   call semi_particle_volumes
endif
! To initialize pressure and density fields
if (.not.restart) call SubCalcPreIdro
! To assess and save the sides with condition "open" (2D) 
#ifdef SPACE_2D     
! Searching for the sides with condition "open" 
   NumOpenSides = 0
   do isi=1,NumBSides
      if ((BoundarySide(isi)%tipo=="leve").or.(BoundarySide(isi)%tipo=="velo") &
         .or.(BoundarySide(isi)%tipo=="flow").or.                              &
         (BoundarySide(isi)%tipo=="crit").or.                                  &
         (BoundarySide(isi)%tipo=="open")) then  
         NumOpenSides = NumOpenSides + 1
         if (NumOpenSides>MAXOPENSIDES) then
            call diagnostic(arg1=10,arg2=6,arg3=nomsub)
         endif
         OpenSide(NumOpenSides) = isi
      endif
   enddo
! To assess and save the faces with condition "open" (3D)
#elif defined SPACE_3D
! Searching for the faces with condition "open" 
      NumOpenFaces = 0
      do nfc=1,NumFacce
         nt = BoundaryFace(nfc)%stretch
         if ((Tratto(nt)%tipo=="leve").or.(Tratto(nt)%tipo=="velo").or.        &
            (Tratto(nt)%tipo=="flow").or.(Tratto(nt)%tipo=="crit").or.         &
            (Tratto(nt)%tipo=="open")) then
            NumOpenFaces = NumOpenFaces + 1
            if (NumOpenFaces>MAXOPENFACES) call                                &
               diagnostic(arg1=10,arg2=7,arg3=nomsub)
            OpenFace(NumOpenFaces) = nfc
         endif
      enddo
#endif
OpCount = 0
EpCount = 0
EpOrdGrid = 0
!------------------------
! Deallocations
!------------------------
return
end subroutine Gest_Input
