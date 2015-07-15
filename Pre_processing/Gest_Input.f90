!cfile gest_input.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name      : gest_input
!
! Last updating : Nov 14, 2012
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Amicarelli/Agate  30Nov11        BSPH approach
! 04  Amicarelli-Agate  13nov12        (AA501b) Body dynamics 
!AA504
! 05  Amicarelli        08Apr14        (v5.04) Several modifications (mainly for granular flows and 3D erosion criterion)
!AA601
! 06  Amicarelli        26Jan15        DBSPH-input (AA601). New DBSPH input management.
!
!************************************************************************************
! Module purpose : Module for input check and management
!
! Calling routine: sphera
!
! Called routines: diagnostic
!                  CompleteBoundaries3D
!                  CreaGrid
!                  DefineBoundaryFaceGeometry3D
!                  DefineBoundarySideGeometry2D
!                  FindBoundaryConvexEdges3D
!                  GeneratePart
!                  Gest_Dealloc
!                  Init_Arrays
!                  ModifyFaces
!                  OrdGrid1
!                  ReadInput
!                  ReadRestartFile
!                  SubCalcPreIdro
!                  Input_Body_Dynamics (AA501b)
!AA601 start
!                  Import_ply_surface_meshes 
!                  DBSPH_IC_surface_elements
!                  DBSPH_find_close_faces
!                  semi_particle_volumes
!AA601 end
!
!************************************************************************************
!
subroutine Gest_Input
!
!.. assign modules
use FILES_ENTITIES
use GLOBAL_MODULE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
integer(4), parameter :: ner0 = 0
!
!AA406!!!test
!real :: dx_w,xmin,xmax,ymin,ymax,zmin,zmax,xmin_obs,xmax_obs,ymin_obs,ymax_obs,zmin_obs,zmax_obs
!integer(4) :: nx,ny,nz,nx_obs,ny_obs,nz_obs   
!
!.. Local Scalars ..
!
!AA406sub
integer(4) :: npi, ier,i,n,isi,nfc,nt !,j
!
!AA504sub
integer(4) :: nrecords,IC_loop
integer(4) :: InputErr
character(len=lencard) :: nomsub = "GEST_INPUT"
character(len=lencard) :: ainp,msg_err
logical :: OnlyTriangle != .FALSE. ! .TRUE. !Opzione trasformazione facce 4 in 2 da 3 lati
!AA504
integer(4) :: machine_Julian_day,machine_hour,machine_minute,machine_second 
!
!.. Local Arrays ..
integer(4), dimension(20)  :: NumberEntities 
!
! External functions and subrotuines
character(10), external :: ltrim
character(80), external :: lcase
!
!AA406 test
!integer(4) :: ParticleCellNumber
!
!.. executable statements
!
  write (nout,'(1x,a)') ">> Input data management starts... "
  write (nscr,'(1x,a)') ">> Input data management starts... "
!
!.. deallocation of arrays
!
  call Gest_Dealloc ( nomsub )
!
!.. initializations
!
  ncord = 0                         ! spatial coordinates
  NumberEntities  = 0
  Domain%istart   = 0               ! initialization for first execution
  Domain%start    = zero            ! initialization for first execution
  Domain%file     = " "             ! initialization for first execution 
  Domain%NormFix  = .FALSE.         !Calcolo SI/NO normali particelle fisse con moto
  Domain%Slip     = .FALSE.
  OnlyTriangle    = .FALSE.
  diffusione      = .FALSE.
  esplosione      = .FALSE.
  erosione        = .FALSE.
  Restart         = .FALSE.
!
!.. time loop parameters are initialized
!
  tempo    = zero
  dt       = zero
  it_start = 0
!
  
!AA504 start
 if (exetype == "linux") then
    call system("date +%j%H%M%S > date_0.txt")
    open (unit_time_elapsed,file='date_0.txt',status="unknown",form="formatted")
    read (unit_time_elapsed,'(i3,i2,i2,i2)') machine_Julian_day,machine_hour,machine_minute,machine_second
    close (unit_time_elapsed)
    Domain%t0 = machine_Julian_day*24*60*60+machine_hour*60*60+machine_minute*60+machine_second
 endif
!AA504 end  
  
!.. allocations of temporary arrays
!
  write (nout,'(1x,a)') ">> Temporary storage allocation in routine "//trim(nomsub)
  allocate ( Vertice(SPACEDIM,1), BoundaryFace(1), Tratto(1), BoundaryVertex(1), stat = ier )
  if (ier /= 0) then
    write (nout,'(1x,a,i2)') "    Arrays VERTICE, BoundaryFace, TRATTO, BOUNDARYVERTEX not allocated. Error code: ",ier
    call diagnostic (arg1=4,arg3=nomsub)
  else
    write (nout,'(1x,a)') "    Arrays VERTICE, BoundaryFace, TRATTO, BOUNDARYVERTEX successfully allocated "
  end if
!
  allocate ( Partz(1), Med(1), Control_Sections(0:1), Control_Points(1), Control_Lines(1), stat = ier )
  if (ier /= 0) then
    write (nout,'(1x,a,i2)') &
    "    Arrays PARTZ, MED, Control_Sections, Control_Points, Control_Lines not allocated. Error code: ",ier
    call diagnostic (arg1=4,arg3=nomsub)
  else
    write (nout,'(1x,a)') "    Arrays PARTZ, MED, Control_Sections, Control_Points, Control_Lines successfully allocated "
  end if
!
!.. input data reading
!.. read for the first time the input file in order to have the parametrs for the arrays dimension
!
  NumVertici   = 1
  NumFacce     = 1
  NumTratti    = 1
  NumBVertices = 1
  NPartZone    = 1
  NMedium      = 1
  NSections    = 1
  NPoints      = 1
  NLines       = 1
!
  call ReadInput( NumberEntities,OnlyTriangle,InputErr,ainp )
!.. an error was detected in the input data deck. Execution fails
!
  msg_err = trim("dimensioning")
  if ( InputErr /= 0 ) then
    InputErr = InputErr + 300
    call diagnostic (arg1=5,arg2=InputErr,arg3=msg_err)
  end if
!
  write (nout,'(1x,a)') " > Data are read from an ASCII input file in the routine ReadInput"
  write (nscr,'(1x,a)') " > Data are read from an ASCII input file in the routine ReadInput"
!
!.. deallocations of temporary arrays
!
  write (nout,'(1x,a)') ">> Deallocation of temporary arrays "
  deallocate ( Vertice, BoundaryFace, Tratto, BoundaryVertex, stat = ier )
  if (ier /= 0) then
    write (nout,'(1x,a,i2)') "    Arrays VERTICE, BoundaryFace, TRATTO, BOUNDARYVERTEX not deallocated. Error code: ",ier
    call diagnostic (arg1=4,arg3=nomsub)
  else
    write (nout,'(1x,a)') "    Arrays VERTICE, BoundaryFace, TRATTO, BOUNDARYVERTEX successfully deallocated "
  end if
!
  deallocate ( Partz , Med , Control_Sections, Control_Points , Control_Lines , stat = ier )
  if (ier /= 0) then
    write (nout,'(1x,a,i2)') &
    "    Arrays PARTZ, MED, Control_Sections, Control_Points, Control_Lines not deallocated. Error code: ",ier
    call diagnostic (arg1=4,arg3=nomsub)
  else
    write (nout,'(1x,a)') "    Arrays PARTZ, MED, Control_Sections, Control_Points, Control_Lines successfully deallocated "
  end if
!
!.. --------------------------
!
!.. --------------------------
!.. a restart procedure has been invoked: restart positioning (step number / step time)
!.. --------------------------
!
  if ( Domain%istart > 0 .or. Domain%start > zero ) then
!
    Restart = .True.
!
!.. open the restart file from which restart data will be restored
!
    open ( unit=nsav, file=trim(Domain%file), form="unformatted", status="old", iostat=ier)
    if ( ier /= 0 ) then
      ainp = Domain%file
      call diagnostic (arg1=5,arg2=201,arg3=trim(ainp))
    else
      write (nout,'(1x,a)') " > Data are read from the restart file"//trim(Domain%file)//" in the routine ReadRestartFile"
      write (nscr,'(1x,a)') " > Data are read from the restart file in the routine ReadRestartFile"
    end if
!
!.. restore the data from the restart file
!
    call ReadRestartFile ( trim("heading"), ier, nrecords)

    msg_err = trim("heading")
    if ( ier /= 0 ) then
      ier = ier + 200
      call diagnostic (arg1=5,arg2=ier,arg3=msg_err)
    end if
!
    NPoints        = NumberEntities(4)
    NPointsl       = NumberEntities(6)
    NPointst       = NumberEntities(4) + NumberEntities(6) + NumberEntities(13)
    NLines         = NumberEntities(5)
    NSections      = NumberEntities(12)
    NPointse       = NumberEntities(13)
!
! .. two fluids case and diffusion are considered
!
!    do i = 1, NMedium
!      if (Med(i)%codif /= zero) diffusione = .TRUE.
!      if (Med(i)%Gamma /= zero) esplosione = .TRUE.
!      if ((index(Med(i)%tipo,"granular") > 0)) then
!        erosione = .TRUE.
!        modelloerosione = Med(i)%modelloerosione
!      end if
!    end do
!
!.. -----------------------------------
!.. no restart. Standard initialization
!.. -----------------------------------
!
  else
!
    ncord          = NumberEntities(1)
    nmedium        = NumberEntities(2)
    NPartZone      = NumberEntities(3)
    NPoints        = NumberEntities(4)
    NPointsl       = NumberEntities(6)
    NPointst       = NumberEntities(4) + NumberEntities(6) + NumberEntities(13)
    NLines         = NumberEntities(5)
    NumVertici     = NumberEntities(7)
    NumTratti      = NumberEntities(8)
    NumBVertices   = NumberEntities(9)
    NumBSides      = NumberEntities(10)
    NumFacce       = NumberEntities(11)
    if ( OnlyTriangle ) NumFacce = NumFacce + NumberEntities(18)
    NSections      = NumberEntities(12)
    NPointse       = NumberEntities(13)
    if ( NumberEntities(19) == 1 )  Domain%Slip     = .TRUE.
    if ( NumberEntities(20) == 1 )  Domain%NormFix = .TRUE.
!
  end if
!
!.. --------------------------
!
!.. allocate the arrays for the calculation depending on the quantities read or restarted
!
  write (nout,'(1x,a)') ">> Final storage allocation in routine "//trim(nomsub)
  allocate ( Vertice       (SPACEDIM,max(1,NumVertici)) , &
             BoundaryFace  (max(1,NumFacce)), &
             Tratto        (max(1,NumTratti)), &
             BoundaryVertex(max(1,NumBVertices)), & 
             BoundarySide  (max(1,NumBSides)), &
             stat = ier )
  if (ier /= 0) then
    write (nout,'(1x,a,i2)') &
    "    Arrays VERTICE, BoundaryFace, TRATTO, BOUNDARYVERTEX, BOUNDARYSIDE not allocated. Error code: ",ier
    call diagnostic (arg1=4,arg3=nomsub)
  else
    write (nout,'(1x,a)') "    Arrays VERTICE, BoundaryFace, TRATTO, BOUNDARYVERTEX, BOUNDARYSIDE successfully allocated "
  end if
!
  allocate ( Partz(NPartZone), &
             Med(NMedium), OpCount(NMedium), SpCount(NMedium), EpCount(NMedium), EpOrdGrid(NMedium), &
             Control_Sections(0:NSections+1), &
             Control_Points(max(1,NPointst)), &
             Control_Lines(max(1,NLines)), &
             Section_Points(1), &
             stat = ier )
  if (ier /= 0) then
    write (nout,'(1x,a,i2)') &
    "    Arrays PARTZ, MED, Control_Sections, Control_Points, Control_Lines, Section_Points not allocated. Error code: ",ier
    call diagnostic (arg1=4,arg3=nomsub)
  else
    write (nout,'(1x,a)') &
    "    Arrays PARTZ, MED, Control_Sections, Control_Points, Control_Lines, Section_Points successfully allocated "
  end if
!
  rewind (ninp)
!
!NumberEntities = 0
!
!.. array initializations
!
  call Init_Arrays
!
!.. ----------------------------
!.. no restart. Data acquisition
!.. ----------------------------
!
  if ( .not. Restart ) then
!
    call ReadInput( NumberEntities,OnlyTriangle,InputErr,ainp )
!
    msg_err = trim("readinput")
    if ( InputErr /= 0 ) then
      InputErr = InputErr + 300
      call diagnostic (arg1=5,arg2=InputErr,arg3=msg_err)
    end if
!
! .. two fluids case and diffusion are considered
!
    do i = 1, NMedium
      if (Med(i)%codif /= zero) diffusione = .TRUE.
      if (Med(i)%Gamma /= zero) esplosione = .TRUE.
      if ((index(Med(i)%tipo,"granular") > 0)) then
        erosione = .TRUE.
        modelloerosione = Med(i)%modelloerosione
      end if
    end do
!
    close (ninp)
!
!.. set the domain parameters
!
    nag    = 0
!
!.. 2D layout
!
    if (ncord == 2)then
!    
!AA406 sub
      Domain%coefke    = 0.682093d0 / squareh  ! 10 / (7 * pigreco) *(3./2.)
!
      Domain%coefkacl = 0.099472d0 / squareh   ! 5 / (16 * pigreco)
! calcolo volume particella
      Domain%PVolume = Domain%dd * Domain%dd
!
!.. 3D layout
!
    else if (ncord == 3)then
      Domain%coefke    = 0.477464d0 / cubich   ! 1 / pigreco
      Domain%coefkacl = 0.074603d0 / cubich    ! 15 / (64 * pigreco)
! calcolo volume particella
      Domain%PVolume = Domain%dd * Domain%dd * Domain%dd
    end if
!
    Control_Sections(NSections+1)%XYZRange(1:3,1) = Domain%coord(1:3,1)
    Control_Sections(NSections+1)%XYZRange(1:3,2) = Domain%coord(1:3,2)
!
!.. an irregular domain is considered (standard or semianalytical)
!
!AA406 sub
    if ( (Domain%tipo == "semi") .or. (Domain%tipo == "bsph") ) then
!
!.. 2D parameter setting
!
      if ( ncord == 2 ) then
!
!.. define the boundary structure
!
        call DefineBoundarySideGeometry2D
!
!.. 3D parameter setting
!
      else if ( ncord == 3 ) then
!
!.. modifies four-sided structures into three-sided structures
!
        if ( OnlyTriangle ) call ModifyFaces ( NumberEntities )
!
        allocate ( BFaceList(NumFacce), stat = ier )
        if (ier /= 0) then
          write (nout,'(1x,a,i2)') "    Array BFACELIST not allocated. Error code: ",ier
          call diagnostic (arg1=4,arg3=nomsub)
        else
          write (nout,'(1x,a)') "    Array BFACELIST successfully allocated "
        end if
!
!.. define the boundary structure
!
        call CompleteBoundaries3D
!
        call DefineBoundaryFaceGeometry3D
!
        
!AA504 start
        allocate (BoundaryConvexEdge(1:Domain%MAXNUMCONVEXEDGES), stat = ier)
        if (ier /= 0) then
           write (nout,'(1x,a,i2)') "   Array BoundaryConvexEdge not allocated. Error code: ",ier
           call diagnostic (arg1=4,arg3=nomsub)
           else
              write (nout,'(1x,a)') "   Array BoundaryConvexEdge successfully allocated "
        end if
!AA504 end

        call FindBoundaryConvexEdges3D
!
      end if
!
!.. creation of the field particles in order to allocate the storage
!
      nagpg = 0

!AA504 start
      if (Granular_flows_options%ID_erosion_criterion > 0) then
          Med(Granular_flows_options%ID_granular)%den0_s = Med(Granular_flows_options%ID_granular)%den0
          Med(Granular_flows_options%ID_granular)%den0 = (1.d0-Med(Granular_flows_options%ID_granular)%gran_vol_frac_max) * &
          Med(Granular_flows_options%ID_main_fluid)%den0 + Med(Granular_flows_options%ID_granular)%gran_vol_frac_max * Med(Granular_flows_options%ID_granular)%den0_s 
          Med(Granular_flows_options%ID_granular)%eps = Med(Granular_flows_options%ID_granular)%eps * &
                                                            Med(Granular_flows_options%ID_granular)%den0/Med(Granular_flows_options%ID_granular)%den0_s
      endif
!AA504 end      
      
!AA504sub
      IC_loop = 1
      call GeneratePart(IC_loop)

      
!AA504
      if (.not.(allocated(pg))) then      

!.. evaluates the number of field particles. Total number of particles is allocated depending on the nag value
!
      if (nag < 100) then
!.. initial domain empty and with source

!AA504 sub
        PARTICLEBUFFER = INIPARTICLEBUFFER * Domain%COEFNMAXPARTI

      else
          
!AA504sub          
        PARTICLEBUFFER = nag * Domain%COEFNMAXPARTI

      end if
!
     
!AA406 sub
      if ( ((Domain%tipo == "semi") .or. (Domain%tipo == "bsph")) ) then
!
        allocate ( pg(PARTICLEBUFFER), stat = ier )
      else
        call diagnostic (arg1=10,arg2=5,arg3=nomsub)
      end if   
      if (ier /= 0) then
        write (nout,'(1x,a,i2)') "    Array PG not allocated. Error code: ",ier
        call diagnostic (arg1=4,arg3=nomsub)
      else
        write (nout,'(1x,a)') "    Array PG successfully allocated "
        Pg(:) = PgZero
      end if
      
!AA504 sub
      endif
      
!
!AA402 start
      if ( Domain%RKscheme > 1 ) then 
        if ( Domain%tipo == "semi" ) then 
          allocate ( ts0_pg(PARTICLEBUFFER), stat = ier )
        else
          call diagnostic (arg1=10,arg2=5,arg3=nomsub)
        end if   
        if (ier /= 0) then
          write (nout,'(1x,a,i2)') "    Array ts0_pg not allocated. Error code: ",ier
          call diagnostic (arg1=4,arg3=nomsub)
        else
          write (nout,'(1x,a)') "    Array ts0_pg successfully allocated "
          ts0_pg(:) = ts_pgZero
        end if
      end if
!AA402 end
!
!.. virtual spatial grid is generated on all the domain
!
      call CreaGrid
!
!.. initial areas are discretized and the particles are created and initialized
!

!AA504sub
      IC_loop = 2
      call GeneratePart(IC_loop)

!AA406test
!
!AA501 sub
!AA601 rm
    else
      call diagnostic (arg1=10,arg2=5,arg3=nomsub)
    end if
!
!.. --------------------------
!.. a restart option is active
!.. --------------------------
!
  else
!
    if (nag < 100) then
!.. initial domain empty and with source
      PARTICLEBUFFER = INIPARTICLEBUFFER * Domain%COEFNMAXPARTI
    else
      PARTICLEBUFFER = nag * Domain%COEFNMAXPARTI
    end if
!
    if ( Domain%tipo == "semi" ) then   
      allocate ( pg(PARTICLEBUFFER), stat = ier )  
    else
      call diagnostic (arg1=10,arg2=5,arg3=nomsub)
    end if   
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "    Array PG not allocated. Error code: ",ier
      call diagnostic (arg1=4,arg3=nomsub)
    else
      write (nout,'(1x,a)') "    Array PG successfully allocated "
      Pg(:) = PgZero
    end if
!
!AA402 start
    if (Domain%RKscheme > 1) then
      if ( Domain%tipo == "semi" ) then   
        allocate ( ts0_pg(PARTICLEBUFFER), stat = ier )  
      else
        call diagnostic (arg1=10,arg2=5,arg3=nomsub)
      end if   
      if (ier /= 0) then
        write (nout,'(1x,a,i2)') "    Array ts0_pg not allocated. Error code: ",ier
        call diagnostic (arg1=4,arg3=nomsub)
      else
        write (nout,'(1x,a)') "    Array ts0_pg successfully allocated "
        ts0_pg(:) = ts_pgZero
      end if
    end if
!AA402 end
!
    allocate ( BFaceList(NumFacce), stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "    Array BFACELIST not allocated. Error code: ",ier
      call diagnostic (arg1=4,arg3=nomsub)
    else
      write (nout,'(1x,a)') "    Array BFACELIST successfully allocated "
    end if
!
    call ReadRestartFile ( trim("reading"), ier, nrecords)
    msg_err = trim("reading")
    if ( ier /= 0 ) then
      ier = ier + 200
      call diagnostic (arg1=5,arg2=ier,arg3=msg_err)
    end if
!
    close (nsav)
!
    call ReadInput( NumberEntities,OnlyTriangle,InputErr,ainp )
!
    msg_err = trim("restart reading?")
    if ( InputErr /= 0 ) then
      InputErr = InputErr + 300
      call diagnostic (arg1=5,arg2=InputErr,arg3=msg_err)
    end if
!
! .. two fluids case and diffusion are considered
!
!$$$$$$$$$$$$$$$
!esplosione = .FALSE.
!do i=1,nag
!  pg(i)%pres = zero
!  pg(i)%IntEn = zero
!  pg(i)%EnVar = zero
!  pg(i)%dEdT = zero
!end do
!$$$$$$$$$$$$$$$
    do i = 1, NMedium
      if (Med(i)%codif /= zero) diffusione = .TRUE.
!      if (Med(i)%Gamma /= zero) esplosione = .TRUE.
      if ((index(Med(i)%tipo,"granular") > 0)) then
        erosione = .TRUE.
        modelloerosione = Med(i)%modelloerosione
      end if
    end do
!
!.. save current time for result_converter
!
    val_time = tempo
!
    close (ninp)
!
  end if
!
!.. --------------------------
!
!
!.. final prints
!
!  write (nout,*) 
!  write (nout,*) "Number of particles          NAG= ",nag
!
! scrittura su file di ouput delle particelle di campo
  if ( Domain%ioutopt < 0 ) then
    write (nout,*) 
    write (nout,*) "======== PARTICLES COORDINATES =========="
    do n = 1, NPartZone
      write (nout,*) 
      write (nout,"(a,i5,3x,a)") "ZONE",n,Partz(n)%label
      do npi = Partz(n)%limit(1), Partz(n)%limit(2)
        write (nout,"(i10,4f14.5)") npi, pg(npi)%coord, pg(npi)%tstop  !, pg(npi)%vel
      end do
    end do
  end if
!

!AA501b start
! Management of body dynamics input
 if (n_bodies > 0) then
    call Input_Body_Dynamics
 endif
!AA501b end

!allocazione della memoria per i vettori di ordinamento delle particelle
!
!AA406 sub
  if ( (Domain%tipo == "semi") .or. (Domain%tipo == "bsph") ) then
!
    allocate ( NPartOrd(PARTICLEBUFFER),Icont(grid%nmax+1), stat = ier ) 
  else
    call diagnostic (arg1=10,arg2=5,arg3=nomsub)
  end if    
!
  if (ier /= 0) then
    write (nout,'(1x,a,i2)') "    Array NPARTORD,ICONT not allocated. Error code: ",ier
    call diagnostic (arg1=4,arg3=nomsub)
  else
    write (nout,'(1x,a)') "    Array NPARTORD,ICONT successfully allocated "
    NPartOrd(:) = 0
    Icont(:) = 0
  end if
!AA601 rm
!AA501b start
  if (n_bodies > 0.) then
     allocate (NPartOrd_bp(n_body_part),Icont_bp(grid%nmax+1),stat = ier) 
     if (ier /= 0) then
       write (nout,'(1x,a,i2)') "    Arrays NPARTORD_bp and ICONT_bp not allocated. Error code: ",ier
       call diagnostic (arg1=4,arg3=nomsub)
       else
          write (nout,'(1x,a)') "    Arrays NPARTORD_bp and ICONT_bp successfully allocated "
          NPartOrd_bp(:) = 0
          Icont_bp(:) = 0
     end if
  endif
!AA501b end
!AA601 start
  if (Domain%tipo == "bsph") then
     call Import_ply_surface_meshes
     call DBSPH_IC_surface_elements
     if (.not.allocated(NPartOrd_w)) then
        allocate(NPartOrd_w(DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet),Icont_w(grid%nmax+1),stat=ier) 
        if (ier/=0) then
           write(nout,*) 'Error! Allocation of NPartOrd_w or Icont_w Gest_Input failed.'           
           call diagnostic (arg1=5,arg2=340)
           else
              write (nout,'(1x,a)') "Arrays NPARTORD_w and ICONT_w successfully allocated."
              NPartOrd_w(:) = 0
              Icont_w(:) = 0
        end if
     endif
  endif
!AA601 end
  call OrdGrid1 ( nout )
!AA601 start
  if (Domain%tipo == "bsph") then
     call DBSPH_find_close_faces 
     call semi_particle_volumes
  endif
!AA601 end
!.. initialize the pressure and density fields
!
  if (.not. Restart) call SubCalcPreIdro
!
!.. evaluates and stores the opened boundary sides for 2D calculation
!
  if ( ncord == 2 ) then     
!
!.. searches for the opened boundary sides
!
    NumOpenSides = 0
!
    do isi = 1, NumBSides
      if (BoundarySide(isi)%tipo == "leve" .OR. BoundarySide(isi)%tipo == "velo" .OR. BoundarySide(isi)%tipo == "flow" .OR. &
          BoundarySide(isi)%tipo == "crit" .OR. BoundarySide(isi)%tipo == "open") then  
        NumOpenSides = NumOpenSides + 1
        if (NumOpenSides > MAXOPENSIDES) call diagnostic (arg1=10,arg2=6,arg3=nomsub)
        OpenSide(NumOpenSides) = isi
      end if
    end do
!
!.. evaluates and stores the opened boundary sides for 3D calculation
!
  else
!
!.. searches for the boundary faces having opened conditions
!
    NumOpenFaces = 0
!
    do nfc = 1, NumFacce
      nt = BoundaryFace(nfc)%stretch
      if (Tratto(nt)%tipo == "leve" .OR. Tratto(nt)%tipo == "velo" .OR. Tratto(nt)%tipo == "flow" .OR. &
          Tratto(nt)%tipo == "crit" .OR. Tratto(nt)%tipo == "open") then
        NumOpenFaces = NumOpenFaces + 1
        if (NumOpenFaces > MAXOPENFACES ) call diagnostic (arg1=10,arg2=7,arg3=nomsub)
        OpenFace(NumOpenFaces) = nfc
      end if
    end do
!
!.. searches for the maximum zone index ("maxzone") among them identifying volumes(?) of particles in the domain
!.. the sources (water material) must have currently the maximum zone index because the particle array is partitioned 
!.. by zones, i.e. the particles belonging to the first zone is loaded at first, after that those of the second one and so on.
!.. Therefore, the source particles can be added only in the last zone (last array section) without restructuring all the 
!.. zone pointers in the array (Part(zone)%limit(1) is the first particle of the zone-th, Part(zone)%limit(2) is the last
!.. one)
!
!    call SearchforParticleZone_3D(maxzone)      
!
  end if
!
  OpCount = 0
  EpCount = 0
  EpOrdGrid = 0

  return
  end subroutine Gest_Input
!---split

