!cfile ReadCheck.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
logical function ReadCheck(IoErr,Ier,Nrighe,ainp,listadati,ninp,nout)
!
!.. implicit declarations
implicit none
!
!.. dummy arguments
integer(4)   :: IoErr,Ier,Nrighe,ninp,nout
character(*) :: ainp, listadati
!
!.. local scalars
!integer(4)   :: MyResult
!character(5) :: Txt
!
  if ( IoErr == 0 ) then
!
    Ier = 0
    ReadCheck = .TRUE.
!
  else
!
    Ier = 4
    ReadCheck = .FALSE.
    write(nout,"(1x,a)")    ">>>>>>>>>>>>>> Warning:"
    write(nout,"(1x,a)")  
    write(nout,"(1x,a,i5)") "Error reading unit:  ",ninp
    write(nout,"(1x,a,a)")  "Expected data:       ",listadati
    write(nout,"(1x,a,i5,a)")  "Last input line read:"//trim(ainp)//"(line number:",Nrighe,")"
!
  end if
!
return
end function ReadCheck
!---split

!cfile readinput.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : readinput
!
! Last updating : May 02, 2014
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Amicarelli/Agate  13nov12        (AA501b) Body dynamics
! 04  Amicarelli/Agate  18apr12        add flag for wall elements reordering and maximum number of neighbours read in input
! 06  Agate             02/05/14       Add more modules check in license
!
!************************************************************************************
! Module purpose : Module for input data reading from ascii input
!
! Calling routine: gest_input
!
! Called routines: readriga
!                  readinputtitle
!                  ReadInputRestart
!                  ReadInputDomain
!                  ReadInputVertices
!                  ReadInputLines
!                  ReadInputFaces
!                  ReadInputExternalFile
!                  ReadInputMedium
!                  ReadInputBoundaries
!                  ReadInputRunParameters
!                  ReadInputGeneralPhysical
!                  ReadInputOutputRegulation
!                  ReadInputControlPoints
!                  ReadInputControlLines
!                  ReadInputControlSections
!                  ReadInputDrawOptions
!                  ReadBodyDynamics
!
!************************************************************************************
!
subroutine ReadInput ( NumberEntities,OnlyTriangle,ier,ainp )

use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE

implicit none

integer(4)    :: ier
character(80) :: ainp
logical       :: OnlyTriangle

integer(4),dimension(20)    :: NumberEntities

integer(4)    :: ioerr, nrighe
character( 1) :: comment = "!"
integer(4)    :: ioutpo2

integer(4)       :: iplot_fr, imemo_fr, irest_fr, icpoi_fr, ipllb_fr, ipllb_md
double precision ::  plot_fr,  memo_fr,  rest_fr,  cpoi_fr,  pllb_fr
integer(4)       :: ioutopt
logical          :: restartOK

character(80), external :: lcase, GetToken
logical,       external :: ReadCheck
!
!.. executable statements
!
!.. parameter initializations
!
  iplot_fr = 0
  plot_fr  = zero
  imemo_fr = 0
  memo_fr  = zero
  irest_fr = 0
  rest_fr  = zero
  icpoi_fr = 0
  cpoi_fr  = zero
  ipllb_fr = 0
  ipllb_md = 0
  pllb_fr  = zero
  ioutopt  = 0
  ioutpo2  = 3
  restartOK= .FALSE.
  ier      = 0
  vtkconv  = .FALSE.
  pesodt   = zero
!AA601
  dt_alfa_Mon = .false.
!
!.. reading loop
!
  ioerr  = 0
  nrighe = 0

  Npoints = NumberEntities(4)
  NumberEntities = 0
!
!.. Reads and checks matching between Program and Input Version
!
  current_version = .TRUE.
  call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
  if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"INPUT VERSION",ninp,nout) ) return
  if ( trim(ainp) /= trim(version) ) then
    ier = 2
    current_version = .FALSE.
    return
  end if
!
!.. Input Sections Loop
!
  SECTION_LOOP:  do while ( ioerr == 0 )

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
!
!.. if EOF has been found, exit, otherwise check the error code
!
    if ( ioerr == -1 ) cycle SECTION_LOOP

    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"INPUT DATA",ninp,nout) ) return

    if ( ncord > 0 ) write(nout,"(//,1x,a,/)") lcase(ainp) 

    select case ( TRIM(lcase(trim(ainp))) )
       

    case ( "##### title #####" )

       call ReadInputTitle ( ainp,comment,nrighe,ier,ninp,nout )

    case ( "##### restart #####" )

       call ReadInputRestart ( ainp,comment,nrighe,ier,ninp,nout )

    case ( "##### domain #####" )

       call ReadInputDomain ( NumberEntities,ainp,comment,nrighe,ier,ninp,nout,nscr )

    case ( "##### vertices #####" )

       call ReadInputVertices ( NumberEntities,Vertice, ainp,comment,nrighe,ier,.TRUE.,ninp,nout )

    case ( "##### lines #####" )

       call ReadInputLines ( NumberEntities,BoundaryVertex,Tratto, ainp,comment,nrighe,ier,ninp,nout )

    case ( "##### faces #####" )
!AA504 sub
       call ReadInputFaces ( NumberEntities,ainp,comment,nrighe,ier,.TRUE.,ninp,nout )

    case ( "##### geometry file #####" )
!AA504 sub
       call ReadInputExternalFile ( NumberEntities,ainp,comment,nrighe,ier,OnlyTriangle,ninp,nout,ninp2 )

!AA501b start
    case ( "##### bed load transport #####" )
        call ReadBedLoadTransport (ainp,comment,nrighe,ier,ninp,nout,nscr)
       
    case ( "##### medium #####" )

       call ReadInputMedium ( NumberEntities,Med,ainp,comment,nrighe,ier,ninp,nout,nscr )

!AA501b start
    case ( "##### body dynamics #####" )

        if (.not. BodyDynamics_Module) then
          write (nscr,"(1x,a)") " "
          write (nout,"(1x,a)") " "
          write (nscr,"(1x,a)") " >>WARNING! - The body dynamics module is not available."
!AA504sub          
          write (nout,"(1x,a)") " >>WARNING! - The body dynamics module is not available."
          ier = 30
          return
        end if
        call ReadBodyDynamics(ainp,comment,nrighe,ier,ninp,nout)

!AA601 start
    case ( "##### dbsph #####" ) ! lower case letters are required
        call ReadDBSPH(ainp,comment,nrighe,ier,ninp,nout)
!AA601 end        
        
    case ( "##### boundaries #####" )

       call ReadInputBoundaries ( NumberEntities,Partz,Tratto,BoundaryVertex,ainp,comment,nrighe,ier, ninp,nout )

    case ( "##### run parameters #####" )

       call ReadInputRunParameters ( ainp,comment,nrighe,ier,ninp,nout,nscr )

    case ( "##### general physical properties #####" )

       call ReadInputGeneralPhysical ( NumberEntities,ainp,comment,nrighe,ier,ninp,nout )

    case ( "##### output regulation #####" )

       call ReadInputOutputRegulation ( Med,ainp,comment,nrighe,ier,ninp,nout )

    case ( "##### control points #####" )

       call ReadInputControlPoints ( NumberEntities,Control_Points, ainp,comment,nrighe,ier,ninp,nout )

    case ( "##### control lines #####" )

       call ReadInputControlLines ( NumberEntities,Control_Points,Control_Lines, ainp,comment,nrighe,ier, ninp,nout )

    case ( "##### control sections #####" )

       call ReadInputControlSections ( NumberEntities,Control_Sections, ainp,comment,nrighe,ier,ninp,nout )

!AA504 start       
    case ( "##### section flow rate #####" )
       call ReadSectionFlowRate (ainp,comment,nrighe,ier,ninp,nout)
!AA504 end
       
    case ( "##### draw options #####" )

       call ReadInputDrawOptions ( ainp,comment,nrighe,ier,ninp,nout )

    case default 

       ier = 1

    end select
!
!.. reading error was detected: return to calling module
!
    if ( ier /= 0 ) then
      write(nscr,"(/,1x,a,i8,//)") ">> Reading error was detected in INPUT FILE.  error code= ",ier
      write(nout,"(/,1x,a,i8,//)") ">> Reading error was detected in INPUT FILE.  error code= ",ier
      return
!
!.. write the end section record
!
    else if ( ncord > 0 ) then
      write(nout,"(1x,a,/)") lcase(ainp)
    end if
!
  end do  SECTION_LOOP 
!
!.. assign the influence sphere parameters of the smoothed particles
!
  if ( ncord > 0 ) then
!
    write(nout,"(/,1x,a,//)") ">> END OF INPUT FILE"
!
    Domain%h = Domain%dd * Domain%trunc
    doubleh = two * Domain%h
    squareh = Domain%h * Domain%h
    doublesquareh = doubleh * doubleh
    cubich = squareh * Domain%h
    unosuh = one / Domain%h
    unosusquareh = one / squareh
    eta = 0.001d0 * Domain%h
    eta2 = 0.01d0 * squareh
!
  end if
!
  return
  end subroutine ReadInput
!---split

!cfile ReadInputBoundaries.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
  subroutine ReadInputBoundaries ( NumberEntities,Partz,Tratto,BoundaryVertex,ainp,comment,  &
                                   nrighe,ier,ninp,nout )

  use GLOBAL_MODULE
  use AdM_USER_TYPE
!
!.. implicit declarations
  implicit none
!
!.. dummy arguments
  integer(4),               dimension(20)           :: NumberEntities
  type (TyZone),            dimension(NPartZone)    :: Partz
  type (TyBoundaryStretch), dimension(NumTratti)    :: Tratto
  integer(4),               dimension(NumBVertices) :: BoundaryVertex
  integer(4)                                        :: nrighe,ier, ninp,nout
  character( 1)                                     :: comment
  character(80)                                     :: ainp
!
!.. local scalars

!AA504 sub start
  integer(4)        :: n,index,numv,indexi,indexf,Izona,ipointer,Medium,icolor,icord,ioerr,npointv,IC_source_type,Car_top_zone,dx_CartTopog,plan_reservoir_points,nag_aux
  integer(4)        :: i,i1,i2,i_point,ID_first_vertex,ID_last_vertex,dam_zone_ID,dam_zone_n_vertices
  double precision  :: pool_value,shear,velocity,trampa,valp,flowrate,H_res
!AA504 sub end
  character(len=1)  :: pool_plane,bends,slip
  character(len=2)  :: pressu
  character(len=3)  :: move
  character(len=4)  :: tipo
  character(len=6)  :: token_color
  character(len=8)  :: label
  character(len=80) :: token
!
!.. local arrays
  double precision, dimension(3) :: values1,values3
  double precision, dimension(0:3,maxpointsvlaw) :: valuev
!AA504
  double precision :: plan_reservoir_pos(4,2),dam_zone_vertices(4,2)

!.. external assignments
  integer(4),    external :: ptcolorrgb
  character(80), external :: lcase, GetToken
  logical,       external :: ReadCheck
!
!.. Executable statements
!
!.. initializations
!
  npointv = 0
  values3 = zero
  valp = zero
!
!.. in case of restart the cards are not read
!
  if (restart) then
    do while ( TRIM(lcase(ainp)) /= "##### end boundaries #####" )
      call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
      if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"BOUNDARIES DATA",ninp,nout) ) return
    end do
    return
  end if
!
!.. input loading of boundary set starts...
!
  call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
  if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"BOUNDARIES DATA",ninp,nout) ) return
!
!.. loads all the cards of the set
!
  do while ( TRIM(lcase(ainp)) /= "##### end boundaries #####" )
!
!.. assign the label of the boundary condition
!
    label = ainp(1:8)
!
!.. read the boundary index card for the first zone having the same condition
!
    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"BOUNDARIES INDEX",ninp,nout) ) return
    token = GetToken(ainp,1,ioerr)
    read ( token,*,iostat=ioerr ) indexi
!
!.. read the boundary index card for the last zone having the same condition
!
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"BOUNDARY INDEX",ninp,nout) ) return
    indexf = indexi
    token = GetToken(ainp,2,ioerr)
    if ( token /= "" )   read ( token,*,iostat=ioerr ) indexf
    NumberEntities(8) = max(indexf,NumberEntities(8))
!
!.. read the boundary type
!
    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BOUNDARY TYPE",ninp,nout) ) return
    tipo = lcase(ainp(1:4))
!
!.. loads the data for the different boundary conditions
!
    numv     = 0
    ipointer = 0
    move     = "   "
    Medium   = 0
    icolor  = Z'FF0000'  
    values1 = zero
    values3 = zero
    pool_plane = " "
    pool_value = zero
    shear      = zero
    velocity   = zero
    flowrate   = zero
    trampa     = zero
    pressu     = "  "
    valp       = zero

!AA504 start
    IC_source_type = 0
    Car_top_zone = 0
    dx_CartTopog = 0.
    H_res = 0.
    ID_first_vertex = 0
    ID_last_vertex = 0
    plan_reservoir_points = 0
    nag_aux = 0
    dam_zone_ID = 0
    dam_zone_n_vertices = 0
    plan_reservoir_pos = 0. 
    dam_zone_vertices = 0.
!AA504 end
    
!    
    select case ( tipo )
!
!.. boundary condition "leve", "crit" or "open"
!
      case ( "leve", "crit", "open" )    
!
        NumberEntities(3) = NumberEntities(3) + 1
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          token = GetToken(ainp,1,ioerr)
          token_color(1:2) = token(5:6)
          token_color(3:4) = token(3:4)
          token_color(5:6) = token(1:2) 
          read ( token_color,'(Z6)',iostat=ioerr) icolor
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FIXED: RRGGBB COLOR",ninp,nout) ) return
!
!.. boundary condition "fixe"
!
     case ( "fixe" )    
!
        NumberEntities(3) = NumberEntities(3) + 1
!

        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        if ( ioerr == 0 ) read (ainp,*,iostat=ioerr) shear
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FIXED: SHEAR STRESS COEFFICIENT",ninp,nout) ) return
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        token = GetToken(ainp,1,ioerr)
        token_color(1:2) = token(5:6)
        token_color(3:4) = token(3:4)
        token_color(5:6) = token(1:2) 
        read ( token_color,'(Z6)',iostat=ioerr) icolor
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FIXED: RRGGBB COLOR",ninp,nout) ) return
!
!.. boundary condition "sour"
!
     case ( "sour" )    
!
        NumberEntities(3) = NumberEntities(3) + 1
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        if ( ioerr == 0 ) read (ainp,*,iostat=ioerr) Medium
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: MEDIUM INDEX",ninp,nout) ) return
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
!        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: NORMAL VELOCITY, TRAMPA ",ninp,nout) ) return
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: FLOW RATE, TRAMPA ",ninp,nout) ) return
!
        token = GetToken(ainp,1,ioerr)
!        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: NORMAL VELOCITY",ninp,nout) ) return
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: FLOW RATE",ninp,nout) ) return
!        read ( token,*,iostat=ioerr ) velocity
!        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: NORMAL VELOCITY",ninp,nout) ) return
        read ( token,*,iostat=ioerr ) flowrate
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: FLOW RATE",ninp,nout) ) return
        token = GetToken(ainp,2,ioerr)
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: TRAMPA",ninp,nout) ) return
        read ( token,*,iostat=ioerr ) trampa
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: TRAMPA",ninp,nout) ) return
        if (trampa /= 0) then
          write (nout,*) ' '
          write (nout,*) 'TRAMPA in SOURCE boundary is not available. TRAMPA is setted to zero; check the VELOCITY boundary.'
          write (nout,*) ' '
          write (*,*) ' '
          write (*,*) 'TRAMPA in SOURCE boundary is not available. TRAMPA is setted to zero;  check the VELOCITY boundary.'
          write (*,*) ' '
          trampa = zero
        end if
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
!        token = GetToken(ainp,1,ioerr)
!        read ( token,*,iostat=ioerr ) pressu
!        pressu = pressu(1:len_trim(pressu))
        pressu = GetToken(ainp,1,ioerr)
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE PRESSURE TYPE",ninp,nout) ) return
        if (pressu == "pa") then  
          token = GetToken(ainp,2,ioerr)
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE PRESSURE VALUES",ninp,nout) ) return
          read ( token,*,iostat=ioerr ) valp
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE PRESSURE VALUES",ninp,nout) ) return
        else if (pressu == "qp") then
          token = GetToken(ainp,2,ioerr)
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE PIEZO LINE",ninp,nout) ) return
          read ( token,*,iostat=ioerr ) valp
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE PIEZO LINE",ninp,nout) ) return
!        else if (pressu == "pl") then      
!          token = GetToken(ainp,2,ioerr)
!          read ( token,*,iostat=ioerr ) valp
        else
          if ( nout > 0 ) write (nout,*) "Unknown option: ",trim(ainp)," in source boundary."
          stop
        end if
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        token = GetToken(ainp,1,ioerr)
        token_color(1:2) = token(5:6)
        token_color(3:4) = token(3:4)
        token_color(5:6) = token(1:2) 
        read ( token_color,'(Z6)',iostat=ioerr) icolor
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FIXED: RRGGBB COLOR",ninp,nout) ) return
!
        move   = "std"
!
!.. boundary condition "velo"
!
     case ( "velo" )    
!
        NumberEntities(3) = NumberEntities(3) + 1
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        if ( ioerr == 0 ) read (ainp,*,iostat=ioerr) velocity,trampa
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VELO: NORMAL VELOCITY, TRAMPA",ninp,nout) ) return
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        token = GetToken(ainp,1,ioerr)
        token_color(1:2) = token(5:6)
        token_color(3:4) = token(3:4)
        token_color(5:6) = token(1:2) 
        read ( token_color,'(Z6)',iostat=ioerr) icolor
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FIXED: RRGGBB COLOR",ninp,nout) ) return
!
        move   = "std"
        pressu = "pa"
        valp = zero
!
!.. boundary condition "flow"
!
     case ( "flow" )    
!
        NumberEntities(3) = NumberEntities(3) + 1
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        if ( ioerr == 0 ) read (ainp,*,iostat=ioerr) flowrate,trampa
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VELO: FLOW RATE, TRAMPA",ninp,nout) ) return
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        token = GetToken(ainp,1,ioerr)
        token_color(1:2) = token(5:6)
        token_color(3:4) = token(3:4)
        token_color(5:6) = token(1:2) 
        read ( token_color,'(Z6)',iostat=ioerr) icolor
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FIXED: RRGGBB COLOR",ninp,nout) ) return
!
        move   = "std"
        pressu = "pa"
        valp = zero
!
!.. boundary condition "peri"
!
      case ( "peri" )    
!
        NumberEntities(3) = NumberEntities(3) + 1
!
       
        call ReadInputParticlesData ( NumberEntities, Medium,icolor,bends,move,slip,npointv,valuev,values3, &
                                      pressu,valp,ainp,comment,nrighe,ier,ninp,nout )

        if ( ier /= 0 ) return  
        
!AA504 start
        call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
        if (ioerr==0) read (ainp,*,iostat=ioerr) IC_source_type,Car_top_zone
        if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"IC_source_type, Car_top_zone",ninp,nout) ) return
        if (IC_source_type==2) then
           call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
           if (ioerr==0) read (ainp,*,iostat=ioerr) dx_CartTopog,H_res
           if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"dx_CartTopog,H_res",ninp,nout) ) return
           call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
           if (ioerr==0) read (ainp,*,iostat=ioerr) ID_first_vertex,ID_last_vertex
           if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ID_first_vertex,ID_last_vertex",ninp,nout) ) return
           call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
           if (ioerr==0) read (ainp,*,iostat=ioerr) plan_reservoir_points,nag_aux
           if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"plan_reservoir_points and nag_aux",ninp,nout) ) return
           do i2=1,plan_reservoir_points
              call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
              if (ioerr==0) read (ainp,*,iostat=ioerr) plan_reservoir_pos(i2,1),plan_reservoir_pos(i2,2)
              if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"plan_reservoir_vertices",ninp,nout) ) return
           end do
           call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
           if (ioerr==0) read (ainp,*,iostat=ioerr) dam_zone_ID,dam_zone_n_vertices
           if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"dam_zone_ID and dam_zone_vertices",ninp,nout) ) return
           if (dam_zone_ID>1) then
              do i2=1,dam_zone_n_vertices
                 call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
                 if (ioerr==0) read (ainp,*,iostat=ioerr) dam_zone_vertices(i2,1),dam_zone_vertices(i2,2)
                 if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"dam zone vertices",ninp,nout) ) return
              end do
           endif
        endif
!AA504 end
        
!
!.. boundary condition "tapi"
!
      case (  "tapi" )    
!
!.. returns an error if the number of vertices is not equal to two in 3D
!     
        if ( numv /= 2 .and. NumberEntities(1) == 3 ) then 
          write (nout,'(a,i15)') "TAPIS boundary type: 2 vertices are requested:",numv
          ier = 103
          return
        end if
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        if ( ioerr == 0 ) read (ainp,*,iostat=ioerr) shear
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"TAPIS: SHEAR STRESS COEFFICIENT",ninp,nout) ) return
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"TAPIS: VELOCITY COMPONENTS",ninp,nout) ) return
        do n = 1, NumberEntities(1)
          icord = icoordp(n,NumberEntities(1)-1)
          token = GetToken(ainp,n,ioerr)
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,xyzlabel(icord)//" VELOCITY COMPONENT (TAPIS)",ninp,nout) ) return
          read ( token,*,iostat=ioerr ) values1(icord)
        end do
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        token = GetToken(ainp,1,ioerr)
        token_color(1:2) = token(5:6)
        token_color(3:4) = token(3:4)
        token_color(5:6) = token(1:2) 
        read ( token_color,'(Z)',iostat=ioerr) icolor
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"TAPIS: RRGGBB COLOR",ninp,nout) ) return
!
!.. boundary condition "pool" : active only for 3D case
!
       case ("pool")
!
         if ( NumberEntities(1) == 3 ) then      
!   
           NumberEntities(3) = NumberEntities(3) + 1
!
           call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
           pool_plane = GetToken(ainp,1,ioerr)
           if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"POOL: X/Y/Z/ PLANE LABEL",ninp,nout) ) return
           token = GetToken(ainp,2,ioerr)
           if ( ioerr == 0 ) read (token,*,iostat=ioerr) pool_value
           if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"POOL: PLANE VALUE",ninp,nout) ) return
!
           call ReadInputParticlesData ( NumberEntities, &
                                         Medium,icolor,bends,move,slip,npointv,valuev,values3, &
                                         pressu,valp,ainp,comment,nrighe,ier,ninp,nout )

           if ( ier /= 0 ) return  
!
         end if 
!
       case default

         write (nout,*) "Unrecognised BOUNDARY type: ",tipo
         ier = 101
         return
!
    end  select
!
!.. load the parameter structure with input data
!
    if ( ncord > 0 ) then
!
       Izona = NumberEntities(3)
!
       Partz(Izona)%label    = label
       Partz(Izona)%tipo     = tipo
       Partz(Izona)%Medium   = Medium

!AA504 start
       Partz(Izona)%IC_source_type = IC_source_type
       Partz(Izona)%Car_top_zone = Car_top_zone
       if (IC_source_type == 2) then
          Partz(Izona)%dx_CartTopog = dx_CartTopog
          Partz(Izona)%H_res = H_res
          Partz(Izona)%ID_first_vertex = ID_first_vertex
          Partz(Izona)%ID_last_vertex = ID_last_vertex
          Partz(Izona)%plan_reservoir_points = plan_reservoir_points
          Partz(Izona)%nag_aux = nag_aux
          Partz(Izona)%plan_reservoir_pos = plan_reservoir_pos
          Partz(Izona)%dam_zone_ID = dam_zone_ID
          Partz(Izona)%dam_zone_n_vertices = dam_zone_n_vertices
          if (dam_zone_ID>1) Partz(Izona)%dam_zone_vertices = dam_zone_vertices
       endif
!AA504 end

       Partz(Izona)%icol     = icolor
       Partz(Izona)%bend     = bends
       Partz(Izona)%move     = move
       Partz(Izona)%slip     = slip
       if (npointv /= 0) then
         Partz(Izona)%npointv = npointv
         Partz(Izona)%vlaw(icoordp(0:ncord,ncord-1),1:npointv)= valuev(0:ncord,1:npointv)
       end if
       Partz(Izona)%vel      = zero
       Partz(Izona)%vel(icoordp(1:ncord,ncord-1)) = values3(1:ncord)
       Partz(Izona)%trampa   = trampa
       Partz(Izona)%pressure = pressu
       Partz(Izona)%valp     = valp
       Partz(Izona)%Indix(1) = indexi
       Partz(Izona)%Indix(2) = indexf
       if ( pool_plane == "X" .OR. pool_plane == "x" ) Partz(Izona)%ipool = 1
       if ( pool_plane == "Y" .OR. pool_plane == "y" ) Partz(Izona)%ipool = 2
       if ( pool_plane == "Z" .OR. pool_plane == "z" ) Partz(Izona)%ipool = 3
       Partz(Izona)%pool = pool_value
!
!.. load the constraints 
!
       MULTI_INDEX_LOOP: do index = indexi, indexf       
!
         Tratto(index)%tipo        = tipo
         if (ncord == 3) then
           Tratto(index)%numvertices = numv
           Tratto(index)%inivertex   = ipointer
         end if
         Tratto(index)%ShearCoeff  = shear
         Tratto(index)%Medium      = Medium
         Tratto(index)%velocity    = values1
         Tratto(index)%NormVelocity= velocity
         Tratto(index)%FlowRate    = flowrate
         Tratto(index)%trampa      = trampa
         Tratto(index)%zone        = Izona
         Tratto(index)%ColorCode   = icolor
!
         if (nout > 0 .and. index == indexi ) then
!
             if ( index > 1 ) write (nout,*)
             if ( indexf == indexi ) write (nout,"(1x,a,i5,1x,a)")    "Boundary        : ",indexi
             if ( indexf /= indexi ) write (nout,"(1x,a,i5,1x,a,i5)") "Boundary        : ",indexi,"   to",indexf
             write (nout,"(1x,a,2x,a)")    "Type            : ",Tratto(index)%tipo 
             if ( tipo == "fixe" ) then
                write (nout,"(1x,a,1pe12.4)") "Shear coeff.    : ",Tratto(index)%ShearCoeff
             else if ( tipo == "peri" ) then
                write (nout,"(1x,a,i3,1x,a)") "Medium Index    : ",Tratto(index)%Medium
             else if ( tipo == "pool" ) then
                write (nout,"(1x,a,i3,1x,a)") "Medium Index    : ",Tratto(index)%Medium
             else if ( tipo == "tapi" ) then
                write (nout,"(1x,a,1pe12.4)") "Shear coeff.    : ",Tratto(index)%ShearCoeff
                do n = 1, ncord
                   icord = icoordp(n,ncord-1)
                   write (nout,"(1x,a,a,1pe12.4)") &
                   xyzlabel(icord)," Velocity      : ",Tratto(index)%velocity(n)
                end do
             end if
             if (ncord == 2) then
               write (nout,"(1x,a)") "Vertices List"
               write (nout,"(1x,10i5)") BoundaryVertex(Tratto(index)%inivertex:Tratto(index)%inivertex+Tratto(index)%numvertices-1)
               write (nout,"(1x,a,z6)") "Color           : ",Tratto(index)%colorCode
             end if
!
          end if
!
          select case (tipo)
!
            case ("fixe")
              if ( ncord == 3) then
                Tratto(index)%ColorCode = icolor
                if ( nout > 0 .AND. index == indexi ) write (nout,"(1x,a,z8)")      "Color           : ",Tratto(index)%colorCode
                Tratto(index)%ColorCode = icolor
                write (nout,"(1x,a,z8)")      "Color           : ",Tratto(index)%colorCode
              end if
!
            case ("tapi")
              if (ncord == 3) then 
                Tratto(index)%ColorCode = icolor
                if ( nout > 0 .AND. index == indexi ) write (nout,"(1x,a,z8)")      "Color           : ",Tratto(index)%colorCode
                Tratto(index)%ColorCode = icolor
                write (nout,"(1x,a,z8)")      "Color           : ",Tratto(index)%colorCode
              else 
                numv = Tratto(index)%numvertices
                if ( numv /= 2  ) then 
                  if ( nout > 0 ) write (nout,'(a,i15)') "TAPIS boundary type: 2 vertices are requested:",numv
                  ier = 103
                  return
                end if
              end if
!
            case ("peri")
              if (ncord == 2) then 
                i1= BoundaryVertex(Tratto(index)%inivertex)
                i = BoundaryVertex(Tratto(index)%inivertex+Tratto(index)%numvertices-1)
                if ( i /= i1 ) then ! errore se primo ed ultimo vertice sono diversi
                  if ( nout > 0 ) write (nout,'(a,2i15)') "PERIMETER boundary type: first and last vertices are different: ",i,i1
                  ier = 102
                  return
                end if

                if ( nout > 0 ) then

                  write (nout,"(1x,a,i3,1x,a)")  "Zone            : ",Izona,Partz(Izona)%label
                  write (nout,"(1x,a,i3)")       "Medium Index    : ",Partz(Izona)%Medium
                  write (nout,"(1x,a,Z6.6)")     "Color           : ",Partz(Izona)%icol
                  write (nout,"(1x,a,2x,a)")     "Bends           : ",Partz(Izona)%bend
                  write (nout,"(1x,a,2x,a)")     "Movement Type   : ",Partz(Izona)%move
                  write (nout,"(1x,a,2x,a)")     "Boundary Cond.  : ",Partz(Izona)%slip
                  if ( Partz(Izona)%move == "law" ) then
                    write (nout,"(1x,a,i3)")      "Velocity Table - Number of Points: ",Partz(Izona)%npointv
                    do i = 1,Partz(Izona)%npointv
                      write (nout,"(a,i3,1p,4(2x,a,e12.4))") &
                            " Point",i,(xyzlabel(icoordp(n,ncord-1)),Partz(Izona)%vlaw(icoordp(n,ncord-1),i),n=0,ncord)
                    end do
                  end if
                  do n = 1, ncord
                    icord = icoordp(n,ncord-1)
                    write (nout,"(1x,a,a,1pe12.4)") &
                    xyzlabel(icord)," velocity       : ",Partz(Izona)%vel(icord) 
                  end do
                  write (nout,"(1x,a,1pe12.4)") "Time Rampa      : ",Partz(Izona)%trampa
                  write (nout,"(1x,a,2x,a)")    "Pressure Type   : ",Partz(Izona)%pressure
                  write (nout,"(1x,a,1pe12.4)") "Pressure Value  : ",Partz(Izona)%valp
!                 write (nout,"(1x,a,1p,e12.4)") "Monaghan Coeff. : ",Med(Tratto(index)%Medium)%alfaMon
!                 do n = 1, ncord
!                   icord = icoordp(n,ncord-1)
!                   write (nout,"(1x,a,a,1pe12.4)") &
!                   xyzlabel(icord)," min        ",Partz(Izona)%coord(icord,1)
!                   write (nout,"(1x,a,a,1pe12.4)") &
!                   xyzlabel(icord)," max        ",Partz(Izona)%coord(icord,2)
!                 end do
          
                end if ! aggiunta       

              else
                Tratto(index)%ColorCode = icolor
                if ( nout > 0 .AND. index == indexi ) write (nout,"(1x,a,z8)")      "Color           : ",Tratto(index)%colorCode
              end if
!AA504 start              
              write (nout,"(1x,a,i12)")        "IC_source_type  : ",Partz(Izona)%IC_source_type
              write (nout,"(1x,a,i12)")        "Car_top_zone    : ",Partz(Izona)%Car_top_zone
              if (IC_source_type == 2) then
              write (nout,"(1x,a,1pe12.4)")    "dx_CartTopog    : ",Partz(Izona)%dx_CartTopog
              write (nout,"(1x,a,1pe12.4)")    "H_res           : ",Partz(Izona)%H_res  
              write (nout,"(1x,a,i12)")        "ID_first_vertex : ",Partz(Izona)%ID_first_vertex
              write (nout,"(1x,a,i12)")        "ID_last_vertex  : ",Partz(Izona)%ID_last_vertex
              write (nout,"(1x,a,i12)")        "plan_reservoir_points: ",Partz(Izona)%plan_reservoir_points
              write (nout,"(1x,a,i12)")        "nag_aux         : ",Partz(Izona)%nag_aux
              do i_point=1,plan_reservoir_points
              write (nout,"(1x,a,3(1pe12.4))") "plan_reservoir_pos   : ",Partz(Izona)%plan_reservoir_pos(i_point,:)                  
              end do
              write (nout,"(1x,a,i12)")        "dam_zone_ID          : ",Partz(Izona)%dam_zone_ID
              write (nout,"(1x,a,i12)")        "dam_zone_n_vertices  : ",Partz(Izona)%dam_zone_n_vertices  
              if (dam_zone_ID>1) then
              do i_point=1,dam_zone_n_vertices
              write (nout,"(1x,a,3(1pe12.4))") "dam_zone_vertices    : ",Partz(Izona)%dam_zone_vertices(i_point,:)                  
              end do
              endif
              endif
!AA504 end

!
            case ("pool")
              Tratto(index)%ColorCode = icolor
              if ( nout > 0 .AND. index == indexi ) write (nout,"(1x,a,z8)")      "Color           : ",Tratto(index)%colorCode

          end select
   
       end do  MULTI_INDEX_LOOP

    end if
!
!.. end of boundary loading
!
    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
!
  end do

  return
  end subroutine ReadInputBoundaries
!---split

!cfile ReadInputControlLines.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputControlLines ( NumberEntities,Control_Points,Control_Lines, &
                                   ainp,comment,nrighe,ier,ninp,nout )

use GLOBAL_MODULE
use AdM_USER_TYPE

implicit none

integer(4),           dimension(20)        :: NumberEntities
type (TyCtlPoint),    dimension(NPointst)  :: Control_Points
type (TyCtlLine),     dimension(NLines)    :: Control_Lines

integer(4)    :: nrighe,ier, ninp,nout
character( 1) :: comment
character(80) :: ainp

integer(4)    :: n,i,ndiv,icord,ioerr,npts
double precision       :: vp
character(5)  :: txt
character(8)  :: label
double precision, dimension(3) :: values1, values2, values3

character(80), external :: lcase
logical,       external :: ReadCheck


 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"CONTROL LINES DATA",ninp,nout) ) return

 npts = npoints

 do while ( TRIM(lcase(ainp)) /= "##### end control lines #####" )
     
    values1 = zero
    values2 = zero
    values3 = zero

    NumberEntities(5) = NumberEntities(5) + 1
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL LINE LABEL",ninp,nout) ) return
    label(1:8) = ainp(1:8)
    write(txt,"(i5)") NumberEntities(5)

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    read ( ainp,*,iostat=ioerr ) values1(1:ncord)
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL LINE"//txt//" - FIRST POINT",ninp,nout) ) return

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    read ( ainp,*,iostat=ioerr ) values2(1:ncord)
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL LINE"//txt//" - SECOND POINT",ninp,nout) ) return

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    read ( ainp,*,iostat=ioerr ) ndiv
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL LINE"//txt//" - POINTS NUMBER",ninp,nout) ) return

    NumberEntities(6) = NumberEntities(6) + ndiv

    if ( ncord > 0 ) then
       control_lines(NumberEntities(5))%label     = label
       control_lines(NumberEntities(5))%icont(1)  = npts + 1
       control_lines(NumberEntities(5))%icont(2)  = npts + ndiv
       npts = npts + ndiv
       values3(:) = (values2(:)-values1(:))/(ndiv-1)
       vp = Dsqrt( values3(1)*values3(1) + values3(2)*values3(2) + values3(3)*values3(3) )
       if ( nout > 0 ) then
          write (nout,"(1x,a,i3,1x,a)") &
          "Control line      ",NumberEntities(5),"("//control_lines(NumberEntities(5))%label//")"
          write (nout,"(1x,a,i12)") &
          "First Point:      ",control_lines(NumberEntities(5))%icont(1)
          write (nout,"(1x,a,i12)") &
          "Last  Point:      ",control_lines(NumberEntities(5))%icont(2)
       end if

       do i = control_lines(NumberEntities(5))%icont(1), control_lines(NumberEntities(5))%icont(2)
          do n = 1, ncord
             icord = icoordp(n,ncord-1)
             control_points(i)%coord(icord) = values1(n)
          end do
          if ( i == control_lines(NumberEntities(5))%icont(1) ) then
             control_points(i)%dist = zero
          else
             control_points(i)%dist = control_points(i-1)%dist + vp
          end if
          values1 = values1 + values3
          if ( nout > 0 ) then
             write (nout,"(1x,a,i5,1pe12.4,3(3x,a,e12.4))") &
             "Point ",i,control_points(i)%dist, &
             (xyzlabel(icoordp(n,ncord-1))//" = ", &
             control_points(i)%coord(icoordp(n,ncord-1)),n=1,ncord)
             
          end if
       end do
       write (nout,"(1x,a)") " "
    end if

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"CONTROL LINES DATA",ninp,nout) ) return

 end do

return
end subroutine ReadInputControlLines
!---split

!cfile ReadInputControlPoints.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputControlPoints ( NumberEntities,Control_Points, &
                                    ainp,comment,nrighe,ier,ninp,nout )

use GLOBAL_MODULE
use AdM_USER_TYPE

implicit none

integer(4),        dimension(20)        :: NumberEntities
type (TyCtlPoint), dimension(NPointst)  :: Control_Points

integer(4)    :: nrighe,ier, ninp,nout
character( 1) :: comment
character(80) :: ainp

integer(4)    :: n,i,icord,ioerr
double precision, dimension(3) :: values1

character(80), external :: lcase
logical,       external :: ReadCheck


 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"CONTROL POINTS DATA",ninp,nout) ) return

 do while ( TRIM(lcase(ainp)) /= "##### end control points #####" )

    NumberEntities(4) = NumberEntities(4) + 1
    read ( ainp,*,iostat=ioerr ) values1(1:ncord)
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL POINT COORDINATES",ninp,nout) ) return

    if ( ncord > 0 ) then
       Control_Points(NumberEntities(4))%coord(1:3)  = zero
       Control_Points(NumberEntities(4))%dist        = zero
       do n = 1, ncord
          icord = icoordp(n,ncord-1)
          Control_Points(NumberEntities(4))%coord(icord) = values1(n)
       end do
       if ( nout > 0 ) then
          i = NumberEntities(4)
          write (nout,"(1x,a,i5,3(3x,a,1p,e12.4))") &
          "Control point ",i,(xyzlabel(icoordp(n,ncord-1))//" = ",Control_Points(i)%coord(icoordp(n,ncord-1)),n=1,ncord)
       end if
    end if

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"CONTROL POINTS DATA",ninp,nout) ) return

 end do

return
end subroutine ReadInputControlPoints
!---split

!cfile ReadInputControlSections.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputControlSections ( NumberEntities,Control_Sections, &
                                      ainp,comment,nrighe,ier,ninp,nout )

use GLOBAL_MODULE
use AdM_USER_TYPE

implicit none

integer(4)    :: nrighe,ier, ninp,nout, npts
integer(4),        dimension(20)             :: NumberEntities
type (TySection),  dimension(0:Nsections+1)  :: Control_Sections

character( 1) :: comment
character(80) :: ainp

integer(4)    :: icord, icor2, icor3, icolor, ndiv, ioerr
character(8)  :: label
character(80) :: token
double precision, dimension(3)   :: vp
double precision, dimension(3,2) :: values

character(1), dimension(3) :: CoordLabel = (/ "x", "y", "z" /)

integer(4),    external :: NumberSectionPoints
character(80), external :: lcase, GetToken
logical,       external :: ReadCheck


 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"CONTROL SECTIONS DATA",ninp,nout) ) return

 npts = npoints+npointsl

 do while ( TRIM(lcase(ainp)) /= "##### end control sections #####" )

    NumberEntities(12) = NumberEntities(12) + 1
!   if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL SECTION LABEL",ninp,nout) ) return
    label(1:8) = ainp(1:8)

    vp = zero
    values(:,1) = -99999999.
    values(:,2) =  99999999.

   !codice coordinata costante e suo valore
    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    token = lcase(GetToken(ainp,1,ioerr))
    if ( ioerr == 0 ) then
       select case ( token(1:1) )
          case ( "x", "y" ,"z" )
             icord = index("xyz",token(1:1))
             token = lcase(GetToken(ainp,2,ioerr))
             if ( ioerr == 0 ) read ( token,*,iostat=ioerr ) vp(icord)
          case default
             ioerr = -1 !forza errore caso opzione NON riconosciuta
       end select
    end if
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL SECTION "//label//" - CONSTANT COORD. DEFINITION",ninp,nout) ) return

   !codice prima coordinata limitata e valori min. max. (esclusa la coord. costante!)
    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    token = lcase(GetToken(ainp,1,ioerr))
    if ( ioerr == 0 ) then
       select case ( token(1:1) )
          case ( "x", "y" ,"z" )
             icor2 = index("xyz",token(1:1))
             if ( icor2 /= icord ) then
                token = lcase(GetToken(ainp,2,ioerr))
                if ( ioerr == 0 ) read ( token,*,iostat=ioerr ) values(icor2,1)
                if ( ioerr == 0 ) token = lcase(GetToken(ainp,3,ioerr))
                if ( ioerr == 0 ) read ( token,*,iostat=ioerr ) values(icor2,2)
             else
                ioerr = -1 !forza errore caso opzione errata
             end if
          case default
             ioerr = -1 !forza errore caso opzione NON riconosciuta
       end select
    end if
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL SECTION "//label//" - FIRST LIMIT DEFINITION",ninp,nout) ) return

   !codice seconda coordinata limitata e valori min. max. (esclusa la prima coord. e quella costante!)
    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    token = lcase(GetToken(ainp,1,ioerr))
    if ( ioerr == 0 ) then
       select case ( token(1:1) )
          case ( "x", "y" ,"z" )
             icor3 = index("xyz",token(1:1))
             if ( icor3 /= icord .AND. icor3 /= icor2 ) then
                token = lcase(GetToken(ainp,2,ioerr))
                if ( ioerr == 0 ) read ( token,*,iostat=ioerr ) values(icor3,1)
                if ( ioerr == 0 ) token = lcase(GetToken(ainp,3,ioerr))
                if ( ioerr == 0 ) read ( token,*,iostat=ioerr ) values(icor3,2)
             else
                ioerr = -1 !forza errore caso opzione errata
             end if
          case default
             ioerr = -1 !forza errore caso opzione NON riconosciuta
       end select
    end if
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL SECTION "//label//" - SECOND LIMIT DEFINITION",ninp,nout) ) return

   !codice colore ptcl
    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    token = lcase(GetToken(ainp,1,ioerr))
    if ( ioerr == 0 ) then
       select case ( token(1:5) )
          case ( "color" )
             token = lcase(GetToken(ainp,2,ioerr))
             if ( ioerr == 0 ) read ( token,*,iostat=ioerr ) icolor
          case default
             ioerr = -1 !forza errore caso opzione NON riconosciuta
       end select
    end if
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL SECTION "//label//" - COLOR INDEX",ninp,nout) ) return
  
   !Calcola il numero di punti da generare
    Ndiv = NumberSectionPoints ( values, CoordLabel(icord) )
    NumberEntities(13) = NumberEntities(13) + Ndiv

    if ( ncord > 0 ) then

       Control_Sections(NumberEntities(12))%Label             = label
       Control_Sections(NumberEntities(12))%Tipo              = CoordLabel(icord)  !//"+"
       Control_Sections(NumberEntities(12))%Constant(:)       = vp(:)
       Control_Sections(NumberEntities(12))%icont(1)          = npts + 1
       Control_Sections(NumberEntities(12))%icont(2)          = npts + Ndiv
       Control_Sections(NumberEntities(12))%XYZRange(:,1)     = -99999999.
       Control_Sections(NumberEntities(12))%XYZRange(:,2)     =  99999999.
       Control_Sections(NumberEntities(12))%XYZRange(icor2,1) = minval(values(icor2,:))
       Control_Sections(NumberEntities(12))%XYZRange(icor2,2) = maxval(values(icor2,:))
       Control_Sections(NumberEntities(12))%XYZRange(icor3,1) = minval(values(icor3,:))
       Control_Sections(NumberEntities(12))%XYZRange(icor3,2) = maxval(values(icor3,:))
       Control_Sections(NumberEntities(12))%TGLsection        = zero
       Control_Sections(NumberEntities(12))%TGLsection(1,1)   = one
       Control_Sections(NumberEntities(12))%TGLsection(2,2)   = one
       Control_Sections(NumberEntities(12))%TGLsection(3,3)   = one
       Control_Sections(NumberEntities(12))%ColorCode         = icolor

       call  CreateSectionPoints ( vp, values, CoordLabel(icord), NumberEntities(12) )
       npts = npts + Ndiv

       if ( nout > 0 ) then
          write (nout,"(1x,a,i3,1x,a)") &
          "Control section   ",NumberEntities(12),"("//Control_Sections(NumberEntities(12))%label//")"
          write (nout,"(1x,a,1x,a,f12.4)") "Constant Coordinate   ", &
          Control_Sections(NumberEntities(12))%Tipo//"=",Control_Sections(NumberEntities(12))%Constant(icord)
          write (nout,"(1x,a,2f12.4)") &
          CoordLabel(icor2)//" Coordinate Limits  ",Control_Sections(NumberEntities(12))%XYZRange(icor2,:)
          write (nout,"(1x,a,2f12.4)") &
          CoordLabel(icor3)//" Coordinate Limits  ",Control_Sections(NumberEntities(12))%XYZRange(icor3,:)
          write (nout,"(1x,a,i12)") &
          "First Point:      ",Control_Sections(NumberEntities(12))%icont(1)
          write (nout,"(1x,a,i12)") &
          "Last  Point:      ",Control_Sections(NumberEntities(12))%icont(2)
       end if

    end if

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"CONTROL SECTIONS DATA",ninp,nout) ) return

 end do

return
end subroutine ReadInputControlSections
!---split

!cfile ReadInputDomain.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputDomain ( NumberEntities,ainp,comment,nrighe,ier,ninp,nout,nscr )

use GLOBAL_MODULE
use AdM_USER_TYPE

implicit none

integer(4)    :: nrighe,ier, ninp,nout,nscr
character( 1) :: comment
character(80) :: ainp

integer(4),dimension(20) :: NumberEntities

integer(4)       :: ioerr
double precision :: dd, trunc
character(80)    :: token

character(80),external :: lcase, GetToken
logical,      external :: ReadCheck

!
!.. in case of restart the cards are not read
!
  if (restart) then
    do while ( TRIM(lcase(ainp)) /= "##### end domain #####" )
      call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
      if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DOMAIN DATA",ninp,nout) ) return
    end do
    return
  end if
!
 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DOMAIN DATA",ninp,nout) ) return

!do while ( lcase(ainp(1:22)) /= "##### end domain #####" )
 do while ( TRIM(lcase(ainp)) /= "##### end domain #####" )


    token = lcase(GetToken(ainp,1,ioerr))
    read ( token,*,iostat=ioerr ) NumberEntities(1)
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DOMAIN COORDINATES NUMBER",ninp,nout) ) return 
    if ( ncord > 0 .AND. nout > 0 ) then
       write (nout,"(1x,a,i3,1x,a)") "Domain Dimension       : ",ncord,ncordlabel(ncord)
    end if

    token = lcase(GetToken(ainp,2,ioerr))
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DOMAIN TYPE",ninp,nout) ) return 

    select case ( token(1:4) )

!
!AA406 sub
       case ( "bsph","semi" ) 
!
         if (.not. DBSPH_Module .and. token(1:4) == "bsph") then
           write (nscr,"(1x,a)") " "
           write (nout,"(1x,a)") " "
           write (nscr,"(1x,a)") " >>WARNING! - The DBSPH module is not available."
           write (nout,"(1x,a)") " >>WARNING! - The DBSPH module is not available."
           ier = 40
           return
         end if
!
         Domain%tipo = token(1:4)
         if ( ncord > 0 .AND. nout > 0 ) then
           write (nout,"(1x,a,1x,a)"   ) "Domain Type            : ",trim(token)
         end if      

       case default
         ier = 3
         return

    end select

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( ioerr == 0 ) read ( ainp,*,iostat=ioerr ) dd, trunc
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DD & TRUNC",ninp,nout) ) return
    Domain%dd    = dd
    Domain%trunc = trunc

    token = lcase(GetToken(ainp,3,ioerr))
!    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RANDOM INITIAL POSITION",ninp,nout) ) return 
    if (token(1:1) == 'r') then
      Domain%RandomPos = 'r'
    else
      Domain%RandomPos = 'n'
    end if

    if ( ncord > 0 .AND. nout > 0 ) then
       write (nout,"(1x,a,1pe12.4)") "Dd                     : ",dd
       write (nout,"(1x,a,1pe12.4)") "Trunc                  : ",trunc
       write (nout,"(1x,a,1x,a)") "Random Initial Position: ",Domain%RandomPos
       write (nout,"(1x,a)")  " "
    end if

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DOMAIN DATA",ninp,nout) ) return

 end do

return
end subroutine ReadInputDomain
!---split

!cfile ReadInputDrawOptions.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputDrawOptions ( ainp,comment,nrighe,ier,ninp,nout )

use GLOBAL_MODULE                            
use AdM_USER_TYPE

implicit none

integer(4)    :: nrighe,ier, ninp,nout
character( 1) :: comment
character(80) :: ainp
character(4)  :: steptime

integer(4)       :: ioerr
character(80)    :: token

character(80), external :: lcase, GetToken
logical,       external :: ReadCheck

 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DRAW OPTIONS DATA",ninp,nout) ) return

!do while ( lcase(ainp(1:28)) /= "##### end draw options #####" )
 do while ( TRIM(lcase(ainp)) /= "##### end draw options #####" )

    select case (lcase(GetToken(ainp,1,ioerr)))

       case ("vtkconverter")
          token = lcase(GetToken(ainp,(2),ioerr))
!          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FREQUENCY OPTION",ninp,nout) ) return
          select case (token)
             case ("any")
                token = lcase(GetToken(ainp,(3),ioerr))
                read ( token,*,iostat=ioerr ) freq_time
                if ( ncord > 0 .AND. nout > 0 ) write (nout,"(1x,a,1pe12.4,a)") "VTKConversion any :",freq_time," seconds."
                val_time  = zero  !freq_time
             case ("at")
                token = lcase(GetToken(ainp,(3),ioerr))
                read ( token,*,iostat=ioerr ) freq_time
                if ( ncord > 0 .AND. nout > 0 ) write (nout,"(1x,a,1pe12.4,a)") "VTKConversion at :",freq_time," second."
                val_time  = freq_time
                freq_time = -freq_time
             case ("all")
                token = lcase(GetToken(ainp,(3),ioerr))
                read ( token,*,iostat=ioerr ) steptime
                if (steptime == 'time') then
                  freq_time = Domain%memo_fr
                  val_time  = zero  !-Domain%memo_fr
                  if ( ncord > 0 .AND. nout > 0 ) write (nout,"(1x,a,1pe12.4,a)") "VTKConversion every :",freq_time," second."
                else if (steptime == 'step') then
                  freq_time = zero
                  val_time  = const_m_9999
                  if ( ncord > 0 .AND. nout > 0 ) write (nout,"(1x,a)") "VTKConversion all steps."
                else
                  freq_time = Domain%memo_fr
                  val_time  = zero  !-Domain%memo_fr
                  if ( ncord > 0 .AND. nout > 0 ) write (nout,"(1x,a,1pe12.4,a)") "VTKConversion every :",freq_time," second."
                end if
             case default
                  freq_time = Domain%memo_fr
                  val_time  = zero  !-Domain%memo_fr
                  if ( ncord > 0 .AND. nout > 0 ) write (nout,"(1x,a,1pe12.4,a)") "VTKConversion every :",freq_time," second."
          end select
          if ( ncord > 0 .AND. nout > 0 ) write (nout,"(1x,a)") " "
          vtkconv = .TRUE.

       case default

          ier = 4
          return

    end select

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DRAW OPTIONS DATA",ninp,nout) ) return

 end do

return
end subroutine ReadInputDrawOptions
!---split

!cfile ReadInputExternalFile.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!AA504 sub
subroutine ReadInputExternalFile (NumberEntities,ainp,comment,nrighe,ier,OnlyTriangle,ninp,nout,ninp2 )

use GLOBAL_MODULE
use AdM_USER_TYPE
!AA504
use ALLOC_Module

implicit none

integer(4),               dimension(20)                    :: NumberEntities
!AA504 rm part

integer(4)    :: nrighe,ier, ninp,nout, ninp2
character( 1) :: comment
character(80) :: ainp
logical       :: OnlyTriangle

integer(4)    :: ioerr

character(80), external :: lcase, GetToken
logical,       external :: ReadCheck

!
!.. in case of restart the cards are not read
!
  if (restart) then
    do while ( TRIM(lcase(ainp)) /= "##### end geometry file #####" )
      call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
      if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"GEOMETRY DATA",ninp,nout) ) return
    end do
    return
  end if
!
 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"GEOMETRY FILE",ninp,nout) ) return

 OnlyTriangle = .TRUE.
 
!do while ( lcase(ainp(1:29)) /= "##### end geometry file #####" )
 do while ( TRIM(lcase(ainp)) /= "##### end geometry file #####" )

    open(ninp2,file=trim(ainp),form="formatted",status="old",iostat=ioerr)

    if ( nout > 0 ) then
       if ( ioerr == 0 ) then
          write (nout,"(1x,3a)") "Geometry File: ",trim(ainp)
       else
          write (nout,"(1x,3a)") "Geometry File: ",trim(ainp)," not found!"
          return
       end if
    end if


   !Legge prima riga file
    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp2 )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"GEOMETRY FILE",ninp2,nout) ) return

    SECTION_LOOP: do while ( ioerr == 0 )

       select case ( TRIM(lcase(trim(ainp))) )

       case ( "##### vertices #####" )

          call ReadInputVertices ( NumberEntities,Vertice, &
                                   ainp,comment,&
                                   nrighe,ier, .FALSE.,ninp2,nout )

       case ( "##### lines #####" )

          call ReadInputLines    ( NumberEntities,BoundaryVertex,Tratto, &
                                   ainp,comment, &
                                   nrighe,ier, ninp2,nout )

       case ( "##### faces #####" )

!AA504 sub           
          call ReadInputFaces    ( NumberEntities,ainp,comment,nrighe,ier,.FALSE.,ninp2,nout)

       case default


       end select

       call ReadRiga ( ainp,comment,nrighe,ioerr,ninp2 )

      !se EOF esce, altrimenti controlla errore
       if ( ioerr == -1 ) cycle SECTION_LOOP

       if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"GEOMETRY FILE",ninp,nout) ) return

    end do  SECTION_LOOP

    close (ninp2)

    if ( nout > 0 ) then
       write (nout,"(1x,3a)") "End Reading Geometry File"
    end if

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"GEOMETRY FILE",ninp,nout) ) return

 end do

return
end subroutine ReadInputExternalFile
!---split

!cfile ReadInputFaces.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!AA504 sub
subroutine ReadInputFaces (NumberEntities,ainp,comment,nrighe,ier,prtopt,ninp,nout)

use GLOBAL_MODULE
use AdM_USER_TYPE
!AA504
use ALLOC_Module

implicit none

integer(4),           dimension(20)        :: NumberEntities
!AA504rm line

integer(4)    :: nrighe,ier, ninp,nout
logical(4)    :: prtopt
character( 1) :: comment
character(80) :: ainp

integer(4)    :: n,i,ioerr, stretch
character(8)  :: label
integer(4), dimension(4) :: ivalues

character(80), external :: lcase, GetToken
logical,       external :: ReadCheck
!
!.. in case of restart the cards are not read
!
  if (restart) then
    do while ( TRIM(lcase(ainp)) /= "##### end faces #####" )
      call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
      if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"FACES DATA",ninp,nout) ) return
    end do
    return
  end if
!
 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"FACES DATA",ninp,nout) ) return

 do while ( TRIM(lcase(ainp)) /= "##### end faces #####" )

    select case ( TRIM(Domain%tipo) )

!
!AA406 sub
       case ( "semi","bsph" ) 
!

          ivalues = 0
          read ( ainp,*,iostat=ioerr ) i, ivalues, stretch  ! ivalues(4) - vertice 4 - deve essere 0 se triangolo

          write (label,"(i8)") i
          if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"FACE n."//label,ninp,nout) ) return

          NumberEntities(11) = max(i,NumberEntities(11))

         !conta i quadrilateri da dividere eventualmente solo su opzione
          if ( ivalues(4)> 0 ) NumberEntities(18) = NumberEntities(18) + 1

          if ( ncord > 0 ) then
             if ( BoundaryFace(i)%Node(1)%name == 0 ) then 
                BoundaryFace(i)%Node(1:MAXFACENODES)%name = ivalues(1:MAXFACENODES)
                BoundaryFace(i)%stretch = stretch
             else
                if ( nout > 0 ) then
                   write (nout,*) "Face definition: ",trim(ainp)
                   write (nout,*) "Face already defined: ",i,BoundaryFace(i)%Node(1:MAXFACENODES)%name
                end if
                ier = 104
                return
             end if
          end if

       case default

          if ( nout > 0 ) then
             write (nout,*) "Unknown Domain Type: ",Domain%tipo
          end if
          ier = 2
          return


    end select

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"FACES DATA",ninp,nout) ) return

 end do


 if ( ncord > 0 .AND. nout > 0 .AND. prtopt ) then

    write (nout,*)
    write (nout,"(1x,a)") "List of faces:"
    do n = 1, NumberEntities(11)
!AA504
       write (nout,"(i10,' - ',4i10,' - ',i8)") n,BoundaryFace(n)%Node(1)%name,BoundaryFace(n)%Node(2)%name,BoundaryFace(n)%Node(3)%name,BoundaryFace(n)%Node(4)%name,BoundaryFace(n)%stretch
    end do

 end if

return
end subroutine ReadInputFaces
!---split

!cfile ReadInputGeneralPhysical.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputGeneralPhysical ( NumberEntities,ainp,comment,nrighe,ier,ninp,nout )

use GLOBAL_MODULE                                     
use AdM_USER_TYPE

implicit none

integer(4)    :: nrighe,ier, ninp,nout
character( 1) :: comment
character(80) :: ainp

integer(4),        dimension(20)     :: NumberEntities

integer(4)    :: n,icord,ioerr
double precision       :: prif
double precision, dimension(3) :: values1

character(80), external :: lcase
logical,       external :: ReadCheck

!
!.. in case of restart the cards are not read
!
  if (restart) then
    do while ( TRIM(lcase(ainp)) /= "##### end general physical properties #####" )
      call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
      if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"GENERAL PHYSICAL PROPERTIES DATA",ninp,nout) ) return
    end do
    return
  end if
!
 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"GENERAL PHYSICAL PROPERTIES DATA",ninp,nout) ) return

!do while ( lcase(ainp(1:43)) /= "##### end general physical properties #####" )
 do while ( TRIM(lcase(ainp)) /= "##### end general physical properties #####" )

    read ( ainp,*,iostat=ioerr ) values1(1:NumberEntities(1))
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"GRAVITAL ACCELERATION VECTOR",ninp,nout) ) return

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    read ( ainp,*,iostat=ioerr ) prif
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"REFERENCE PRESSURE",ninp,nout) ) return

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"GENERAL PHYSICAL PROPERTIES DATA",ninp,nout) ) return

 end do

 if ( ncord > 0 ) then
    Domain%grav(:) = zero            
    do n = 1, NumberEntities(1)
       icord = icoordp(n,ncord-1)
       Domain%grav(icord) = values1(n)
       if ( nout > 0 ) write (nout,"(1x,a,a,1p,e12.4)") &
       xyzlabel(icord),"gravity acceler. :",Domain%grav(icord)
    end do
    Domain%prif = prif
    if ( nout > 0 ) write (nout,"(1x,a,1p,e12.4)") "P rif:            :",Domain%prif
    if ( nout > 0 ) write (nout,"(1x,a)") " "
 end if

return
end subroutine ReadInputGeneralPhysical
!---split

!cfile ReadInputLines.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputLines ( NumberEntities,BoundaryVertex,Tratto, &
                            ainp,comment,nrighe,ier,ninp,nout )

use GLOBAL_MODULE
use AdM_USER_TYPE

implicit none

integer(4),               dimension(20)            :: NumberEntities
integer(4),               dimension(NumBVertices)  :: BoundaryVertex
type (TyBoundaryStretch), dimension(NumTratti)     :: Tratto

integer(4)    :: nrighe,ier, ninp,nout
character( 1) :: comment
character(80) :: ainp

integer(4)    :: n,i,ioerr
integer(4)    :: i1,index,numv,numv_line,ipointer
character(5)  :: txt
character(80) :: token

integer(4), parameter :: MAXLINENODES = 20
!integer(4), dimension(MAXLINENODES) :: ivalues

character(80), external :: lcase, GetToken
logical,       external :: ReadCheck

!
!.. in case of restart the cards are not read
!
  if (restart) then
    do while ( TRIM(lcase(ainp)) /= "##### end lines #####" )
      call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
      if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"LINES DATA",ninp,nout) ) return
    end do
    return
  end if
!
 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"LINES DATA",ninp,nout) ) return

 do while ( TRIM(lcase(ainp)) /= "##### end lines #####" )

    select case ( TRIM(Domain%tipo) )

!
!AA406 sub
       case ( "semi","bsph" ) 
!

         !Lettura dei vertici che definiscono il boundary
          numv = 0

         !Legge l'indice della linea 
          numv_line = 1
          token = GetToken(ainp,numv_line,ioerr)
          if ( ioerr == 0 ) read ( token,*,iostat=ioerr ) index
          NumberEntities(8) = max(NumberEntities(8),index)  !NumTratti

          VERTEX_LOOP: do while ( ioerr == 0 )
             numv_line = numv_line + 1
             token = GetToken(ainp,numv_line,ioerr)
            !esce quando trova inizio commento o EOR (NON ci sono piu' dati sulla linea di input)
             if ( ioerr /= 0 .OR. &
                  trim(token) == "" .OR. &
                  ichar(trim(token(1:1))) == 9 .OR. &
                  token(1:1) == "!" ) exit VERTEX_LOOP
            !controlla caso codice continuazione linea
             if ( token(1:1) == "&" .OR. token(1:1) == ">" ) then
                call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
                if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VERTICES LIST (continue...)",ninp,nout) ) return
                numv_line = 0
                cycle VERTEX_LOOP
             end if
             numv  = numv + 1
             if ( numv > MAXLINENODES ) then
               stop 'ERRORE in ReadInputLines numv > MAXLINENODES'
             end if
             NumberEntities(9) = NumberEntities(9) + 1  !contatore vertici dei boundaries
             read (token,*,iostat=ioerr) i
             write(txt,"(i5)") i
             if ( numv == 1 ) i1 = i
             if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VERTEX n."//txt,ninp,nout) ) return
             if ( ncord > 0 ) then
                if ( numv == 1 ) ipointer = NumberEntities(9)
                BoundaryVertex(NumberEntities(9)) = i
             end if
          end do VERTEX_LOOP

         !Conta numero di BoundarySide
          NumberEntities(10) = NumberEntities(10) + numv - 1

          if ( ncord > 0 ) then
             Tratto(index)%numvertices = numv
             Tratto(index)%inivertex   = ipointer
          end if

       case default

          if ( nout > 0 ) then
             write (nout,*) "Unknown Domain Type: ",Domain%tipo
          end if
          ier = 2
          return


    end select

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"LINES DATA",ninp,nout) ) return

 end do


 if ( ncord > 0 .AND. nout > 0 ) then

    write (nout,"(1x,a)") "List of lines"
    write (nout,*)
    do n = 1, NumberEntities(8)
       write (nout,"(1x,a,i3,1x,a)") "Line: ",n
       write (nout,"(1x,a,i3,1x,a)") "Number of Vertices:  ",Tratto(n)%numvertices
       write (nout,"(1x,a,i3,1x,a)") "Vertices Pointer:    ",Tratto(n)%inivertex
       write (nout,"(1x,a,i3,1x,a)") "Vertices List"
       write (nout,"(1x,10i5)") BoundaryVertex(Tratto(n)%inivertex:Tratto(n)%inivertex+Tratto(n)%numvertices-1)
       write (nout,"(1x,a)") " "
    end do

 end if

return
end subroutine ReadInputLines
!---split

!cfile ReadInputMedium.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ReadInputMedium
!
! Last updating : Jul 03, 2012
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Agate             2009-2011      varie
! 04  Agate             03/07/12       Module exclusion (Diffusion, Esplosion, Temporal schemes)
!
!************************************************************************************
! Module purpose : Read Medium from input file
!
! Calling routine: readinput
!
! Called routines: none
!
!************************************************************************************
!
subroutine ReadInputMedium ( NumberEntities,Med, &
                             ainp,comment,nrighe,ier, ninp,nout,nscr )

use GLOBAL_MODULE                            
use AdM_USER_TYPE

implicit none

integer(4),     dimension(20)      :: NumberEntities
type (TyMedium),dimension(NMedium) :: Med

integer(4)    :: nrighe,ier, ninp,nout,nscr
character( 1) :: comment
character(80) :: ainp

integer(4)       :: index,nitersol
integer(4)       :: ioerr
!AA504 sub
double precision :: den0, eps, alfaMon, betaMon, visc, viscmx, taucri, cuin, phi, Cs, &
                    cons, codif, Settling, coes, Rough, D50, Gamma, InitialIntEn,d_90,porosity
character(8)     :: tipo, erosionmodel
character(80)    :: token

character(80), external :: GetToken, lcase
logical,       external :: ReadCheck

!
!.. in case of restart the cards are not read
!
  if (restart) then
    do while ( TRIM(lcase(ainp)) /= "##### end medium #####" )
      call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
      if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"MEDIUM DATA",ninp,nout) ) return
    end do
    return
  end if
!
 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"MEDIUM DATA",ninp,nout) ) return

 do while ( TRIM(lcase(ainp)) /= "##### end medium #####" )

    read ( ainp,*,iostat=ioerr ) tipo
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"MEDIUM TYPE",ninp,nout) ) return

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    read ( ainp,*,iostat=ioerr ) index
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"MEDIUM INDEX",ninp,nout) ) return
    NumberEntities(2) = max(NumberEntities(2),index)
!
    den0         = zero
    eps          = zero
    alfaMon      = zero
    betaMon      = zero
    codif        = zero
    Settling     = zero
    Gamma        = zero
    InitialIntEn = zero
    visc         = zero
    Rough        = zero
    Cs           = zero
    cons         = zero
    viscmx       = zero
    taucri       = zero
    cuin         = zero
    coes         = zero
    phi          = zero
    D50          = zero
    nitersol     = 0
!AA504 start
    porosity     = 0.d0
    d_90         = 0.d0
!AA504 end
    erosionmodel = "-"
!
    tipo = lcase(tipo)
    select case ( tipo )
!
! gas per esplosioni con viscosita' dinamica assegnata
!
       case ( "gas     " )
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) den0, eps
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"WATER DENSITY & COMPRIMIBILITY",ninp,nout) ) return
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) alfaMon, betaMon
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ALPHA e BETA MONAGHAN",ninp,nout) ) return
!.. modello diffusione
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) codif, Settling
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DIFFUSION COEFF",ninp,nout) ) return
!.. modello diffusione
!.. modello esplosioni
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) gamma, InitialIntEn
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"EXPLOSION COEFF",ninp,nout) ) return
!.. modello esplosioni
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) visc
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DYNAMIC VISCOSITY",ninp,nout) ) return            
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) Rough
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ROUGH COEFFICIENT",ninp,nout) ) return
!
! fluidi Newtoniani con viscosita' dinamica assegnata
!
       case ( "liquid  " )
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) den0, eps
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"WATER DENSITY & COMPRIMIBILITY",ninp,nout) ) return
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) alfaMon, betaMon
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ALPHA e BETA MONAGHAN",ninp,nout) ) return
!.. modello diffusione
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) codif, Settling
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DIFFUSION COEFF",ninp,nout) ) return
!.. modello diffusione
!.. modello esplosioni
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) gamma, InitialIntEn
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"EXPLOSION COEFF",ninp,nout) ) return
!.. modello esplosioni
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) visc
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DYNAMIC VISCOSITY",ninp,nout) ) return            
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) Rough
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ROUGH COEFFICIENT",ninp,nout) ) return
!
! fluidi Newtoniani con viscosita' dinamica calcolata con modello di turbolenza
!        mixing-length legata alla dimensione della particella vedi Smagorinsky) 
!
       case ( "smagorin" )
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) den0, eps
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"WATER DENSITY & COMPRIMIBILITY",ninp,nout) ) return
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) alfaMon, betaMon
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ALPHA e BETA MONAGHAN",ninp,nout) ) return
!.. modello diffusione
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) codif, Settling
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DIFFUSION COEFF",ninp,nout) ) return
!.. modello diffusione
!.. modello esplosioni
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) gamma, InitialIntEn
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"EXPLOSION COEFF",ninp,nout) ) return
!.. modello esplosioni
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) visc
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DYNAMIC VISCOSITY",ninp,nout) ) return            
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) Cs
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SMAGORINSKY CONST",ninp,nout) ) return
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) Rough
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ROUGH COEFFICIENT",ninp,nout) ) return
!
! fluidi non-Newtoniani  con formula di viscosita' apparente di Chen
!
       case ( "general " )
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) den0, eps
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"WATER DENSITY & COMPRIMIBILITY",ninp,nout) ) return
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) alfaMon, betaMon
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ALPHA e BETA MONAGHAN",ninp,nout) ) return
!.. modello diffusione
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) codif, Settling
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DIFFUSION COEFF",ninp,nout) ) return
!.. modello diffusione
!.. modello esplosioni
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) gamma, InitialIntEn
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"EXPLOSION COEFF",ninp,nout) ) return
!.. modello esplosioni
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) cons, viscmx
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONSISTENCY & VISCO MAX",ninp,nout) ) return            
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) taucri, cuin
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"TAUCRI & CUIN",ninp,nout) ) return
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) Rough
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ROUGH COEFFICIENT",ninp,nout) ) return
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) nitersol
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"N° ITERATION SOLID",ninp,nout) ) return
          visc     = viscmx
!
! fluidi non-Newtoniani (materiale granulare) con formula di viscosita' apparente
!
       case ( "granular" )
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) den0, eps
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"WATER DENSITY & COMPRIMIBILITY",ninp,nout) ) return
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) alfaMon, betaMon
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ALPHA e BETA MONAGHAN",ninp,nout) ) return
!.. modello diffusione
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) codif, Settling
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DIFFUSION COEFF",ninp,nout) ) return
!.. modello diffusione
!.. modello esplosioni
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          read ( ainp,*,iostat=ioerr ) gamma, InitialIntEn
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"EXPLOSION COEFF",ninp,nout) ) return
!.. modello esplosioni
!AA504 start
          if (Granular_flows_options%ID_erosion_criterion==1) then
             call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
             read ( ainp,*,iostat=ioerr ) phi
             if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"PHI",ninp,nout) ) return
             call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
             read ( ainp,*,iostat=ioerr ) porosity,D50,d_90
             if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"POROSITY, D50 and D_90",ninp,nout) ) return
          else
!AA504 end              
            call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
            read ( ainp,*,iostat=ioerr ) coes, viscmx, visc
            if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"COHESION, VISCO MAX & VISCO",ninp,nout) ) return            
            call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
            read ( ainp,*,iostat=ioerr ) phi
            if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"PHI",ninp,nout) ) return
!!            call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
!!            read ( ainp,*,iostat=ioerr ) Rough
!!            if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ROUGH COEFFICIENT",ninp,nout) ) return
            call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
            token = lcase(GetToken(ainp,1,ioerr))
            read ( token,*,iostat=ioerr ) Rough
            if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ROUGH COEFFICIENT",ninp,nout) ) return
            token = lcase(GetToken(ainp,2,ioerr))
            if (token == "") then
              D50 = zero
              erosionmodel = "mohr    "
              write (nout,"(a)") " !!!!!!!!!!!!!!!!!!!!!"
              write (nscr,"(a)") " !!!!!!!!!!!!!!!!!!!!!"
              write (nout,"(a,a)") " ATTENTION!!! The erosion model has not been declared,", &
                                   " the Mohr-Coulomb model will be used."
              write (nscr,"(a,a)") " ATTENTION!!! The erosion model has not been declared,", &
                                   " the Mohr-Coulomb model will be used."
              write (nout,"(a)") " !!!!!!!!!!!!!!!!!!!!!"
              write (nscr,"(a)") " !!!!!!!!!!!!!!!!!!!!!"
            else
              read ( token,*,iostat=ioerr ) D50
              if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"D50 granular dimension",ninp,nout) ) return
              token = lcase(GetToken(ainp,3,ioerr))
              if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"EROSION MODEL",ninp,nout) ) return
              if ( token(1:7) == "shields" ) then
                erosionmodel = "shields "
              else if ( token(1:4) == "mohr" ) then
                erosionmodel = "mohr    "
              else
                ier = 6
                return
              end if
            end if
            call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
            read ( ainp,*,iostat=ioerr ) nitersol
            if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"N° ITERATION SOLID",ninp,nout) ) return
!AA504
          endif
!          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
!          read ( ainp,*,iostat=ioerr ) coes, cuin
!          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"COHESION & CUIN",ninp,nout) ) return                      
     case default
    end select
!
!.. check for the limited modules
!
    if (.not. Erosion_Module_Shields_Mohr .and. Granular_flows_options%ID_erosion_criterion > 1 .and. tipo == "granular") then
      tipo = " "
      erosionmodel = "  "
      write (nscr,"(1x,a)") " "
      write (nout,"(1x,a)") " "
      write (nscr,"(1x,a)") " >>WARNING! - The erosion module Shields/Mohr is not available."
      write (nout,"(1x,a)") " >>WARNING! - The erosion module Shields/Mohr is not available."
      ier = 6
      return
    end if
    if (.not. Diffusion_Module .and. codif /= zero) then
      codif = zero
      Settling = zero
      write (nscr,"(1x,a)") " "
      write (nout,"(1x,a)") " "
      write (nscr,"(1x,a)") " >>WARNING! - The diffusion module is not available."
      write (nout,"(1x,a)") " >>WARNING! - The diffusion module is not available."
      ier = 7
      return
    end if
    if (.not. Explosion_Module .and. gamma /= zero) then
      gamma = zero
      InitialIntEn = zero
      write (nscr,"(1x,a)") " "
      write (nout,"(1x,a)") " "
      write (nscr,"(1x,a)") " >>WARNING! - The esplosion module is not available. Gamma and Initial Internal Energy are setted to zero."
      write (nout,"(1x,a)") " >>WARNING! - The esplosion module is not available. Gamma and Initial Internal Energy are setted to zero."
      ier = 8
      return
    end if
    if (.not. MultiFluid_Module .and.  index > 1 .and. tipo /= "liquid  ") then
      write (nscr,"(1x,a)") " "
      write (nout,"(1x,a)") " "
      write (nscr,"(1x,a)") " >>WARNING! - Only single fluid simulation is available."
      write (nout,"(1x,a)") " >>WARNING! - Only single fluid simulation is available."
      ier = 9
      return
    end if
    if (.not. MoreFluids_Module .and. index > 1) then
      write (nscr,"(1x,a)") " "
      write (nout,"(1x,a)") " "
      write (nscr,"(1x,a)") " >>WARNING! - Only one fluid in the simulation is available."
      write (nout,"(1x,a)") " >>WARNING! - Only one fluid in the simulation is available."
      ier = 10
      return
    end if
!
!.. assign the values read
!
    if ( ncord > 0 ) then
       Med(index)%tipo         = tipo
       Med(index)%index        = index
       Med(index)%NIterSol     = nitersol
       Med(index)%den0         = den0
       Med(index)%eps          = eps
       Med(index)%celerita     = Dsqrt(eps/den0)
       Med(index)%alfaMon      = alfaMon
       Med(index)%betaMon      = betaMon
       Med(index)%visc         = visc
       Med(index)%mumx         = viscmx
       Med(index)%cons         = cons
       Med(index)%taucri       = taucri
       Med(index)%cuin         = cuin
       Med(index)%phi          = Dabs(phi)
       Med(index)%coes         = coes
       Med(index)%Cs           = Cs
       Med(index)%RoughCoef    = Rough
       Med(index)%D50          = D50
!AA504 start
       Med(index)%gran_vol_frac_max = (1.d0-porosity)
       Med(index)%d_90 = d_90
!AA504 end       
       Med(index)%modelloerosione = erosionmodel
!.. modello diffusione
       Med(index)%codif        = codif
       Med(index)%SettlingCoef = Settling
!.. modello diffusione
!.. modello esplosioni
       Med(index)%gamma        = gamma
       Med(index)%InitialIntEn = InitialIntEn
!.. modello esplosioni
       if ( nout > 0 ) then
          write (nout,"(1x,a,i3,1x,a)")  "Medium:.....................",Med(index)%index,"("//Med(index)%tipo//")"
          write (nout,"(1x,a,1p,e12.4)") "Density:....................",Med(index)%den0
          write (nout,"(1x,a,1p,e12.4)") "Comprimibility:.............",Med(index)%eps
          write (nout,"(1x,a,1p,e12.4)") "Celerity:...................",Med(index)%celerita
          write (nout,"(1x,a,1p,e12.4)") "Alpha Monaghan Coeff.:......",Med(index)%alfaMon
          write (nout,"(1x,a,1p,e12.4)") "Beta Monaghan Coeff.:.......",Med(index)%betaMon
          write (nout,"(1x,a,1p,e12.4)") "Dynamic Viscosity:..........",Med(index)%visc
          write (nout,"(1x,a,1p,e12.4)") "Max. Dynamic Viscosity:.....",Med(index)%mumx
          write (nout,"(1x,a,1p,e12.4)") "Tau Critical:...............",Med(index)%taucri
          write (nout,"(1x,a,1p,e12.4)") "Index curve:................",Med(index)%cuin
          write (nout,"(1x,a,1p,e12.4)") "Friction Angle:.............",Med(index)%phi
          write (nout,"(1x,a,1p,e12.4)") "Cohesion:...................",Med(index)%coes
          write (nout,"(1x,a,1p,e12.4)") "Consistency:................",Med(index)%cons
          write (nout,"(1x,a,1p,e12.4)") "Smagorinsky Constant:.......",Med(index)%Cs
          write (nout,"(1x,a,1p,e12.4)") "Roughness Coeff.:...........",Med(index)%RoughCoef
          write (nout,"(1x,a,1p,e12.4)") "D50:........................",Med(index)%D50
          write (nout,"(1x,a,1p,e12.4)") "D90:........................",Med(index)%d_90
          write (nout,"(1x,a,1p,e12.4)") "(1-prosity):................",Med(index)%gran_vol_frac_max
          write (nout,"(1x,a,1x,a)")     "Erosion Model:..............",Med(index)%modelloerosione
          write (nout,"(1x,a,1p,e12.4)") "Diffusion Coeff.:...........",Med(index)%codif       
          write (nout,"(1x,a,1p,e12.4)") "Settling Velocity Coeff.:...",Med(index)%SettlingCoef       
          write (nout,"(1x,a,1p,e12.4)") "Explosion Gamma Coeff.:.....",Med(index)%Gamma       
          write (nout,"(1x,a,1p,e12.4)") "Initial Internal Energy.:...",Med(index)%InitialIntEn       
          write (nout,"(1x,a)")  " "
       end if
       Med(index)%visc  = visc / den0
       Med(index)%numx  = viscmx / den0
!???       Med(index)%taucri= taucri / den0
       Med(index)%phi   = Med(index)%phi * PIGRECO / 180.  !/ 57.29583 ! da gradi sessagesimali a radianti
!???       Med(index)%cons  = cons / den0
!???       Med(index)%coes  = coes / den0
!???       Med(index)%codif = codif  / den0
    end if           

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"MEDIUM DATA",ninp,nout) ) return

 end do

return
end subroutine ReadInputMedium
!---split

!AA501b (whole subroutine)
!cfile ReadBodyDynamics.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ReadInputFile.f90
!
! Creation      : 13nov12 Amicarelli-Agate
!
!************************************************************************************
! Module purpose : Reading input for body dynamics
!
! Calling routines: ReadInput
!
! Called routines: /
!
!************************************************************************************
subroutine ReadBodyDynamics (ainp,comment,nrighe,ier,ninp,nout)

! modules
use GLOBAL_MODULE                            
use AdM_USER_TYPE
use ALLOC_Module

!Declarations
implicit none

integer(4)    :: nrighe,ier,ninp,nout,ioerr,i,Id_body,n_elem,j,Id_elem,imposed_kinematics,n_records,Ic_imposed !,k
character( 1) :: comment
character(80) :: ainp,lcase !,token,GetToken
double precision :: mass                          
double precision :: L_geom(3),x_CM(3),alfa(3),u_CM(3),omega(3),mass_deact(6),Ic(3,3),x_rotC(3)
integer(4) :: normal_act(6)

logical,       external :: ReadCheck

! in case of restart the cards are not read
if (restart) then
  do while (TRIM(lcase(ainp)) /= "##### end body dynamics #####")
    call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
    if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BODY DYNAMICS DATA",ninp,nout)) return
  end do
  return
end if

call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BODY DYNAMICS DATA",ninp,nout)) return

 do while (TRIM(lcase(ainp)) /= "##### end body dynamics #####")
!Reading the number of bodies and the ratio between fluid and body particle size
    read (ainp,*,iostat=ioerr) n_bodies,dx_dxbodies,imping_body_grav
    if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BODY DYNAMICS GENERAL INPUT",ninp,nout)) return
! Writing the number of bodies and dx_dxbodies on the log file
!AA504 sub
    if ((ncord>0).and.(nout > 0)) then
       write (nout,"(1x,a,1p,i12)")   "n_bodies:.....................",n_bodies
       write (nout,"(1x,a,1p,e12.4)") "dx_dxbodies:..................",dx_dxbodies   
       write (nout,"(1x,a,1p,i12)")   "imping_body_grav:.............",imping_body_grav
       write (nout,"(1x,a)")  " "
    end if
! Allocation of the array of the bodies
    if (allocated(body_arr)) then
    else
    allocate(body_arr(n_bodies))  
    endif
! Loop over the transported bodies
    do i=1,n_bodies
!Reading the body parameters
       call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) Id_body,n_elem
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ID_BODY-N_ELEM",ninp,nout) ) return
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) mass
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"MASS",ninp,nout)) return
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) x_CM(1),x_CM(2),x_CM(3)
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"X_CM",ninp,nout)) return
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) Ic_imposed
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"IC_IMPOSED",ninp,nout)) return
       if (Ic_imposed == 1) then
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) Ic(1,1),Ic(1,2),Ic(1,3)
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"IC(1,1-3)",ninp,nout)) return
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) Ic(2,1),Ic(2,2),Ic(2,3)
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"IC(2,1-3)",ninp,nout)) return
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) Ic(3,1),Ic(3,2),Ic(3,3)
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"IC(3,1-3)",ninp,nout)) return
       endif
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) alfa(1),alfa(2),alfa(3)
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ALFA",ninp,nout)) return
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) x_rotC(1),x_rotC(2),x_rotC(3)
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"X_ROTC",ninp,nout)) return
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) u_CM(1),u_CM(2),u_CM(3)
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"U_CM",ninp,nout)) return
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) omega(1),omega(2),omega(3)
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"OMEGA",ninp,nout)) return
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) imposed_kinematics,n_records
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BODY_KINEMATICS",ninp,nout)) return       
!Assignation to the body parameters 
       body_arr(Id_body)%n_elem = n_elem
       body_arr(Id_body)%mass = mass
       body_arr(Id_body)%x_CM = x_CM
       body_arr(Id_body)%Ic_imposed = Ic_imposed
       if (Ic_imposed == 1) then
          body_arr(Id_body)%Ic = Ic
          else
             body_arr(Id_body)%Ic = 0.
       endif   
       body_arr(Id_body)%alfa = alfa
       body_arr(Id_body)%x_rotC = x_rotC
       body_arr(Id_body)%u_CM = u_CM
       body_arr(Id_body)%omega = omega
       body_arr(Id_body)%imposed_kinematics = imposed_kinematics
       body_arr(Id_body)%n_records = n_records
!Writing on the log file
       if ( ncord > 0 ) then
          if ( nout > 0 ) then
             write (nout,"(1x,a,1p,i12)")    "body:.......................",Id_body
             write (nout,"(1x,a,1p,e12.4)")  "mass:.......................",mass
             write (nout,"(1x,a,1p,3e12.4)") "x_CM:.......................",x_CM
             write (nout,"(1x,a,1p,i12)")    "IC_imposed:.................",Ic_imposed
             if (Ic_imposed == 1) then
             write (nout,"(1x,a,1p,3e12.4)") "Ic(1,1-3):..................",Ic(1,1),Ic(1,2),Ic(1,3)
             write (nout,"(1x,a,1p,3e12.4)") "Ic(2,1-3):..................",Ic(2,1),Ic(2,2),Ic(2,3)
             write (nout,"(1x,a,1p,3e12.4)") "Ic(3,1-3):..................",Ic(3,1),Ic(3,2),Ic(3,3)
             endif
!AA504             
             write (nout,"(1x,a,1p,3e12.4)") "alfa:.......................",alfa(1),alfa(2),alfa(3) 
             write (nout,"(1x,a,1p,3e12.4)") "x_rotC:.....................",x_rotC 
             write (nout,"(1x,a,1p,3e12.4)") "u_CM:.......................",u_CM
             write (nout,"(1x,a,1p,3e12.4)") "omega:......................",omega
             write (nout,"(1x,a,1p,i12)")    "imposed_kinematics:.........",imposed_kinematics    
             write (nout,"(1x,a,1p,i12)")    "n_records:..................",n_records
             write (nout,"(1x,a)")  " "
          end if
       endif
!Allocating body elements
       if (ncord>0) then
          else
          allocate(body_arr(Id_body)%elem(body_arr(Id_body)%n_elem))
          allocate(body_arr(Id_body)%body_kinematics(n_records,7))
       endif
!Reading the eventual imposed kinematics
       do j=1,n_records
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) body_arr(Id_body)%body_kinematics(j,:) 
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BODY_KINEMATICS_RECORDS",ninp,nout)) return            
       enddo       
! Reading element parameters                     
       do j=1,n_elem
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) Id_elem
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ID_ELEM",ninp,nout)) return
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) L_geom(1),L_geom(2),L_geom(3)
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"L_GEOM",ninp,nout)) return
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) x_CM(1),x_CM(2),x_CM(3)
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"X_CM_ELEM",ninp,nout)) return
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) alfa(1),alfa(2),alfa(3)
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ALFA_ELEM",ninp,nout)) return
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) normal_act(1),normal_act(2),normal_act(3),normal_act(4),normal_act(5),normal_act(6)
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"NORMAL_ACT",ninp,nout)) return
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) mass_deact(1),mass_deact(2),mass_deact(3),mass_deact(4),mass_deact(5),mass_deact(6)
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"MASS_DEACT",ninp,nout)) return
!Assignation of the element parameters
          body_arr(Id_body)%elem(Id_elem)%L_geom = L_geom
          body_arr(Id_body)%elem(Id_elem)%x_CM = x_CM
          body_arr(Id_body)%elem(Id_elem)%alfa = alfa
          body_arr(Id_body)%elem(Id_elem)%normal_act = normal_act
          body_arr(Id_body)%elem(Id_elem)%mass_deact = mass_deact          
!Writing on the log file
          if ( ncord > 0 ) then
             if ( nout > 0 ) then
                write (nout,"(1x,a,1p,i12)")    "element:....................",Id_elem
                write (nout,"(1x,a,1p,3e12.4)") "L_geom:.....................",L_geom
                write (nout,"(1x,a,1p,3e12.4)") "x_CM_elem:..................",x_CM
!AA504                
                write (nout,"(1x,a,1p,3e12.4)") "alfa_elem:..................",alfa(1),alfa(2),alfa(3)    
                write (nout,"(1x,a,1p,6i12)")   "normal_act:.................",normal_act    
                write (nout,"(1x,a,1p,6e12.4)") "mass_deact:.................",mass_deact                                          
                write (nout,"(1x,a)")  " "
             end if
          endif          
       end do
    enddo         
    call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
    if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BODY DYNAMICS DATA",ninp,nout)) return
 end do

return
end subroutine ReadBodyDynamics
!---split

!AA504 (whole subroutine)
!cfile ReadBedLoadTransport.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Subroutine name : ReadBedLoadTransport
!
! Creation        : 20Aug13 Amicarelli-Agate
!
!************************************************************************************
! Module purpose : Reading input for sediment transport
!
! Calling routines: ReadInput
!
! Called routines: /
!
!************************************************************************************
subroutine ReadBedLoadTransport (ainp,comment,nrighe,ier,ninp,nout,nscr)

! modules
use GLOBAL_MODULE                            
use AdM_USER_TYPE
use ALLOC_Module

!Declarations
implicit none
!AA504 sub
integer(4) :: nrighe,ier,ninp,nout,nscr,ioerr,ID_erosion_criterion,ID_main_fluid,ID_granular,monitoring_lines,i,line_ID,n_max_iterations,erosion_flag,viscosity_blt_formula
!AA504
integer(4) :: deposition_at_frontiers,Gamma_slope_flag
!AA504 sub
double precision :: dt_out,x_fixed,y_fixed,conv_crit_erosion,velocity_fixed_bed,Chezy_friction_coeff,x_min_dt,x_max_dt,y_min_dt,y_max_dt,z_min_dt,z_max_dt
character( 1) :: comment
character(80) :: ainp,lcase 

!External functions
logical, external :: ReadCheck

! in case of restart the cards are not read
if (restart) then
  do while (TRIM(lcase(ainp)) /= "##### end bed load transport #####")
    call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
    if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BED LOAD TRANSPORT DATA",ninp,nout)) return
  end do
  return
end if

call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BED LOAD TRANSPORT DATA",ninp,nout)) return

do while (TRIM(lcase(ainp)) /= "##### end bed load transport #####")
!Reading input parameters (first part)
   read (ainp,*,iostat=ioerr) ID_erosion_criterion,ID_main_fluid,ID_granular
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BED LOAD TRANSPORT GENERAL INPUT",ninp,nout)) return
   if (ID_erosion_criterion>0) then
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
      read (ainp,*,iostat=ioerr) velocity_fixed_bed,erosion_flag
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VELOCITY FIXED BED, EROSION FLAG",ninp,nout)) return
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)   
      read (ainp,*,iostat=ioerr) viscosity_blt_formula,deposition_at_frontiers,Gamma_slope_flag
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VISCOSITY BLT FORMULA",ninp,nout)) return
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)   
      read (ainp,*,iostat=ioerr) monitoring_lines,dt_out,conv_crit_erosion,n_max_iterations
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BED LOAD TRANSPORT MONITORING LINES",ninp,nout)) return
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)   
      read (ainp,*,iostat=ioerr) Chezy_friction_coeff
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CHEZY FRICTION COEFFICIENT",ninp,nout)) return 
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)   
      read (ainp,*,iostat=ioerr) x_min_dt,x_max_dt
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"X_MIN_DT,X_MAX_DT",ninp,nout)) return      
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)   
      read (ainp,*,iostat=ioerr) y_min_dt,y_max_dt
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"Y_MIN_DT,Y_MAX_DT",ninp,nout)) return 
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)   
      read (ainp,*,iostat=ioerr) z_min_dt,z_max_dt
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"Z_MIN_DT,Z_MAX_DT",ninp,nout)) return       
   endif
!
!.. check for the limited modules
!
    if (.not. Erosion_Module_Shields_Mohr .and. ID_erosion_criterion > 1) then
      write (nscr,"(1x,a)") " "
      write (nout,"(1x,a)") " "
      write (nscr,"(1x,a)") " >>WARNING! - The erosion module Shields/Mohr is not available."
      write (nout,"(1x,a)") " >>WARNING! - The erosion module Shields/Mohr is not available."
      ier = 6
      return
    end if
    if (.not. Granular_flux .and. erosion_flag /= 1) then
      write (nscr,"(1x,a)") " "
      write (nout,"(1x,a)") " "
      write (nscr,"(1x,a)") " >>WARNING! - The granular flux module in the simulation are not available."
      write (nout,"(1x,a)") " >>WARNING! - The granular flux module in the simulation are not available."
      ier = 11
      return
    end if
    if (.not. Erosion_Module_Shields_Seminara .and. ID_erosion_criterion == 1) then
      write (nscr,"(1x,a)") " "
      write (nout,"(1x,a)") " "
      write (nscr,"(1x,a)") " >>WARNING! - The erosion module Shields-VanRijn-Seminara in the simulation are not available."
      write (nout,"(1x,a)") " >>WARNING! - The erosion module Shields-VanRijn-Seminara in the simulation are not available."
      ier = 12
      return
    end if
!
! Writing input parameters (first part)
!
   if ((ncord>0).and.(nout > 0)) then
      select case (ID_erosion_criterion)
      case (0)
          write (nout,"(1x,a,1p,i12,a)")   "ID_erosion_criterion:.........",ID_erosion_criterion," no bed load transport."
      case (1)
          write (nout,"(1x,a,1p,i12,a)")   "ID_erosion_criterion:.........",ID_erosion_criterion," Shields-Seminara bed load transport."
      case (2)
          write (nout,"(1x,a,1p,i12,a)")   "ID_erosion_criterion:.........",ID_erosion_criterion," Shields bed load transport."
      case (3)
          write (nout,"(1x,a,1p,i12,a)")   "ID_erosion_criterion:.........",ID_erosion_criterion," Mohr-Coulomb bed load transport."
      end select
      if (ID_erosion_criterion>0) then      
         write (nout,"(1x,a,1p,i12)")   "ID_main_fluid:................",ID_main_fluid
         write (nout,"(1x,a,1p,i12)")   "ID_granular:..................",ID_granular
         write (nout,"(1x,a,1p,g12.5)") "velocity_fixed_bed:...........",velocity_fixed_bed    
         write (nout,"(1x,a,1p,i12)")   "erosion_flag:.................",erosion_flag
         write (nout,"(1x,a,1p,i12)")   "viscosity_blt_formula:........",viscosity_blt_formula
         write (nout,"(1x,a,1p,i12)")   "deposition_at_frontiers:......",deposition_at_frontiers
         write (nout,"(1x,a,1p,i12)")   "deposition_at_frontiers:......",Gamma_slope_flag         
         write (nout,"(1x,a,1p,i12)")   "monitoring_lines:.............",monitoring_lines
         write (nout,"(1x,a,1p,g12.5)") "dt_out:.......................",dt_out
         write (nout,"(1x,a,1p,g12.5)") "conv_crit_erosion:............",conv_crit_erosion
         write (nout,"(1x,a,1p,i12)")   "n_max_iterations:.............",n_max_iterations   
         if (viscosity_blt_formula==2) then
         write (nout,"(1x,a,1p,g12.5)")   "Chezy_friction_coeff:.........",Chezy_friction_coeff
         write (nout,"(1x,a,1p,g12.5)")   "x_min_dt:.....................",x_min_dt
         write (nout,"(1x,a,1p,g12.5)")   "x_max_dt:.....................",x_max_dt
         write (nout,"(1x,a,1p,g12.5)")   "y_min_dt:.....................",y_min_dt
         write (nout,"(1x,a,1p,g12.5)")   "y_max_dt:.....................",y_max_dt
         write (nout,"(1x,a,1p,g12.5)")   "z_min_dt:.....................",z_min_dt
         write (nout,"(1x,a,1p,g12.5)")   "z_max_dt:.....................",z_max_dt         
         endif
         write (nout,"(1x,a)")  " "
      endif
   end if
!Assignation to the body parameters (first part)
   Granular_flows_options%ID_erosion_criterion = ID_erosion_criterion
   if (ID_erosion_criterion>0) then      
      Granular_flows_options%ID_main_fluid = ID_main_fluid
      Granular_flows_options%ID_granular = ID_granular
      Granular_flows_options%velocity_fixed_bed = velocity_fixed_bed
      Granular_flows_options%erosion_flag = erosion_flag
      Granular_flows_options%viscosity_blt_formula = viscosity_blt_formula
      Granular_flows_options%deposition_at_frontiers = deposition_at_frontiers
      Granular_flows_options%Gamma_slope_flag = Gamma_slope_flag     
      Granular_flows_options%monitoring_lines = monitoring_lines 
      Granular_flows_options%dt_out = dt_out 
      Granular_flows_options%conv_crit_erosion = conv_crit_erosion
      if (viscosity_blt_formula==2) then
         Granular_flows_options%Chezy_friction_coeff = Chezy_friction_coeff  
         else
            Granular_flows_options%Chezy_friction_coeff = 0.d0
      endif
      Granular_flows_options%n_max_iterations = n_max_iterations
      Granular_flows_options%x_min_dt = x_min_dt
      Granular_flows_options%x_max_dt = x_max_dt
      Granular_flows_options%y_min_dt = y_min_dt
      Granular_flows_options%y_max_dt = y_max_dt
      Granular_flows_options%z_min_dt = z_min_dt
      Granular_flows_options%z_max_dt = z_max_dt      
! Allocation of the array of the monitoring lines
      if (allocated(Granular_flows_options%lines)) then
         else
            allocate(Granular_flows_options%lines(monitoring_lines,2)) 
!Initializing the auxiliary variable to print results
            Granular_flows_options%it_out_last = 0
      endif
! Loop over the monitoring lines
      do i=1,monitoring_lines
         call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
         read (ainp,*,iostat=ioerr) line_ID
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BED LOAD TRANSPORT MONITORING LINES",ninp,nout)) return      
         call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
         read (ainp,*,iostat=ioerr) x_fixed,y_fixed
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BED LOAD TRANSPORT MONITORING LINES",ninp,nout)) return
!Assignation to the section parameters 
         Granular_flows_options%lines(i,1) = x_fixed
         Granular_flows_options%lines(i,2) = y_fixed
!Writing on the log file
         if ( ncord > 0 ) then
            if ( nout > 0 ) then
               write (nout,"(1x,a,i12)")       "ID_line:....................",i
               write (nout,"(1x,a,1p,2e12.4)") "x_fixed,y_fixed:............",Granular_flows_options%lines(i,1),Granular_flows_options%lines(i,2)
               write (nout,"(1x,a)")  " "
            end if
         endif
      enddo  
   endif   
!Reading last line 
   call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BED LOAD TRANSPORT DATA",ninp,nout)) return
end do
 
return
end subroutine ReadBedLoadTransport
!---split

!cfile ReadInputOutputRegulation.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputOutputRegulation ( Med,ainp,comment,nrighe,ier,ninp,nout )

use GLOBAL_MODULE                                      
use AdM_USER_TYPE

implicit none

type (TyMedium), dimension(NMedium) :: Med

integer(4)    :: nrighe,ier, ninp,nout
character( 1) :: comment
character(80) :: ainp

integer(4)       :: iplot_fr, imemo_fr, irest_fr, icpoi_fr, ipllb_fr, ipllb_md
double precision ::  plot_fr,  memo_fr,  rest_fr,  cpoi_fr,  pllb_fr
integer(4)       :: n,ioutopt,ioutpo2,ioerr
character(80)    :: token

character( 7), dimension(3) :: outopt = (/ "full   ","partial","null   " /)

character(80), external :: GetToken,lcase
logical,       external :: ReadCheck

!AA504_B4
double precision :: depth_dt_out

 iplot_fr = 0
 plot_fr  = zero
 imemo_fr = 0
 memo_fr  = zero
 irest_fr = 0
 rest_fr  = zero
 icpoi_fr = 0
 cpoi_fr  = zero
 ipllb_fr = 0
 ipllb_md = 0
 pllb_fr  = zero
 ioutopt  = 0
 ioutpo2  = 3
 depth_dt_out = 0.0d0

 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"OUTPUT REGULATION DATA",ninp,nout) ) return

!do while ( lcase(ainp(1:33)) /= "##### end output regulation #####" )
 do while ( TRIM(lcase(ainp)) /= "##### end output regulation #####" )

    select case ( TRIM(lcase(GetToken(ainp,1,ioerr))) )

       case ("display")
          token = lcase(GetToken(ainp,2,ioerr))
          if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DISPLAY FREQUENCY STEP/TIME",ninp,nout) ) return
          if      ( token(1:4) == "step" ) then 
             token = lcase(GetToken(ainp,3,ioerr))
             if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DISPLAY FREQUENCY STEP value",ninp,nout) ) return
             read ( token,*,iostat=ioerr ) iplot_fr
             if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DISPLAY FREQUENCY STEP value",ninp,nout) ) return
          else if ( token(1:4) == "time" ) then 
             token = lcase(GetToken(ainp,3,ioerr))
             if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DISPLAY FREQUENCY TIME value",ninp,nout) ) return
             read ( token,*,iostat=ioerr )  plot_fr
             if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DISPLAY FREQUENCY TIME value",ninp,nout) ) return
          end if

       case ("results")
          token = lcase(GetToken(ainp,2,ioerr))
          if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESULTS SAVING FREQUENCY STEP/TIME",ninp,nout) ) return
          if      ( token(1:4) == "step" ) then 
             token = lcase(GetToken(ainp,3,ioerr))
             if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESULTS SAVING FREQUENCY STEP value",ninp,nout) ) return
             read ( token,*,iostat=ioerr ) imemo_fr
             if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESULTS SAVING FREQUENCY STEP value",ninp,nout) ) return
          else if ( token(1:4) == "time" ) then 
             token = lcase(GetToken(ainp,3,ioerr))
             if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESULTS SAVING FREQUENCY TIME value",ninp,nout) ) return
             read ( token,*,iostat=ioerr )  memo_fr
             if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESULTS SAVING FREQUENCY TIME value",ninp,nout) ) return
          end if

       case ("restart")
          token = lcase(GetToken(ainp,2,ioerr))
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"RESTART SAVING FREQUENCY STEP/TIME",ninp,nout) ) return
          if      ( token(1:4) == "step" ) then 
             token = lcase(GetToken(ainp,3,ioerr))
             if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART SAVING FREQUENCY STEP value",ninp,nout) ) return
             read ( token,*,iostat=ioerr ) irest_fr
             if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART SAVING FREQUENCY STEP value",ninp,nout) ) return
          else if ( token(1:4) == "time" ) then 
             token = lcase(GetToken(ainp,3,ioerr))
             if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART SAVING FREQUENCY TIME value",ninp,nout) ) return
             read ( token,*,iostat=ioerr )  rest_fr
             if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART SAVING FREQUENCY TIME value",ninp,nout) ) return
          end if

       case ("control")
          token = lcase(GetToken(ainp,2,ioerr))
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL POINTS SAVING FREQUENCY STEP/TIME",ninp,nout) ) return
          if      ( token(1:4) == "step" ) then 
             token = lcase(GetToken(ainp,3,ioerr))
             if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL POINTS SAVING FREQUENCY STEP value",ninp,nout) ) return
             read ( token,*,iostat=ioerr ) icpoi_fr
             if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL POINTS SAVING FREQUENCY STEP value",ninp,nout) ) return
          else if ( token(1:4) == "time" ) then 
             token = lcase(GetToken(ainp,3,ioerr))
             if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL POINTS SAVING FREQUENCY TIME value",ninp,nout) ) return
             read ( token,*,iostat=ioerr )  cpoi_fr
             if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL POINTS SAVING FREQUENCY TIME value",ninp,nout) ) return
          end if

       case ("level")
          do n = 1,2
             token = lcase(GetToken(ainp,2*n,ioerr))
             if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"LEVEL OPTION",ninp,nout) ) return
             if ( token(1:4) == "step" ) then 
                token = lcase(GetToken(ainp,2*n+1,ioerr))
                if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"LEVEL STEP FREQUENCY SAVING",ninp,nout) ) return
                read ( token,*,iostat=ioerr ) ipllb_fr
                if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"LEVEL STEP FREQUENCY SAVING",ninp,nout) ) return
             else if ( token(1:4) == "time" ) then 
                token = lcase(GetToken(ainp,2*n+1,ioerr))
                if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"LEVEL TIME FREQUENCY SAVING",ninp,nout) ) return
                read ( token,*,iostat=ioerr )  pllb_fr
                if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"LEVEL TIME FREQUENCY SAVING",ninp,nout) ) return
             else if ( token(1:6) == "medium" ) then 
                token = lcase(GetToken(ainp,2*n+1,ioerr))
                if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"LEVEL MEDIUM INDEX",ninp,nout) ) return
                read ( token,*,iostat=ioerr ) ipllb_md
                if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"LEVEL MEDIUM INDEX",ninp,nout) ) return
                if (ipllb_md == 0) then
                  ioerr = 5
                  if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"LEVEL MEDIUM INDEX",ninp,nout) ) return
                end if
             else 
                ier = 5
                return
             end if
          end do

!AA504_B4 start
       case ("depth")
          token = lcase(GetToken(ainp,2,ioerr))
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DEPTH OPTION",ninp,nout) ) return
          if ( token(1:4) == "dt_o" ) then 
             token = lcase(GetToken(ainp,2+1,ioerr))
             if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DEPTH DT_OUT SAVING",ninp,nout) ) return
             read ( token,*,iostat=ioerr ) depth_dt_out
             if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DEPTH DT_OUT SAVING",ninp,nout) ) return
          end if
!A504_B4 end          
          
       case ("print")
          token = lcase(GetToken(ainp,2,ioerr))
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"PRINT MODE AND FREQUENCY FULL/MINIMUM/NULL",ninp,nout) ) return
          if      ( token(1:4) == "full" ) then 
             token = lcase(GetToken(ainp,3,ioerr))
             if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"PRINT MINIMUM FREQUENCY step",ninp,nout) ) return
             read ( token,*,iostat=ioerr )  ioutopt
             if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"PRINT MINIMUM FREQUENCY step",ninp,nout) ) return
             ioutopt =-1*abs(ioutopt)
             ioutpo2 = 1
          else if ( token(1:7) == "partial" ) then 
             token = lcase(GetToken(ainp,3,ioerr))
             if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"PRINT MINIMUM FREQUENCY step",ninp,nout) ) return
             read ( token,*,iostat=ioerr )  ioutopt
             if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"PRINT MINIMUM FREQUENCY step",ninp,nout) ) return
             ioutpo2 = 2
          else if ( token(1:4) == "null" ) then 
             ioutopt = 0
             ioutpo2 = 3
          end if

       case default

          ier = 4
          return

    end select

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"OUTPUT REGULATION DATA",ninp,nout) ) return

 end do

 if ( ncord > 0 ) then
    Domain%iplot_fr = iplot_fr
    Domain%plot_fr  = plot_fr
    Domain%imemo_fr = imemo_fr
    Domain%memo_fr  = memo_fr
    Domain%irest_fr = irest_fr
    Domain%rest_fr  = rest_fr
    Domain%icpoi_fr = icpoi_fr
    Domain%cpoi_fr  = cpoi_fr
    Domain%ipllb_md = ipllb_md
    Domain%ipllb_fr = ipllb_fr
    Domain%pllb_fr  = pllb_fr
!AA504_B4 start
    Domain%depth_dt_out = depth_dt_out
    Domain%depth_it_out_last = 0.d0    
!AA504_B4 end    
    Domain%ioutopt  = ioutopt
    if ( nout > 0 ) then 
       write (nout,"(1x,a,i12)")   "Displaying     step frequency : ",Domain%iplot_fr
       write (nout,"(1x,a,e12.4)") "Displaying     time frequency : ",Domain%plot_fr
       write (nout,"(1x,a,i12)")   "Results saving step frequency : ",Domain%imemo_fr
       write (nout,"(1x,a,e12.4)") "Results saving time frequency : ",Domain%memo_fr
       write (nout,"(1x,a,i12)")   "Restart saving step frequency : ",Domain%irest_fr
       write (nout,"(1x,a,e12.4)") "Restart saving time frequency : ",Domain%rest_fr
       write (nout,"(1x,a,i12)")   "Ctrl.Points    step frequency : ",Domain%icpoi_fr
       write (nout,"(1x,a,e12.4)") "Ctrl.Points    time frequency : ",Domain%cpoi_fr
       write (nout,"(1x,a,i12)")   "Level Medium Index            : ",Domain%ipllb_md
       if ( Domain%ipllb_md > 0 ) &
         write (nout,"(1x,a,e12.4)") "Level Medium Density Limit    : ",med(Domain%ipllb_md)%den0 * half
       write (nout,"(1x,a,i12)")   "Level          step frequency : ",Domain%ipllb_fr
       write (nout,"(1x,a,e12.4)") "Level          time frequency : ",Domain%pllb_fr
!AA504_B4 
       write (nout,"(1x,a,e12.4)") "Depth          dt_out(s)      : ",Domain%depth_dt_out
       write (nout,"(1x,a,a)")     "Printing Option               : ",outopt(ioutpo2)
       if ( Domain%ioutopt > 0 ) &
         write (nout,"(1x,a,i12)")   "Printing       step frequency : ",iabs(Domain%ioutopt)
       write (nout,"(1x,a)") " "
    end if
 end if

return
end subroutine ReadInputOutputRegulation
!---split

!cfile ReadInputParticlesData.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
  subroutine ReadInputParticlesData ( NumberEntities, &
                                      Medium,icolor,bends,move,slip,npointv,valuev,values3, &
                                      pressu,valp,ainp,comment,nrighe,ier,ninp,nout )
!
  use GLOBAL_MODULE
  use AdM_USER_TYPE
!
  implicit none
!
  integer(4),        dimension(20)          :: NumberEntities
  integer(4)    :: nrighe,ier, ninp,nout
  character( 1) :: comment
!  character( 8) :: label
  character(80) :: ainp
!  character(4)  :: tipo
!
  integer(4)    :: ioerr
  character(80) :: token
!
  integer(4)    :: Medium
  integer(4)    :: i,n, icord
  integer(4)    :: npointv, icolor
!  double precision       :: trampa
  double precision       :: valp
  character(1)  :: bends
  character(1)  :: slip
  character(3)  :: move
  character(2)  :: pressu
  character(6)  :: token_color
  double precision, dimension(3) :: values3
  double precision, dimension(0:3,maxpointsvlaw) :: valuev
!
  character(80), external :: lcase, GetToken
  logical,       external :: ReadCheck

!
  call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
  read ( ainp,*,iostat=ioerr ) Medium
  if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"MEDIUM INDEX",ninp,nout) ) return
!
  call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
  token = GetToken(ainp,1,ioerr)
  if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"COLOURING TYPE",ninp,nout) ) return
  bends = lcase(token(1:1))
  token = GetToken(ainp,2,ioerr)
!read ( token,*,iostat=ioerr ) icolor
!
!salva colore per successiva inversione codice
  token_color(1:2) = token(5:6)
  token_color(3:4) = token(3:4)
  token_color(5:6) = token(1:2) 
  read ( token_color,'(Z6)',iostat=ioerr) icolor
  if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"COLOUR",ninp,nout) ) return
!
  call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
  move = GetToken(ainp,1,ioerr)
  if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL MOVE STATUS TYPE",ninp,nout) ) return
  values3(:) = zero
!
  slip = " "                             
  npointv = 0
!
  select case (lcase(move))

    case ("std")
       npointv = 0
!       do n = 1, NumberEntities(1)
       do n = 1, 3
          token = GetToken(ainp,(n+1),ioerr)
          read ( token,*,iostat=ioerr ) values3(n)
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL VELOCITY",ninp,nout) ) return
       end do
!       if (tipo == 'sour' .or. tipo == 'velo' .or. tipo == 'flow') then
!         token = GetToken(ainp,(3+2),ioerr)
!         read ( token,*,iostat=ioerr ) trampa
!         if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"TEMPO di RAMPA",ninp,nout) ) return
!       end if

    case ("fix")
      !NumberEntities(17) = NumberEntities(17) + 1
       npointv = 1
       valuev  = zero  
       do n = 1, NumberEntities(1)
          token = GetToken(ainp,(n+1),ioerr)
          read ( token,*,iostat=ioerr ) values3(n)
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL VELOCITY",ninp,nout) ) return
       end do
      !Aderenza/Free Slip
       token = GetToken(ainp,(NumberEntities(1)+2),ioerr)
       if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"NO_SLIP/FREE_SLIP",ninp,nout) ) return
       if      ( lcase(token(1:7)) == "no_slip"   ) then
          slip = "n"
       else if ( lcase(token(1:9)) == "free_slip" ) then
          slip = "f"
       else if ( lcase(token(1:9)) == "cont_slip" ) then  ! 20051212 +2
          slip = "c"
          NumberEntities(19) = 1
       else
          slip = "?"
          ioerr= 1
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"NO_SLIP/FREE_SLIP/CONT_SLIP",ninp,nout) ) return
       end if

    case ("law")
       ioerr   = 0
       npointv = 0
       valuev  = zero  
       token = GetToken(ainp,2,ioerr)
       read ( token,*,iostat=ioerr ) npointv
       if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"POINTS OF VELOCITY LAW",ninp,nout) ) return
      
       slip = "n"

       do i = 1, npointv
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VALUES OF VELOCITY LAW",ninp,nout) ) return
          if ( ioerr /= 0 .OR. ncord <= 0 ) cycle
          do n = 0, NumberEntities(1)
             icord = icoordp(n,ncord-1)
             token = GetToken(ainp,(n+1),ioerr)
             if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VALUES OF VELOCITY LAW",ninp,nout) ) return
             read ( token,*,iostat=ioerr ) valuev(n,i)
          end do
          if ( i == 1 ) values3(1:3) = valuev(1:3,1)
       end do

    case default

       if ( nout > 0 ) write (nout,*) "Unknown option: ",trim(ainp)
       stop

  end select
!
!.. loads the pressure type and values for all the boundaries
!
  call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
  pressu = GetToken(ainp,1,ioerr)
  if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL PRESSURE TYPE",ninp,nout) ) return
!
!.. the pressure value is assigned ("pa" type)
!
  if ( pressu == "pa" ) then  
!    
    token = lcase(GetToken(ainp,2,ioerr))
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"PRESSURE VALUES",ninp,nout) ) return
    read ( token,*,iostat=ioerr ) valp
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"PRESSURE VALUES",ninp,nout) ) return
!
!.. the pressure is evaluated as hydrostatic with respect a reference piezo line ("qp" type)
!
  else if ( pressu == "qp" ) then 
!     
    token = lcase(GetToken(ainp,2,ioerr))
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL PIEZO LINE",ninp,nout) ) return
    read ( token,*,iostat=ioerr ) valp
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL PIEZO LINE",ninp,nout) ) return
!
!.. the pressure is evaluated as hydrostatic with respect the maximum level of an assigned medium ("pl" type)
!
  else if ( pressu == "pl" ) then      
    token = lcase(GetToken(ainp,2,ioerr))
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FREE LEVEL LINE",ninp,nout) ) return
    read ( token,*,iostat=ioerr ) valp
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FREE LEVEL LINE",ninp,nout) ) return
  else
    if ( nout > 0 ) write (nout,*) "Unknown option: ",trim(ainp)
    stop
  end if
!
  if ( ncord > 0 ) then
   !controllo bends
    if ( bends == "b" .AND. icolor > 5 ) then
       icolor = 5
       if ( nout > 0 ) write(nout,*) "Maximum number of bends is 5!"
    end if
  end if
!
  return
!
  end subroutine ReadInputParticlesData
!---split

!cfile ReadInputRestart.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputRestart ( ainp,comment,nrighe,ier,ninp,nout )

use GLOBAL_MODULE                             
use AdM_USER_TYPE

implicit none

integer(4)    :: nrighe,ier, ninp,nout
logical       :: restartOK
character( 1) :: comment
character(80) :: ainp

integer(4)    :: ioerr
character(80) :: token

character(80), external :: lcase, GetToken
logical,       external :: ReadCheck

!
!.. in case of restart the cards are not read
!
!  if (restart) then
!    do while ( TRIM(lcase(ainp)) /= "##### end restart #####" )
!      call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
!      if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART DATA",ninp,nout) ) return
!    end do
!    return
!  end if
!
 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART DATA",ninp,nout) ) return

 do while ( TRIM(lcase(ainp)) /= "##### end restart #####" )

    select case ( TRIM(lcase(GetToken(ainp,1,ioerr))) )

       case ("step")
          token = lcase(GetToken(ainp,2,ioerr))
          if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART DATA STEP value",ninp,nout) ) return
          read ( token,*,iostat=ioerr ) Domain%istart
          if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART DATA STEP value",ninp,nout) ) return
          if ( ncord > 0 .AND. nout > 0 ) then
             write (nout,"(1x,a,i12)") "Restart from step: ",Domain%istart
             if ( Domain%istart < 0 ) write (nout,"(1x,a)") "Negative restart step!"
            !resta attiva solo l'ultima opzione letta
             if ( Domain%start > zero ) then
                write (nout,"(1x,a,f20.12,a)") "Restart from time: ",Domain%start," option ignored!"
                Domain%start = zero
             end if
          end if

       case ("time")
          token = lcase(GetToken(ainp,2,ioerr))
          if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART DATA TIME value",ninp,nout) ) return
          read ( token,*,iostat=ioerr ) Domain%start
          if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART DATA TIME value",ninp,nout) ) return
          if ( ncord > 0 .AND. nout > 0 ) then
             write (nout,"(1x,a,f20.12)") "Restart from time: ",Domain%start
             if ( Domain%start < zero ) write (nout,"(1x,a)") "Negative restart time!"
            !resta attiva solo l'ultima opzione letta
             if ( Domain%istart > 0 ) then
                write (nout,"(1x,a,i12,a)") "Restart from step: ",Domain%istart," option ignored!"
                Domain%istart = 0 
             end if
          end if

       case default
          Domain%file = ainp
          if ( ncord > 0 .AND. nout > 0 ) then
             inquire ( file = Domain%file, exist = restartOK )
             if ( restartOK ) then
                write (nout,"(1x,3a)") "Restart file: ",trim(Domain%file)
             else
                write (nout,"(1x,3a)") "Restart file: ",trim(Domain%file)," not found!"
             end if
          end if

    end select

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART DATA",ninp,nout) ) return

 end do

return
end subroutine ReadInputRestart
!---split

!cfile ReadInputRunParameters.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ReadInputRunParameters
!
! Last updating : April 18, 2013
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Agate             2009-2011      varie
! 04  Agate             03/07/12       Module exclusion (Diffusion, Esplosion, Temporal schemes)
! 05  Amicarelli/Agate  18apr12        add flag for wall elements reordering and maximum number of neighbours read in input
!
!************************************************************************************
! Module purpose : Read parameters from input file
!
! Calling routine: readinput
!
! Called routines: none
!
!************************************************************************************
!
subroutine ReadInputRunParameters ( ainp,comment,nrighe,ier,ninp,nout,nscr )

use GLOBAL_MODULE                                   
use AdM_USER_TYPE

implicit none

integer(4)    :: nrighe,ier, ninp,nout,nscr
character( 1) :: comment
character(80) :: ainp
!
!.. Local Scalars ..
integer(4)       :: itmax
double precision :: tmax
!AA401 sub start
!double precision :: cote, TetaP, TetaV !, TetaX
!integer(4)       :: ioerr
!AA504sub start
double precision :: CFL, TetaP, TetaV, COEFNMAXPARTJ,COEFNMAXPARTI
integer(4)       :: ioerr, time_split, RKscheme, body_part_reorder,MAXCLOSEBOUNDFACES,MAXNUMCONVEXEDGES,GCBFVecDim_loc,density_thresholds
!AA504 sub end
character(1)     :: Psurf
character(80)    :: token
!
character(80), external :: lcase, GetToken
logical,       external :: ReadCheck

 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RUN PARAMETERS DATA",ninp,nout) ) return

 do while ( TRIM(lcase(ainp)) /= "##### end run parameters #####" )

    read ( ainp,*,iostat=ioerr ) tmax
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"MAX. TRANSIENT TIME & ITERATIONS",ninp,nout) ) return
    if ( ioerr == 0 ) then
      token = GetToken(ainp,2,ioerr)
      if ( ioerr == 0 ) then
        read (token,*,iostat=ioerr) itmax
      else
        itmax = 1000000000
        ioerr = 0
      end if
    end if

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
!AA401 sub
!    read ( ainp,*,iostat=ioerr ) cote, pesodt
!AA601 sub
    read ( ainp,*,iostat=ioerr ) CFL, time_split, RKscheme, pesodt, dt_alfa_Mon
!AA401 sub
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CFL",ninp,nout) ) return

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    read ( ainp,*,iostat=ioerr ) TetaP, TetaV
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"TETAP & TETAV",ninp,nout) ) return
    token = GetToken(ainp,3,ioerr)
    if ( ioerr == 0 ) then
      read (token,*,iostat=ioerr) Psurf
      Psurf = lcase(Psurf)
      if (Psurf /= 'o' .and. Psurf /= 's' .and. Psurf /= 'a') then
        write (nscr,"(1x,a)") "Error setting run parameters. SMOOTHING Pressure Surface not setted."
        write (nout,"(1x,a)") "Error setting run parameters. SMOOTHING Pressure Surface not setted."
        stop
      end if
    else
      if ( nout > 0 ) write (nout,*) "Unknown option: ",trim(ainp)," in run parameters."
      write (nscr,*) "Unknown option: ",trim(ainp)," in run parameters."
      stop
    end if

!AA504sub start    
    call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
    read (ainp,*,iostat=ioerr) COEFNMAXPARTI,COEFNMAXPARTJ,body_part_reorder
    if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"COEFNMAXPARTI and COEFNMAXPARTJ ",ninp,nout) ) return
    call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
    read (ainp,*,iostat=ioerr) MAXCLOSEBOUNDFACES,MAXNUMCONVEXEDGES
    if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"MAXCLOSEBOUNDFACES and MAXNUMCONVEXEDGES ",ninp,nout) ) return
    call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
    read (ainp,*,iostat=ioerr) GCBFVecDim_loc
    if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"GCBFVecDim_loc ",ninp,nout) ) return
    call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
    read (ainp,*,iostat=ioerr) density_thresholds
    if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DENSITY_THRESHOLDS ",ninp,nout) ) return
!AA504sub end

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RUN PARAMETERS DATA",ninp,nout) ) return

 end do
!
!.. check for the limited modules
!
    if (.not. TemporalScheme_Module .and. (time_split /= 1 .or. RKscheme /= 1)) then
      time_split = 1
      RKscheme = 1
      write (nscr,"(1x,a)") " "
      write (nout,"(1x,a)") " "
      write (nscr,"(1x,a)") " >>WARNING! - The temporal schemes module is not available. Only time_split=1 and RKscheme=1 is available."
      write (nout,"(1x,a)") " >>WARNING! - The temporal schemes module is not available. Only time_split=1 and RKscheme=1 is available."
      ier = 20
      return
    end if
!
!.. assign the values read
!
 if ( ncord > 0 ) then
    Domain%tmax  = tmax
    Domain%itmax = itmax
!AA401 sub start
!    Domain%cote  = cote
    Domain%CFL  = CFL
    Domain%time_split  = time_split
    if (time_split == 1) then
      Domain%RKscheme  = 1
    else
      if (RKscheme == 0) then
        write (nscr,"(1x,a)") " "
        write (nout,"(1x,a)") " "
        write (nscr,"(1x,a)") "Error setting run parameters. RKscheme must be greather than 0."
        write (nout,"(1x,a)") "Error setting run parameters. RKscheme must be greather than 0."
        stop
      end if
      Domain%RKscheme  = RKscheme
    end if
!AA401 sub end
    Domain%TetaP = TetaP
    Domain%TetaV = TetaV
    Domain%Psurf = Psurf

!AA504
    Domain%COEFNMAXPARTI = COEFNMAXPARTI

!AA503 start   
    Domain%COEFNMAXPARTJ = COEFNMAXPARTJ
    Domain%body_part_reorder = body_part_reorder
!AA503 end

!AA504 start
    Domain%MAXCLOSEBOUNDFACES = MAXCLOSEBOUNDFACES
    Domain%MAXNUMCONVEXEDGES = MAXNUMCONVEXEDGES 
    GCBFVecDim = GCBFVecDim_loc
    Domain%density_thresholds = density_thresholds
!AA504 end
    if ( nout > 0 ) then
       write (nout,"(1x,a,1p,e12.4)") "TMAX                       : ",Domain%tmax
       write (nout,"(1x,a,i12)")      "ITMAX                      : ",Domain%itmax
!AA401 sub start
!       write (nout,"(1x,a,1p,e12.4)") "COTE                       : ",Domain%cote
       write (nout,"(1x,a,1p,e12.4)") "CFL                        : ",Domain%CFL
       write (nout,"(1x,a,1p,i1)")    "staggering option          : ",Domain%time_split
       write (nout,"(1x,a,1p,i1)")    "RKscheme                   : ",Domain%RKscheme
!AA401 sub end
       write (nout,"(1x,a,1p,e12.4)") "SMOOTHING PRES             : ",Domain%TetaP
       write (nout,"(1x,a,1p,e12.4)") "SMOOTHING VEL              : ",Domain%TetaV
       if (Domain%Psurf == 'o') then
         write (nout,"(1x,a,a)") "SMOOTHING Pressure Surface : ","Original"
       else if (Domain%Psurf == 's') then
         write (nout,"(1x,a,a)") "SMOOTHING Pressure Surface : ","Delta Pressure from hydrostatic"
       else if (Domain%Psurf == 'a') then
         write (nout,"(1x,a,a)") "SMOOTHING Pressure Surface : ","Weight calculation with atmospheric pressure"
       end if

!AA504
       write (nout,"(1x,a,1p,e12.4)") "COEFNMAXPARTI              : ",Domain%COEFNMAXPARTI
       
!AA503       
       write (nout,"(1x,a,1p,e12.4)") "COEFNMAXPARTJ              : ",Domain%COEFNMAXPARTJ
       write (nout,"(1x,a,1p,i1)")    "body_part_reorder          : ",Domain%body_part_reorder
       
!AA504 start
       write (nout,"(1x,a,1p,i12)")   "MAXCLOSEBOUNDFACES         : ",Domain%MAXCLOSEBOUNDFACES
       write (nout,"(1x,a,1p,i12)")   "MAXNUMCONVEXEDGES          : ",Domain%MAXNUMCONVEXEDGES
       write (nout,"(1x,a,1p,i12)")   "GCBFVecDim_loc             : ",GCBFVecDim
       write (nout,"(1x,a,1p,i1)")    "density_thresholds         : ",Domain%density_thresholds       
!AA504 end       
       write (nout,"(1x,a)") " "
    end if
 end if

return
end subroutine ReadInputRunParameters
!---split

!cfile ReadInputTitle.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputTitle ( ainp,comment,nrighe,ier,ninp,nout )

use GLOBAL_MODULE                            

implicit none

integer(4)    :: nrighe,ier, ninp,nout
character( 1) :: comment
character(80) :: ainp

integer(4)    :: n,ioerr

character(80), external :: lcase
logical,       external :: ReadCheck

 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"TITLE DATA",ninp,nout) ) return 
 n = 0

 do while ( TRIM(lcase(ainp)) /= "##### end title #####" )
 
    n = n + 1
    if ( n <= maxtit ) then
       title(n) = ainp
       if ( ncord > 0 .AND. nout > 0 ) write (nout,"(1x,a)") title(n)
    end if

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"TITLE DATA",ninp,nout) ) return 

 end do
 if ( ncord > 0 .AND. nout > 0 ) write (nout,"(1x,a)") " "

return
end subroutine ReadInputTitle
!---split

!cfile ReadInputVertices.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputVertices ( NumberEntities,Vertice, &
                               ainp,comment,nrighe,ier,prtopt,ninp,nout )

use GLOBAL_MODULE                              
use AdM_USER_TYPE

implicit none

integer(4),      dimension(20)                    :: NumberEntities
double precision,dimension(1:SPACEDIM,NumVertici) :: Vertice

integer(4)    :: nrighe,ier, ninp,nout
logical(4)    :: prtopt
character( 1) :: comment
character(80) :: ainp

integer(4)    :: n,i,icord,ioerr
character(8)  :: label
double precision, dimension(3) :: values1

character(80), external :: lcase, GetToken
logical,       external :: ReadCheck

!
!.. in case of restart the cards are not read
!
  if (restart) then
    do while ( TRIM(lcase(ainp)) /= "##### end vertices #####" )
      call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
      if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"VERTICES DATA",ninp,nout) ) return
    end do
    return
  end if
!
 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"VERTICES DATA",ninp,nout) ) return

 if ( ncord > 0 .AND. nout > 0 .AND. prtopt ) then
    write (nout,"(1x,a)") "List of vertices:"
 end if

 do while ( TRIM(lcase(ainp)) /= "##### end vertices #####" )

    select case ( TRIM(Domain%tipo) )
!
!AA406 sub
       case ( "semi","bsph" ) 
!
          read ( ainp,*,iostat=ioerr ) i, values1(1:NumberEntities(1))

! arrotondamento a 1.0d-5
!if ( ncord > 0 ) then
!write (*,'(i10,2f15.10)') i,values1(1:NumberEntities(1))
!values1(1:NumberEntities(1)) = real((nint(values1(1:NumberEntities(1)) * 1.d5)) / 1e5)
!values1(1:NumberEntities(1)) = (anint(values1(1:NumberEntities(1)) * 1.0d5)) / 1.0d5
!write (*,'(i10,2f15.10)') NumberEntities(1),values1(1:NumberEntities(1))
!end if
!

          write (label,"(i8)") i
          if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"VERTEX n."//label,ninp,nout) ) return
!if ( i > 8 ) i = i - 3420
          NumberEntities(7) = max(i,NumberEntities(7))

          if ( ncord > 0 ) then
             do n = 1, NumberEntities(1)
                icord = icoordp(n,ncord-1)
                if ( NumberEntities(7) == 1 ) then
                   Domain%coord(icord,1) = values1(n)
                   Domain%coord(icord,2) = values1(n)
                end if
                Vertice(icord,i) = values1(n)
                Domain%coord(icord,1) = min(values1(n),Domain%coord(icord,1))
                Domain%coord(icord,2) = max(values1(n),Domain%coord(icord,2))
             end do
          end if

       case default

          if ( nout > 0 ) then
             write (nout,*) "Unknown Domain Type: ",Domain%tipo
          end if
          ier = 2
          return


    end select

    if ( ncord > 0 .AND. nout > 0 .AND. prtopt ) then
       write (nout,"(i6,1p,3(2x,a,e12.4))") &
       i,(xyzlabel(icoordp(n,ncord-1)),Vertice(icoordp(n,ncord-1),i),n=1,ncord)
    end if

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"VERTICES DATA",ninp,nout) ) return

 end do

 if ( ncord > 0 .AND. nout > 0 ) then

!   write (nout,*)
!   write (nout,"(1x,a)") "List of vertices:"
!   do n = 1, NumberEntities(7)
!      write (nout,"(i6,1p,3(2x,a,e12.4))") &
!      n,(xyzlabel(icoordp(i,ncord-1)),Vertice(icoordp(i,ncord-1),n),i=1,ncord)
!   end do

    do n = 1, NumberEntities(1)
       icord = icoordp(n,ncord-1)
       write (nout,"(1x,a,a,1p,e12.4)") xyzlabel(icord)," coordinate min. ",Domain%coord(icord,1)
       write (nout,"(1x,a,a,1p,e12.4)") xyzlabel(icord)," coordinate max. ",Domain%coord(icord,2)
    end do
    write (nout,"(1x,a)") " "

 end if

return
end subroutine ReadInputVertices
!---split

!cfile ReadRestartFile.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : readrestartfile
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
! Module purpose : Module to read results file for rstart purpose
!
! Calling routine: gest_input
!
! Called routines: diagnostic
!
!************************************************************************************

subroutine ReadRestartFile ( option, ier, nrecords )
!
!.. assign modules
use files_entities
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
character(7),intent(IN)  :: option
integer(4),intent(INOUT) :: ier
integer(4),intent(INOUT) :: nrecords 
!
!.. Local Scalars ..
!
integer(4)       :: restartcode, save_istart
integer(4)       :: ioerr
double precision :: save_start
character(12)    :: ainp = "Restart File"
character(len=8) :: versionerest
!
!.. External Routines ..
character(80), external :: lcase
logical,       external :: ReadCheck
!
!.. Executable Statements ..
!
 ier = 0
!
!.. Restart heading 
!
 if ( TRIM(lcase(option)) == TRIM(lcase("heading")) ) then
!
   rewind (nsav)
!
   write(nout,'(a)')    "-------------------"
   write(nout,"(1x,a)") ">> Restart heading."
   write(nout,'(a)')    "-------------------"
   write(nscr,'(a)')    "-------------------"
   write(nscr,"(1x,a)") ">> Restart heading."
   write(nscr,'(a)')    "-------------------"
!
   read (nsav,iostat=ioerr) versionerest,nrecords
   if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"versionerest,nrecords",nsav,nout) ) return
! check sulla versione
   if (TRIM(lcase(version)) /= TRIM(lcase(versionerest))) then
     write(nout,'(a)')    "---------------------------------------------------------------"
     write(nout,"(1x,a)") ">> ERROR! The Restart version is not equal the current version."
     write(nout,"(1x,a)") ">>        The Run is stopped."
     write(nout,'(a)')    "---------------------------------------------------------------"
     flush(nout)
     write(nscr,'(a)')    "---------------------------------------------------------------"
     write(nscr,"(1x,a)") ">> ERROR! The Restart version is not equal the current version."
     write(nscr,"(1x,a)") ">>        The Run is stopped."
     write(nscr,'(a)')    "---------------------------------------------------------------"
     flush(nscr)
     stop
   end if
!
   read (nsav,iostat=ioerr) Ncord, Nag, NMedium, NPartZone, &
                            NumVertici, NumFacce, NumTratti, NumBVertices, NumBSides, &
                            NPointst,NPoints,NPointsl,NPointse, NLines, NSections, GCBFVecDim, &
                            doubleh
   if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"ncord, nag, ...",nsav,nout) ) return
!
 else if ( TRIM(lcase(option)) == "reading" ) then
!
   write(nout,'(a)')    "-----------------------------------------------------------------------"
   write(nout,"(1x,a)") ">> Restart reading:  step          time      interval    num.particles"
   write(nscr,'(a)')    "-----------------------------------------------------------------------"
   write(nscr,"(1x,a)") ">> Restart reading:  step          time      interval    num.particles"
!
   save_istart = Domain%istart
   save_start  = Domain%start
!
   read (nsav,iostat=ioerr) domain
   if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"domain",nsav,nout) ) return
!
   read (nsav,iostat=ioerr) grid
   if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"grid",nsav,nout) ) return
!.. allocazione matrice 2d per calcolo pelolibero (caso erosione)
!AA504 sub (removed everywhere the fifth element of the array
  allocate (ind_interfaces(Grid%ncd(1),Grid%ncd(2),4), stat = ioerr)
  if (ioerr /= 0) then
!AA504 sub      
    write (nout,'(1x,a,i2)') "    Array ind_interfaces not allocated. Error code: ",ioerr
    stop ' routine ReadRestartFile'
  else
!AA504 sub      
    write (nout,'(1x,a)') "    Array ind_interfaces successfully allocated "
  end if
!
   read (nsav,iostat=ioerr) Med(1:NMedium)
   if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"Med",nsav,nout) ) return
!
   if ( NumVertici   > 0 ) then
     read (nsav,iostat=ioerr) Vertice(1:SPACEDIM,1:NumVertici)
     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"Vertice",nsav,nout) ) return
   end if
!
   if ( NumFacce     > 0 ) then 
     read (nsav,iostat=ioerr) BoundaryFace(1:NumFacce)
     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"BoundaryFace",nsav,nout) ) return
   end if
!
   if ( NumFacce     > 0 ) then
     read (nsav,iostat=ioerr) BFaceList(1:NumFacce)
     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"BFaceList",nsav,nout) ) return
   end if
!
   if ( NumTratti    > 0 ) then
     read (nsav,iostat=ioerr) Tratto(1:NumTratti)
     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"Tratto",nsav,nout) ) return
   end if
!
   if ( NPartZone    > 0 ) then
     read (nsav,iostat=ioerr) Partz(1:NPartZone)
     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"Partz",nsav,nout) ) return
   end if
!
   if ( NumBVertices > 0 ) then
     read (nsav,iostat=ioerr) BoundaryVertex(1:NumBVertices)
     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"BoundaryVertex",nsav,nout) ) return
   end if
!
   if ( NumBSides    > 0 ) then
     read (nsav,iostat=ioerr) BoundarySide(1:NumBSides)
     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"BoundarySide",nsav,nout) ) return
   end if
!
!   if ( NPointst     > 0 ) then
!     read (nsav,iostat=ioerr) control_points(1:NPointst)
!     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"control_points",nsav,nout) ) return
!   end if
!
!   if ( NLines       > 0 ) then
!     read (nsav,iostat=ioerr) control_lines(1:NLines)
!     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"control_lines",nsav,nout) ) return
!   end if
!
!   if ( NSections    > 0 ) then
!     read (nsav,iostat=ioerr) Control_Sections(0:NSections+1)
!     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"control_sections",nsav,nout) ) return
!   end if
!
!
!.. Restart positioning is based on the step number
!
   it_start = 0 
!
!   if ( Domain%istart > 0 ) then
   if ( save_istart > 0 ) then
!
!     do while ( Domain%istart > it_start )
     do while ( save_istart > it_start )
!
       read(nsav,iostat=ioerr) it_start,tempo,dt,nag,ncord,restartcode
       if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"it_start,tempo,dt,nag,ncord,restartcode",nsav,nout) ) return
!
       write(nout,"(16x,i10,2(2x,g12.5),7x,i10)") it_start,tempo,dt,nag; flush(nout)
       write(nscr,"(16x,i10,2(2x,g12.5),7x,i10)") it_start,tempo,dt,nag; flush(nscr)
!
!       if ( it_start < Domain%istart ) then
       if ( it_start < save_istart ) then
         read (nsav,iostat=ioerr) 
         if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"...",nsav,nout) ) return
!
       else
! leggo per restart
         if (restartcode == 1) then
           read (nsav,iostat=ioerr) pg(1:nag)
           if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"pg",nsav,nout) ) return
!
           write(nout,'(a)') " "
           write(nout,'(a,i10,a,g12.5)') "   Located Restart Step :",it_start,"   Time :",tempo; flush(nout)
           write(nscr,'(a)') " "
           write(nscr,'(a,i10,a,g12.5)') "   Located Restart Step :",it_start,"   Time :",tempo; flush(nscr)
! leggo per risultati
         else if (restartcode == 0) then
           read (nsav,iostat=ioerr) pg(1:nag)%coord(1),pg(1:nag)%coord(2),pg(1:nag)%coord(3), &
                                    pg(1:nag)%vel(1),  pg(1:nag)%vel(2),  pg(1:nag)%vel(3), &
                                    pg(1:nag)%pres,   &
                                    pg(1:nag)%dens,   &
                                    pg(1:nag)%mass,   &
                                    pg(1:nag)%visc,   &
                                    pg(1:nag)%IntEn,  &
                                    pg(1:nag)%VolFra, &
                                    pg(1:nag)%imed,   &
                                    pg(1:nag)%icol
!
           if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"pg",nsav,nout) ) return
!
           write(nout,'(a)') " "
           write(nout,'(a,i10,a,g12.5)') "   Located Result Step :",it_start,"   Time :",tempo; flush(nout)
           write(nout,'(a)') "       But this step is not a restart step. Check the correct step for restart in the restart file."; flush(nout)
           write(nout,'(a)') " The program is terminated."; flush(nout)
           write(nscr,'(a)') " "
           write(nscr,'(a,i10,a,g12.5)') "   Located Result Step :",it_start,"   Time :",tempo; flush(nscr)
           write(nscr,'(a)') "       But this step is not a restart step. Check the correct step for restart in the restart file."; flush(nscr)
           write(nscr,'(a)') " The program is terminated."; flush(nscr)
           stop
         end if
         return
       end if
     end do
     write(nout,'(a,i10,a)') "   Restart Step Number:",it_start," has not been found"
     write(nscr,'(a,i10,a)') "   Restart Step Number:",it_start," has not been found"
     ier = 3
!
!.. Restart positioning is based on the time step
!
   else if ( save_start > zero ) then
!
     tempo = zero
!
     do while ( save_start > tempo )
!
       read(nsav,iostat=ioerr) it_start,tempo,dt,nag,ncord,restartcode
       if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"it_start,tempo,dt,nag,ncord,restartcode",nsav,nout) ) return
!
       write(nout,"(16x,i10,2(2x,g12.5),7x,i10)") it_start,tempo,dt,nag; flush(nout)
       write(nscr,"(16x,i10,2(2x,g12.5),7x,i10)") it_start,tempo,dt,nag; flush(nscr)
!
       if ( tempo < Domain%start ) then
         read (nsav,iostat=ioerr) 
         if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"...",nsav,nout) ) return
       else
! leggo per restart
         if (restartcode == 1) then
           read (nsav,iostat=ioerr) pg(1:nag)
           if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"pg",nsav,nout) ) return
           write(nout,'(a)') 
           write(nout,'(a,i10,a,g12.5)') "   Located Restart Step :",it_start,"   Time :",tempo; flush(nout)
           write(nscr,'(a)') 
           write(nscr,'(a,i10,a,g12.5)') "   Located Restart Step :",it_start,"   Time :",tempo; flush(nscr)
! leggo per risultati
         else if (restartcode == 0) then
           read (nsav,iostat=ioerr) pg(1:nag)%coord(1),pg(1:nag)%coord(2),pg(1:nag)%coord(3), &
                                    pg(1:nag)%vel(1),  pg(1:nag)%vel(2),  pg(1:nag)%vel  (3), &
                                    pg(1:nag)%pres,   &
                                    pg(1:nag)%dens,   &
                                    pg(1:nag)%mass,   &
                                    pg(1:nag)%visc,   &
                                    pg(1:nag)%IntEn,  &
                                    pg(1:nag)%VolFra, &
                                    pg(1:nag)%imed,   &
                                    pg(1:nag)%icol
           if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"pg",nsav,nout) ) return
           write(nout,'(a)') 
           write(nout,'(a,i10,a,g12.5)') "   Located Result Time :",it_start,"   Time :",tempo; flush(nout)
           write(nout,'(a)') "       But this time is not a restart time. Check the correct time for restart in the restart file."; flush(nout)
           write(nout,'(a)') " The program is terminated."; flush(nout)
           write(nscr,'(a)') 
           write(nscr,'(a,i10,a,g12.5)') "   Located Result Time :",it_start,"   Time :",tempo; flush(nscr)
           write(nscr,'(a)') "       But this time is not a restart time. Check the correct time for restart in the restart file."; flush(nscr)
           write(nscr,'(a)') " The program is terminated."; flush(nscr)
           stop
         end if
         return
       end if
     end do
     write(nout,'(a,i10,a)') "   Restart Time Step:",Domain%start," has not been found"
     write(nscr,'(a,i10,a)') "   Restart Time Step:",Domain%start," has not been found"
     ier = 3
!
   else
     write (nout,'(a)' ) "  > Restart cannot be read at step:",it_start,"  time:",tempo
     write (nscr,'(a)' ) "  > Restart cannot be read at step:",it_start,"  time:",tempo
     ier = 4
   end if
!
   write (nout,'(a)' ) "  > Restart read successfully at step:",it_start,"  time:",tempo
   write (nscr,'(a)' ) "  > Restart read successfully at step:",it_start,"  time:",tempo
!
 else
!
   ier = 5
!
 end if
!
return
endsubroutine ReadRestartFile
!---split
!
!cfile ReadRiga.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadRiga ( ainp,comment,nrighe,ier,ninp )
implicit none
!
integer(4)   :: ier, ninp
character(1) :: comment
character(*) :: ainp
!
integer(4)   :: ioerr, n, l
integer(4)   :: nrighe
!
  ioerr = 0
!
  READ_LOOP: do while ( ioerr == 0 )
!
    read  ( ninp, "(a)", iostat=ioerr ) ainp
    nrighe = nrighe + 1
!
    if ( ioerr == 0 .AND. trim(ainp) /= "" ) then
!
!.. sostituisce i tabulatori con blank
      if ( ainp(1:1) /= comment ) then
        do n = 1,len(trim(ainp))
          if ( iachar(ainp(n:n)) == 9 ) ainp(n:n) = " "
        end do
        l = index(ainp,comment)
!
        if ( l > 0 ) then
          do n = l,len(trim(ainp))
            ainp(n:n)=" "
          end do
        end if
        exit READ_LOOP 
      end if
!
    end if
!
  end do  READ_LOOP
!
  ier = ioerr
!
return
end subroutine ReadRiga
!---split


!AA504 the whole subroutine
!cfile ReadSectionFlowRate.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: ReadSectionFlowRate
!
!************************************************************************************
! Module purpose : Input management for flow rate mnitoring sections 
!
! Calling routines: Gest_Input
!
! Called subroutines: ReadRiga
!
! Creation: Amicarelli-Agate 5July13
!
!************************************************************************************

subroutine ReadSectionFlowRate (ainp,comment,nrighe,ier,ninp,nout)

! modules
use GLOBAL_MODULE                            
use AdM_USER_TYPE
use ALLOC_Module

!Declarations
implicit none

integer(4)    :: nrighe,ier,ninp,nout,ioerr,i,n_sect,n_vertices,section_ID,n_fluid_types
character(1)  :: comment
character(80) :: ainp,lcase 
double precision :: dt_out,aux_dis,area
double precision :: plane_normal(3),vec_aux_1(3),vec_aux_2(3),vec_aux_3(3),vec_aux_4(3)
double precision :: vertex(4,3)
logical,external :: ReadCheck

! in case of restart the cards are not read
if (restart) then
  do while (TRIM(lcase(ainp)) /= "##### end Section_flow_rate #####")
    call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
    if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"Section_flow_rate DATA",ninp,nout)) return
  end do
  return
end if

call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"Section_flow_rate DATA",ninp,nout)) return

 do while (TRIM(lcase(ainp)) /= "##### end section flow rate #####")
!Reading the number of monitoring sections for the flow rate and their writing time step
    read (ainp,*,iostat=ioerr) n_sect,dt_out,n_fluid_types
    if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"Section_flow_rate GENERAL INPUT",ninp,nout)) return
! Writing the number of sections and the writing time step on the log file
    if (nout > 0) then
       write (nout,"(1x,a,1p,i12)")   "n_sect:.............................",n_sect
       write (nout,"(1x,a,1p,e12.4)") "dt_out:.............................",dt_out
       write (nout,"(1x,a,1p,i12)")   "n_sect:.............................",n_fluid_types
       write (nout,"(1x,a)")  " "
    end if
! Allocation of the array of the flow rate monitoring sections
    if (allocated(Q_sections%section)) then
       else
          allocate(Q_sections%section(n_sect)) 
          Q_sections%n_sect = n_sect
          Q_sections%dt_out = dt_out
          Q_sections%n_fluid_types = n_fluid_types
!Initializinf the auxiliary variable to print results
          Q_sections%it_out_last = 0
    endif
! Loop over the sections
    do i=1,n_sect
!Reading the section parameters
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) section_ID
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"section_ID",ninp,nout)) return
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) n_vertices
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"n_vertices",ninp,nout)) return
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) vertex(1,1),vertex(1,2),vertex(1,3)
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"vertex_1",ninp,nout)) return
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) vertex(2,1),vertex(2,2),vertex(2,3)
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"vertex_2",ninp,nout)) return
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) vertex(3,1),vertex(3,2),vertex(3,3)
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"vertex_3",ninp,nout)) return
       if (n_vertices==4) then
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) vertex(4,1),vertex(4,2),vertex(4,3)
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"vertex_4",ninp,nout)) return           
       endif
!Assignation to the section parameters 
       Q_sections%section(i)%n_vertices = n_vertices
       Q_sections%section(i)%vertex(:,:) = vertex(:,:)
! Computation of the section area       
       vec_aux_1(:) = vertex(1,:)
       vec_aux_2(:) = vertex(2,:)
       vec_aux_3(:) = vertex(3,:)
       vec_aux_4(:) = vertex(4,:)
       call area_quadrilateral(vec_aux_1,vec_aux_2,vec_aux_3,vec_aux_4,area)
       Q_sections%section(i)%area = area
       call dis_point_plane(vec_aux_1,vec_aux_1,vec_aux_2,vec_aux_3,aux_dis,plane_normal)
       Q_sections%section(i)%normal(:) = plane_normal(:)
!Writing on the log file
       if ( ncord > 0 ) then
          if ( nout > 0 ) then
             write (nout,"(1x,a,i12)")       "n_vertices:.................",n_vertices
             write (nout,"(1x,a,1p,e12.4)")  "area:.......................",Q_sections%section(i)%area
             write (nout,"(1x,a,1p,3e12.4)") "normal:.....................",Q_sections%section(i)%normal(1),Q_sections%section(i)%normal(2),Q_sections%section(i)%normal(3)
             write (nout,"(1x,a,1p,3e12.4)") "vertex(1,1-3):..............",vertex(1,1),vertex(1,2),vertex(1,3)
             write (nout,"(1x,a,1p,3e12.4)") "vertex(2,1-3):..............",vertex(2,1),vertex(2,2),vertex(2,3)
             write (nout,"(1x,a,1p,3e12.4)") "vertex(3,1-3):..............",vertex(3,1),vertex(3,2),vertex(3,3)
             if (n_vertices==4) then
             write (nout,"(1x,a,1p,3e12.4)") "vertex(4,1-3):..............",vertex(4,1),vertex(4,2),vertex(4,3)       
             endif
             write (nout,"(1x,a)")  " "
          end if
       endif
    enddo         
    call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
    if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"Section_flow_rate DATA",ninp,nout)) return
 end do

 return
end subroutine ReadSectionFlowRate
!---split



!AA601 (whole subroutine)
!cfile ReadInputFile.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name   : ReadDBSPH
!
! Creation      : Amicarelli A., 26Jan15
!
!************************************************************************************
! Module purpose : Reading input for DBSPH boundary treatment
!
! Calling routines: ReadInput
!
! Called routines: /
!
!************************************************************************************
subroutine ReadDBSPH (ainp,comment,nrighe,ier,ninp,nout)

! modules
use GLOBAL_MODULE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE

!Declarations
implicit none
character(80),intent(inout) :: ainp 
character(1),intent(inout) :: comment
integer(4),intent(inout) :: nrighe,ier,ninp,nout
character(80) :: lcase 
logical :: MUSCL_boundary_flag,in_built_monitors
integer(4) :: ioerr,n_monitor_points,n_monitor_regions,i,alloc_stat,dealloc_stat,n_kinematics_records,j,n_inlet,n_outlet
double precision :: dx_dxw,k_w
double precision,dimension(:) :: monitor_region(6)           
integer(4),allocatable,dimension(:) :: monitor_IDs
logical,external :: ReadCheck

! in case of restart the cards are not read
if (restart) then
  do while (TRIM(lcase(ainp)) /= "##### end dbsph #####") ! lower case letters are required
    call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
    if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH DATA",ninp,nout)) return
  end do
  return
end if
call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH DATA",ninp,nout)) return
do while (TRIM(lcase(ainp)) /= "##### end dbsph #####")
!Reading the ratio between the fluid and the semi-particle sizes (dx/dx_w)
   read (ainp,*,iostat=ioerr) dx_dxw,MUSCL_boundary_flag,k_w
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH GENERAL INPUT",ninp,nout)) return
   call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
   read (ainp,*,iostat=ioerr) n_monitor_points,n_monitor_regions
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH_monitor_numbers",ninp,nout)) return
   if (n_monitor_points>0) then
      if (.not.allocated(monitor_IDs)) allocate (monitor_IDs(n_monitor_points),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(nout,*) 'Allocation of monitor_IDs in ReadDBSPH failed; the program terminates here'
         stop ! Stop the main program
      endif
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
      read (ainp,*,iostat=ioerr) monitor_IDs(:)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH_monitor_IDs",ninp,nout)) return
   endif
   if (n_monitor_regions==1) then
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
      read (ainp,*,iostat=ioerr) monitor_region(:)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH_monitor_region",ninp,nout)) return
   endif
   call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
   read (ainp,*,iostat=ioerr) n_kinematics_records,in_built_monitors
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH_KINEMATICS",ninp,nout)) return  
   if (.not.(allocated(DBSPH%kinematics))) then
      allocate (DBSPH%kinematics(n_kinematics_records,4),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(nout,*) 'Error! Allocation of DBSPH%kinematics in ReadDBSPH failed; the program terminates here.'
         call diagnostic (arg1=5,arg2=340)
         stop ! Stop the main program
         else
            write (nout,'(1x,a)') "Array DBSPH%kinematics successfully allocated in subrouitne ReadDBSPH."
      endif
   endif  
   do j=1,n_kinematics_records
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
      read (ainp,*,iostat=ioerr) DBSPH%kinematics(j,:)  
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH_KINEMATICS_RECORDS",ninp,nout)) return            
   enddo
   call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
   read (ainp,*,iostat=ioerr) n_inlet,n_outlet
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH_INLET_OUTLET",ninp,nout)) return 
   if (n_inlet>0) then
      if (.not.allocated(DBSPH%inlet_sections)) then
         allocate (DBSPH%inlet_sections(n_inlet,10),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(nout,*) 'Allocation of DBSPH%inlet_sections in ReadDBSPH failed; the program terminates here'
            stop ! Stop the main program
         endif
      endif
   endif
   do j=1,n_inlet
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
! Reading position, normal and velocity of an inlet surface element      
      read (ainp,*,iostat=ioerr) DBSPH%inlet_sections(j,:)  
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH_INLET_SECTIONS",ninp,nout)) return            
   enddo 
   if (n_outlet>0) then
! Reading position and normal of an outlet surface element       
      if (.not.allocated(DBSPH%outlet_sections)) then
         allocate (DBSPH%outlet_sections(n_outlet,8),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(nout,*) 'Allocation of DBSPH_outlet_sections in ReadDBSPH failed; the program terminates here'
            stop ! Stop the main program
         endif
      endif   
   endif
   do j=1,n_outlet
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
      read (ainp,*,iostat=ioerr) DBSPH%outlet_sections(j,:)  
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH_OUTLET_SECTIONS",ninp,nout)) return            
   enddo
! Writing the DBSPH input parameters on the log file
   if ((ncord>0).and.(nout > 0)) then
      write (nout,"(1x,a,1p,e12.4)")   "dx/dx_w:........................",dx_dxw
      write (nout,"(1x,a,1p,l12)")     "MUSCL_boundary_flag:............",MUSCL_boundary_flag
      write (nout,"(1x,a,1p,e12.4)")   "k_w(semi-particle coefficient)..",k_w
      write (nout,"(1x,a,1p,i12)")     "n_monitor_points................",n_monitor_points       
      if (n_monitor_points>0) then
      do i=1,n_monitor_points
      write (nout,"(1x,a,1p,e12.4)")   "ID_monitor......................",monitor_IDs(i)        
      end do    
      endif
      write (nout,"(1x,a,1p,i12)")     "n_monitor_regions...............",n_monitor_regions        
      if (n_monitor_regions==1) then
      write (nout,"(1x,a,1p,g12.5)")   "monitor_region_x_min: ..........",monitor_region(1)
      write (nout,"(1x,a,1p,g12.5)")   "monitor_region_x_max: ..........",monitor_region(2)
      write (nout,"(1x,a,1p,g12.5)")   "monitor_region_y_min: ..........",monitor_region(3)
      write (nout,"(1x,a,1p,g12.5)")   "monitor_region_y_max: ..........",monitor_region(4)
      write (nout,"(1x,a,1p,g12.5)")   "monitor_region_z_min: ..........",monitor_region(5)
      write (nout,"(1x,a,1p,g12.5)")   "monitor_region_z_max: ..........",monitor_region(6)
      endif
      write (nout,"(1x,a,1p,i12)")     "n_kinematics_records............",n_kinematics_records 
      write (nout,"(1x,a,1p,l12)")     "in-built_monitor_flag:..........",in_built_monitors      
      do i=1,n_kinematics_records
      write (nout,"(1x,a,1p,4(g12.4))") "time(s),u(m/s),v(m/s),w(m/s):...",DBSPH%kinematics(i,:)        
      end do 
      write (nout,"(1x,a,i12)")         "n_inlet:........................",n_inlet
      do i=1,n_inlet
      write (nout,"(1x,a,1p,9(g12.4))") "x(m),y(m),z(m),n_x,n_y,n_z,u(m/s),v(m/s),w(m/s),length(m): ",DBSPH%inlet_sections(i,:)        
      end do 
      write (nout,"(1x,a,i12)")         "n_outlet:.......................",n_outlet
      do i=1,n_outlet
      write (nout,"(1x,a,1p,6(g12.4))") "x(m),y(m),z(m),n_x,n_y,n_z,length(m),pres(Pa)............: ",DBSPH%outlet_sections(i,:)        
      end do       
      write (nout,"(1x,a)")  " "
! Assignation to the DBSPH parameters 
      DBSPH%dx_dxw = dx_dxw
      DBSPH%MUSCL_boundary_flag = MUSCL_boundary_flag
      DBSPH%k_w = k_w
      DBSPH%n_monitor_points = n_monitor_points 
      DBSPH%n_monitor_regions = n_monitor_regions
      DBSPH%monitor_region(:) = monitor_region(:)  
      if (n_monitor_points>0) then
         if (.not.(allocated(DBSPH%monitor_IDs))) then
            allocate (DBSPH%monitor_IDs(n_monitor_points),STAT=alloc_stat)
            if (alloc_stat/=0) then
               write(nout,*) 'Allocation of DBSPH%n_monitor_points in ReadDBSPH failed; the program terminates here.'
               stop ! Stop the main program
            endif   
         endif       
         DBSPH%monitor_IDs(:) = monitor_IDs(:)
      endif
      DBSPH%n_kinematics_records = n_kinematics_records 
      DBSPH%in_built_monitors = in_built_monitors
      DBSPH%n_inlet = n_inlet   
      DBSPH%n_outlet = n_outlet
   end if
   if (allocated(monitor_IDs)) then
      deallocate (monitor_IDs,STAT=dealloc_stat)
      if (dealloc_stat/=0) then
         write(nout,*) 'Deallocation of monitor_IDs in ReadDBSPH failed; the program terminates here.'
         stop ! Stop the main program
      endif   
   endif   
   call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH DATA",ninp,nout)) return
end do
return
end subroutine ReadDBSPH
!---split