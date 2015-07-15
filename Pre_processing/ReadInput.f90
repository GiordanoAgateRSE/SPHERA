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

