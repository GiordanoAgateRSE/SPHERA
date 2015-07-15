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

