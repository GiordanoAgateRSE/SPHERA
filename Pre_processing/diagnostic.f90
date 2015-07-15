!cfile diagnostic.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : diagnostic
!
! Last updating : May 02, 2014
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Agate             08/05/2012     Limited and licensed version
! 04  Agate             27/03/2013     Limited and licensed version for single and multi users
! 05  Agate             02/05/14       Add more modules check in license
!
!************************************************************************************
! Module purpose : Diagnostic management
!
! Calling routine: all
!
! Called routines: memo_results
!                  print_results
!                  result_converter
!
!************************************************************************************
!
  subroutine diagnostic (ierr,ivalue,avalue)
!
!.. assign modules
  use FILES_ENTITIES
  use GLOBAL_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
  integer(4), intent(in)                      :: ierr
  integer(4), intent(in), optional            :: ivalue
  character(LEN=lencard), intent(in),optional :: avalue
!
!.. Local Scalars ..
  character(LEN=2)       :: error
  integer(4)             :: dato,la,it1,it2,it3
  double precision       :: dtvel
  character(len=lencard) :: stringa
  logical                :: print_out
!
!.. Executable Statements ..
!
  print_out = .false.
  if (present(avalue)) then
    stringa = adjustl(avalue)
    la = len_trim(stringa)
  else
    stringa = ' '
    la = 0
  end if
  if (present(ivalue)) then
    dato = ivalue
  else
    dato = -9999
  end if
!
  write (error,'(i2)') ierr
  write (nscr,'(a,a2,a)') ' >> Diagnostic ',error,' detected * Run ended unsuccessfully'
  if (ierr /= 1) then
    write (nout,'(1x,100(''=''))')
    write (nout,'(31x,5a)') ' *  ',acode,'  version  ',version,'  * '
    write (nout,*)
    write (nout,'(a,a2,a)') 'Diagnostic ',error,' detected * Run ended unsuccessfully'
    write (nout,*)
    flush(nout)
  end if
!
  select case (ierr)
    case ( 1)    ! Errors in arguments
      write (nscr,*) "-----------------------------------"
      write (nscr,*) "ERROR in command line argument !!!."
      write (nscr,*) "-----------------------------------"
      write (nscr,*) " Correct command sequence-> <executablename> <casename> [euristic]|[original]"
!
    case ( 2)    ! Errors in files
      write (nout,'(1x,a)') 'Checkfile function detected an error:'
      select case (dato)
        case (1); write (nout,'(1x,a)') '  1 - Output file cannot be opened.'
        case (2); write (nout,'(1x,a)') '  2 - Input file does not exist.'
        case (3); write (nout,'(1x,a)') '  3 - Input file cannot be opened.'
      end select
!
    case ( 3)    ! Errors in arrays
      write (nout,'(1x,2a)') 'Array deallocation failed.  Routine -> ',stringa(1:la)
      write (nout,'(1x,a)') 'Run is stopped.'
!
    case ( 4)    ! Errors in arrays
      write (nout,'(1x,2a)') 'Array allocation/deallocation failed.  Routine -> ',stringa(1:la)
      write (nout,'(1x,a)') 'Run is stopped.'
!
    case ( 5)    ! Errors in input
      write (nout,'(1x,a)') 'An error was detected analyzing the input data deck for the line:'
      write (nout,'(1x,a)') '->'//stringa(1:la)//'<-'
      write (nout,'(1x,a)') 'Specific code error is:'
      select case (dato)
        case (1); write (nout,'(1x,a)')   '    1 - Unknown Input Option'
        case (2); write (nout,'(1x,a)')   "    2 - Input file doesn't match program versions"
        case (3); write (nout,'(1x,a)')   '    3 - Unknown Domain Type. Routine -> '
        case (4); write (nout,'(1x,a)')   '    4 - Wrong type data or unsufficient number of data'
        case (5); write (nout,'(1x,a)')   '    5 - Type of Domain = rect is possible only with 2D problems '
        case (6); write (nout,'(1x,a)')   '    6 - Type of Erosion model is not correct (shields, mohr) '
        case (101); write (nout,'(1x,a)') '  101 - Unrecognised boundary type'
        case (102); write (nout,'(1x,a)') '  102 - PERIMETER boundary type: different first and last vertex'
        case (103); write (nout,'(1x,a)') '  103 - TAPIS boundary type: 2 vertices are requested'
        case (104); write (nout,'(1x,a)') '  104 - FACE definition: face already defined'
        case (201); write (nout,'(1x,a)') '  201 - Restart file cannot be opened'
        case (203); write (nout,'(1x,a)') '  203 - Restart step location failed'
        case (204); write (nout,'(1x,a)') '  204 - Restart '//stringa(1:la)//' fails at the record reported above'
        case (205); write (nout,'(1x,a)') '  205 - Restart option is unknown'
        case (301); write (nout,'(1x,a)') '  301 - Ascii input option is unknown'
        case (302); write (nout,'(1x,a)') '  302 - Ascii input version does not match code version'
        case (304); write (nout,'(1x,a)') '  304 - Ascii input '//stringa(1:la)//' fails at the record reported above'
        case (305); write (nout,'(1x,a)') '  305 - Ascii input '//stringa(1:la)//' fails at the record for level option'
        case (306); write (nout,'(1x,a)') '  306 - Ascii input '//stringa(1:la)//' fails at the record for erosion Shields Mohr module'
        case (307); write (nout,'(1x,a)') '  307 - Ascii input '//stringa(1:la)//' fails at the record for diffusion module'
        case (308); write (nout,'(1x,a)') '  308 - Ascii input '//stringa(1:la)//' fails at the record for esplosion module'
        case (309); write (nout,'(1x,a)') '  309 - Ascii input '//stringa(1:la)//' fails at the record for multi fluid module'
        case (310); write (nout,'(1x,a)') '  310 - Ascii input '//stringa(1:la)//' fails at the record for more fluid module'
        case (311); write (nout,'(1x,a)') '  311 - Ascii input '//stringa(1:la)//' fails at the record for Granular flux module'
        case (312); write (nout,'(1x,a)') '  312 - Ascii input '//stringa(1:la)//' fails at the record for erosion Shields Van Rijn Seminara module'
        case (320); write (nout,'(1x,a)') '  320 - Ascii input '//stringa(1:la)//' fails at the record for temporal scheme module'
        case (330); write (nout,'(1x,a)') '  330 - Ascii input '//stringa(1:la)//' fails at the record for body dynamics module'
        case (340); write (nout,'(1x,a)') '  340 - Ascii input '//stringa(1:la)//' fails at the record for DB-SPH module'
        case default; write (nout,'(1x,a)') '  ??? - Unknown error type'
      end select
!
    case ( 6)    ! Errors in Source
      write (nout,'(1x,a)') 'An error was detected analyzing the source of particles'
      write (nout,'(1x,a)') 'Specific code error is:'
      select case (dato)
        case (1); write (nout,'(1x,2a)') ' 1 - Number of particles nag > PARTICLEBUFFER (a).  Routine ->',stringa(1:la)
        case (2); write (nout,'(1x,2a)') ' 2 - Number of particles nag > PARTICLEBUFFER (b).  Routine ->',stringa(1:la)
      end select
!
    case ( 7); write (nout,'(1x,2a)') 'Number of cels in grid: NumCellmax < Grid%nmax.  Routine ->',stringa(1:la)
!
    case ( 8)    ! Errors in boundaries
      write (nout,'(1x,a)') 'An error was detected analyzing the boundaries'
      write (nout,'(1x,a)') 'Specific code error is:'
      select case (dato)
        case (1);  write (nout,'(1x,a,i10,2a)') ' 1 - New Max num particles*BoundaryCloseSides the new value is: MaxNcbs = ', &
                                                MaxNcbs,'.  Routine ->',stringa(1:la)
        case (2);  write (nout,'(1x,a,i10,2a)') ' 2 - New Max num particles*BoundaryCloseFaces the new value is: MaxNcbf = ', &
                                                MaxNcbf,'.  Routine ->',stringa(1:la)
        case (3);  write (nout,'(1x,2a)')       ' 3 - The kernel table is not implemented for 2D case.  Routine ->',stringa(1:la)
        case (4);  write (nout,'(1x,2a)')       ' 4 - Sides are not consecutive.  Routine ->',stringa(1:la)
        case (5);  write (nout,'(1x,3a)')       ' 5 - Number of cell calculated is different from the number stored for the ', &
                                               'current particle.  Routine ->',stringa(1:la)
        case (6);  write (nout,'(1x,2a)')       ' 6 - Number of close boundary faces ncbf > MAXCLOSEBOUNDFACES (a).  Routine ->', &
                                                stringa(1:la)
        case (7);  write (nout,'(1x,2a)')       ' 7 - Number of close boundary faces ncbf > MAXCLOSEBOUNDFACES (b).  Routine ->', &
                                                stringa(1:la)
        case (8);  write (nout,'(1x,2a)')       ' 8 - Number of close boundary sides Ncbs > MAXCLOSEBOUNDSIDES.  Routine ->', &
                                                stringa(1:la)
        case (9);  write (nout,'(1x,2a)')       ' 9 - Number of close boundary sides Ncbs > LIMCLOSEBOUNDSIDES.  Routine ->', &
                                                stringa(1:la)
        case (10); write (nout,'(1x,2a)')       '10 - Number of convex edges NumBEdges > MAXNUMCONVEXEDGES.  Routine ->', &
                                                stringa(1:la)
      end select
!
    case ( 9)    ! Errors in Loop
      write (nout,'(1x,a)') 'An error was detected during the loop'
      write (nout,'(1x,a)') 'Specific code error is:'
      select case (dato)
        case (1); write (nout,'(1x,a,i10,2a)') ' 1 - Number of fluid particles array Array_Flu greater than max ',PARTICLEBUFFER, &
                                               ' (a).  Routine ->',stringa(1:la)
        case (2); write (nout,'(1x,a,i10,2a)') ' 2 - Number of fluid particles array Array_Flu greater than max ',PARTICLEBUFFER, &
                                               ' (b).  Routine ->',stringa(1:la)
        case (3); write (nout,'(1x,2a)')       ' 3 - Increase parameter MAXCLOSEBOUNDFACES.  Routine ->',stringa(1:la)
      end select
      print_out = .true.
!
    case (10)    ! Errors in Tools
      write (nout,'(1x,a)') 'An error was detected during the tools execution'
      write (nout,'(1x,a)') 'Specific code error is:'
      select case (dato)
        case (1); write (nout,'(1x,a,i10,2a)') ' 1 - Number of surrounding particles of current particle is greater than max ', &
                                                NMAXPARTJ,' (a).  Routine ->',stringa(1:la)
        case (2); write (nout,'(1x,2a)')       ' 2 - Array BoundaryVertex bounds exceed.  Routine ->',stringa(1:la)
        case (3); write (nout,'(1x,2a)')       ' 3 - Array Vertice bounds exceed.  Routine ->',stringa(1:la)
        case (4); write (nout,'(1x,2a)')       ' 4 - Number of particles nag > PARTICLEBUFFER.  Routine ->',stringa(1:la)
        case (5); write (nout,'(1x,2a)')       ' 5 - Unknown Domain Type. Routine -> ',stringa(1:la)
        case (6); write (nout,'(1x,2a)')       ' 6 - Number of open sides NumOpenSides > MAXOPENSIDES. Routine -> ',stringa(1:la)
        case (7); write (nout,'(1x,2a)')       ' 7 - Number of open faces NumOpenFaces > MAXOPENFACES. Routine -> ',stringa(1:la)
        case (8); write (nout,'(1x,2a)')       ' 8 - Free level condition is not found. Routine -> ',stringa(1:la)
        case (88); write (nout,'(1x,2a)')      '88 - Error deltat law calculation for "law" particles. Routine -> ',stringa(1:la)
      end select
      print_out = .true.
!
    case (11)    ! Errors in Erosion Model
      write (nout,'(1x,a)') 'An error was detected during execution of erosion model'
      write (nout,'(1x,a)') 'Specific code error is:'
      select case (dato)
        case (1); write (nout,'(1x,2a)') ' 1 - The Ustar in Shields erosion model is NaN.  Routine ->',stringa(1:la)
        case (2); write (nout,'(1x,3a)') ' 2 - It is Impossible to find liquid interface and solid interface for the ', &
                                         'current particle.  Routine ->',stringa(1:la)
      end select
      print_out = .true.
!
    case (73)    ! Errors in License file
      write (nout,'(1x,a)')'License is not correct. Code error:'
      select case (dato)
        case (1) ; write (nout,'(1x,a)')  '1  - The license file has not been provided with the installation'
        case (2) ; write (nout,'(1x,3a)') '2  - Running machine name ',stringa(1:la),' does not match the license machine name'
        case (3) ; write (nout,'(1x,3a)') '3  - Running machine address ',stringa(1:la),' does not match the license machine address'
        case (4) ; write (nout,'(1x,3a)') '4  - Code release ',stringa(1:la),' does not match the license release'
        case (5) ; write (nout,'(1x,3a)') '5  - The user ',stringa(1:la),' does not match the license user'
        case (6) ; write (nout,'(1x,3a)') '6  - No modules are licensed: at least one module must be licensed'
        case (7) ; write (nout,'(1x,a)')  '7  - License has been expired',stringa(1:la),' days ago'
        case (8) ; write (nout,'(1x,a)')  '8  - License file is corrupted'
        case (9) ; write (nout,'(1x,a)')  '9  - The master node of the Linux cluster has not been setted in the proper environmental variable (MASTER_NAME)'
        case (10); write (nout,'(1x,a)')  '10 - Error in <whoami> call. License key cannot be verified'
        case (11); write (nout,'(1x,a)')  '11 - Error in <',stringa(1:la),'> call. License key cannot be verified'
        case (12); write (nout,'(1x,a)')  '12 - Error in the name of the executable file: the extension must be "exe" (Windows) or "x" (Linux/Unix)'
        case default 
          write (nout,*) 'Generic error reading the license. Please contact the assistance.'
       end select
!
    case default 
      write (nout,*) 'Unpredictable error'
      print_out = .true.
  end select
  if (ierr /= 1) then
    write (nout,'(1x,100(''=''))')
    flush(nout)
  end if
!
  if (print_out) then
    write (nout,'(1x,(a))') ' ---------------------------'
    write (nout,'(1x,(a))') ' Print results at last step.'
    write (nout,'(1x,(a))') ' ---------------------------'
    write (nscr,'(1x,(a))') ' ---------------------------'
    write (nscr,'(1x,(a))') ' Print results at last step.'
    write (nscr,'(1x,(a))') ' ---------------------------'
!
!=== SCRITTURA SU FILE RISULTATI ######################
!
    it1 = 0
    it2 = 0
    it3 = 0
    dtvel = 99999.0
!Salvataggio risultati e restart ULTIMO STEP (sempre?)
    if (nout > 0) then
      call print_results ( it1, it2, 'fine__' )
    end if
    if ( nres > 0 ) then
      call memo_results ( it1, it2, it3, dtvel, 'fine__' )
    end if
!
!=== SCRITTURA SU FILE FORMATO VTK ######################
!
    if ( vtkconv ) then
      call result_converter ('fine__')
    end if
!
  end if
!
  stop
!
  end subroutine diagnostic
!---split

