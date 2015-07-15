!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : time_usertype
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! 00  Agate/Guandalini  18/12/07       INITIAL
! 01  Agate/Guandalini  2008           Check and review entire code
!AA504
! 02  Amicarelli        08Apr14        Minor modification for the mixture viscosity (granular flows)
!
!************************************************************************************
! Module purpose : Assignment of the common data
!
! Calling routine: 
!
! Called routines: start_and_stop
!
!************************************************************************************
!
module time_usertype

use GLOBAL_MODULE
!
  integer(4), parameter :: redcard   = 20      ! reduced length of strings
!
!AA504
  integer(4), parameter :: numb_subr = 21      ! number of subroutine to count + 1
!
! common /strings   / including the significative character information
!  
  character(len=redcard) :: routine
!  character(len=lencard) :: prefix,date_exec,nomecaso,nomefile_tempi,nomefile_sta,  &
  character(len=lencard) :: prefix,date_exec,nomefile_sta,  &
                            nomecaso_rest,nomefile_inp,nomefile_msh,nomefile_cad,   &
                            nomefile_sec,nomefile_spy,nomefile_res,nomefile_txt,    &
                            nomefile_chk,nomefile_neu,case_title,case_subtitle

  character(len=32)      :: host_name
  character(LEN=16)      :: case_data  
  character(LEN=8)       :: case_hour
!
! common /timetable / including all the time statistic information
!
  type iiar                                 
    sequence
    integer(4) ihr
    integer(4) imin
    integer(4) isec
  end type
  type(iiar) itime_struct
!
  integer(4), dimension(numb_subr)              :: tot_call
  double precision, dimension(numb_subr,2)      :: tot_times
  character(LEN=40), dimension(numb_subr)       :: tot_routines
!
! .. common assignments ..
!
  common /strings   / tot_routines,routine,date_exec,prefix,nomefile_sta,nomecaso_rest,    &
                      nomefile_inp,nomefile_msh,nomefile_cad,nomefile_sec,                 &
                      nomefile_spy,host_name,nomefile_res,nomefile_txt,           &
                      nomefile_chk,nomefile_neu,case_data,case_hour,case_title,            &
                      case_subtitle
!
  common /timetable / tot_times,tot_call
!
!..
! file uscita tempi cpu
  integer(4)     ::  ifiout

  data tot_routines(        1) / "Initialization                          " /
  data tot_routines(        2) / "Gest_Input                              " /
  data tot_routines(        3) / "Gest_Trans                              " /
  data tot_routines(        4) / "   Loop_ghost                           " /
  data tot_routines(        5) / "   Loop_Irre_2D/3D                      " /
  data tot_routines(        6) / "     Motion Equation & Update velocity  " /
  data tot_routines(        7) / "     Velocity smoothing                 " /
  data tot_routines(        8) / "     Update position                    " /
  data tot_routines(        9) / "     Outgone/Source/Ordgrid             " /
  data tot_routines(       10) / "     Search for neighbourhood particles " /
  data tot_routines(       11) / "     Boundary integrals                 " /
  data tot_routines(       12) / "     Continuity equation                " /
  data tot_routines(       13) / "     State equation                     " /
  data tot_routines(       14) / "     Pressure smoothing                 " /
  data tot_routines(       15) / "     Apparent viscosity                 " /
  data tot_routines(       16) / "     Diffusion model                    " /
  data tot_routines(       17) / "     RK time integration                " /
!
!AA406
  data tot_routines(       18) / "     Wall parameter update              " /
!
!AA501b
  data tot_routines(       19) / "     Rigid body transport               " /
!
!AA504  
  data tot_routines(       20) / "     Mixture viscosity                  " /
  data tot_routines(numb_subr) / "Totale                                  " /

 end module time_usertype
 