!cfile start_and_stop.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : start_and_stop
!
! Last updating : May 26, 2014
!
! Improvement traceback:
!
! 00  Agate/Guandalini  18/12/07       INITIAL
! 01  Agate             07/05/12       Limited and licensed version
! 02  Agate/Guandalini  27/03/13       modify the key system
!
!************************************************************************************
! Module purpose : management of the time statistics and machine checks
!
! Calling modules: sphera, gest_trans
!
! Called modules : getarg       :acquire the command line parameters
!                  system_clock :check the time independently on the machine system
!                  iargc        :detect the number of command line parameter
!                  hostname     :check the operating machine
!                  getcwd       :check the current directory
!                  getenv       : check the contents of the SWEET_DIR env variable
!      
!************************************************************************************
!
  subroutine start_and_stop (iset,itot)
!
!.. directive to compiler on Windows in order to get the machine identifier
!.. it must be removed for Linux or Unix compilers
!
!....... MS$ ATTRIBUTES ALIAS:'_HOSTNAM@8' :: hostnm  ! compilazione windows fortran 6.x
!
use time_usertype
use FILES_ENTITIES
!
!.. Implicits ..
  implicit none
!
!.. Formal Arguments ..
  integer(4), intent(in) :: iset,itot
!
!.. Local variables
  integer(4)          :: is_cwd,i,hours,days,minutes,seconds,ishost
  character(LEN=1024) :: pwd_name,sphera_dir
  character(LEN=256)  :: exe_name  !text, 
  character(LEN=8)    :: dat
  character(LEN=5)    :: zone
  character(LEN=10)   :: ct
  double precision    :: appo
!
!.. Local Arrays ..
  integer(4),      dimension(8)  :: dat_array
  character(LEN=3),dimension(12) :: mesi
!
!.. External Functions ..
  integer :: HOSTNM
  integer :: GETCWD
  double precision,external :: omp_get_wtime
!
!.. Data Assignments ..
  data mesi/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/
!
!.. Executable Statements ...
!
  routine='start_and_stop      '
!
!AA406
 ishost = 0 
! 
!.. start detecting the executable module and the project prefix. if not
!.. assigned, the excution aborts ..
!
  select case (iset)
    case (0)
!
!..detect the data and time when execution starts ..
!
!
!.. initialize the CPU time statistics ..
!
!.. tot_times(1) = total initialization time ???
!.. tot_times(2) = total Gest_Input time
!.. tot_times(3) = total Gest_Trans time
!.. tot_times(4) = total Loop_ghost time
!.. tot_times(5) = total Loop_Irre_2D/3D time
!.. tot_times(6) = total Storage array & Boundary integrals time
!.. tot_times(7) = total Motion Equation time
!.. tot_times(8) = total Velocity smoothing time
!.. tot_times(9) = total Outgone/Source/Ordgrid
!.. tot_times(10) = total Search for neighbourhood particles
!.. tot_times(11) = total Boundary integrals
!.. tot_times(12) = total Continuity equation time
!.. tot_times(13) = total State equation time
!.. tot_times(14) = total Pressure smoothing time
!.. tot_times(15) = total Apparent viscosity time
!.. tot_times(16) = total Diffusion model time
!.. tot_times(17) = RK time integration
!AA406
!.. tot_times(18) = wall parameter updates
!AA501b comment
!.. tot_times(19) = rigid body transport 
!..
!.. tot_times(.) = total ... time
!..
!.. tot_times(numb_subr) = total time
!
!
      tot_times = zero
      tot_call = 0
      call cpu_time(tot_times(numb_subr,1))
!$ tot_times(numb_subr,1) = omp_get_wtime()

      call DATE_AND_TIME(dat,ct,zone,dat_array)
      date_exec = mesi(dat_array(2))//" "//dat(7:8)//", "//dat(1:4)// &
                " at "//ct(1:2)//":"//ct(3:4)//":"//ct(5:10)//" "//zone//" GMT"
      case_data(1:12) = date_exec(1:12)
      case_hour(1:8) = date_exec(17:24)
      read (ct(1:2),'(i2)') itime_struct%ihr
      read (ct(3:4),'(i2)') itime_struct%imin
      read (ct(5:6),'(i2)') itime_struct%isec
!
!.. read the name of the machine from environment ..
!
      ishost = hostnm(host_name)
!
!.. read the pathname of the working directory ..
!
      is_cwd = getcwd(pwd_name)
!
!.. read the SPHERA_DIR environmental variable 
!
!      is_env = getenv ("SPHERA_DIR",sphera_dir)
      call getenv ("SPHERA_DIR",sphera_dir)
!      sphera_dir = "./"
      is_cwd = len_trim(sphera_dir)
!
!.. read the name of executable file
      call getarg (0,exe_name)
      exe_name = adjustl(exe_name)
!      exe_name = text(is_cwd+1:len_trim(text))
!
! .. Open the output listing file (file ASCII) and assign the other file names
!
!      nomefile_tempi =  TRIM(nomecaso)//"_tempi.sta"
!
!      open (unit=nout,file=nomefile_tempi,form='formatted',recl=124,status='REPLACE',err=997,iostat=ierr)
!
! 997 stop 'errore open file: stoppato programma'
!
! ..Write the calculation heading on the output listing ..
!
      write (nout,*)
      write (nout,'(a,a)') " >>> SPHERA execution started the ",TRIM(date_exec)
      write (nout,*)
      write (nout,'(a,a)') " >>> SPHERA execution module :",TRIM(exe_name)
      write (nout,*)
      write (nout,'(a,a)') " >>> SPHERA execution machine :",TRIM(host_name)
!
      write (nout,*)
      write (nout,'(a,a)') " >>> Project working directory :",TRIM(pwd_name)
      write (nout,*)
!      write (nout,'(a,a)') " >>> Project identifier :",TRIM(prefix)
!      write (nout,*)
      write (nout,'(a,a)') " >>> Case identifier :",TRIM(nomecaso)
!
!.. check for the license correctness
!

!AA601 test sub start (with no license)
!call KeyDecoderCheck (sphera_dir,exe_name,dat)
Erosion_Module_Shields_Seminara = .true.
Erosion_Module_Shields_Mohr = .true.
Diffusion_Module = .true.
Explosion_Module = .true.
TemporalScheme_Module = .true.
BodyDynamics_Module = .true.
DBSPH_Module  = .true.
MultiFluid_Module = .true.
MoreFluids_Module = .true.
Granular_flux = .true.
!AA601 test sub end (with no license)

!
!.. detect the current time and check the progressive time ..
!
    case (1)
!
!.. end the execution evaluating the CPU and elapsed time statistics
!
      call cpu_time(appo)
!$ appo = omp_get_wtime()
      tot_times(numb_subr,2) = appo - tot_times(numb_subr,1)


      call DATE_AND_TIME(dat,ct,zone,dat_array)
      date_exec = mesi(dat_array(2))//" "//dat(7:8)//", "//dat(1:4)// &
                " at "//ct(1:2)//":"//ct(3:4)//":"//ct(5:10)//" "//zone//" GMT"
!
      write (nout,'(a)')   "----------------------------------------------------------------------------------------"
      write (nout,'(a,a)') " SPHERA execution ended the ",TRIM(date_exec)
      write (nout,'(a)')   "----------------------------------------------------------------------------------------"
      write (nout,'(a)')   " "
      write (nout,'(a)')   " Summary of the execution times (s) and internal iteration statistics:"
      write (nout,'(a)')   "----------------------------------------------------------------------------------------"
      write (nout,'(a)')   "   Task                                                  Total          %      number of"
      write (nout,'(a)')   "                                                         Elapsed                 calls" 
      write (nout,'(a)')   "----------------------------------------------------------------------------------------"
      do i = 1,numb_subr
        if (tot_call(i) == 0) cycle
        write (nout,"(4x,i3,2x,a,(f15.5,3x,f6.2,1x),i14)") &
        i,tot_routines(i),tot_times(i,2),tot_times(i,2)*100./tot_times(numb_subr,2),tot_call(i)
      end do
      write (nout,'(a)')   "----------------------------------------------------------------------------------------"
      days = int(tot_times(numb_subr,2)/86400)
      hours = int((tot_times(numb_subr,2)-days*86400)/3600)
      minutes = int((tot_times(numb_subr,2)-days*86400-hours*3600)/60)
      seconds = ceiling(tot_times(numb_subr,2)-days*86400-hours*3600-minutes*60)
      write (nout,"(4x,a,f10.2,a,4(i2,a3))")  &
        "Total elapsed time    : ",tot_times(numb_subr,2)," s equal to ",days," d ",hours," h ",  &
        minutes," m ",seconds," s " 
      write (nout,'(a)')   "----------------------------------------------------------------------------------------"
!
! .. close the files
!
!      close (unit=nout)
!
! ..Detect the time for the initial call to subroutines ..
!
    case (2)
      tot_call(itot) = tot_call(itot) + 1
      tot_call(numb_subr) = tot_call(numb_subr) + 1

      call cpu_time(appo)
!$ appo = omp_get_wtime()
      tot_times(itot,1) = appo
!
!.. evaluates the incremental time for the different subroutines
!
    case (3)
      call cpu_time(appo)
!$ appo = omp_get_wtime()
      tot_times(itot,2) = tot_times(itot,2) + (appo - tot_times(itot,1))
!
!.. ends the execution of the binary restart file converter
!
  end select
!
!.. I/O Formats ..
!
  return
end subroutine start_and_stop
!---split

