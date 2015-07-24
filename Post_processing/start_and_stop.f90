!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-; SPHERA has been authored for RSE SpA by 
!    Andrea Amicarelli, Antonio Di Monaco, Sauro Manenti, Elia Bon, Daria Gatti, Giordano Agate, Stefano Falappi, 
!    Barbara Flamini, Roberto Guandalini, David Zuccalà).
! Main numerical developments of SPHERA: 
!    Amicarelli et al. (2015,CAF), Amicarelli et al. (2013,IJNME), Manenti et al. (2012,JHE), Di Monaco et al. (2011,EACFM). 
! Email contact: andrea.amicarelli@rse-web.it

! This file is part of SPHERA.
! SPHERA is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with SPHERA. If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: start_and_stop                 
! Description: Time recording.          
!----------------------------------------------------------------------------------------------------------------------------------

subroutine start_and_stop (iset,itot)
!------------------------
! Modules
!------------------------ 
use Time_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4), intent(in) :: iset,itot
integer(4) :: is_cwd,i,hours,days,minutes,seconds,ishost
double precision :: appo
integer(4), dimension(8) :: dat_array
character(LEN=1024) :: pwd_name,sphera_dir
character(LEN=256) :: exe_name  !text, 
character(LEN=8) :: dat
character(LEN=5) :: zone
character(LEN=10) :: ct
character(LEN=3),dimension(12) :: mesi
integer :: HOSTNM,GETCWD
double precision,external :: omp_get_wtime
data mesi/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov",   &
   "Dec"/
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
routine='start_and_stop      '
ishost = 0 
!------------------------
! Statements
!------------------------
! Start detecting SPHERA executable file and the input file prefix.
! If this is not present, then the execution is stopped.
select case(iset)
   case (0)
      tot_times = zero
      tot_call = 0
      call cpu_time(tot_times(numb_subr,1))
      call DATE_AND_TIME(dat,ct,zone,dat_array)
      date_exec = mesi(dat_array(2))//" "//dat(7:8)//", "//dat(1:4)//          &
                " at "//ct(1:2)//":"//ct(3:4)//":"//ct(5:10)//" "//zone//" GMT"
      case_data(1:12) = date_exec(1:12)
      case_hour(1:8) = date_exec(17:24)
      read (ct(1:2),'(i2)') itime_struct%ihr
      read (ct(3:4),'(i2)') itime_struct%imin
      read (ct(5:6),'(i2)') itime_struct%isec
! Getting the machine name 
      ishost = hostnm(host_name)
! Getting the path of the working directory
      is_cwd = getcwd(pwd_name)
! Getting SPHERA_DIR environmental variable 
      call getenv ("SPHERA_DIR",sphera_dir)
      is_cwd = len_trim(sphera_dir)
! Reading the name of SPHERA executable file
      call getarg (0,exe_name)
      exe_name = adjustl(exe_name)
! Headings in log file 
      write(nout,*)
      write(nout,'(a,a)') " >>> SPHERA execution started the ",TRIM(date_exec)
      write(nout,*)
      write(nout,'(a,a)') " >>> SPHERA execution module :",TRIM(exe_name)
      write(nout,*)
      write(nout,'(a,a)') " >>> SPHERA execution machine :",TRIM(host_name)
      write(nout,*)
      write(nout,'(a,a)') " >>> Project working directory :",TRIM(pwd_name)
      write(nout,*)
      write(nout,'(a,a)') " >>> Case identifier :",TRIM(nomecaso)
! Detection of the current time
   case (1)
! End the execution by assessing the computational time &
! and the phisical/elapsed time statistics.
      call cpu_time(appo)
      tot_times(numb_subr,2) = appo - tot_times(numb_subr,1)
      call DATE_AND_TIME(dat,ct,zone,dat_array)
      date_exec = mesi(dat_array(2))//" "//dat(7:8)//", "//dat(1:4)//          &
                " at "//ct(1:2)//":"//ct(3:4)//":"//ct(5:10)//" "//zone//" GMT"
      write(nout,'(a)')                                                        &
"----------------------------------------------------------------------------------------"
      write(nout,'(a,a)') " SPHERA execution ended the ",TRIM(date_exec)
      write(nout,'(a)')                                                        &
"----------------------------------------------------------------------------------------"
      write(nout,'(a)')   " "
      write(nout,'(a)')                                                        &
" Summary of the execution times (s) and internal iteration statistics:"
      write(nout,'(a)')                                                        &
"----------------------------------------------------------------------------------------"
      write(nout,'(a)')                                                        &
"   Task                                                  Total          %      number of"
      write(nout,'(a)')                                                        &
"                                                         Elapsed                 calls" 
      write(nout,'(a)')                                                        &
"----------------------------------------------------------------------------------------"
      do i=1,numb_subr
         if (tot_call(i)==0) cycle
         write(nout,"(4x,i3,2x,a,(f15.5,3x,f6.2,1x),i14)") i,tot_routines(i),  &
            tot_times(i,2),tot_times(i,2)*100./tot_times(numb_subr,2),         &
            tot_call(i)
      enddo
      write(nout,'(a)')                                                        &
"----------------------------------------------------------------------------------------"
      days = int(tot_times(numb_subr,2) / 86400)
      hours = int((tot_times(numb_subr,2) - days * 86400) / 3600)
      minutes = int((tot_times(numb_subr,2) - days * 86400 - hours * 3600) / 60)
      seconds = ceiling(tot_times(numb_subr,2) - days * 86400 - hours * 3600 - &
                minutes * 60)
      write(nout,"(4x,a,f10.2,a,4(i2,a3))") "Total elapsed time    : ",        &
         tot_times(numb_subr,2)," s equal to ",days," d ",hours," h ",minutes, &
         " m ",seconds," s " 
      write(nout,'(a)')                                                        &
"----------------------------------------------------------------------------------------"
! To detect the time for the initial call to subroutines
   case (2)
      tot_call(itot) = tot_call(itot) + 1
      tot_call(numb_subr) = tot_call(numb_subr) + 1
      call cpu_time(appo)
      tot_times(itot,1) = appo
! To assess the incremental time for the different subroutines
   case (3)
      call cpu_time(appo)
      tot_times(itot,2) = tot_times(itot,2) + (appo - tot_times(itot,1))
! To end the execution of the binary restart file converter
end select
!------------------------
! Deallocations
!------------------------
return
end subroutine start_and_stop

