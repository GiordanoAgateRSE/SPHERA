!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.8.0.
! SPHERA v.8.0 is free software: you can redistribute it and/or modify
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
! Program unit: sphera           
! Description: Main program unit.                   
!-------------------------------------------------------------------------------
program sphera
!------------------------
! Modules
!------------------------ 
use I_O_file_module
use Static_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: ier,i,n,narg
double precision :: starttime,endtime
character(len=255) :: nomearg
character(len=lencard) :: nomsub = "Sphera"
character(100),external :: lcase
double precision,external :: omp_get_wtime
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
! Initializations icoordp for compatibility with xlf90 (IBM Fortran)
icoordp(0,1) = 0
icoordp(1,1) = 1
icoordp(2,1) = 3
icoordp(3,1) = 0
icoordp(0,2) = 0
icoordp(1,2) = 1
icoordp(2,2) = 2
icoordp(3,2) = 3
pinttimeratio = 0
NumOpenSides = 0
SourceSide = 0
Openside = 0
NumOpenFaces = 0
SourceFace = 0
OpenFace = 0
mat = 0
itime_jet = 0
do i=1,MAXOPENSIDES
   RowVelocity(i) = zero
   NumPartperLine(i) = 0
   NumPartFace(i) = 0
enddo
RowPeriod = zero
yfila = zero
zfila = zero
ParticleVolume = zero
kill_flag = .False.
PgZero%vel_type = "   "
PgZero%slip = " "
PgZero%state = "   "
PgZero%kodvel = 0
PgZero%koddens = 0
PgZero%CloseBcOut = 0
PgZero%cella = 0
PgZero%izona = 0
PgZero%icol = 0
PgZero%imed = 0
PgZero%coord(1) = zero
PgZero%coord(2) = zero
PgZero%coord(3) = zero
PgZero%CoordOld(1) = zero
PgZero%CoordOld(2) = zero
PgZero%CoordOld(3) = zero
PgZero%vel(1) = zero
PgZero%vel(2) = zero
PgZero%vel(3) = zero
PgZero%velass(1) = zero
PgZero%velass(2) = zero
PgZero%velass(3) = zero
PgZero%densass = zero
PgZero%acc(1) = zero
PgZero%acc(2) = zero
PgZero%acc(3) = zero
PgZero%mass = zero
PgZero%dens = zero
PgZero%pres = zero
PgZero%dden = zero
PgZero%uni = zero
PgZero%zer(1) = zero
PgZero%zer(2) = zero
PgZero%zer(3) = zero
PgZero%var(1) = zero
PgZero%var(2) = zero
PgZero%var(3) = zero
PgZero%vpres = zero
PgZero%vden = zero
PgZero%secinv = zero
PgZero%dudx = zero
PgZero%dudy = zero
PgZero%dvdx = zero
PgZero%dvdy = zero
PgZero%visc = zero
PgZero%mu = zero
PgZero%vstart(1) = zero
PgZero%vstart(2) = zero
PgZero%vstart(3) = zero
PgZero%velmorr(1) = zero
PgZero%velmorr(2) = zero
PgZero%velmorr(3) = zero
PgZero%tstop = zero
PgZero%mno = zero
PgZero%ang = zero
PgZero%VolFra = zero
PgZero%rhoc = zero
PgZero%rhow = zero
PgZero%tiroc = zero
PgZero%cden = zero
PgZero%wden = zero
PgZero%diffu = zero
PgZero%coefdif = zero
PgZero%veldif(1) = zero
PgZero%veldif(2) = zero
PgZero%veldif(3) = zero
PgZero%IntEn = zero
PgZero%Envar = zero
PgZero%dEdT  = zero
PgZero%Csound = zero
ts_pgZero%ts_coord(1) = zero
ts_pgZero%ts_coord(2) = zero
ts_pgZero%ts_coord(3) = zero
ts_pgZero%ts_vel(1) = zero
ts_pgZero%ts_vel(2) = zero
ts_pgZero%ts_vel(3) = zero
ts_pgZero%ts_acc(1) = zero
ts_pgZero%ts_acc(2) = zero
ts_pgZero%ts_acc(3) = zero
ts_pgZero%ts_dden = zero
ts_pgZero%ts_dens = zero
ts_pgZero%ts_var(1) = zero
ts_pgZero%ts_var(2) = zero
ts_pgZero%ts_var(3) = zero
ts_pgZero%ts_IntEn = zero
ts_pgZero%ts_dEdT = zero
! Max number of particles
PARTICLEBUFFER = INIPARTICLEBUFFER
! Initializing the array of the cell indices 
indicecelle(1,1) = 0
indicecelle(1,2) = 0
indicecelle(1,3) = 0
indicecelle(2,1) = 1
indicecelle(2,2) = 0
indicecelle(2,3) = 0
indicecelle(3,1) = 1
indicecelle(3,2) = 0
indicecelle(3,3) = 1
indicecelle(4,1) = 0
indicecelle(4,2) = 0
indicecelle(4,3) = 1
indicecelle(5,1) = - 1
indicecelle(5,2) =  0
indicecelle(5,3) =  1
indicecelle(6,1) =  1
indicecelle(6,2) =  1
indicecelle(6,3) =  0
indicecelle(7,1) =  0
indicecelle(7,2) =  1 
indicecelle(7,3) =  0
indicecelle(8,1) = - 1
indicecelle(8,2) =  1
indicecelle(8,3) =  0
indicecelle(9,1) =  1
indicecelle(9,2) =  1
indicecelle(9,3) =  1
indicecelle(10,1) =  0
indicecelle(10,2) =  1 
indicecelle(10,3) =  1
indicecelle(11,1) = - 1
indicecelle(11,2) =  1
indicecelle(11,3) =  1
indicecelle(12,1) =  1
indicecelle(12,2) =  1
indicecelle(12,3) = - 1
indicecelle(13,1) =  0
indicecelle(13,2) =  1
indicecelle(13,3) = - 1
indicecelle(14,1) = - 1
indicecelle(14,2) =  1
indicecelle(14,3) = - 1
nomecaso = "??????"
narg = iargc()
!------------------------
! Statements
!------------------------
! Check the command line arguments
if (narg<1) call diagnostic(arg1=1)
do n=1,narg
   do i=1, len(nomearg)
      nomearg(i:i) = " "
   enddo
   call getarg(n,nomearg)
   if (nomecaso=="??????") then
      nomecaso = nomearg
      nomefile(0) = trim(nomecaso)//".inp"
      nomefile(1) = trim(nomecaso)//".out"
      nomefile(2) = trim(nomecaso)//".ris"
      nomefile(3) = trim(nomecaso)//".rst"
      nomefile(4) = trim(nomecaso)//".plb"
      nomefile(5) = trim(nomecaso)//".fro"
      nomefilekill = trim(nomecaso)//".kll"
      nomefileerr = trim(nomecaso)//".err"
      else
! error in the case name: execution is stopped
         call diagnostic(arg1=1)
   endif
enddo
! I/O assignments
call check_files
write(nout,1000) 
1000 format (                                                                  &
1x,'------------------------------------------------------------------------'/ &
1x,'This output file is generated by SPHERA (Smoothed Particle Hydrodynamics'/ &
1x,'research software; mesh-less Computational Fluid Dynamics code).        '/ &
1x,'Copyright 2005-2017 (RSE SpA-formerly ERSE SpA, formerly CESI RICERCA,  '/ &
1x,'formerly CESI-Ricerca di Sistema; SPHERA has been authored for RSE SpA by:'/&
1x,'Andrea Amicarelli, Antonio Di Monaco, Sauro Manenti, Elia Bon,          '/ &
1x,'Daria Gatti, Giordano Agate, Stefano Falappi, Barbara Flamini,          '/ &
1x,'Roberto Guandalini, David Zuccalà).                                     '/ &
1x,'Main numerical developments of SPHERA:                                  '/ &
1x,'Amicarelli et al. (2015,CAF), Amicarelli et al. (2013,IJNME),           '/ &
1x,'Manenti et al. (2012,JHE), Di Monaco et al. (2011,EACFM).               '/ & 
1x,'Email contact: andrea.amicarelli@rse-web.it                             '/ &
1x,'SPHERA is released under the terms of GNU General Public License as     '/ &
1x,'published by the Free Software Foundation,                              '/ &
1x,'either version 3 of the License, or (at your option) any later version. '/ &
1x,'SPHERA is distributed in the hope that it will be                       '/ &
1x,'useful, but WITHOUT ANY WARRANTY; without even the implied warranty of  '/ &
1x,'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    '/ &
1x,' ',//)
call start_and_stop(0,0)
write(nout,*)
write(nout,*) " >> The following files have been assigned and checked:"
write(nout,*) "    Input  file      : ",trim(nomefile(0))
write(nout,*) "    Output file      : ",trim(nomefile(1))
write(nout,*) "    Result file      : ",trim(nomefile(2))
write(nout,*) "    Restart file     : ",trim(nomefile(3))
write(nout,*) "    Free surface file: ",trim(nomefile(4))
write(nout,*) "    Front file       : ",trim(nomefile(5))
write(nout,*)
! Input 
call start_and_stop(2,2)
call Gest_Input
call start_and_stop(3,2)
starttime = zero
! Main algorithm 
call start_and_stop(2,3)
call Gest_Trans
call start_and_stop(3,3)
!------------------------
! Deallocations
!------------------------
call Gest_Dealloc(nomsub)
endtime = zero
if (endtime/=zero) then
   write(nout,*)
   write(nout,*)
   write(nout,*) "Parallel work took ",endtime-starttime," seconds"
   write(nout,*)
   write(nout,*)
endif
write(nout,*)
call start_and_stop(1,0)
write(nout,*)
write(nout,*)
call s_ctime(nout )
write(nout,*)
write(nout,*)
! If the killer file exists (.kll file), then the execution is stopped and 
! results are saved.
if (kill_flag) then
   open(unit=unitkill,file=nomefilekill,form='formatted',status='old')
   close (unit=unitkill,status='delete')
   write(nout,*) " Execution Terminated : kill file found by user request."
   else
      write(nout,*) " Execution Terminated."
endif
write(nout,*)
! Close the files
if (nout>0) close(nout)
if (nres>0) close(nres)
if (nplb>0) close(nplb)
if (nfro>0) close(nfro)
if (ncpt>0) close(ncpt)
stop
!------------------------
! Contains
!------------------------
contains

subroutine check_files
implicit none
logical fileexist
! To assign the output file
open(nout,file=nomefile(1),status="unknown",form="formatted",iostat=ier)
if (ier/=0 ) call diagnostic(arg1=2,arg2=1)
fileexist = .false.
! To check the input file
inquire(file=nomefile(0),exist=fileexist)
if (.not.fileexist) then
   call diagnostic(arg1=2,arg2=2)
   else
      open(ninp,file=nomefile(0),status="old",form="formatted",iostat=ier)
      if (ier/=0) call diagnostic(arg1=2,arg2=3)
endif
return
end subroutine check_files

end program sphera

