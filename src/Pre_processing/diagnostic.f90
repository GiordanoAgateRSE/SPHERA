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
! Program unit: diagnostic                   
! Description: Diagnostic (error) messages.                
!-------------------------------------------------------------------------------
subroutine diagnostic(ierr,ivalue,avalue)
!------------------------
! Modules
!------------------------ 
use I_O_file_module
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: ierr
integer(4),intent(in),optional :: ivalue
character(LEN=lencard),intent(in),optional :: avalue
logical :: print_out
integer(4) :: dato,la,it1,it2,it3
double precision :: dtvel
character(len=lencard) :: stringa
character(LEN=2) :: error
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
print_out = .false.
if (present(avalue)) then
   stringa = adjustl(avalue)
   la = len_trim(stringa)
   else
      stringa = ' '
      la = 0
endif
if (present(ivalue)) then
   dato = ivalue
   else
      dato = -9999
endif
!------------------------
! Statements
!------------------------
write(error,'(i2)') ierr
write(nscr,'(a,a2,a)') ' >> Diagnostic ',error,                                &
   ' detected * Run ended unsuccessfully'
if (ierr/=1) then
   write(nout,'(1x,100(''=''))')
   write(nout,*)
   write(nout,'(a,a2,a)') 'Diagnostic ',error,                                 &
      ' detected * Run ended unsuccessfully'
   write(nout,*)
   flush(nout)
endif
select case (ierr)
   case (1)    
! Errors in arguments
      write(nscr,*) "-----------------------------------"
      write(nscr,*) "ERROR in command line argument !!!."
      write(nscr,*) "-----------------------------------"
      write(nscr,*)                                                            &
" Correct command sequence-> <executablename> <casename> [euristic]|[original]"
   case (2)
! Error in files
      write(nout,'(1x,a)') 'Checkfile function detected an error:'
      select case (dato)
         case (1); write(nout,'(1x,a)') '  1 - Output file cannot be opened.'
         case (2); write(nout,'(1x,a)') '  2 - Input file does not exist.'
         case (3); write(nout,'(1x,a)') '  3 - Input file cannot be opened.'
      endselect
   case (3)
! Error in arrays
      write(nout,'(1x,2a)') 'Array deallocation failed.  Routine -> ',         &
         stringa(1:la)
      write(nout,'(1x,a)') 'Run is stopped.'
   case (4)
! Error in arrays
      write(nout,'(1x,2a)')                                                    &
         'Array allocation/deallocation failed.  Routine -> ',stringa(1:la)
      write(nout,'(1x,a)') 'Run is stopped.'
   case (5)
! Error in input
      write(nout,'(1x,a)')                                                     &
         'An error was detected analyzing the input data deck for the line:'
      write(nout,'(1x,a)') '->'//stringa(1:la)//'<-'
      write(nout,'(1x,a)') 'Specific code error is:'
      select case (dato)
         case (1); write(nout,'(1x,a)') '    1 - Unknown Input Option'
         case (2); write(nout,'(1x,a)')                                        &
            "    2 - Input file doesn't match program versions"
         case (3); write(nout,'(1x,a)')                                        &
            '    3 - Unknown Domain Type. Routine -> '
         case (4); write(nout,'(1x,a)')                                        &
            '    4 - Wrong type data or unsufficient number of data'
         case (5); write(nout,'(1x,a)')                                        &
            '    5 - Type of Domain = rect is possible only with 2D problems '
         case (6); write(nout,'(1x,a)')                                        &
            '    6 - Type of Erosion model is not correct (shields, mohr) '
         case (101); write(nout,'(1x,a)') '  101 - Unrecognised boundary type'
         case (102); write(nout,'(1x,a)')                                      &
            '  102 - PERIMETER boundary type: different first and last vertex'
         case (103); write(nout,'(1x,a)')                                      &
            '  103 - TAPIS boundary type: 2 vertices are requested'
         case (104); write(nout,'(1x,a)')                                      &
            '  104 - FACE definition: face already defined'
         case (201); write(nout,'(1x,a)')                                      &
            '  201 - Restart file cannot be opened'
         case (203); write(nout,'(1x,a)')                                      &
            '  203 - Restart step location failed'
         case (204); write(nout,'(1x,a)')                                      &
'  204 - Restart '//stringa(1:la)//' fails at the record reported above'
         case (205); write(nout,'(1x,a)') '  205 - Restart option is unknown'
         case (301); write(nout,'(1x,a)')                                      &
            '  301 - Ascii input option is unknown' 
         case (302); write(nout,'(1x,a)')                                      &
            '  302 - Ascii input version does not match code version'
         case (304); write(nout,'(1x,a)')                                      &
'  304 - Ascii input '//stringa(1:la)//' fails at the record reported above'
         case (305); write(nout,'(1x,a)')                                      &
'  305 - Ascii input '//stringa(1:la)//' fails at the record for level option'
         case (306); write(nout,'(1x,a)')                                      &
'  306 - Ascii input '//stringa(1:la)//' fails at the record for erosion Shields Mohr module'
         case (307); write(nout,'(1x,a)')                                      &
'  307 - Ascii input '//stringa(1:la)//' fails at the record for diffusion module'
         case (308); write(nout,'(1x,a)')                                      &
'  308 - Ascii input '//stringa(1:la)//' fails at the record for esplosion module'
         case (309); write(nout,'(1x,a)')                                      &
'  309 - Ascii input '//stringa(1:la)//' fails at the record for multi fluid module'
         case (310); write(nout,'(1x,a)')                                      &
'  310 - Ascii input '//stringa(1:la)//' fails at the record for more fluid module'
         case (311); write(nout,'(1x,a)')                                      &
'  311 - Ascii input '//stringa(1:la)//' fails at the record for Granular flux module'
         case (312); write(nout,'(1x,a)')                                      &
'  312 - Ascii input '//stringa(1:la)//' fails at the record for erosion Shields Van Rijn Seminara module'
         case (320); write(nout,'(1x,a)')                                      &
'  320 - Ascii input '//stringa(1:la)//' fails at the record for temporal scheme module'
         case (330); write(nout,'(1x,a)')                                      &
'  330 - Ascii input '//stringa(1:la)//' fails at the record for body dynamics module'
         case (340); write(nout,'(1x,a)')                                      &
'  340 - Ascii input '//stringa(1:la)//' fails at the record for DB-SPH module'
         case default; write(nout,'(1x,a)') '  ??? - Unknown error type'
      endselect
   case (6)
! Error in source
      write(nout,'(1x,a)')                                                     &
         'An error was detected analyzing the source of particles'
      write(nout,'(1x,a)') 'Specific code error is:'
      select case (dato)
         case (1); write(nout,'(1x,2a)')                                       &
' 1 - Number of particles nag > PARTICLEBUFFER (a).  Routine ->',stringa(1:la)
         case (2); write(nout,'(1x,2a)')                                       &
' 2 - Number of particles nag > PARTICLEBUFFER (b).  Routine ->',stringa(1:la)
      endselect
   case (7); write(nout,'(1x,2a)')                                             &
      'Number of cels in grid: NumCellmax < Grid%nmax.  Routine ->',           &
      stringa(1:la)
   case (8)
! Error at boundaries
      write(nout,'(1x,a)') 'An error was detected analyzing the boundaries'
      write(nout,'(1x,a)') 'Specific code error is:'
      select case (dato)
         case (1); write(nout,'(1x,a,i10,2a)')                                 &
' 1 - New Max num particles*BoundaryCloseSides the new value is: MaxNcbs = ',  &
            MaxNcbs,'.  Routine ->',stringa(1:la)
         case (2); write(nout,'(1x,a,i10,2a)')                                 &
' 2 - New Max num particles*BoundaryCloseFaces the new value is: MaxNcbf = ',  &
            MaxNcbf,'.  Routine ->',stringa(1:la)
         case (3); write(nout,'(1x,2a)')                                       &
' 3 - The kernel table is not implemented for 2D case.  Routine ->',           &
            stringa(1:la)
         case (4); write(nout,'(1x,2a)')                                       &
' 4 - Sides are not consecutive.  Routine ->',stringa(1:la)
         case (5); write(nout,'(1x,3a)')                                       &
' 5 - Number of cell calculated is different from the number stored for the ', &
            'current particle.  Routine ->',stringa(1:la)
         case (6); write(nout,'(1x,2a)')                                       &
' 6 - Number of close boundary faces ncbf > MAXCLOSEBOUNDFACES (a).  Routine ->'&
            ,stringa(1:la)
         case (7); write(nout,'(1x,2a)')                                       &
' 7 - Number of close boundary faces ncbf > MAXCLOSEBOUNDFACES (b).  Routine ->'&
            ,stringa(1:la)
         case (8); write(nout,'(1x,2a)')                                       &
' 8 - Number of close boundary sides Ncbs > MAXCLOSEBOUNDSIDES.  Routine ->',  &
            stringa(1:la)
         case (9); write(nout,'(1x,2a)')                                       &
' 9 - Number of close boundary sides Ncbs > LIMCLOSEBOUNDSIDES.  Routine ->',  &
            stringa(1:la)
         case (10); write(nout,'(1x,2a)')                                      &
'10 - Number of convex edges NumBEdges > MAXNUMCONVEXEDGES.  Routine ->',      &
            stringa(1:la)
      endselect
   case (9)
! Error during the simulation
      write(nout,'(1x,a)') 'An error was detected during the loop'
      write(nout,'(1x,a)') 'Specific code error is:'
      select case (dato)
         case (1); write(nout,'(1x,a,i10,2a)')                                 &
' 1 - Number of fluid particles array Array_Flu greater than max ',            &
            PARTICLEBUFFER,' (a).  Routine ->',stringa(1:la)
         case (2); write(nout,'(1x,a,i10,2a)')                                 &
' 2 - Number of fluid particles array Array_Flu greater than max ',            &
            PARTICLEBUFFER,' (b).  Routine ->',stringa(1:la)
         case (3); write(nout,'(1x,2a)')                                       &
' 3 - Increase parameter MAXCLOSEBOUNDFACES.  Routine ->',stringa(1:la)
      endselect
      print_out = .true.
   case (10)    
! Error during the execution of auxiliary subroutines
      write(nout,'(1x,a)') 'An error was detected during the tools execution'
      write(nout,'(1x,a)') 'Specific code error is:'
      select case (dato)
         case (1); write(nout,'(1x,a,i10,2a)')                                 &
' 1 - Too many fluid particle neighbours. Maximum number allowed is: ',        &
                      NMAXPARTJ,' (a).  Routine ->',stringa(1:la)
         case (2); write(nout,'(1x,2a)')                                       &
            ' 2 - Array BoundaryVertex bounds exceed.  Routine ->',            &
            stringa(1:la)
         case (3); write(nout,'(1x,2a)')                                       &
            ' 3 - Array Vertice bounds exceed.  Routine ->',                   &
            stringa(1:la)
         case (4); write(nout,'(1x,2a)')                                       &
' 4 - Number of particles nag > PARTICLEBUFFER.  Routine ->',                  &
            stringa(1:la)
         case (5); write(nout,'(1x,2a)')                                       &
            ' 5 - Unknown Domain Type. Routine -> ',stringa(1:la)
         case (6); write(nout,'(1x,2a)')                                       &
' 6 - Number of open sides NumOpenSides > MAXOPENSIDES. Routine -> ',          &
            stringa(1:la)
         case (7); write(nout,'(1x,2a)')                                       &
' 7 - Number of open faces NumOpenFaces > MAXOPENFACES. Routine -> ',          &
            stringa(1:la)
         case (8); write(nout,'(1x,2a)')                                       &
' 8 - Free level condition is not found. Routine -> ',                         &
            stringa(1:la)
         case (9); write(nout,'(1x,a,i10,2a)')                                 &
' 9 - Too many wall element neighbours. Maximum number allowed is: ',          &
                      NMAXPARTJ,' (a).  Routine ->',stringa(1:la)
         case (88); write(nout,'(1x,2a)')                                      &
'88 - Error deltat law calculation for "law" particles. Routine -> '           &
            ,stringa(1:la)
      endselect
      print_out = .true.
   case (11)
! Error in erosion scheme
      write(nout,'(1x,a)')                                                     &
         'An error was detected during execution of erosion model'
      write(nout,'(1x,a)') 'Specific code error is:'
      select case (dato)
         case (1); write(nout,'(1x,2a)')                                       &
' 1 - The Ustar in Shields erosion model is NaN.  Routine ->'                  &
            ,stringa(1:la)
         case (2); write(nout,'(1x,3a)')                                       &
' 2 - It is Impossible to find liquid interface and solid interface for the ', &
            'current particle.  Routine ->',stringa(1:la)
      endselect
      print_out = .true.
   case default 
      write(nout,*) 'Unpredictable error'
      print_out = .true.
endselect
if (ierr/=1) then
   write(nout,'(1x,100(''=''))')
   flush(nout)
endif
if (print_out) then
   write(nout,'(1x,(a))') ' ---------------------------'
   write(nout,'(1x,(a))') ' Print results at last step.'
   write(nout,'(1x,(a))') ' ---------------------------'
   write(nscr,'(1x,(a))') ' ---------------------------'
   write(nscr,'(1x,(a))') ' Print results at last step.'
   write(nscr,'(1x,(a))') ' ---------------------------'
   it1 = 0
   it2 = 0
   it3 = 0
   dtvel = 99999.0d0
! Result saving and restart file
   if (nout>0) then
      call Print_Results(it1,it2,'fine__')
   endif
   if (nres>0) then
      call memo_results(it1,it2,it3,dtvel,'fine__')
   endif
   if (vtkconv) then
      call result_converter('fine__')
   endif
endif
stop
!------------------------
! Deallocations
!------------------------
end subroutine diagnostic

