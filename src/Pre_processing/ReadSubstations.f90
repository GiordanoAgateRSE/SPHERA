!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2018 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.9.0.0
! SPHERA v.9.0.0 is free software: you can redistribute it and/or modify
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
! Program unit: ReadSubstations                               
! Description: Input management for monitoring the electrical substations (ref. 
!              template input file)                  
!-------------------------------------------------------------------------------
subroutine ReadSubstations(ainp,comment,nrighe,ier)
!------------------------
! Modules
!------------------------
use I_O_file_module
use Static_allocation_module                            
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nrighe,ier,alloc_stat,ioerr,i,n_sub,n_vertices,substation_ID
integer(4) :: type_ID
double precision :: dt_out,area
character(1) :: comment
character(100) :: ainp,lcase
double precision :: vec_aux_4(3)
double precision :: vertex(6,2)
logical,external :: ReadCheck
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"Substations DATA",ninp,ulog)) return
do while (TRIM(lcase(ainp))/="##### end substations #####")
! Reading the number of substations and their writing time step
   read (ainp,*,iostat=ioerr) n_sub,dt_out
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"Substations GENERAL INPUT",ninp,  &
      ulog)) return
! Writing the number of substations and the writing time step on the log file
   if (ulog>0) then
      write(ulog,"(1x,a,1p,i12)")   "n_substations:......................",    &
         n_sub
      write(ulog,"(1x,a,1p,e12.4)") "dt_out_substations:.................",    &
         dt_out
      write(ulog,"(1x,a)")  " "
   endif
! Allocation of the array of the substations
   if ((n_sub>0).and.(.not.allocated(substations%sub))) then
      allocate(substations%sub(n_sub),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(uerr,*) 'Allocation of "substations%sub" in the ',              &
            'subroutine "ReadSubstations" failed; the execution terminates ',  &
            'here.'
         stop
         else
            write(ulog,*) 'Allocation of "substations%sub" in the subroutine ',&
               '"ReadSubstations" is successfully completed.'
      endif
      substations%sub(:)%area = 0.d0
      substations%sub(:)%POS_fsum(1) = 0.d0
      substations%sub(:)%POS_fsum(2) = 0.d0
      substations%sub(:)%Ymax(1) = -1.d9
      substations%sub(:)%Ymax(2) = -1.d9
      substations%sub(:)%EOT(1) = 0.d0
      substations%sub(:)%EOT(2) = 0.d0
      substations%sub(:)%Val = 0.d0
      substations%sub(:)%n_DEM_vertices = 0
   endif
   substations%n_sub = n_sub
   substations%dt_out = dt_out
   substations%it_out_last = 0
! Loop over the substations
   do i=1,n_sub
! Reading the substation variables
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      read (ainp,*,iostat=ioerr) substation_ID
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"substation_ID",ninp,ulog))     &
         return
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      read (ainp,*,iostat=ioerr) type_ID,n_vertices
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"type_ID,n_vertices",ninp,ulog))&
         return
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      read (ainp,*,iostat=ioerr) vertex(1,1),vertex(1,2)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"vertex_1",ninp,ulog)) return
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      read (ainp,*,iostat=ioerr) vertex(2,1),vertex(2,2)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"vertex_2",ninp,ulog)) return
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      read (ainp,*,iostat=ioerr) vertex(3,1),vertex(3,2)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"vertex_3",ninp,ulog)) return
      if (n_vertices>3) then
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read (ainp,*,iostat=ioerr) vertex(4,1),vertex(4,2)
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"vertex_4",ninp,ulog)) return           
      endif
      if (n_vertices>4) then
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read (ainp,*,iostat=ioerr) vertex(5,1),vertex(5,2)
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"vertex_5",ninp,ulog)) return           
      endif
      if (n_vertices>5) then
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read (ainp,*,iostat=ioerr) vertex(6,1),vertex(6,2)
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"vertex_6",ninp,ulog)) return           
      endif
! Assignation to the substation variables
      substations%sub(i)%type_ID = type_ID
      substations%sub(i)%n_vertices = n_vertices
      substations%sub(i)%vert(:,:) = vertex(:,:)
! Writing on the log file
      if (ncord>0) then
         if (ulog>0) then
            write(ulog,"(1x,a,i12)")       "substation_ID:..............",     &
               substation_ID
            write(ulog,"(1x,a,i12)")       "type_ID:....................",     &
               type_ID
            write(ulog,"(1x,a,i12)")       "n_vertices:.................",     &
               n_vertices
            write(ulog,"(1x,a,1p,3e12.4)") "vertex(1,1:2):..............",     &
               vertex(1,1),vertex(1,2)
            write(ulog,"(1x,a,1p,3e12.4)") "vertex(2,1:2):..............",     &
               vertex(2,1),vertex(2,2)
            write(ulog,"(1x,a,1p,3e12.4)") "vertex(3,1:2):..............",     &
               vertex(3,1),vertex(3,2)
            if (n_vertices>3) then
               write(ulog,"(1x,a,1p,3e12.4)") "vertex(4,1:2):..............",  &
                  vertex(4,1),vertex(4,2)
            endif
            if (n_vertices>4) then
               write(ulog,"(1x,a,1p,3e12.4)") "vertex(5,1:2):..............",  &
                  vertex(5,1),vertex(5,2)
            endif
            if (n_vertices>5) then
               write(ulog,"(1x,a,1p,3e12.4)") "vertex(6,1:2):..............",  &
                  vertex(6,1),vertex(6,2)
            endif
            write(ulog,"(1x,a)")  " "
         endif
      endif
   enddo
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"Substations DATA",ninp,ulog))     &
      return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadSubstations

