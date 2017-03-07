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
! Program unit: ReadSectionFlowRate                                
! Description: Input management for the flow rate monitoring sections.                         
!-------------------------------------------------------------------------------
subroutine ReadSectionFlowRate(ainp,comment,nrighe,ier,ninp,nout)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module                            
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nrighe,ier,ninp,nout
character(1) :: comment
character(100) :: ainp,lcase
integer(4) :: n_fluid_types,ioerr,i,n_sect,n_vertices,section_ID
double precision :: dt_out,aux_dis,area
double precision :: plane_normal(3),vec_aux_1(3),vec_aux_2(3),vec_aux_3(3)
double precision :: vec_aux_4(3)
double precision :: vertex(4,3)
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
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"Section_flow_rate DATA",ninp,nout))  &
   return
do while (TRIM(lcase(ainp))/="##### end section flow rate #####")
! Reading the number of monitoring sections for the flow rate and their writing 
! time step
   read (ainp,*,iostat=ioerr) n_sect,dt_out,n_fluid_types
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"Section_flow_rate GENERAL INPUT", &
      ninp,nout)) return
! Writing the number of sections and the writing time step on the log file
   if (nout>0) then
      write (nout,"(1x,a,1p,i12)")    "n_sect:.............................",  &
         n_sect
      write (nout,"(1x,a,1p,e12.4)") "dt_out:.............................",   &
         dt_out
      write (nout,"(1x,a,1p,i12)")   "n_sect:.............................",   &
         n_fluid_types
      write (nout,"(1x,a)")  " "
   endif
! Allocation of the array of the flow rate monitoring sections
   if (allocated(Q_sections%section)) then
      else
         allocate(Q_sections%section(n_sect)) 
         Q_sections%n_sect = n_sect
         Q_sections%dt_out = dt_out
         Q_sections%n_fluid_types = n_fluid_types
! Initializing the auxiliary variable to print results
         Q_sections%it_out_last = 0
   endif
! Loop over the monitoring sections for flow rate
   do i=1,n_sect
! Reading the section parameters
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      read (ainp,*,iostat=ioerr) section_ID
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"section_ID",ninp,nout)) return
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      read (ainp,*,iostat=ioerr) n_vertices
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"n_vertices",ninp,nout)) return
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      read (ainp,*,iostat=ioerr) vertex(1,1),vertex(1,2),vertex(1,3)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"vertex_1",ninp,nout)) return
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      read (ainp,*,iostat=ioerr) vertex(2,1),vertex(2,2),vertex(2,3)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"vertex_2",ninp,nout)) return
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      read (ainp,*,iostat=ioerr) vertex(3,1),vertex(3,2),vertex(3,3)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"vertex_3",ninp,nout)) return
      if (n_vertices==4) then
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read (ainp,*,iostat=ioerr) vertex(4,1),vertex(4,2),vertex(4,3)
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"vertex_4",ninp,nout)) return           
      endif
! Assignation to the section parameters 
      Q_sections%section(i)%n_vertices = n_vertices
      Q_sections%section(i)%vertex(:,:) = vertex(:,:)
! Computation of the section areas       
      vec_aux_1(:) = vertex(1,:)
      vec_aux_2(:) = vertex(2,:)
      vec_aux_3(:) = vertex(3,:)
      vec_aux_4(:) = vertex(4,:)
      call area_quadrilateral(vec_aux_1,vec_aux_2,vec_aux_3,vec_aux_4,area)
      Q_sections%section(i)%area = area
      call dis_point_plane(vec_aux_1,vec_aux_1,vec_aux_2,vec_aux_3,aux_dis,    &
         plane_normal)
      Q_sections%section(i)%normal(:) = plane_normal(:)
! Writing on the log file
      if (ncord>0) then
         if (nout>0) then
            write (nout,"(1x,a,i12)")       "n_vertices:.................",    &
               n_vertices
            write (nout,"(1x,a,1p,e12.4)")  "area:.......................",    &
               Q_sections%section(i)%area
            write (nout,"(1x,a,1p,3e12.4)") "normal:.....................",    &
               Q_sections%section(i)%normal(1),Q_sections%section(i)%normal(2) &
               ,Q_sections%section(i)%normal(3)
            write (nout,"(1x,a,1p,3e12.4)") "vertex(1,1-3):..............",    &
               vertex(1,1),vertex(1,2),vertex(1,3)
            write (nout,"(1x,a,1p,3e12.4)") "vertex(2,1-3):..............",    &
               vertex(2,1),vertex(2,2),vertex(2,3)
            write (nout,"(1x,a,1p,3e12.4)") "vertex(3,1-3):..............",    &
               vertex(3,1),vertex(3,2),vertex(3,3)
            if (n_vertices==4) then
               write (nout,"(1x,a,1p,3e12.4)") "vertex(4,1-3):..............", &
                  vertex(4,1),vertex(4,2),vertex(4,3)       
            endif
            write (nout,"(1x,a)")  " "
         endif
      endif
   enddo         
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"Section_flow_rate DATA",ninp,     &
      nout)) return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadSectionFlowRate

