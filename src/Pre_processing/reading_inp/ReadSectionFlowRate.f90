!-------------------------------------------------------------------------------
! SPHERA v.10.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.10.0.0
! SPHERA v.10.0.0 is free software: you can redistribute it and/or modify
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
#ifdef SPACE_3D
subroutine ReadSectionFlowRate(ainp,comment,nrighe,ier)
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
integer(4) :: nrighe,ier,alloc_stat
character(1) :: comment
character(100) :: lcase
character(len=lencard) :: ainp
integer(4) :: n_fluid_types,ioerr,i,n_sect,n_vertices,section_ID
double precision :: dt_out,aux_dis,area
double precision :: plane_normal(3),vec_aux_1(3),vec_aux_2(3),vec_aux_3(3)
double precision :: vec_aux_4(3)
double precision :: vertex(4,3)
logical,external :: ReadCheck
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine ReadRiga(ninp,ainp,io_err,comment_sym,lines_treated)
      implicit none
      integer(4),intent(in) :: ninp
      character(*),intent(inout) :: ainp
      integer(4),intent(out) :: io_err
      character(1),intent(in),optional :: comment_sym
      integer(4),intent(inout),optional :: lines_treated
   end subroutine ReadRiga
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"Section_flow_rate DATA",ninp,ulog))  &
   return
do while (trim(lcase(ainp))/="##### end section flow rate #####")
! Reading the number of monitoring sections for the flow rate and their writing 
! time step
   read(ainp,*,iostat=ioerr) n_sect,dt_out,n_fluid_types
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"Section_flow_rate GENERAL INPUT", &
      ninp,ulog)) return
! Writing the number of sections and the writing time step on the log file
   if (ulog>0) then
      write(ulog,"(1x,a,1p,i12)")   "n_sect:.............................",    &
         n_sect
      write(ulog,"(1x,a,1p,e12.4)") "dt_out:.............................",    &
         dt_out
      write(ulog,"(1x,a,1p,i12)")   "n_fluid_types:......................",    &
         n_fluid_types
      write(ulog,"(1x,a)")  " "
   endif
! Allocation of the array of the flow rate monitoring sections
   if (.not.allocated(Q_sections%section)) then
      allocate(Q_sections%section(n_sect),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(uerr,*) 'Allocation of Q_sections%section failed; the ',        &
            'execution terminates here.'
         stop
         else
            write(ulog,*) 'Allocation of Q_sections%section is successfully ', &
               'completed.'
      endif
      Q_sections%n_sect = n_sect
      Q_sections%dt_out = dt_out
      Q_sections%n_fluid_types = n_fluid_types
! Initializing the auxiliary variable to print results
      Q_sections%it_out_last = 0
   endif
! Loop over the monitoring sections for flow rate
   do i=1,n_sect
! Reading the section variables
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
      read(ainp,*,iostat=ioerr) section_ID
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"section_ID",ninp,ulog)) return
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
      read(ainp,*,iostat=ioerr) n_vertices
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"n_vertices",ninp,ulog)) return
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
      read(ainp,*,iostat=ioerr) vertex(1,1),vertex(1,2),vertex(1,3)
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"vertex_1",ninp,ulog)) return
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
      read(ainp,*,iostat=ioerr) vertex(2,1),vertex(2,2),vertex(2,3)
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"vertex_2",ninp,ulog)) return
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
      read(ainp,*,iostat=ioerr) vertex(3,1),vertex(3,2),vertex(3,3)
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"vertex_3",ninp,ulog)) return
      if (n_vertices==4) then
         call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
         read(ainp,*,iostat=ioerr) vertex(4,1),vertex(4,2),vertex(4,3)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"vertex_4",ninp,ulog)) return           
      endif
! Assignation to the section variables 
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
      if ((input_second_read.eqv..true.).and.(ulog>0)) then
         write(ulog,"(1x,a,i12)")       "n_vertices:.................",        &
            n_vertices
         write(ulog,"(1x,a,1p,e12.4)")  "area:.......................",        &
            Q_sections%section(i)%area
         write(ulog,"(1x,a,1p,3e12.4)") "normal:.....................",        &
            Q_sections%section(i)%normal(1),Q_sections%section(i)%normal(2)    &
            ,Q_sections%section(i)%normal(3)
         write(ulog,"(1x,a,1p,3e12.4)") "vertex(1,1-3):..............",        &
            vertex(1,1),vertex(1,2),vertex(1,3)
         write(ulog,"(1x,a,1p,3e12.4)") "vertex(2,1-3):..............",        &
            vertex(2,1),vertex(2,2),vertex(2,3)
         write(ulog,"(1x,a,1p,3e12.4)") "vertex(3,1-3):..............",        &
            vertex(3,1),vertex(3,2),vertex(3,3)
         if (n_vertices==4) then
            write(ulog,"(1x,a,1p,3e12.4)") "vertex(4,1-3):..............",     &
               vertex(4,1),vertex(4,2),vertex(4,3)       
         endif
         write(ulog,"(1x,a)")  " "
      endif
   enddo         
   call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"Section_flow_rate DATA",ninp,     &
      ulog)) return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadSectionFlowRate
#endif
