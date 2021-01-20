!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2021 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: ReadDBSPH                    
! Description: Reading input data for the DB-SPH boundary treatment scheme 
!              (Amicarelli et al., 2013, IJNME).                   
!-------------------------------------------------------------------------------
subroutine ReadDBSPH(ainp,comment,nrighe,ier,ninp,ulog)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(inout) :: nrighe,ier,ninp,ulog
character(1),intent(inout) :: comment
character(len=lencard),intent(inout) :: ainp
logical :: MUSCL_boundary_flag,in_built_monitors,Gamma_limiter_flag
logical :: negative_wall_p_allowed,FS_allowed
integer(4) :: ioerr,n_monitor_points,n_monitor_regions,i,alloc_stat   
integer(4) :: dealloc_stat,j,n_inlet,n_outlet,slip_ID
integer(4) :: ply_n_face_vert,surface_mesh_files
double precision :: dx_dxw,k_w
character(100) :: lcase
integer(4) :: n_kinematics_records(100)
integer(4),allocatable,dimension(:) :: monitor_IDs
double precision :: monitor_region(6)
double precision :: rotation_centre(100,3)
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
n_kinematics_records(:) = 0
!------------------------
! Statements
!------------------------
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH DATA",ninp,ulog)) return
do while (trim(lcase(ainp))/="##### end dbsph #####")
! Reading the ratio between the fluid and the semi-particle sizes (dx/dx_w)
   read(ainp,*,iostat=ioerr) dx_dxw,MUSCL_boundary_flag,k_w,slip_ID,           &
      Gamma_limiter_flag
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH GENERAL INPUT",ninp,ulog))  &
      return
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   read(ainp,*,iostat=ioerr) negative_wall_p_allowed,FS_allowed
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH GENERAL INPUT 2",ninp,ulog))&
      return
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   read(ainp,*,iostat=ioerr) n_monitor_points,n_monitor_regions
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH_monitor_numbers",ninp,ulog))&
      return
   if (n_monitor_points>0) then
      if (.not.allocated(monitor_IDs)) allocate(monitor_IDs(n_monitor_points), &
         STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(ulog,*)                                                         &
'Allocation of monitor_IDs in ReadDBSPH failed; the program terminates here'
! Stop the main program
         stop 
      endif
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      read(ainp,*,iostat=ioerr) monitor_IDs(:)
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH_monitor_IDs",ninp,ulog)) &
         return
      endif
      if (n_monitor_regions==1) then
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) monitor_region(:)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH_monitor_region",ninp, &
            ulog)) return
      endif
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      read(ainp,*,iostat=ioerr) surface_mesh_files,in_built_monitors
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"SURFACE_MESH_FILES",ninp,ulog))&
         return  
      do i=1,surface_mesh_files
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) n_kinematics_records(i),                    &
            rotation_centre(i,1:3)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH_KINEMATICS",ninp,ulog &
            )) return
         if ((.not.(allocated(DBSPH%kinematics))).and.                         &
            (n_kinematics_records(i)/=0)) then
            allocate(DBSPH%kinematics(surface_mesh_files,                      &
               n_kinematics_records(i),7),STAT=alloc_stat)
            if (alloc_stat/=0) then
               write(ulog,*)                                                   &
'Error! Allocation of DBSPH%kinematics in ReadDBSPH failed; the program terminates here.'
            call diagnostic(arg1=5,arg2=340)
! Stop the main program
               stop
               else
                  write(ulog,'(1x,a)')                                         &
"Array DBSPH%kinematics successfully allocated in subroutine ReadDBSPH."
            endif
         endif
         do j=1,n_kinematics_records(i)
            call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
            read(ainp,*,iostat=ioerr) DBSPH%kinematics(i,j,1:7)
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH_KINEMATICS_RECORDS"&
               ,ninp,ulog)) return
         enddo
      enddo
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      read(ainp,*,iostat=ioerr) n_inlet,n_outlet,ply_n_face_vert
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                                &
         "DBSPH_INLET_OUTLET_PLY_N_FACE_VERT",ninp,ulog)) return
      if (n_inlet>0) then
         if (.not.allocated(DBSPH%inlet_sections)) then
            allocate(DBSPH%inlet_sections(n_inlet,10),STAT=alloc_stat)
            if (alloc_stat/=0) then
               write(ulog,*)                                                   &
'Allocation of DBSPH%inlet_sections in ReadDBSPH failed; the program terminates here'
! Stop the main program
               stop 
         endif
      endif
   endif
   do j=1,n_inlet
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
! Reading position, normal and velocity of an inlet surface element      
      read(ainp,*,iostat=ioerr) DBSPH%inlet_sections(j,:)  
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH_INLET_SECTIONS",ninp,    &
         ulog)) return            
   enddo 
   if (n_outlet>0) then
! Reading position and normal of an outlet surface element       
      if (.not.allocated(DBSPH%outlet_sections)) then
         allocate(DBSPH%outlet_sections(n_outlet,8),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(ulog,*)                                                      &
'Allocation of DBSPH_outlet_sections in ReadDBSPH failed; the program terminates here'
! Stop the main program
            stop 
         endif
      endif   
   endif
   do j=1,n_outlet
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      read(ainp,*,iostat=ioerr) DBSPH%outlet_sections(j,:)  
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH_OUTLET_SECTIONS",ninp,   &
         ulog)) return            
   enddo
! Writing the DB-SPH input parameters on the log file
   if ((input_second_read.eqv..true.).and.(ulog>0)) then
      write(ulog,"(1x,a,1p,e12.4)")                                            &
"dx/dx_w:........................",dx_dxw
      write(ulog,"(1x,a,1p,l12)")                                              &
"MUSCL_boundary_flag:............",MUSCL_boundary_flag
      write(ulog,"(1x,a,1p,e12.4)")                                            &
"k_w(semi-particle coefficient)..",k_w
      write(ulog,"(1x,a,1p,i12)")                                              &
"DB-SPH slip_ID..................",slip_ID       
      write(ulog,"(1x,a,1p,l12)")                                              &
"Gamma_limiter_flag:.............",Gamma_limiter_flag
      write(ulog,"(1x,a,1p,l12)")                                              &
"negative_wall_p_allowed:........",negative_wall_p_allowed
      write(ulog,"(1x,a,1p,l12)")                                              &
"FS_allowed:.....................",FS_allowed
      write(ulog,"(1x,a,1p,i12)")                                              &
"n_monitor_points................",n_monitor_points       
      if (n_monitor_points>0) then
         do i=1,n_monitor_points
            write(ulog,"(1x,a,1p,i12)")                                        &
"ID_monitor......................",monitor_IDs(i)        
         enddo    
      endif
      write(ulog,"(1x,a,1p,i12)")                                              &
"n_monitor_regions...............",n_monitor_regions        
      if (n_monitor_regions==1) then
         write(ulog,"(1x,a,1p,g12.5)")                                         &
"monitor_region_x_min: ..........",monitor_region(1)
         write(ulog,"(1x,a,1p,g12.5)")                                         &
"monitor_region_x_max: ..........",monitor_region(2)
         write(ulog,"(1x,a,1p,g12.5)")                                         &
"monitor_region_y_min: ..........",monitor_region(3)
         write(ulog,"(1x,a,1p,g12.5)")                                         &
"monitor_region_y_max: ..........",monitor_region(4)
         write(ulog,"(1x,a,1p,g12.5)")                                         &
"monitor_region_z_min: ..........",monitor_region(5)
         write(ulog,"(1x,a,1p,g12.5)")                                         &
"monitor_region_z_max: ..........",monitor_region(6)
      endif
      write(ulog,"(1x,a,1p,i12)")                                              &
"surface_mesh_files..............",surface_mesh_files
      write(ulog,"(1x,a,1p,l12)")                                              &
"in-built_monitor_flag:..........",in_built_monitors
      do j=1,surface_mesh_files
         write(ulog,"(1x,a,1p,i12)")                                           &
"n_kinematics_records............",n_kinematics_records(j)
         write(ulog,"(1x,a,1p,3(g12.4))")                                      &
"rotation_centre:................",rotation_centre(j,1:3)
         do i=1,n_kinematics_records(j)
            write(ulog,"(1x,a,1p,7(g12.4))")                                   &
"time(s),u(m/s),v(m/s),w(m/s),omega_x(rad/s),omega_y(rad/s),omega_z(rad/s):..."&
               ,DBSPH%kinematics(j,i,1:7)
         enddo
      enddo 
      write(ulog,"(1x,a,i12)")                                                 &
"n_inlet:........................",n_inlet
      do i=1,n_inlet
         write(ulog,"(1x,a,1p,9(g12.4))")                                      &
"x(m),y(m),z(m),n_x,n_y,n_z,u(m/s),v(m/s),w(m/s),length(m): ",                 &
            DBSPH%inlet_sections(i,:)        
      enddo 
      write(ulog,"(1x,a,i12)")                                                 &
"n_outlet:.......................",n_outlet
      do i=1,n_outlet
         write(ulog,"(1x,a,1p,6(g12.4))")                                      &
"x(m),y(m),z(m),n_x,n_y,n_z,length(m),pres(Pa)............: ",                 &
            DBSPH%outlet_sections(i,:)        
      enddo
      write(ulog,"(1x,a,i12)")                                                 &
"ply_n_face_vert:................",ply_n_face_vert
      write(ulog,"(1x,a)")  " "
! Assignment of the DB-SPH parameters 
      DBSPH%dx_dxw = dx_dxw
      DBSPH%MUSCL_boundary_flag = MUSCL_boundary_flag
      DBSPH%k_w = k_w
      DBSPH%slip_ID = slip_ID
      DBSPH%Gamma_limiter_flag = Gamma_limiter_flag
      DBSPH%negative_wall_p_allowed = negative_wall_p_allowed
      DBSPH%FS_allowed = FS_allowed
      DBSPH%n_monitor_points = n_monitor_points 
      DBSPH%n_monitor_regions = n_monitor_regions
      DBSPH%monitor_region(:) = monitor_region(:)  
      if (n_monitor_points>0) then
         if (.not.(allocated(DBSPH%monitor_IDs))) then
            allocate(DBSPH%monitor_IDs(n_monitor_points),STAT=alloc_stat)
            if (alloc_stat/=0) then
               write(ulog,*)                                                   &
'Allocation of DBSPH%n_monitor_points in ReadDBSPH failed; the program terminates here.'
! Stop the main program
               stop 
            endif   
         endif       
         DBSPH%monitor_IDs(:) = monitor_IDs(:)
      endif
      DBSPH%surface_mesh_files = surface_mesh_files
      DBSPH%in_built_monitors = in_built_monitors
      if (.not.(allocated(DBSPH%n_kinematics_records))) then
         allocate(DBSPH%n_kinematics_records(surface_mesh_files),              &
            STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(ulog,*)                                                      &
'Error! Allocation of DBSPH%n_kinematics_records in ReadDBSPH failed; the program terminates here.'
            call diagnostic(arg1=5,arg2=340)
! Stop the main program
            stop
            else
               write(ulog,'(1x,a)')                                            &
"Array DBSPH%n_kinematics_records successfully allocated in subrouitne ReadDBSPH."
         endif
      endif
      if (.not.(allocated(DBSPH%rotation_centre))) then
         allocate(DBSPH%rotation_centre(surface_mesh_files,3),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(ulog,*)                                                      &
'Error! Allocation of DBSPH%rotation_centre in ReadDBSPH failed; the program terminates here.'
            call diagnostic(arg1=5,arg2=340)
! Stop the main program
            stop
            else
               write(ulog,'(1x,a)')                                            &
"Array DBSPH%rotation_centre successfully allocated in subrouitne ReadDBSPH."
         endif
      endif
      do i=1,surface_mesh_files
         DBSPH%n_kinematics_records(i) = n_kinematics_records(i)
         DBSPH%rotation_centre(i,:) = rotation_centre(i,:)
      enddo
      DBSPH%n_inlet = n_inlet   
      DBSPH%n_outlet = n_outlet
      DBSPH%ply_n_face_vert = ply_n_face_vert
   endif
   if (allocated(monitor_IDs)) then
      deallocate(monitor_IDs,STAT=dealloc_stat)
      if (dealloc_stat/=0) then
         write(ulog,*)                                                         &
'Deallocation of monitor_IDs in ReadDBSPH failed; the program terminates here.'
! Stop the main program
         stop 
      endif   
   endif   
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DBSPH DATA",ninp,ulog)) return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadDBSPH
