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
! Program unit: ReadRestartFile
! Description: To read the restart file.                       
!-------------------------------------------------------------------------------
subroutine ReadRestartFile(option,ier,nrecords)
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
character(7),intent(IN) :: option
integer(4),intent(INOUT) :: ier,nrecords
integer(4) :: restartcode,save_istart,ioerr,i,alloc_stat,i_aux,aux_integer
double precision :: save_start
character(12) :: ainp = "Restart File"
character(len=8) :: versionerest
logical,external :: ReadCheck
character(100),external :: lcase
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
ier = 0
aux_integer = 0
!------------------------
! Statements
!------------------------
! Restart heading 
if (TRIM(lcase(option))==TRIM(lcase("heading"))) then
   rewind(nsav)
   write(nout,'(a)')    "-------------------"
   write(nout,"(1x,a)") ">> Restart heading."
   write(nout,'(a)')    "-------------------"
   read(nsav,iostat=ioerr) versionerest,nrecords
   if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"versionerest,nrecords",nsav,    &
      nout)) return
! Check the program version
   if (TRIM(lcase(version))/=TRIM(lcase(versionerest))) then
      write(nscr,'(a)')                                                        &
         "---------------------------------------------------------------"
      write(nscr,"(1x,a)")                                                     &
         ">> ERROR! The Restart version is not equal the current version."
      write(nscr,"(1x,a)") ">>        The Run is stopped."
      write(nscr,'(a)')                                                        &
         "---------------------------------------------------------------"
      flush(nscr)
      stop
   endif
   read(nsav,iostat=ioerr) ncord,nag,NMedium,NPartZone,NumVertici,NumFacce,    &
      NumTratti,NumBVertices,NumBSides,GCBFVecDim,Grid%nmax,NPointst,NPoints,  &
      NPointsl,NPointse,NLines,NSections,doubleh
   if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"ncord, nag, ...",nsav,nout))    &
      return
! The parameter is read from the restart file and not from the input file
   write(nout,"(1x,a,1p,i12)")   "GCBFVecDim (from restart file): ",GCBFVecDim
! Allocation of the array "GCBFVector"
   if ((Domain%tipo=="semi").and.(GCBFVecDim>0).and.                           &
      (.not.allocated(GCBFVector))) then
      allocate(GCBFVector(GCBFVecDim),stat=ier)    
      if (ier/=0) then
         write(nout,'(1x,2a)') "Allocation of GCBFVector in ",                 &
            "ReadRestartFile failed. "
         else
            write(nout,'(1x,2a)') "Allocation of GCBFVector in ",              &
            "ReadRestartFile successully completed. "
      endif
   endif
   elseif (TRIM(lcase(option))=="reading") then
      write(nout,'(a)')                                                        &
         "---------------------------------------------------------------------"
      write(nout,"(1x,a)")                                                     &
         ">> Restart reading:  step         time      interval    num.particles"
      save_istart = Domain%istart
      save_start = Domain%start
! Since here, Domain%istart and Domain%start are zeroed, until next reading of 
! the main input file.
      read(nsav,iostat=ioerr) Domain
      if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"domain",nsav,nout)) return
      read(nsav,iostat=ioerr) Grid
      if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"Grid",nsav,nout)) return
! Allocating the 2D matrix to detect free surface
      if (.not.allocated(ind_interfaces)) then
         allocate(ind_interfaces(Grid%ncd(1),Grid%ncd(2),6),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(nout,*)                                                      &
            'Allocation of ind_interfaces in ReadRestartFile failed;',         &
            ' the program terminates here.'
            stop
            else
               write (nout,*)                                                  &
                  'Allocation of ind_interfaces in ReadRestartFile ',          &
                  'successfully completed.'
         endif
      endif
! Allocation of the array "GCBFPointers"
      if (Domain%tipo=="semi") then
         allocate(GCBFPointers(Grid%nmax,2),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(nout,*) "Allocation of GCBFPointers in ReadRestartFile ",    &
               "failed; the program stops here."
            stop
            else
               write (nout,*) "Allocation of GCBFPointers in ReadRestartFile ",&
                  "is successfully completed."
         endif
      endif
      read(nsav,iostat=ioerr) Med(1:NMedium)
      if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"Med",nsav,nout)) return
      if (NumVertici>0) then
         read(nsav,iostat=ioerr) Vertice(1:SPACEDIM,1:NumVertici)
         if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"Vertice",nsav,nout)) return
      endif
      if (NumFacce>0) then 
         read(nsav,iostat=ioerr) BoundaryFace(1:NumFacce)
         if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"BoundaryFace",nsav,nout)) &
            return
      endif
      if (NumFacce>0) then
         read(nsav,iostat=ioerr) BFaceList(1:NumFacce)
         if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"BFaceList",nsav,nout))    &
            return
      endif
      if (NumTratti>0) then
         read(nsav,iostat=ioerr) Tratto(1:NumTratti)
         if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"Tratto",nsav,nout)) return
      endif
      if (NPartZone>0) then
         read(nsav,iostat=ioerr) Partz(1:NPartZone)
         if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"Partz",nsav,nout)) return
      endif
      if (NumBVertices>0) then
        read(nsav,iostat=ioerr) BoundaryVertex(1:NumBVertices)
        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"BoundaryVertex",nsav,nout))&
           return
      endif
      if (NumBSides>0) then
         read(nsav,iostat=ioerr) BoundarySide(1:NumBSides)     
         if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"BoundarySide",nsav,nout)) &
            return
      endif
      if (Domain%tipo=="semi") then
         if (GCBFVecDim>0) then
            read(nsav,iostat=ioerr) GCBFVector(1:GCBFVecDim)
            if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"GCBFVector",nsav,nout))&
               return
         endif
         if (Grid%nmax>0) then
            read(nsav,iostat=ioerr) GCBFPointers(1:Grid%nmax,1:2)
            if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"GCBFPointers",nsav,    &
               nout)) return
         endif
      endif
! Allocation of the array of the maximum water depth
      do i_aux=1,NPartZone
         if(Partz(i_aux)%IC_source_type==2) then
            aux_integer = Partz(i_aux)%ID_last_vertex -                        &
                          Partz(i_aux)%ID_first_vertex + 1
            exit
         endif
      enddo
      if ((aux_integer>0).and.(.not.allocated(Z_fluid_max))) then
         allocate(Z_fluid_max(Grid%ncd(1)*Grid%ncd(2)),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(nout,*)                                                      &
            'Allocation of Z_fluid_max in ReadRestartFile failed;',            &
            ' the program terminates here.'
            stop
            else
               write (nout,*)                                                  &
                  'Allocation of Z_fluid_max in ReadRestartFile successfully', &
                  ' completed.'
         endif
      endif
! Allocation of the array of the maximum specific flow rate
      if ((aux_integer>0).and.(.not.allocated(q_max))) then
         allocate(q_max(aux_integer),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(nout,*)                                                      &
            'Allocation of q_max in ReadRestartFile failed;',                  &
            ' the program terminates here.'
            stop
            else
               write (nout,*)                                                  &
                  'Allocation of q_max in ReadRestartFile successfully ',      &
                  'completed.'
         endif
      endif
! Allocation of the 2D array of the minimum saturation flag (bed-load transport)
      if ((Granular_flows_options%ID_erosion_criterion>0).and.                 &
         (.not.allocated(Granular_flows_options%minimum_saturation_flag))) then
         allocate(Granular_flows_options%minimum_saturation_flag(Grid%ncd(1),  &
            Grid%ncd(2)),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(nout,*)                                                      &
            'Allocation of Granular_flows_options%minimum_saturation_flag ',   &
            ' in ReadRestartFile failed; the program stops here.'
            stop 
            else
               write (nout,*)                                                  &
                  'Allocation of ',                                            &
                  'Granular_flows_options%minimum_saturation_flag in ',        &
                  'ReadRestartFile is successfully completed.'
         endif
      endif
! Allocation of the 2D array of the maximum saturation flag (bed-load transport)
      if ((Granular_flows_options%ID_erosion_criterion>0).and.                 &
         (.not.allocated(Granular_flows_options%maximum_saturation_flag))) then
         allocate(Granular_flows_options%maximum_saturation_flag(Grid%ncd(1),  &
            Grid%ncd(2)),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(nout,*) "Allocation of ",                                    &
               "Granular_flows_options%maximum_saturation_flag in ",           &
               "ReadRestartFile failed; the program stops here."
            stop 
            else
               write (nout,*) "Allocation of ",                                &
                  "Granular_flows_options%maximum_saturation_flag in ",        &
                  "ReadRestartFile is successfully completed."
         endif
      endif      
! Restart positions are based on the step number
      it_start = 0 
      if (save_istart>0) then
! It reads all the saved steps and overwrite the restart values until the
! restart time is reached.      
         do while (save_istart>it_start)
            read(nsav,iostat=ioerr) it_start,simulation_time,dt,nag,ncord,     &
               restartcode
            if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,                        &
               "it_start,simulation_time,dt,nag,ncord,restartcode",nsav,nout)) &
               return
            write(nout,"(16x,i10,2(2x,g12.5),7x,i10)") it_start,               &
               simulation_time,dt,nag
            flush(nout)
            if (it_start<save_istart) then
               read(nsav,iostat=ioerr) 
               if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"...",nsav,nout))    &
                  return  
               if (restartcode==1) then                 
                  if (allocated(pg_w)) then
                     read(nsav,iostat=ioerr) 
                     if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"pg_w",        &
                        nsav,nout)) return
                  endif
                  if (n_bodies>0) then
                     do i=1,n_bodies
                        read(nsav,iostat=ioerr)
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,            &
                           "body_arr_1_of_2",nsav,nout)) return
                        if (body_arr(i)%n_records>0) then
                           read(nsav,iostat=ioerr)
                           if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,         &
                              "body_arr_2_of_2",nsav,nout)) return
                        endif                  
                     enddo
                     read(nsav,iostat=ioerr) 
                     if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"bp_arr",      &
                              nsav,nout)) return
                     read(nsav,iostat=ioerr) 
                     if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,               &
                              "surf_body_part",nsav,nout)) return
                  endif
                  if (allocated(Z_fluid_max)) then
                     read(nsav,iostat=ioerr)
                     if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,               &
                        "Z_fluid_max",nsav,nout)) return
                  endif
                  if (allocated(q_max)) then
                     read(nsav,iostat=ioerr) 
                     if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"q_max",       &
                        nsav,nout)) return
                  endif  
                  if (allocated(Granular_flows_options%minimum_saturation_flag)&
                     ) then
                     read(nsav,iostat=ioerr) 
                     if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,               &
                        "minimum_saturation_flag",nsav,nout)) return
                  endif
                  if (allocated(Granular_flows_options%maximum_saturation_flag)&
                     ) then
                     read(nsav,iostat=ioerr) 
                     if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,               &
                        "maximum_saturation_flag",nsav,nout)) return
                  endif  
               endif              
               else
! Actual array reading for restart
                  if (restartcode==1) then
                     read(nsav,iostat=ioerr) pg(1:nag)
                     if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"pg",nsav,     &
                        nout)) return
                     if (allocated(pg_w)) then
                        read(nsav,iostat=ioerr) pg_w(1:DBSPH%n_w+DBSPH%n_inlet+&
                           DBSPH%n_outlet)
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"pg_w",nsav,&
                           nout)) return
                     endif
                     if (n_bodies>0) then
                        do i=1,n_bodies
                           read(nsav,iostat=ioerr) body_arr(i)%npart,          &
                              body_arr(i)%Ic_imposed,                          &
                              body_arr(i)%imposed_kinematics,                  &
                              body_arr(i)%n_records,body_arr%mass,             &
                              body_arr(i)%umax,body_arr(i)%pmax,               &
                              body_arr(i)%x_CM(1:3),body_arr(i)%alfa(1:3),     &
                              body_arr(i)%u_CM(1:3),body_arr(i)%omega(1:3),    &
                              body_arr(i)%Force(1:3),body_arr(i)%Moment(1:3),  &
                              body_arr(i)%Ic(1:3,1:3),                         &
                              body_arr(i)%Ic_inv(1:3,1:3)
                           if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,         &
                              "body_arr_1_of_2",nsav,nout)) return                              
                           if (body_arr(i)%n_records>0) then
                              if (.not.allocated(body_arr(i)%body_kinematics)) &
allocate(body_arr(i)%body_kinematics(body_arr(i)%n_records,7))
                              read(nsav,iostat=ioerr)                          &
body_arr(i)%body_kinematics(1:body_arr(i)%n_records,1:7)
                              if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,      &
                                 "body_arr_2_of_2",nsav,nout)) return
                           endif
                        enddo
                        read(nsav,iostat=ioerr) bp_arr(1:n_body_part)
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"bp_arr",   &
                           nsav,nout)) return
                        read(nsav,iostat=ioerr) surf_body_part(1:              &
                           n_surf_body_part)
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,            &
                          "surf_body_part",nsav,nout)) return
                     endif
                     if (allocated(Z_fluid_max)) then
                        read(nsav,iostat=ioerr)                                &
                           Z_fluid_max(1:Grid%ncd(1)*Grid%ncd(2))
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,            &
                           "Z_fluid_max",nsav,nout)) return
                     endif
                     if (allocated(q_max)) then
                        read(nsav,iostat=ioerr) q_max(1:size(q_max))
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"q_max",nsav&
                           ,nout)) return
                     endif   
                     if (allocated                                             &
                        (Granular_flows_options%minimum_saturation_flag)) then
                        read(nsav,iostat=ioerr)                                &
                           Granular_flows_options%minimum_saturation_flag(     &
                           1:Grid%ncd(1),1:Grid%ncd(2))
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,            &
                           "minimum_saturation_flag",nsav,nout)) return
                     endif 
                     if (allocated                                             &
                        (Granular_flows_options%maximum_saturation_flag)) then
                        read(nsav,iostat=ioerr)                                &
                           Granular_flows_options%maximum_saturation_flag(     &
                           1:Grid%ncd(1),1:Grid%ncd(2))
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,            &
                           "maximum_saturation_flag",nsav,nout)) return
                     endif                                           
                     write(nout,'(a)') " "
                     write(nout,'(a,i10,a,g12.5)') "   Located Restart Step :",&
                        it_start,"   Time :",simulation_time
                     flush(nout)
! Reading for post-processing
                     elseif (restartcode==0) then
                        read(nsav,iostat=ioerr) pg(1:nag)%coord(1),            &
                           pg(1:nag)%coord(2),pg(1:nag)%coord(3),              &
                           pg(1:nag)%vel(1),pg(1:nag)%vel(2),pg(1:nag)%vel(3), &
                           pg(1:nag)%pres,pg(1:nag)%dens,pg(1:nag)%mass,       &
                           pg(1:nag)%visc,pg(1:nag)%IntEn,pg(1:nag)%VolFra,    &
                           pg(1:nag)%imed,pg(1:nag)%icol
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"pg",nsav,  &
                           nout)) return
                        write(nscr,'(a)') " "
                        write(nscr,'(a,i10,a,g12.5)')                          &
                           "   Located Result Step :",it_start,"   Time :",    &
                           simulation_time
                        flush(nscr)
                        write(nscr,'(a)')                                      &
"       But this step is not a restart step. Check the correct step for restart in the restart file."
                        flush(nscr)
                        write(nscr,'(a)') " The program is terminated."
                        flush(nscr)
                        stop
                  endif
                  return
            endif
         enddo
         write(nscr,'(a,i10,a)') "   Restart Step Number:",it_start,           &
            " has not been found"
         ier = 3
! Restart positions are based on the step number
         elseif (save_start>zero) then
            simulation_time = zero
            do while (save_start>simulation_time)
               read(nsav,iostat=ioerr) it_start,simulation_time,dt,nag,ncord,  &
                  restartcode
               if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,                     &
                  "it_start,simulation_time,dt,nag,ncord,restartcode",nsav,    &
                  nout)) return
               write(nout,"(16x,i10,2(2x,g12.5),7x,i10)") it_start,            &
                  simulation_time,dt,nag
               flush(nout)
               if (simulation_time<save_start) then
                  read(nsav,iostat=ioerr)
                  if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"...",nsav,nout)) &
                     return
                  if (restartcode==1) then
                     if (allocated(pg_w)) then
                        read(nsav,iostat=ioerr) 
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"pg_w",     &
                           nsav,nout)) return
                     endif
                     if (n_bodies>0) then
                        do i=1,n_bodies
                           read(nsav,iostat=ioerr)
                           if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,         &
                              "body_arr_1_of_2",nsav,nout)) return
                           if (body_arr(i)%n_records>0) then   
                              read(nsav,iostat=ioerr)
                              if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,      &
                                 "body_arr_2_of_2",nsav,nout)) return
                           endif
                        enddo
                        read(nsav,iostat=ioerr) 
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"bp_arr",   &
                                 nsav,nout)) return
                        read(nsav,iostat=ioerr)
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,            &
                                 "surf_body_part",nsav,nout)) return
                     endif
                     if (allocated(Z_fluid_max)) then
                        read(nsav,iostat=ioerr)
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,            &
                           "Z_fluid_max",nsav,nout)) return
                     endif
                     if (allocated(q_max)) then
                        read(nsav,iostat=ioerr) 
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"q_max",    &
                           nsav,nout)) return
                     endif
                     if (allocated                                             &
                        (Granular_flows_options%minimum_saturation_flag)) then
                        read(nsav,iostat=ioerr) 
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,            &
                        "minimum_saturation_flag",nsav,nout)) return
                     endif
                     if (allocated                                             &
                        (Granular_flows_options%maximum_saturation_flag)) then
                        read(nsav,iostat=ioerr) 
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,            &
                        "maximum_saturation_flag",nsav,nout)) return
                     endif
                  endif 
                  else
! Actual array reading for restart
                     if (restartcode==1) then
                        read(nsav,iostat=ioerr) pg(1:nag)
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"pg",nsav,  &
                           nout)) return
                        if (allocated(pg_w)) then
                           read(nsav,iostat=ioerr) pg_w(1:                     &
                              DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet)
                           if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"pg_w",  &
                              nsav,nout)) return
                        endif
                        if (n_bodies>0) then
                           do i=1,n_bodies
                              read(nsav,iostat=ioerr) body_arr(i)%npart,       &
                                 body_arr(i)%Ic_imposed,                       &
                                 body_arr(i)%imposed_kinematics,               &
                                 body_arr(i)%n_records,body_arr%mass,          &
                                 body_arr(i)%umax,body_arr(i)%pmax,            &
                                 body_arr(i)%x_CM(1:3),body_arr(i)%alfa(1:3),  &
                                 body_arr(i)%u_CM(1:3),body_arr(i)%omega(1:3), &
                                 body_arr(i)%Force(1:3),                       &
                                 body_arr(i)%Moment(1:3),                      &
                                 body_arr(i)%Ic(1:3,1:3),                      &
                                 body_arr(i)%Ic_inv(1:3,1:3)
                              if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,      &
                                 "body_arr_1_of_2",nsav,nout)) return
                              if (body_arr(i)%n_records>0) then
                               if (.not.allocated(body_arr(i)%body_kinematics))&
allocate(body_arr(i)%body_kinematics(body_arr(i)%n_records,7))
                                 read(nsav,iostat=ioerr)                       &
body_arr(i)%body_kinematics(1:body_arr(i)%n_records,1:7)
                                 if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,   &
                                 "body_arr_2_of_2",nsav,nout)) return
                              endif
                           enddo
                           read(nsav,iostat=ioerr) bp_arr(1:n_body_part)
                           if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"bp_arr",&
                              nsav,nout)) return
                           read(nsav,iostat=ioerr) surf_body_part(1:           &
                              n_surf_body_part)
                           if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,         &
                              "surf_body_part",nsav,nout)) return
                        endif
                        if (allocated(Z_fluid_max)) then
                           read(nsav,iostat=ioerr) Z_fluid_max(1:Grid%ncd(1)*  &
                              Grid%ncd(2))
                           if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,         &
                              "Z_fluid_max",nsav,nout)) return
                        endif
                        if (allocated(q_max)) then
                           read(nsav,iostat=ioerr) q_max(1:size(q_max))
                           if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"q_max", &
                              nsav,nout)) return
                        endif
if (allocated(Granular_flows_options%minimum_saturation_flag)) then
                           read(nsav,iostat=ioerr)                             &
                              Granular_flows_options%minimum_saturation_flag(  &
                              1:Grid%ncd(1),1:Grid%ncd(2))
                           if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,         &
                              "minimum_saturation_flag",nsav,nout)) return
                        endif 
                        if (allocated                                          &
                           (Granular_flows_options%maximum_saturation_flag)) &
                           then
                           read(nsav,iostat=ioerr)                             &
                              Granular_flows_options%maximum_saturation_flag(  &
                              1:Grid%ncd(1),1:Grid%ncd(2))
                           if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,         &
                              "maximum_saturation_flag",nsav,nout)) return
                        endif
                        write(nout,'(a)') 
                        write(nout,'(a,i10,a,g12.5)')                          &
                           "   Located Restart Step :",it_start,"   Time :",   &
                           simulation_time
                        flush(nout)
! Reading for post-processing
                        elseif (restartcode==0) then
                           read(nsav,iostat=ioerr) pg(1:nag)%coord(1),         &
                              pg(1:nag)%coord(2),pg(1:nag)%coord(3),           &
                              pg(1:nag)%vel(1),pg(1:nag)%vel(2),               &
                              pg(1:nag)%vel(3),pg(1:nag)%pres,pg(1:nag)%dens,  &
                              pg(1:nag)%mass,pg(1:nag)%visc,pg(1:nag)%IntEn,   &
                              pg(1:nag)%VolFra,pg(1:nag)%imed,pg(1:nag)%icol
                           if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"pg",    &
                              nsav,nout)) return
                           write(nout,'(a)') 
                           write(nout,'(a,i10,a,g12.5)')                       &
                              "   Located Result Time :",it_start,"   Time :", &
                              simulation_time
                           flush(nout)
                           write(nscr,'(a)')                                   &
"       But this time is not a restart time. Check the correct time for restart in the restart file."
                           flush(nscr)
                           write(nscr,'(a)') " The program is terminated."
                           flush(nscr)
                           stop
                  endif 
                  return
               endif
            enddo
            write(nscr,'(a,i10,a)') "   Restart Time Step:",save_start,        &
               " has not been found"
            ier = 3
            else
               write (nscr,'(a)') "  > Restart cannot be read at step:",       &
                  it_start,"  time:",simulation_time
               ier = 4
      endif
      write (nout,'(a)') "  > Restart read successfully at step:",it_start,    &
         "  time:",simulation_time
      else
         ier = 5
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadRestartFile

