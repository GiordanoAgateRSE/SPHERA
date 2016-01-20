!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-)



! SPHERA authors and email contact are provided on SPHERA documentation.

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
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: ReadRestartFile
! Description: To read the restart file.                       
!----------------------------------------------------------------------------------------------------------------------------------

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
integer(4) :: restartcode,save_istart,ioerr,i
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
!------------------------
! Statements
!------------------------
ier = 0
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
      NumTratti,NumBVertices,NumBSides,NPointst,NPoints,NPointsl,NPointse,     &
      NLines,NSections,GCBFVecDim,doubleh
   if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"ncord, nag, ...",nsav,nout))    &
      return
   elseif (TRIM(lcase(option))=="reading") then
      write(nout,'(a)')                                                        &
         "---------------------------------------------------------------------"
      write(nout,"(1x,a)")                                                     &
         ">> Restart reading:  step         time      interval    num.particles"
      save_istart = Domain%istart
      save_start = Domain%start
      read(nsav,iostat=ioerr) domain
      if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"domain",nsav,nout)) return
      read(nsav,iostat=ioerr) grid
      if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"grid",nsav,nout)) return
! Allocating the 2D matrix to detect free surface (erosion criterion)
      allocate(ind_interfaces(Grid%ncd(1),Grid%ncd(2),4),stat=ioerr)
      if (ioerr/=0) then
         write (nout,'(1x,a,i2)')                                              &
            "    Array ind_interfaces not allocated. Error code: ",ioerr
         stop ' routine ReadRestartFile'
         else
            write (nout,'(1x,a)')                                              &
               "    Array ind_interfaces successfully allocated "
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
      if (NumBSides>1) then
         read(nsav,iostat=ioerr) BoundarySide(1:NumBSides)     
         if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"BoundarySide",nsav,nout)) &
            return
         else
            read(nsav,iostat=ioerr) BoundarySide(1)
            if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"BoundarySide",nsav, &
               nout)) return
      endif
! Restart positions are based on the step number
      it_start = 0 
      if (save_istart>0) then
! It reads all the saved steps and overwrite the restart values until the
! restart time is reached.      
         do while (save_istart>it_start)
            read(nsav,iostat=ioerr) it_start,tempo,dt,nag,ncord,restartcode
            if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,                        &
               "it_start,tempo,dt,nag,ncord,restartcode",nsav,nout)) return
            write(nout,"(16x,i10,2(2x,g12.5),7x,i10)") it_start,tempo,dt,nag
            flush(nout)
            if (it_start<save_istart) then
               read(nsav,iostat=ioerr) 
               if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"...",nsav,nout))    &
                  return
               else
! Reading for restart
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
                              body_arr(i)%Ic_imposed,body_arr(i)%n_elem,       &
                              body_arr(i)%imposed_kinematics,                  &
                              body_arr(i)%n_records,body_arr%mass,             &
                              body_arr(i)%umax,body_arr(i)%pmax,               &
                              body_arr(i)%x_CM,body_arr(i)%alfa,               &
                              body_arr(i)%x_rotC,body_arr(i)%u_CM,             &
                              body_arr(i)%omega,body_arr(i)%Force,             &
                              body_arr(i)%Moment,body_arr(i)%Ic,               &
                              body_arr(i)%Ic_inv,body_arr(i)%body_kinematics,  &
                              body_arr(i)%elem
                        enddo
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"body_arr", &
                           nsav,nout)) return
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
                        read(nsav,iostat=ioerr)                                &
                           q_max(1:Grid%ncd(1)*Grid%ncd(2))
                        if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"q_max",nsav&
                           ,nout)) return
                     endif                        
                     write(nout,'(a)') " "
                     write(nout,'(a,i10,a,g12.5)') "   Located Restart Step :",&
                        it_start,"   Time :",tempo; flush(nout)
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
                           "   Located Result Step :",it_start,"   Time :",tempo
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
            tempo = zero
            do while (save_start>tempo)
               read(nsav,iostat=ioerr) it_start,tempo,dt,nag,ncord,restartcode
               if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,                     &
                  "it_start,tempo,dt,nag,ncord,restartcode",nsav,nout)) return
               write(nout,"(16x,i10,2(2x,g12.5),7x,i10)") it_start,tempo,dt,nag
               flush(nout)
               if (tempo<Domain%start) then
                  read(nsav,iostat=ioerr) 
                  if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"...",nsav,nout)) &
                     return
                  else
! Reading for restart
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
                                 body_arr(i)%Ic_imposed,body_arr(i)%n_elem,    &
                                 body_arr(i)%imposed_kinematics,               &
                                 body_arr(i)%n_records,body_arr%mass,          &
                                 body_arr(i)%umax,body_arr(i)%pmax,            &
                                 body_arr(i)%x_CM,body_arr(i)%alfa,            &
                                 body_arr(i)%x_rotC,body_arr(i)%u_CM,          &
                                 body_arr(i)%omega,body_arr(i)%Force,          &
                                 body_arr(i)%Moment,body_arr(i)%Ic,            &
                                 body_arr(i)%Ic_inv,body_arr(i)%body_kinematics&
                                 ,body_arr(i)%elem
                           enddo
                           if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,         &
                              "body_arr",nsav,nout)) return
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
                           read(nsav,iostat=ioerr) q_max(1:                    &
                              Grid%ncd(1)*Grid%ncd(2))
                           if (.NOT.ReadCheck(ioerr,ier,it_start,ainp,"q_max", &
                              nsav,nout)) return
                        endif                          
                        write(nout,'(a)') 
                        write(nout,'(a,i10,a,g12.5)')                          &
                           "   Located Restart Step :",it_start,"   Time :",   &
                           tempo
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
                              tempo
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
            write(nscr,'(a,i10,a)') "   Restart Time Step:",Domain%start,      &
               " has not been found"
            ier = 3
            else
               write (nscr,'(a)') "  > Restart cannot be read at step:",       &
                  it_start,"  time:",tempo
               ier = 4
      endif
      write (nout,'(a)') "  > Restart read successfully at step:",it_start,    &
         "  time:",tempo
      else
         ier = 5
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadRestartFile

