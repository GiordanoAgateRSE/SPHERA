!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: ReadBodyDynamics                    
! Description: Reading input data for body transport in fluid flows (Amicarelli
!              et al., 2015, CAF).                  
!-------------------------------------------------------------------------------
#ifdef SOLID_BODIES
subroutine ReadBodyDynamics(ainp,comment,nrighe,ier,ninp,ulog)
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
logical :: CAE
integer(4) :: nrighe,ier,ninp,ulog,ioerr,i,Id_body,n_elem,j,Id_elem,alloc_stat
integer(4) :: imposed_kinematics,n_records,Ic_imposed,surface_detection
double precision :: mass,teta_R_IO
integer(4) :: normal_act(6)
double precision :: L_geom(3),x_CM(3),n_R_IO(3),u_CM(3),omega(3),x_rotC(3)
#ifdef SPACE_3D
double precision :: vec_bp_CAE_trans(3)
#endif
double precision :: mass_deact(6)
double precision :: Ic(3,3)
character(1) :: comment
character(100) :: lcase
character(len=lencard) :: ainp
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
#ifdef SPACE_3D
vec_bp_CAE_trans(1:3) = 0.d0
surface_detection = 0
#endif
!------------------------
! Statements
!------------------------
call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"BODY DYNAMICS DATA",ninp,ulog))      &
   return
do while (trim(lcase(ainp)) /= "##### end body dynamics #####")
! Reading the number of bodies and the ratio between fluid and body particle 
! size
   read(ainp,*,iostat=ioerr) n_bodies,dx_dxbodies,friction_angle,              &
      time_max_no_body_gravity_force,time_max_no_body_frontier_impingements,   &
      body_minimum_pressure_limiter,body_maximum_pressure_limiter,             &
      FSI_slip_conditions,thin_walls
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"BODY DYNAMICS GENERAL INPUT",ninp,&
      ulog)) return
! Writing the number of bodies and "dx_dxbodies" on the log file
! In case of restart, Domain%dx>0 since the first (and only) reading of the main
! input file
   if ((input_second_read.eqv..true.).and.(ulog>0)) then
      write(ulog,"(1x,a,1p,i12)") "n_bodies:.....................",n_bodies
      write(ulog,"(1x,a,1p,e12.4)") "dx_dxbodies:..................",          &
         dx_dxbodies
      write(ulog,"(1x,a,1p,e12.4)") "friction_angle................",          &
         friction_angle
      write(ulog,"(1x,a,1p,e12.4)") "time_max_no_body_gravity_force",          &
         time_max_no_body_gravity_force
      write(ulog,"(1x,a,1p,e12.4)") "time_max_no_body_frontier_imp.",          &
         time_max_no_body_frontier_impingements
      write(ulog,"(1x,a,1p,l12)") "body_minimum_pressure_limiter:",            &
         body_minimum_pressure_limiter
      write(ulog,"(1x,a,1p,l12)") "body_maximum_pressure_limiter:",            &
         body_maximum_pressure_limiter
      write(ulog,"(1x,a,1p,i12)") "FSI_slip_conditions:..........",            &
         FSI_slip_conditions
      write(ulog,"(1x,a,1p,l12)") "thin_walls:...................",            &
         thin_walls
      write(ulog,"(1x,a)")  " "
   endif
! Allocation of the array of the bodies
   if (.not.allocated(body_arr)) then
      allocate(body_arr(n_bodies),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(ulog,*)                                                         &
            'Allocation of body_arr in ReadBodyDynamics failed;',              &
            ' the program terminates here.'
         stop ! Stop the main program
         else
            write(ulog,*)                                                      &
               "Allocation of body_arr in ReadBodyDynamics ",                  &
               "successfully completed."
      endif
   endif
! Loop over the transported bodies
   do i=1,n_bodies
! Reading the body parameters
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
      read(ainp,*,iostat=ioerr) Id_body,n_elem,CAE
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"ID_BODY-N_ELEM-CAE",ninp,      &
         ulog)) return
#ifdef SPACE_3D
      if (CAE) then
         call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
         read(ainp,*,iostat=ioerr) vec_bp_CAE_trans(1:3)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"vec_bp_CAE_trans",ninp,ulog)&
            ) return
         call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
         read(ainp,*,iostat=ioerr) surface_detection
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"surface_detection",ninp,    &
            ulog)) return
      endif
#endif
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
      read(ainp,*,iostat=ioerr) mass
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"MASS",ninp,ulog)) return
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
      read(ainp,*,iostat=ioerr) x_CM(1),x_CM(2),x_CM(3)
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"X_CM",ninp,ulog)) return
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
      read(ainp,*,iostat=ioerr) Ic_imposed
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"IC_IMPOSED",ninp,ulog)) return
      if (Ic_imposed==1) then
         call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
         read(ainp,*,iostat=ioerr) Ic(1,1),Ic(1,2),Ic(1,3)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"IC(1,1-3)",ninp,ulog))      &
            return
         call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
         read(ainp,*,iostat=ioerr) Ic(2,1),Ic(2,2),Ic(2,3)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"IC(2,1-3)",ninp,ulog))      &
            return
         call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
         read(ainp,*,iostat=ioerr) Ic(3,1),Ic(3,2),Ic(3,3)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"IC(3,1-3)",ninp,ulog))      &
            return
      endif
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
      read(ainp,*,iostat=ioerr) n_R_IO(1),n_R_IO(2),n_R_IO(3),teta_R_IO
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                                &
         'BODY ROTATION AXIS AND ANGLE FOR IC',ninp,ulog)) return
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
      read(ainp,*,iostat=ioerr) x_rotC(1),x_rotC(2),x_rotC(3)
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"X_ROTC",ninp,ulog)) return
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
      read(ainp,*,iostat=ioerr) u_CM(1),u_CM(2),u_CM(3)
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"U_CM",ninp,ulog)) return
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
      read(ainp,*,iostat=ioerr) omega(1),omega(2),omega(3)
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"OMEGA",ninp,ulog)) return
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
      read(ainp,*,iostat=ioerr) imposed_kinematics,n_records
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"BODY_KINEMATICS",ninp,ulog))   &
         return
! Assignment of the body parameters 
      body_arr(Id_body)%CAE = CAE
      body_arr(Id_body)%n_elem = n_elem
      body_arr(Id_body)%mass = mass
#ifdef SPACE_3D
      body_arr(Id_body)%vec_bp_CAE_trans(1:3) = vec_bp_CAE_trans(1:3)
      body_arr(Id_body)%surface_detection = surface_detection
#endif
      body_arr(Id_body)%x_CM = x_CM
      body_arr(Id_body)%Ic_imposed = Ic_imposed
      if (Ic_imposed==1) then
         body_arr(Id_body)%Ic = Ic
         else
            body_arr(Id_body)%Ic = 0.d0
      endif   
      body_arr(Id_body)%n_R_IO(:) = n_R_IO(:)
      body_arr(Id_body)%teta_R_IO = teta_R_IO   
      body_arr(Id_body)%x_rotC = x_rotC
      body_arr(Id_body)%u_CM = u_CM
      body_arr(Id_body)%omega = omega
      body_arr(Id_body)%imposed_kinematics = imposed_kinematics
      body_arr(Id_body)%n_records = n_records
! Initialization of other non-read variables
      body_arr(Id_body)%p_max_limiter = 0.d0
! Writing on the log file
      if ((input_second_read.eqv..true.).and.(ulog>0)) then
         write(ulog,"(1x,a,1p,i12)") "body:.......................",Id_body
         write(ulog,"(1x,a,1p,l12)") "CAE:........................",CAE
         write(ulog,"(1x,a,1p,e12.4)") "mass:.......................",mass
#ifdef SPACE_3D
         if (CAE) then
            write(ulog,"(1x,a,1p,3e12.4)") "vec_bp_CAE_trans(1:3):......",     &
               vec_bp_CAE_trans(1:3)
            write(ulog,"(1x,a,1p,i12)") "surface_detection:..........",        &
               surface_detection
         endif
#endif
         write(ulog,"(1x,a,1p,3e12.4)") "x_CM:.......................",x_CM
         write(ulog,"(1x,a,1p,i12)") "IC_imposed:.................",           &
            Ic_imposed
         if (Ic_imposed==1) then
            write(ulog,"(1x,a,1p,3e12.4)") "Ic(1,1-3):..................",     &
               Ic(1,1),Ic(1,2),Ic(1,3)
            write(ulog,"(1x,a,1p,3e12.4)") "Ic(2,1-3):..................",     &
               Ic(2,1),Ic(2,2),Ic(2,3)
            write(ulog,"(1x,a,1p,3e12.4)") "Ic(3,1-3):..................",     &
               Ic(3,1),Ic(3,2),Ic(3,3)
         endif
         write(ulog,"(1x,a,1p,3e12.4)") "n_R_IO:.....................",        &
            n_R_IO(1),n_R_IO(2),n_R_IO(3)
         write(ulog,"(1x,a,1p,e12.4)") "teta_R_IO:..................",         &
            teta_R_IO
         write(ulog,"(1x,a,1p,3e12.4)") "x_rotC:.....................",x_rotC(:)
         write(ulog,"(1x,a,1p,3e12.4)") "u_CM:.......................",u_CM
         write(ulog,"(1x,a,1p,3e12.4)") "omega:......................",omega
         write(ulog,"(1x,a,1p,i12)") "imposed_kinematics:.........",           &
            imposed_kinematics    
         write(ulog,"(1x,a,1p,i12)") "n_records:..................",           &
            n_records
         write(ulog,"(1x,a)")  " "
      endif
! Allocating body elements
      if (input_second_read.eqv..true.) then
         else
            if ((.not.CAE).and.(.not.allocated(body_arr(Id_body)%elem))) then
               allocate(body_arr(Id_body)%elem(body_arr(Id_body)%n_elem),      &
                  STAT=alloc_stat)
               if (alloc_stat/=0) then
                  write(ulog,*) 'Allocation of body_arr(Id_body)%elem in ',    &
                     'ReadBodyDynamics failed; the program terminates here.'
                  stop ! Stop the main program
                  else
                     write(ulog,*) 'Allocation of body_arr(Id_body)%elem in ',&
                        'ReadBodyDynamics successfully completed.'
               endif
            endif
            if ((n_records>0) .and.                                            &
               (.not.allocated(body_arr(Id_body)%body_kinematics))) then
               allocate(body_arr(Id_body)%body_kinematics(n_records,7),        &
                  STAT=alloc_stat)
               if (alloc_stat/=0) then
                  write(ulog,*) 'Allocation of ',                              &
                     'body_arr(Id_body)%body_kinematics in ReadBodyDynamics ', &
                     'failed; the program terminates here.'
                  stop ! Stop the main program
                  else
                     write(ulog,*) 'Allocation of ',                           &
                        'body_arr(Id_body)%body_kinematics in ',               &
                        'ReadBodyDynamics successfully completed.'
               endif
            endif
      endif
! Reading the eventual imposed kinematics
      do j=1,n_records
         call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
         read(ainp,*,iostat=ioerr) body_arr(Id_body)%body_kinematics(j,:) 
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"BODY_KINEMATICS_RECORDS",   &
            ninp,ulog)) return            
      enddo
! Reading element parameters
      do j=1,n_elem
         call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,                    &
            lines_treated=nrighe)
         read(ainp,*,iostat=ioerr) Id_elem
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"ID_ELEM",ninp,ulog))        &
            return
         call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,                    &
            lines_treated=nrighe)
         read(ainp,*,iostat=ioerr) L_geom(1),L_geom(2),L_geom(3)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"L_GEOM",ninp,ulog)) return
         call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,                    &
            lines_treated=nrighe)
         read(ainp,*,iostat=ioerr) x_CM(1),x_CM(2),x_CM(3)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"X_CM_ELEM",ninp,ulog))      &
            return
         call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,                    &
            lines_treated=nrighe)
         read(ainp,*,iostat=ioerr) n_R_IO(1),n_R_IO(2),n_R_IO(3),teta_R_IO
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "ELEMENT ROTATION AXIS AND ANGLE FOR IC",ninp,ulog)) return         
         call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,                    &
            lines_treated=nrighe)
         read(ainp,*,iostat=ioerr) normal_act(1),normal_act(2),                &
            normal_act(3),normal_act(4),normal_act(5),normal_act(6)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"NORMAL_ACT",ninp,ulog))     &
            return
         call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,                    &
            lines_treated=nrighe)
         read(ainp,*,iostat=ioerr) mass_deact(1),mass_deact(2),                &
            mass_deact(3),mass_deact(4),mass_deact(5),mass_deact(6)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"MASS_DEACT",ninp,ulog))     &
            return
! Assignment of the element parameters
         body_arr(Id_body)%elem(Id_elem)%L_geom = L_geom
         body_arr(Id_body)%elem(Id_elem)%x_CM = x_CM
         body_arr(Id_body)%elem(Id_elem)%n_R_IO(:) = n_R_IO(:)
         body_arr(Id_body)%elem(Id_elem)%teta_R_IO = teta_R_IO
         body_arr(Id_body)%elem(Id_elem)%normal_act = normal_act
         body_arr(Id_body)%elem(Id_elem)%mass_deact = mass_deact
! Writing on the log file
         if ((input_second_read.eqv..true.).and.(ulog>0)) then
            write(ulog,"(1x,a,1p,i12)") "element:....................",        &
               Id_elem
            write(ulog,"(1x,a,1p,3e12.4)") "L_geom:.....................",     &
               L_geom
            write(ulog,"(1x,a,1p,3e12.4)") "x_CM_elem:..................",     &
               x_CM(:)
            write(ulog,"(1x,a,1p,3e12.4)") "n_R_IO_elem:................",     &
               n_R_IO(1),n_R_IO(2),n_R_IO(3)
            write(ulog,"(1x,a,1p,e12.4)") "teta_R_IO_elem:.............",      &
               teta_R_IO
            write(ulog,"(1x,a,1p,6i12)") "normal_act:.................",       &
               normal_act    
            write(ulog,"(1x,a,1p,6e12.4)") "mass_deact:.................",     &
               mass_deact
            write(ulog,"(1x,a)")  " "
         endif
      enddo
   enddo
   call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"BODY DYNAMICS DATA",ninp,ulog))   &
      return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadBodyDynamics
#endif
