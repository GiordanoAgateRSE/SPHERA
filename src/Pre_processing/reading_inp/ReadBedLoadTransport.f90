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
! Program unit: ReadBedLoadTransport                   
! Description: Reading input data for bed-load transport.                 
!-------------------------------------------------------------------------------
subroutine ReadBedLoadTransport(ainp,comment,nrighe,ier,ninp,ulog,uerr)
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
integer(4) :: nrighe,ier,ninp,ulog,uerr,ioerr,KTGF_config,alloc_stat
integer(4) :: ID_main_fluid,monitoring_lines,i,line_ID
integer(4) :: n_max_iterations,erosion_flag
integer(4) :: deposition_at_frontiers,saturation_scheme
#ifdef SPACE_3D
integer(4) :: Gamma_slope_flag
#endif
double precision :: dt_out,x_fixed,y_fixed,conv_crit_erosion,velocity_fixed_bed
double precision :: x_min_dt,x_max_dt,y_min_dt,y_max_dt,time_minimum_saturation
double precision :: z_min_dt,z_max_dt,t_q0,t_liq,time_maximum_saturation
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
!------------------------
! Statements
!------------------------
call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"BED LOAD TRANSPORT DATA",ninp,ulog)) &
   return
do while (trim(lcase(ainp)) /= "##### end bed load transport #####")
! Reading input parameters (first part)
   read(ainp,*,iostat=ioerr) KTGF_config,ID_main_fluid
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"BED-LOAD TRANSPORT INPUT LINE 1", &
      ninp,ulog)) return
   call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
   read(ainp,*,iostat=ioerr) saturation_scheme,time_minimum_saturation,        &
      time_maximum_saturation
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"BED-LOAD TRANSPORT INPUT LINE 2", &
      ninp,ulog)) return
   if (KTGF_config>0) then
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
      read(ainp,*,iostat=ioerr) velocity_fixed_bed,erosion_flag
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                                &
         "VELOCITY FIXED BED, EROSION FLAG",ninp,ulog)) return
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
#ifdef SPACE_3D
         read(ainp,*,iostat=ioerr) deposition_at_frontiers,Gamma_slope_flag
#elif defined SPACE_2D
            read(ainp,*,iostat=ioerr) deposition_at_frontiers
#endif
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"OTHER EROSION PARAMETERS",ninp,&
         ulog)) return
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)   
      read(ainp,*,iostat=ioerr) monitoring_lines,dt_out,conv_crit_erosion,     &
         n_max_iterations
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                                &
         "BED LOAD TRANSPORT MONITORING LINES",ninp,ulog)) return
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)   
      read(ainp,*,iostat=ioerr) x_min_dt,x_max_dt
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"X_MIN_DT,X_MAX_DT",ninp,ulog)) &
         return      
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)   
      read(ainp,*,iostat=ioerr) y_min_dt,y_max_dt
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"Y_MIN_DT,Y_MAX_DT",ninp,ulog)) &
         return 
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)   
      read(ainp,*,iostat=ioerr) z_min_dt,z_max_dt
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"Z_MIN_DT,Z_MAX_DT",ninp,ulog)) &
         return
      call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
      read(ainp,*,iostat=ioerr) t_q0,t_liq
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"t_q0,t_liq",ninp,ulog)) return          
   endif
! Writing input parameters (first part)
   if ((input_second_read.eqv..true.).and.(ulog>0)) then
      select case (KTGF_config)
      case (0)
         write(ulog,"(1x,a,1p,i12,a)") "KTGF_config:...................",      &
            KTGF_config," KTGF scheme deactivated."
      case (1)
         write(ulog,"(1x,a,1p,i12,a)") "KTGF_config:...................",      &
            KTGF_config," KTGF scheme with possible 3D erosion criterion."
      case (2)
         write(ulog,"(1x,a,1p,i12,a)") "KTGF_config:...................",      &
            KTGF_config," KTGF scheme with possible 2D erosion criterion."
      endselect
      if (KTGF_config>0) then      
         write(ulog,"(1x,a,1p,i12)") "ID_main_fluid:................",         &
            ID_main_fluid
         write(ulog,"(1x,a,1p,i12)") "saturation_scheme:............",         &
            saturation_scheme
         write(ulog,"(1x,a,1p,g12.5)") "time_minimum_saturation:......",       &
            time_minimum_saturation
         write(ulog,"(1x,a,1p,g12.5)") "time_maximum_saturation:......",       &
            time_maximum_saturation
         write(ulog,"(1x,a,1p,g12.5)") "velocity_fixed_bed:...........",       &
            velocity_fixed_bed    
         write(ulog,"(1x,a,1p,i12)") "erosion_flag:.................",         &
            erosion_flag
         write(ulog,"(1x,a,1p,i12)") "deposition_at_frontiers:......",         &
            deposition_at_frontiers
#ifdef SPACE_3D
         write(ulog,"(1x,a,1p,i12)") "Gamma_slope_flag:.............",         &
            Gamma_slope_flag
#endif         
         write(ulog,"(1x,a,1p,i12)") "monitoring_lines:.............",         &
            monitoring_lines
         write(ulog,"(1x,a,1p,g12.5)") "dt_out:.......................",       &
            dt_out
         write(ulog,"(1x,a,1p,g12.5)") "conv_crit_erosion:............",       &
            conv_crit_erosion
         write(ulog,"(1x,a,1p,i12)") "n_max_iterations:.............",         &
            n_max_iterations   
         write(ulog,"(1x,a,1p,g12.5)") "x_min_dt:.....................",x_min_dt
         write(ulog,"(1x,a,1p,g12.5)") "x_max_dt:.....................",x_max_dt
         write(ulog,"(1x,a,1p,g12.5)") "y_min_dt:.....................",y_min_dt
         write(ulog,"(1x,a,1p,g12.5)") "y_max_dt:.....................",y_max_dt
         write(ulog,"(1x,a,1p,g12.5)") "z_min_dt:.....................",z_min_dt
         write(ulog,"(1x,a,1p,g12.5)") "z_max_dt:.....................",z_max_dt
         write(ulog,"(1x,a,1p,g12.5)") "t_q0:.........................",t_q0
         write(ulog,"(1x,a,1p,g12.5)") "t_liq:........................",t_liq               
         write(ulog,"(1x,a)")  " "
      endif
   endif
! Assignment to the KTGF variables (first part)
   Granular_flows_options%KTGF_config = KTGF_config
   if (KTGF_config>0) then      
      Granular_flows_options%ID_main_fluid = ID_main_fluid
      Granular_flows_options%saturation_scheme = saturation_scheme
      Granular_flows_options%time_minimum_saturation = time_minimum_saturation
      Granular_flows_options%time_maximum_saturation = time_maximum_saturation
      Granular_flows_options%velocity_fixed_bed = velocity_fixed_bed
      Granular_flows_options%erosion_flag = erosion_flag
      Granular_flows_options%deposition_at_frontiers = deposition_at_frontiers
#ifdef SPACE_3D
      Granular_flows_options%Gamma_slope_flag = Gamma_slope_flag     
#endif
      Granular_flows_options%monitoring_lines = monitoring_lines 
      Granular_flows_options%dt_out = dt_out 
      Granular_flows_options%conv_crit_erosion = conv_crit_erosion
      Granular_flows_options%n_max_iterations = n_max_iterations
      Granular_flows_options%x_min_dt = x_min_dt
      Granular_flows_options%x_max_dt = x_max_dt
      Granular_flows_options%y_min_dt = y_min_dt
      Granular_flows_options%y_max_dt = y_max_dt
      Granular_flows_options%z_min_dt = z_min_dt
      Granular_flows_options%z_max_dt = z_max_dt
      Granular_flows_options%t_q0 = t_q0
      Granular_flows_options%t_liq = t_liq
! Allocation of the array of the monitoring lines
      if (allocated(Granular_flows_options%lines)) then
         else
            if (.not.allocated(Granular_flows_options%lines)) then
               allocate(Granular_flows_options%lines(monitoring_lines,2),      &
                  stat=alloc_stat)
               if (alloc_stat/=0) then
                  write(uerr,'(a)')                                            &
                     "Allocation of Granular_flows_options%lines failed. "
                  write(uerr,'(a)') "The execution stops here. "
                  stop
               endif
            endif
! Initializing the auxiliary variable to print results
            Granular_flows_options%it_out_last = 0
      endif
! Loop over the monitoring lines
      do i=1,monitoring_lines
         call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
         read(ainp,*,iostat=ioerr) line_ID
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "BED LOAD TRANSPORT MONITORING LINES",ninp,ulog)) return      
         call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
         read(ainp,*,iostat=ioerr) x_fixed,y_fixed
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "BED LOAD TRANSPORT MONITORING LINES",ninp,ulog)) return
! Assignation of the section parameters 
         Granular_flows_options%lines(i,1) = x_fixed
         Granular_flows_options%lines(i,2) = y_fixed
! Writing on the log file
         if ((input_second_read.eqv..true.).and.(ulog>0)) then
            write(ulog,"(1x,a,i12)") "ID_line:....................",i
            write(ulog,"(1x,a,1p,2e12.4)") "x_fixed,y_fixed:............",     &
               Granular_flows_options%lines(i,1),                              &
               Granular_flows_options%lines(i,2)
            write(ulog,"(1x,a)")  " "
         endif
      enddo  
   endif   
! Reading the last line 
   call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"BED LOAD TRANSPORT DATA",ninp,ulog&
      )) return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadBedLoadTransport
