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
! Program unit: ReadBedLoadTransport                   
! Description: Reading input data for bed-load transport.                 
!-------------------------------------------------------------------------------
subroutine ReadBedLoadTransport(ainp,comment,nrighe,ier,ninp,nout,nscr)
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
integer(4) :: nrighe,ier,ninp,nout,nscr,ioerr,ID_erosion_criterion,alloc_stat
integer(4) :: ID_main_fluid,monitoring_lines,i,line_ID
integer(4) :: n_max_iterations,erosion_flag
integer(4) :: deposition_at_frontiers,Gamma_slope_flag,saturation_scheme
double precision :: dt_out,x_fixed,y_fixed,conv_crit_erosion,velocity_fixed_bed
double precision :: x_min_dt,x_max_dt,y_min_dt,y_max_dt,time_minimum_saturation
double precision :: z_min_dt,z_max_dt,t_q0,t_liq,time_maximum_saturation
character(1) :: comment
character(100) :: ainp,lcase
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
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BED LOAD TRANSPORT DATA",ninp,nout)) &
   return
do while (TRIM(lcase(ainp)) /= "##### end bed load transport #####")
! Reading input parameters (first part)
   read(ainp,*,iostat=ioerr) ID_erosion_criterion,ID_main_fluid
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BED-LOAD TRANSPORT INPUT LINE 1", &
      ninp,nout)) return
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   read(ainp,*,iostat=ioerr) saturation_scheme,time_minimum_saturation,        &
      time_maximum_saturation
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BED-LOAD TRANSPORT INPUT LINE 2", &
      ninp,nout)) return
   if (ID_erosion_criterion>0) then
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      read(ainp,*,iostat=ioerr) velocity_fixed_bed,erosion_flag
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                                &
         "VELOCITY FIXED BED, EROSION FLAG",ninp,nout)) return
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)   
      read(ainp,*,iostat=ioerr) deposition_at_frontiers,Gamma_slope_flag
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"OTHER EROSION PARAMETERS",ninp,&
         nout)) return
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)   
      read(ainp,*,iostat=ioerr) monitoring_lines,dt_out,conv_crit_erosion,     &
         n_max_iterations
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                                &
         "BED LOAD TRANSPORT MONITORING LINES",ninp,nout)) return
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)   
      read(ainp,*,iostat=ioerr) x_min_dt,x_max_dt
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"X_MIN_DT,X_MAX_DT",ninp,nout)) &
         return      
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)   
      read(ainp,*,iostat=ioerr) y_min_dt,y_max_dt
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"Y_MIN_DT,Y_MAX_DT",ninp,nout)) &
         return 
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)   
      read(ainp,*,iostat=ioerr) z_min_dt,z_max_dt
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"Z_MIN_DT,Z_MAX_DT",ninp,nout)) &
         return
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      read(ainp,*,iostat=ioerr) t_q0,t_liq
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"t_q0,t_liq",ninp,nout)) return          
   endif
! Writing input parameters (first part)
   if ((ncord>0).and.(nout > 0)) then
      select case (ID_erosion_criterion)
      case (0)
         write(nout,"(1x,a,1p,i12,a)") "ID_erosion_criterion:.........",       &
            ID_erosion_criterion," no-bed load transport."
      case (1)
         write(nout,"(1x,a,1p,i12,a)") "ID_erosion_criterion:.........",       &
            ID_erosion_criterion," Shields-Seminara."
      case (2)
         write(nout,"(1x,a,1p,i12,a)") "ID_erosion_criterion:.........",       &
            ID_erosion_criterion," Shields."
      case (3)
         write(nout,"(1x,a,1p,i12,a)") "ID_erosion_criterion:.........",       &
            ID_erosion_criterion," Mohr-Coulomb."
      endselect
      if (ID_erosion_criterion>0) then      
         write(nout,"(1x,a,1p,i12)") "ID_main_fluid:................",         &
            ID_main_fluid
         write(nout,"(1x,a,1p,i12)") "saturation_scheme:............",         &
            saturation_scheme
         write(nout,"(1x,a,1p,g12.5)") "time_minimum_saturation:......",       &
            time_minimum_saturation
         write(nout,"(1x,a,1p,g12.5)") "time_maximum_saturation:......",       &
            time_maximum_saturation
         write(nout,"(1x,a,1p,g12.5)") "velocity_fixed_bed:...........",       &
            velocity_fixed_bed    
         write(nout,"(1x,a,1p,i12)") "erosion_flag:.................",         &
            erosion_flag
         write(nout,"(1x,a,1p,i12)") "deposition_at_frontiers:......",         &
            deposition_at_frontiers
         write(nout,"(1x,a,1p,i12)") "deposition_at_frontiers:......",         &
            Gamma_slope_flag         
         write(nout,"(1x,a,1p,i12)") "monitoring_lines:.............",         &
            monitoring_lines
         write(nout,"(1x,a,1p,g12.5)") "dt_out:.......................",       &
            dt_out
         write(nout,"(1x,a,1p,g12.5)") "conv_crit_erosion:............",       &
            conv_crit_erosion
         write(nout,"(1x,a,1p,i12)") "n_max_iterations:.............",         &
            n_max_iterations   
         write(nout,"(1x,a,1p,g12.5)") "x_min_dt:.....................",x_min_dt
         write(nout,"(1x,a,1p,g12.5)") "x_max_dt:.....................",x_max_dt
         write(nout,"(1x,a,1p,g12.5)") "y_min_dt:.....................",y_min_dt
         write(nout,"(1x,a,1p,g12.5)") "y_max_dt:.....................",y_max_dt
         write(nout,"(1x,a,1p,g12.5)") "z_min_dt:.....................",z_min_dt
         write(nout,"(1x,a,1p,g12.5)") "z_max_dt:.....................",z_max_dt
         write(nout,"(1x,a,1p,g12.5)") "t_q0:.........................",t_q0
         write(nout,"(1x,a,1p,g12.5)") "t_liq:.........................",t_liq               
         write(nout,"(1x,a)")  " "
      endif
   endif
! Assignment to the body parameters (first part)
   Granular_flows_options%ID_erosion_criterion = ID_erosion_criterion
   if (ID_erosion_criterion>0) then      
      Granular_flows_options%ID_main_fluid = ID_main_fluid
      Granular_flows_options%saturation_scheme = saturation_scheme
      Granular_flows_options%time_minimum_saturation = time_minimum_saturation
      Granular_flows_options%time_maximum_saturation = time_maximum_saturation
      Granular_flows_options%velocity_fixed_bed = velocity_fixed_bed
      Granular_flows_options%erosion_flag = erosion_flag
      Granular_flows_options%deposition_at_frontiers = deposition_at_frontiers
      Granular_flows_options%Gamma_slope_flag = Gamma_slope_flag     
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
                  write(nscr,'(a)')                                            &
                     "Allocation of Granular_flows_options%lines failed. "
                  write(nscr,'(a)') "The execution stops here. "
                  stop
               endif
            endif
! Initializing the auxiliary variable to print results
            Granular_flows_options%it_out_last = 0
      endif
! Loop over the monitoring lines
      do i=1,monitoring_lines
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) line_ID
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "BED LOAD TRANSPORT MONITORING LINES",ninp,nout)) return      
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) x_fixed,y_fixed
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "BED LOAD TRANSPORT MONITORING LINES",ninp,nout)) return
! Assignation of the section parameters 
         Granular_flows_options%lines(i,1) = x_fixed
         Granular_flows_options%lines(i,2) = y_fixed
! Writing on the log file
         if (ncord>0) then
            if (nout>0) then
               write(nout,"(1x,a,i12)") "ID_line:....................",i
               write(nout,"(1x,a,1p,2e12.4)") "x_fixed,y_fixed:............",  &
                  Granular_flows_options%lines(i,1),                           &
                  Granular_flows_options%lines(i,2)
               write(nout,"(1x,a)")  " "
            endif
         endif
      enddo  
   endif   
! Reading the last line 
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BED LOAD TRANSPORT DATA",ninp,nout&
      )) return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadBedLoadTransport

