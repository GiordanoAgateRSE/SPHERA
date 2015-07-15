!AA504 (whole subroutine)
!cfile ReadBedLoadTransport.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Subroutine name : ReadBedLoadTransport
!
! Creation        : 20Aug13 Amicarelli-Agate
!
!************************************************************************************
! Module purpose : Reading input for sediment transport
!
! Calling routines: ReadInput
!
! Called routines: /
!
!************************************************************************************
subroutine ReadBedLoadTransport (ainp,comment,nrighe,ier,ninp,nout,nscr)

! modules
use GLOBAL_MODULE                            
use AdM_USER_TYPE
use ALLOC_Module

!Declarations
implicit none
!AA504 sub
integer(4) :: nrighe,ier,ninp,nout,nscr,ioerr,ID_erosion_criterion,ID_main_fluid,ID_granular,monitoring_lines,i,line_ID,n_max_iterations,erosion_flag,viscosity_blt_formula
!AA504
integer(4) :: deposition_at_frontiers,Gamma_slope_flag
!AA504 sub
double precision :: dt_out,x_fixed,y_fixed,conv_crit_erosion,velocity_fixed_bed,Chezy_friction_coeff,x_min_dt,x_max_dt,y_min_dt,y_max_dt,z_min_dt,z_max_dt
character( 1) :: comment
character(80) :: ainp,lcase 

!External functions
logical, external :: ReadCheck

! in case of restart the cards are not read
if (restart) then
  do while (TRIM(lcase(ainp)) /= "##### end bed load transport #####")
    call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
    if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BED LOAD TRANSPORT DATA",ninp,nout)) return
  end do
  return
end if

call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BED LOAD TRANSPORT DATA",ninp,nout)) return

do while (TRIM(lcase(ainp)) /= "##### end bed load transport #####")
!Reading input parameters (first part)
   read (ainp,*,iostat=ioerr) ID_erosion_criterion,ID_main_fluid,ID_granular
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BED LOAD TRANSPORT GENERAL INPUT",ninp,nout)) return
   if (ID_erosion_criterion>0) then
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
      read (ainp,*,iostat=ioerr) velocity_fixed_bed,erosion_flag
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VELOCITY FIXED BED, EROSION FLAG",ninp,nout)) return
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)   
      read (ainp,*,iostat=ioerr) viscosity_blt_formula,deposition_at_frontiers,Gamma_slope_flag
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VISCOSITY BLT FORMULA",ninp,nout)) return
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)   
      read (ainp,*,iostat=ioerr) monitoring_lines,dt_out,conv_crit_erosion,n_max_iterations
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BED LOAD TRANSPORT MONITORING LINES",ninp,nout)) return
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)   
      read (ainp,*,iostat=ioerr) Chezy_friction_coeff
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CHEZY FRICTION COEFFICIENT",ninp,nout)) return 
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)   
      read (ainp,*,iostat=ioerr) x_min_dt,x_max_dt
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"X_MIN_DT,X_MAX_DT",ninp,nout)) return      
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)   
      read (ainp,*,iostat=ioerr) y_min_dt,y_max_dt
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"Y_MIN_DT,Y_MAX_DT",ninp,nout)) return 
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)   
      read (ainp,*,iostat=ioerr) z_min_dt,z_max_dt
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"Z_MIN_DT,Z_MAX_DT",ninp,nout)) return       
   endif
!
!.. check for the limited modules
!
    if (.not. Erosion_Module_Shields_Mohr .and. ID_erosion_criterion > 1) then
      write (nscr,"(1x,a)") " "
      write (nout,"(1x,a)") " "
      write (nscr,"(1x,a)") " >>WARNING! - The erosion module Shields/Mohr is not available."
      write (nout,"(1x,a)") " >>WARNING! - The erosion module Shields/Mohr is not available."
      ier = 6
      return
    end if
    if (.not. Granular_flux .and. erosion_flag /= 1) then
      write (nscr,"(1x,a)") " "
      write (nout,"(1x,a)") " "
      write (nscr,"(1x,a)") " >>WARNING! - The granular flux module in the simulation are not available."
      write (nout,"(1x,a)") " >>WARNING! - The granular flux module in the simulation are not available."
      ier = 11
      return
    end if
    if (.not. Erosion_Module_Shields_Seminara .and. ID_erosion_criterion == 1) then
      write (nscr,"(1x,a)") " "
      write (nout,"(1x,a)") " "
      write (nscr,"(1x,a)") " >>WARNING! - The erosion module Shields-VanRijn-Seminara in the simulation are not available."
      write (nout,"(1x,a)") " >>WARNING! - The erosion module Shields-VanRijn-Seminara in the simulation are not available."
      ier = 12
      return
    end if
!
! Writing input parameters (first part)
!
   if ((ncord>0).and.(nout > 0)) then
      select case (ID_erosion_criterion)
      case (0)
          write (nout,"(1x,a,1p,i12,a)")   "ID_erosion_criterion:.........",ID_erosion_criterion," no bed load transport."
      case (1)
          write (nout,"(1x,a,1p,i12,a)")   "ID_erosion_criterion:.........",ID_erosion_criterion," Shields-Seminara bed load transport."
      case (2)
          write (nout,"(1x,a,1p,i12,a)")   "ID_erosion_criterion:.........",ID_erosion_criterion," Shields bed load transport."
      case (3)
          write (nout,"(1x,a,1p,i12,a)")   "ID_erosion_criterion:.........",ID_erosion_criterion," Mohr-Coulomb bed load transport."
      end select
      if (ID_erosion_criterion>0) then      
         write (nout,"(1x,a,1p,i12)")   "ID_main_fluid:................",ID_main_fluid
         write (nout,"(1x,a,1p,i12)")   "ID_granular:..................",ID_granular
         write (nout,"(1x,a,1p,g12.5)") "velocity_fixed_bed:...........",velocity_fixed_bed    
         write (nout,"(1x,a,1p,i12)")   "erosion_flag:.................",erosion_flag
         write (nout,"(1x,a,1p,i12)")   "viscosity_blt_formula:........",viscosity_blt_formula
         write (nout,"(1x,a,1p,i12)")   "deposition_at_frontiers:......",deposition_at_frontiers
         write (nout,"(1x,a,1p,i12)")   "deposition_at_frontiers:......",Gamma_slope_flag         
         write (nout,"(1x,a,1p,i12)")   "monitoring_lines:.............",monitoring_lines
         write (nout,"(1x,a,1p,g12.5)") "dt_out:.......................",dt_out
         write (nout,"(1x,a,1p,g12.5)") "conv_crit_erosion:............",conv_crit_erosion
         write (nout,"(1x,a,1p,i12)")   "n_max_iterations:.............",n_max_iterations   
         if (viscosity_blt_formula==2) then
         write (nout,"(1x,a,1p,g12.5)")   "Chezy_friction_coeff:.........",Chezy_friction_coeff
         write (nout,"(1x,a,1p,g12.5)")   "x_min_dt:.....................",x_min_dt
         write (nout,"(1x,a,1p,g12.5)")   "x_max_dt:.....................",x_max_dt
         write (nout,"(1x,a,1p,g12.5)")   "y_min_dt:.....................",y_min_dt
         write (nout,"(1x,a,1p,g12.5)")   "y_max_dt:.....................",y_max_dt
         write (nout,"(1x,a,1p,g12.5)")   "z_min_dt:.....................",z_min_dt
         write (nout,"(1x,a,1p,g12.5)")   "z_max_dt:.....................",z_max_dt         
         endif
         write (nout,"(1x,a)")  " "
      endif
   end if
!Assignation to the body parameters (first part)
   Granular_flows_options%ID_erosion_criterion = ID_erosion_criterion
   if (ID_erosion_criterion>0) then      
      Granular_flows_options%ID_main_fluid = ID_main_fluid
      Granular_flows_options%ID_granular = ID_granular
      Granular_flows_options%velocity_fixed_bed = velocity_fixed_bed
      Granular_flows_options%erosion_flag = erosion_flag
      Granular_flows_options%viscosity_blt_formula = viscosity_blt_formula
      Granular_flows_options%deposition_at_frontiers = deposition_at_frontiers
      Granular_flows_options%Gamma_slope_flag = Gamma_slope_flag     
      Granular_flows_options%monitoring_lines = monitoring_lines 
      Granular_flows_options%dt_out = dt_out 
      Granular_flows_options%conv_crit_erosion = conv_crit_erosion
      if (viscosity_blt_formula==2) then
         Granular_flows_options%Chezy_friction_coeff = Chezy_friction_coeff  
         else
            Granular_flows_options%Chezy_friction_coeff = 0.d0
      endif
      Granular_flows_options%n_max_iterations = n_max_iterations
      Granular_flows_options%x_min_dt = x_min_dt
      Granular_flows_options%x_max_dt = x_max_dt
      Granular_flows_options%y_min_dt = y_min_dt
      Granular_flows_options%y_max_dt = y_max_dt
      Granular_flows_options%z_min_dt = z_min_dt
      Granular_flows_options%z_max_dt = z_max_dt      
! Allocation of the array of the monitoring lines
      if (allocated(Granular_flows_options%lines)) then
         else
            allocate(Granular_flows_options%lines(monitoring_lines,2)) 
!Initializing the auxiliary variable to print results
            Granular_flows_options%it_out_last = 0
      endif
! Loop over the monitoring lines
      do i=1,monitoring_lines
         call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
         read (ainp,*,iostat=ioerr) line_ID
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BED LOAD TRANSPORT MONITORING LINES",ninp,nout)) return      
         call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
         read (ainp,*,iostat=ioerr) x_fixed,y_fixed
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BED LOAD TRANSPORT MONITORING LINES",ninp,nout)) return
!Assignation to the section parameters 
         Granular_flows_options%lines(i,1) = x_fixed
         Granular_flows_options%lines(i,2) = y_fixed
!Writing on the log file
         if ( ncord > 0 ) then
            if ( nout > 0 ) then
               write (nout,"(1x,a,i12)")       "ID_line:....................",i
               write (nout,"(1x,a,1p,2e12.4)") "x_fixed,y_fixed:............",Granular_flows_options%lines(i,1),Granular_flows_options%lines(i,2)
               write (nout,"(1x,a)")  " "
            end if
         endif
      enddo  
   endif   
!Reading last line 
   call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BED LOAD TRANSPORT DATA",ninp,nout)) return
end do
 
return
end subroutine ReadBedLoadTransport
!---split

