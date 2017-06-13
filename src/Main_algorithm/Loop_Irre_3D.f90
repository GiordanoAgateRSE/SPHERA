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
! Program unit: Loop_Irre_3D         
! Description: 3D main algorithm.                    
!-------------------------------------------------------------------------------
subroutine Loop_Irre_3D
!------------------------
! Modules
!------------------------ 
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
logical :: done_flag
integer(4) :: Ncbf,i,npi,it,it_print,it_memo,it_rest,ir,ii,num_out,Ncbf_Max
integer(4) :: OpCountot,SpCountot,EpCountot,EpOrdGridtot,ncel,aux,igridi
integer(4) :: jgridi,kgridi,machine_Julian_day,machine_hour,machine_minute
integer(4) :: machine_second,alloc_stat
real :: time_aux_2
double precision :: pretot,BCtorodivV,dt_previous_step,TetaV1,xmax,ymax,zmax
double precision :: appo1,appo2,appo3,dtvel
real :: time_aux(2)
double precision,dimension(1:SPACEDIM) :: tpres,tdiss,tvisc,BoundReaction
character(len=lencard) :: nomsub = "Loop_Irre_3D"
integer(4),external :: ParticleCellNumber,CellIndices,CellNumber
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
if (esplosione) then
   do npi=1,nag
      if (index(Med(pg(npi)%imed)%tipo,"gas")>0) then
         pg(npi)%pres = (Med(pg(npi)%imed)%gamma - one) * pg(npi)%IntEn *      &
                        pg(npi)%dens
         pg(npi)%Csound = Dsqrt(Med(pg(npi)%imed)%gamma *                      &
                          (Med(pg(npi)%imed)%gamma - one) * pg(npi)%IntEn)
         else
            pg(npi)%Csound = Med(pg(npi)%imed)%Celerita
            pg(npi)%IntEn = pg(npi)%pres / ((Med(pg(npi)%imed)%gamma - one) *  &
                            pg(npi)%dens)
      endif
      pg(npi)%state = "flu"
   enddo
endif
SpCount = 0
OpCount = 0
EpCount = 0
EpOrdGrid = 0
num_out = 0
SourceFace = 0
it_eff = it_start
it_print = it_start
it_memo = it_start
it_rest = it_start
! Variable to count the particles, which are not "sol"
indarrayFlu = 0
do npi=1,nag
   if ((pg(npi)%cella==0).or.(pg(npi)%vel_type/="std")) cycle
   if (pg(npi)%state=="flu") then
      indarrayFlu = indarrayFlu + 1
! To check the maximum dimension of the array and possible resizing
      if (indarrayFlu>PARTICLEBUFFER) then
         call diagnostic (arg1=9,arg2=1,arg3=nomsub)
      endif
      Array_Flu(indarrayFlu) = npi
   endif
enddo
! Introductory procedure for inlet conditions
call PreSourceParticles_3D
! Initializing the time stage for time integration
if (Domain%time_split==0) Domain%time_stage = 1
!------------------------
! Statements
!------------------------
! SPH parameters 
call start_and_stop(2,10)
if ((on_going_time_step==it_start).and.(Domain%tipo=="bsph"))                  &
   on_going_time_step = -2
call CalcVarLength
if ((on_going_time_step==-2).and.(Domain%tipo=="bsph")) on_going_time_step =   &
   it_start
call start_and_stop(3,10)
! Removing fictitious reservoirs used for DB-SPH initialization
if ((on_going_time_step==it_start).and.(Domain%tipo=="bsph")) then
   call start_and_stop(2,9) 
!$omp parallel do default(none) shared(nag,pg,OpCount,Partz) private(npi)
   do npi=1,nag
! Fictitious air reservoirs
      if (Partz(pg(npi)%izona)%DBSPH_fictitious_reservoir_flag.eqv.(.true.))  &
         then   
         OpCount(pg(npi)%imed) = OpCount(pg(npi)%imed) + 1    
         pg(npi)%cella = -1
      endif
! Fictitious fluid reservoir top
      if (Partz(pg(npi)%izona)%IC_source_type==2) then
         if (pg(npi)%coord(3)>Partz(pg(npi)%izona)%H_res) then
            OpCount(pg(npi)%imed) = OpCount(pg(npi)%imed) + 1
            pg(npi)%cella = -1
         endif
      endif
   enddo
!$omp end parallel do
   call OrdGrid1 (nout)
! Variable to count the particles, which are not "sol"
   indarrayFlu = 0
   do npi=1,nag
      if ((pg(npi)%cella==0).or.(pg(npi)%vel_type/="std")) cycle
      if (pg(npi)%state=="flu") then
         indarrayFlu = indarrayFlu + 1
! To check the maximum dimension of the array and possible resizing
         if (indarrayFlu>PARTICLEBUFFER) then
            call diagnostic(arg1=9,arg2=1,arg3=nomsub)
         endif
         Array_Flu(indarrayFlu) = npi
      endif
   enddo
   call start_and_stop(3,9)
   call start_and_stop(2,10)
! This fictitious value avoid CalcVarlength computing Gamma, sigma and density 
! for this very call     
   on_going_time_step = - 1
   call CalcVarLength
! The correct value of "on_going_time_step" is restored
   on_going_time_step = it_start
   call start_and_stop(3,10)
endif
! Subroutine for wall BC treatment (DB-SPH)
! Density and pressure updates for wall elements: MUSCL + LPRS + 
! equation of state 
call start_and_stop(2,18)
if ((Domain%tipo=="bsph").and.(nag>0).and.(DBSPH%n_w>0)) then
   call Gradients_to_MUSCL   
   call BC_wall_elements
endif
call start_and_stop(3,18)
! Pressure initialization for body particles
if (n_bodies>0) then
   call start_and_stop(2,19)
   call body_pressure_mirror
   call body_pressure_postpro
   call start_and_stop(3,19)
endif
! To evaluate the close boundaries and integrals for the current particle in 
! every loop and storing them in the general storage array. Computation and 
! storage of the intersections between the kernel support and the frontier and 
! the corresponding boundary integrals (SA-SPH).
if (Domain%tipo=="semi") then
   call start_and_stop(2,11)
   call ComputeBoundaryDataTab
   call start_and_stop(3,11)
endif
! To evaluate the properties that must be attributed to the fixed particles
if (Domain%NormFix) call NormFix
if (nscr>0) write (nscr,"(a,1x,a)") " Running case:",trim(nomecas2)
if (nout>0) then
   it_print = it_eff
   call Print_Results(it_eff,it_print,'inizio')
endif
if (nres>0) then
   it_memo = it_eff
   it_rest = it_eff
   call Memo_Results(it_eff,it_memo,it_rest,dtvel,'inizio')
endif
if (vtkconv) then
   call result_converter('inizio')
endif
! To assess the initial time step
if (it_start==0) call inidt2
it = it_start
if (exetype=="linux") then
   if (Domain%tmax>0.d0) then
      call system("date +%j%H%M%S>date_pre_iterations.txt")
      open(unit_time_elapsed,file='date_pre_iterations.txt',status="unknown",  &
         form="formatted")
      read(unit_time_elapsed,'(i3,i2,i2,i2)') machine_Julian_day,machine_hour, &
         machine_minute,machine_second
      close(unit_time_elapsed)
      Domain%t_pre_iter = machine_Julian_day * 24 * 60 * 60 + machine_hour * 60&
                          * 60 + machine_minute * 60 + machine_second
   endif
endif
ITERATION_LOOP: do while (it<=Domain%itmax)
   done_flag = .false.
   it = it + 1
   on_going_time_step = it
   it_eff = it
! To store the old time step duration, to evaluate the new value of time step 
! duration and the total time value
   dt_previous_step = dt
! Stability criteria
   if (nag>0) call rundt2
   simulation_time = simulation_time + dt
   if (nscr>0) write (nscr,"(a,i8,a,g13.6,a,g12.4,a,i8,a,i5)") " it= ",it,     &
      "   time= ",simulation_time,"  dt= ",dt,"    npart= ",nag,"   out= ",    &
      num_out
   do while ((done_flag.eqv.(.false.)).or.((Domain%RKscheme>1).and.            &
      (Domain%time_stage>1)))
      if (Domain%time_split==0) then
! Explicit RK schemes 
! SPH parameters 
         call start_and_stop(2,10)
         if ((Domain%tipo=="semi").and.(Domain%time_stage<2)) call CalcVarLength
         call start_and_stop(3,10)
! To evaluate the close boundaries and integrals
! for the current particle in every loop and storing them in the general 
! storage array.
! Computation and storage of the boundary integrals
         if (Domain%tipo=="semi") then
            call start_and_stop(2,11)
            call ComputeBoundaryDataTab
            call start_and_stop(3,11)
            endif
      endif
! Possible imposed velocities
! Stop particles belonging to a rigid block and moving with a fixed 
! velocity.
! One may loop over particles instead of over zones
      if (Domain%tipo=="semi") then
         do ir = 1,NPartZone 
            if (Partz(ir)%move/="fix") cycle
! It assigns the movement with an imposed kinematics ("npointv" velocity data)
            if (Partz(ir)%npointv>1) then
               call vellaw(Partz(ir)%vlaw,Partz(ir)%vel,Partz(ir)%npointv)
               write(nout,"(f12.4,a,i2,a,3e12.4)") simulation_time,"  zona",ir,&
                  "  vel.",Partz(ir)%vel
! As an alternative, one may loop over particles 
               do npi=Partz(ir)%limit(1),Partz(ir)%limit(2) 
                  if (pg(npi)%cella==0) cycle
                  pg(npi)%var(:) = Partz(ir)%vel(:)
                  if (simulation_time>=pg(npi)%tstop) then
                     pg(npi)%vstart(:) = zero
                     pg(npi)%vel(:) = zero
                     else
                        pg(npi)%vstart(:) = Partz(ir)%vel(:)
                        pg(npi)%vel(:) = Partz(ir)%vel(:)
                  endif
               enddo
! Fixed value of velocity (npointv = 1)
               elseif (Partz(ir)%npointv==1) then
! As an alternative, one may loop over particles 
                  do npi=Partz(ir)%limit(1),Partz(ir)%limit(2)
                     if (pg(npi)%cella==0) cycle
                     pg(npi)%var(:) = Partz(ir)%vel(:)
                     if (simulation_time>=pg(npi)%tstop) then
                        pg(npi)%vstart(:) = zero
                        pg(npi)%vel(:)    = zero
                        else
                           pg(npi)%vel(:) = Partz(ir)%vel(:)
                     endif
                  enddo
            endif
         enddo
      endif
      if ((Domain%time_split==0).and.(Domain%time_stage==1)) then               
! Erosion criterium + continuity equation RHS  
         call start_and_stop(2,12)
         if ((Granular_flows_options%ID_erosion_criterion>0).and.              &
            (.not.esplosione)) then
            select case (Granular_flows_options%ID_erosion_criterion)
               case(1)
!$omp parallel do default(none) shared(pg,nag) private(npi,ncel)
                  do npi=1,nag
                     pg(npi)%vel_old(:) = pg(npi)%vel(:)
                     pg(npi)%normal_int_old(:) = pg(npi)%normal_int(:)
                     call initialization_fixed_granular_particle(npi)             
                  enddo
!$omp end parallel do 
!$omp parallel do default(none) shared(pg,nag) private(npi) 
                  do npi=1,nag
                     call Shields(npi) 
                  enddo
!$omp end parallel do
! Initializing viscosity for fixed particles
!$omp parallel do default(none) shared(pg,nag,Granular_flows_options,Med)      &
!$omp private(npi,ncel,aux,igridi,jgridi,kgridi)
                  do npi=1,nag
                     ncel = ParticleCellNumber(pg(npi)%coord)
                     aux = CellIndices(ncel,igridi,jgridi,kgridi)
                     if (pg(npi)%state=="sol") then
                        pg(npi)%mu = Med(pg(npi)%imed)%mumx
                        pg(npi)%visc = pg(npi)%mu / pg(npi)%dens
                     endif
                  enddo
!$omp end parallel do
               case(2)
!$omp parallel do default(none) shared(pg,nag) private(npi)
                  do npi=1,nag
                     call Shields(npi) 
                  enddo
!$omp end parallel do 
               case(3)
! To compute the second invariant of the rate-strain tensor and density 
! derivatives
                  call inter_EqCont_3D 
                  call MohrC
               case default
            endselect
! Update auxiliary vector for counting particles, whose status is not "sol"
            indarrayFlu = 0
            do npi=1,nag
               if ((pg(npi)%cella==0).or.(pg(npi)%vel_type/="std")) cycle
               if (pg(npi)%state=="flu") then
                  indarrayFlu = indarrayFlu + 1
! Check the boundary sizes and possible resizing
                  if (indarrayFlu>PARTICLEBUFFER) then
                     call diagnostic(arg1=9,arg2=2,arg3=nomsub)
                  endif
                  Array_Flu(indarrayFlu) = npi
               endif
            enddo
            else
! No erosion criterion
               indarrayFlu = 0
               do npi=1,nag
                  if ((pg(npi)%cella==0).or.(pg(npi)%vel_type/="std")) cycle
                  indarrayFlu = indarrayFlu + 1
! Checking not to overpass array sizes. Possible resizing.
                  if (indarrayFlu>PARTICLEBUFFER) then
                     call diagnostic(arg1=9,arg2=2,arg3=nomsub)
                  endif
                  Array_Flu(indarrayFlu) = npi
               enddo
         endif
      endif
! Momentum equation 
      call start_and_stop(2,6)
      Ncbf_Max = 0
!$omp parallel do default(none)                                                &
!$omp private(npi,ii,tpres,tdiss,tvisc,ncbf,boundreaction)                     &
!$omp shared(nag,pg,Domain,BoundaryDataPointer,Ncbf_Max,indarrayFlu,Array_Flu) &
!$omp shared(it,Med,Granular_flows_options)
! Loop over particles
      do ii = 1,indarrayFlu
         npi = Array_Flu(ii)
! The mixture particles in the elastic-plastic strain regime are held fixed.
         if (pg(npi)%mu==Med(pg(npi)%imed)%mumx) then
            pg(npi)%acc(:) = zero
            cycle
         endif
         call inter_EqMoto(npi,tpres,tdiss,tvisc)
! Searching for the boundary faces, which are the nearest the npi-th current 
! particle
         if ((Domain%time_stage==1).or.(Domain%time_split==1)) then 
            pg(npi)%kodvel = 0
            pg(npi)%velass = zero
         endif
         if (Domain%tipo=="semi") then
            Ncbf = BoundaryDataPointer(1,npi)
            else
               Ncbf = 0
         endif
         if ((Domain%tipo=="semi").and.(Ncbf>0)) then
               Ncbf_Max = max(Ncbf_Max, Ncbf)
               call AddBoundaryContributions_to_ME3D(npi,Ncbf,tpres,tdiss,tvisc)
               if (pg(npi)%kodvel==0) then
                  BoundReaction = zero
                  call AddElasticBoundaryReaction_3D(npi,Ncbf,BoundReaction)
                  pg(npi)%acc(:) = tpres(:) + tdiss(:) + tvisc(:) +            &
                                   Domain%grav(:) + BoundReaction(:)
                  else
                     pg(npi)%acc(:) = zero
               endif
               else
                  if (Domain%tipo=="semi") then
                     pg(npi)%acc(:) = tpres(:) + tdiss(:) + tvisc(:) +         &
                                      Domain%grav(:)
                     else
                        if (Domain%tipo=="bsph") then
                           pg(npi)%acc(:) = (tpres(:) + tdiss(:) + tvisc(:)) / &
                                            pg(npi)%Gamma + Domain%grav(:)
                        endif         
                  endif
         endif
      enddo
!$omp end parallel do
      if (Ncbf_Max>Domain%MAXCLOSEBOUNDFACES) then
         write (nout,"(a,i5,a,i5)")                                            &
            "Increase parameter MAXCLOSEBOUNDFACES from",                      &
            Domain%MAXCLOSEBOUNDFACES," to ",Ncbf_Max
         call diagnostic(arg1=9,arg2=3,arg3=nomsub)
      endif
! Time integration for body dynamics
      if (n_bodies>0) then
         call start_and_stop(3,6)
         call start_and_stop(2,19)
         call RHS_body_dynamics(dtvel)
         call start_and_stop(3,19)
         call start_and_stop(2,6)
      endif
! Time integration scheme for momentum and energy equations 
      if (Domain%time_split==0) then   
! Explicit RK schemes
         call start_and_stop(3,6)
! Velocity smoothing, trajectory equation, BC, neighboring parameters (start)
         elseif (Domain%time_split==1) then   
! dt computation 
            dtvel = half * (dt + dt_previous_step) 
!$omp parallel do default(none)                                                &
!$omp private(npi,ii)                                                          &
!$omp shared(nag,Pg,dtvel,indarrayFlu,Array_Flu)
            do ii=1,indarrayFlu
               npi = Array_Flu(ii)
! kodvel = 0: the particle is internal to the domain. 
               if (pg(npi)%kodvel==0) then 
                  pg(npi)%vel(:) = pg(npi)%vel(:) + dtvel * pg(npi)%acc(:)
! kodvel = 1: the particle has a critical flux condition. The vertical 
! velocity component is assigned.
                  elseif (pg(npi)%kodvel==1) then                                                    
                     pg(npi)%vel(:) = pg(npi)%vel(:) + dtvel * pg(npi)%acc(:)                            
                     pg(npi)%vel(3) = pg(npi)%velass(3)
! kodvel = 2: the particle has an assigned normal velocity or source 
! condition. All the velocity components are assigned.
                     elseif (pg(npi)%kodvel==2) then                                                    
                        pg(npi)%vel(:) = pg(npi)%velass(:)                                                
               endif 
            enddo
!$omp end parallel do
            call start_and_stop(3,6)
! Energy equation: start
            if (esplosione) then
               do ii = 1,indarrayFlu
                 npi = Array_Flu(ii)
                 pg(npi)%IntEn = pg(npi)%IntEn + dtvel * pg(npi)%dEdT
               enddo
            endif
! Energy equation: end
! Time integration for body dynamics
            if (n_bodies>0) then
               call start_and_stop(2,19)
               call time_integration_body_dynamics(dtvel)
               call start_and_stop(3,19)
            endif
! Partial smoothing for velocity: start 
            call start_and_stop(2,7)
            call inter_SmoothVelo_3D
!$omp parallel do default(none) private(npi,ii,TetaV1)                         &
!$omp shared(nag,Pg,Med,Domain,dt,indarrayFlu,Array_Flu,esplosione)
! Loop over all the active particles
            do ii=1,indarrayFlu
               npi = Array_Flu(ii)
               if (esplosione) then
                  TetaV1 = Domain%TetaV * pg(npi)%Csound * dt / Domain%h
                  else
! TetaV depending on the time step
                     TetaV1 = Domain%TetaV * Med(pg(npi)%imed)%Celerita * dt / &
                              Domain%h
               endif
               if (esplosione) pg(npi)%IntEn = pg(npi)%IntEn + TetaV1 *        &
                                               pg(npi)%Envar
               if (pg(npi)%kodvel==0) then                                            
! The particle is inside the domain and far from boundaries
                  pg(npi)%var(:) = pg(npi)%vel(:) + TetaV1 * pg(npi)%var(:)     
                  pg(npi)%vel(:) = pg(npi)%var(:)                                        
                  else
! The particle is close to a "source", "level" or "normal velocity boundary
! (kodvel = 1 or = 2): the final velocity is kept unmodified
                     pg(npi)%var(:) = pg(npi)%vel(:)                                         
               endif                                                                 
            enddo
!$omp end parallel do
            call start_and_stop(3,7)
! Partial smoothing for velocity: end
! Diffusion coefficient: start (input option not recommended)
            if (diffusione) then
               call start_and_stop(2,15)
!$omp parallel do default(none) private(npi,ii,appo1,appo2,appo3)              &
!$omp shared(nag,pg,Med,indarrayFlu,Array_Flu)
               do ii=1,indarrayFlu
                  npi = Array_Flu(ii)
                  if ((pg(npi)%VolFra==VFmx).and.                              &
                     (pg(npi)%visc==Med(pg(npi)%imed)%mumx/pg(npi)%dens)) then
                     pg(npi)%coefdif = zero
                     else
                        call inter_CoefDif(npi)
                        if (pg(npi)%uni>zero) pg(npi)%veldif = pg(npi)%veldif  &
                                                               / pg(npi)%uni 
                        appo1 = (pg(npi)%veldif(1) - pg(npi)%var(1)) *         &
                                (pg(npi)%veldif(1) - pg(npi)%var(1))
                        appo2 = (pg(npi)%veldif(2) - pg(npi)%var(2)) *         &
                                (pg(npi)%veldif(2) - pg(npi)%var(2))
                        appo3 = (pg(npi)%veldif(3) - pg(npi)%var(3)) *         &
                                (pg(npi)%veldif(3) - pg(npi)%var(3))
                        pg(npi)%coefdif = pg(npi)%coefdif * Dsqrt(appo1 + appo2&
                                          + appo3)
                  endif
               enddo
!$omp end parallel do
               call start_and_stop(3,15)
            endif
! Diffusion coefficient: end
! Update the particle positions
            call start_and_stop(2,8)
!$omp parallel do default(none) private(npi) shared(nag,pg,dt)
! Loop over the active particles
            do npi=1,nag
               if (pg(npi)%cella==0) cycle
! To save the old coordinates
               pg(npi)%CoordOld(:) = pg(npi)%coord(:)
               if (pg(npi)%vel_type/="std") then
! If the motion type is not "std", velocities are assigned
                  pg(npi)%coord(:) = pg(npi)%coord(:) + dt * pg(npi)%vstart(:)
                  else
! Otherwise, the partial smoothed velocity field is integrated in time
                     pg(npi)%coord(:) = pg(npi)%coord(:) + dt * pg(npi)%vel(:)  
               endif
            enddo
!$omp end parallel do
            call start_and_stop(3,8)
! Check on the particles gone out of the domain throughout the opened faces
            call start_and_stop(2,9)
            if (NumOpenFaces>0) call CancelOutgoneParticles_3D
! Adding new particles from the inlet sections
            if (SourceFace/=0) call GenerateSourceParticles_3D
! Particle reordering
            call OrdGrid1(nout)
            call start_and_stop(3,9)
! Set the parameters for the fixed particles 
            if (Domain%NormFix) call NormFix
! SPH parameters
            call start_and_stop(2,10)
            call CalcVarLength
            call start_and_stop(3,10)
! Assessing the close boundaries and the integrals
! for the current particle in every loop and storing them in the general 
! storage array.
! Computation and storage of the boundary integrals
            if (Domain%tipo=="semi") then
               call start_and_stop(2,11)
               call ComputeBoundaryDataTab
               call start_and_stop(3,11)
            endif
      endif
! Continuity equation 
! Erosion criterion + continuity equation RHS 
      call start_and_stop(2,12)
      Ncbf_Max = 0
      if ((Granular_flows_options%ID_erosion_criterion>0).and.                 &
         (.not.esplosione)) then  
         if (Domain%time_split==1) then 
! Assessing particle status ("flu" or "sol") of the mixture particles
! Calling the proper subroutine for the erosion criterion 
            select case (Granular_flows_options%ID_erosion_criterion)
               case(1)
!$omp parallel do default(none) shared(pg,nag) private(npi,ncel)
                  do npi=1,nag
                     pg(npi)%vel_old(:) = pg(npi)%vel(:)
                     pg(npi)%normal_int_old(:) = pg(npi)%normal_int(:)        
                     call initialization_fixed_granular_particle(npi)             
                  enddo
!$omp end parallel do
!$omp parallel do default(none) shared(pg,nag) private(npi)
                  do npi=1,nag
                     call Shields(npi) 
                  enddo
!$omp end parallel do 
!$omp parallel do default(none) shared(pg,nag,Granular_flows_options,Med)      &
!$omp private(npi,ncel,aux,igridi,jgridi,kgridi)
                  do npi=1,nag
                     ncel = ParticleCellNumber(pg(npi)%coord)
                     aux = CellIndices(ncel,igridi,jgridi,kgridi)
                     if (pg(npi)%state=="sol") then
                        pg(npi)%mu = Med(pg(npi)%imed)%mumx
                        pg(npi)%visc = pg(npi)%mu / pg(npi)%dens
                     endif
                  enddo
!$omp end parallel do
               case(2)
!$omp parallel do default(none) shared(pg,nag) private(npi)
                  do npi=1,nag
                     call Shields(npi) 
                  enddo
!$omp end parallel do 
               case(3)
! To compute the second invariant of the rate-strain tensor and density 
! derivatives
                  call inter_EqCont_2D 
                  call MohrC
               case default
            endselect
! Update auxiliary vector for counting particles, whose status is not "sol"
            indarrayFlu = 0
            do npi=1,nag
               if (pg(npi)%cella==0.or.pg(npi)%vel_type/="std") cycle
               if (pg(npi)%state=="flu") then
                  indarrayFlu = indarrayFlu + 1
! Check the boundary sizes and possible resizing
                  if (indarrayFlu>PARTICLEBUFFER) then
                     call diagnostic(arg1=9,arg2=2,arg3=nomsub)
                  endif
                  Array_Flu(indarrayFlu) = npi
               endif
            enddo
         endif
! To compute the second invariant of the strain-rate tensor and density 
! derivatives
         call inter_EqCont_3D
         else 
! No erosion criterion  
! To compute the second invariant of the strain-rate tensor and density 
! derivatives
            call inter_EqCont_3D
            if (Domain%time_split==1) then 
               indarrayFlu = 0
               do npi=1,nag
                  if (pg(npi)%cella==0.or.pg(npi)%vel_type/="std") cycle
                  indarrayFlu = indarrayFlu + 1
! Check array sizes and possible resizing 
                  if (indarrayFlu>PARTICLEBUFFER) then
                     call diagnostic(arg1=9,arg2=2,arg3=nomsub)
                  endif
                  Array_Flu(indarrayFlu) = npi
               enddo
            endif
      endif 
      if ((Domain%time_stage==1).or.(Domain%time_split==1)) then 
!$omp parallel do default(none) private(npi) shared(nag,pg)
         do npi=1,nag
            pg(npi)%koddens = 0
            pg(npi)%densass = pg(npi)%dens
         enddo
!$omp end parallel do
      endif
      if (Domain%tipo=="semi") then
!$omp parallel do default(none) private(npi,ii,BCtorodivV,Ncbf)                &
!$omp shared(nag,pg,BoundaryDataPointer,Ncbf_Max,indarrayFlu,Array_Flu,it)
! Loop over all the active particles
! Density update 
! Boundary contributions to the continuity equation (SA-SPH)
         do ii=1,indarrayFlu
            npi = Array_Flu(ii)
! Seaching for the neighbouring faces of the particle "npi"
            BCtorodivV = zero
            Ncbf = BoundaryDataPointer(1,npi)
! Detecting the faces with actual contributions
            if (Ncbf>0) then        
               Ncbf_Max = max(Ncbf_Max, Ncbf)
               Call AddBoundaryContribution_to_CE3D (npi, Ncbf, BCtorodivV)
            endif
            if (pg(npi)%koddens==0) then
               pg(npi)%dden = pg(npi)%dden - BCtorodivV
               elseif (pg(npi)%koddens==1) then
                  pg(npi)%dden = zero  
! Boundary type is velocity or source
                   elseif (pg(npi)%koddens==2) then
                      pg(npi)%dden = zero
             endif
          enddo
!$omp end parallel do
      endif  
      if (Ncbf_Max>Domain%MAXCLOSEBOUNDFACES) then
         write (nout,"(a,i5,a,i5)")                                            &
            "Increase parameter MAXCLOSEBOUNDFACES from "                      &
            ,Domain%MAXCLOSEBOUNDFACES," to ",Ncbf_Max
         call diagnostic(arg1=9,arg2=4,arg3=nomsub)
      endif
! Loop over all the active particles
! Time integration of the continuity equation 
      if (Domain%time_split==0) then   
! Explicit RK schemes
         call start_and_stop(3,12)
         elseif (Domain%time_split==1) then 
! Leapfrog scheme
!$omp parallel do default(none) private(npi,ii)                                &
!$omp shared(nag,pg,dt,indarrayFlu,Array_Flu,Domain,med)
            do ii=1,indarrayFlu
               npi = Array_Flu(ii)
               if (pg(npi)%cella==0.or.pg(npi)%vel_type/="std") cycle
               if (Domain%tipo=="bsph") pg(npi)%dden=pg(npi)%dden/pg(npi)%uni
! Boundary type is "fixe" or "tapis" or "level"
               if (pg(npi)%koddens==0) then
                  if (Domain%tipo=="semi") pg(npi)%dens = pg(npi)%dens + dt *  &
                     pg(npi)%dden
                  if (Domain%density_thresholds==1) then        
                     if (pg(npi)%dens<(0.9d0*med(pg(npi)%imed)%den0))          &
                        pg(npi)%dens = 0.9d0*med(pg(npi)%imed)%den0
                     if (pg(npi)%dens>(1.1d0*med(pg(npi)%imed)%den0))          &
                        pg(npi)%dens = 1.1d0 * med(pg(npi)%imed)%den0
                  endif
                  pg(npi)%densass = zero
! Boundary type is "velocity" or "source"
                  elseif (pg(npi)%koddens==2) then
                     pg(npi)%dens = pg(npi)%densass  
! Density is kept constant
               endif
            enddo
!$omp end parallel do
            call start_and_stop(3,12)
! Continuity equation: end
            if (diffusione) then
               call start_and_stop(2,16)
               call aggdens
               call start_and_stop(3,16)
            endif
! Equation of state 
            call start_and_stop(2,13)
            call calcpre 
            call start_and_stop(3,13)
            if (n_bodies>0) then
               call start_and_stop(2,19)
               call body_pressure_mirror
               call start_and_stop(3,19)
            endif
      endif
      if (Domain%time_split==0) call time_integration
! Explicit RK schemes
! Partial smoothing for pressure and density update
      if (Domain%TetaP>zero) then
         call start_and_stop(2,14)
         if (Domain%Psurf=='s') then
            call inter_SmoothPres
            elseif (Domain%Psurf=='a') then
               call PressureSmoothing_3D
         endif
         call start_and_stop(3,14)
      endif
      if (n_bodies>0) then
         call start_and_stop(2,19)
         call body_pressure_postpro
         call start_and_stop(3,19)
      endif
      if (diffusione) then
         call start_and_stop(2,16)
!$omp parallel do default(none) private(npi,ii)                                &
!$omp shared(nag,Pg,Med,indarrayFlu,Array_Flu)
         do ii = 1,indarrayFlu
            npi = Array_Flu(ii)
            if (pg(npi)%koddens/=0) cycle
            if (pg(npi)%imed==1) then
               pg(npi)%dens = pg(npi)%pres / (Med(1)%celerita *                &
                              Med(1)%celerita) + (Med(2)%den0 * VFmn +         &
                              Med(1)%den0 * (one - VFmn))
               elseif (pg(npi)%imed==2) then
                  pg(npi)%dens = pg(npi)%pres / (Med(2)%celerita *             &
                                 Med(2)%celerita) + (Med(2)%den0 * VFmx +      &
                                 Med(1)%den0 * (one - VFmx))
            endif
            Pg(npi)%rhoc = pg(npi)%pres / (med(2)%celerita * med(2)%celerita) +&
                           med(2)%den0
            Pg(npi)%rhow = pg(npi)%pres / (med(1)%celerita * med(1)%celerita) +&
                           med(1)%den0
         enddo
!$omp end parallel do
         call start_and_stop(3,16)
      endif
      call start_and_stop(2,20)
      if (Granular_flows_options%ID_erosion_criterion==1) call mixture_viscosity 
      call start_and_stop(3,20)
! Apparent viscosity (input option not recommended)  
      if (diffusione.or.esplosione) then
         if ((Domain%time_split==1).or.(Domain%time_stage==Domain%RKscheme))   &
            then
            call start_and_stop(2,15)
            call viscapp 
            call start_and_stop(3,15)
         endif
      endif
      if (Domain%tipo=="semi") then
! Boundary Conditions: start
! BC: checks for the particles gone out of the domain throughout the opened 
! faces

         if ((Domain%time_split==0).and.(Domain%time_stage==Domain%RKscheme))  &
            then
            call start_and_stop(2,9)
            if (NumOpenFaces>0) call CancelOutgoneParticles_3D
! Adding new particles at the inlet sections
            if (SourceFace/=0) call GenerateSourceParticles_3D
! Particle reordering on the background positioning grid
            call OrdGrid1 (nout)
            call start_and_stop(3,9)
! Set the parameters for the fixed particles 
            if (Domain%NormFix) call NormFix
            indarrayFlu = 0
            do npi=1,nag
               if ((pg(npi)%cella==0).or.(pg(npi)%vel_type/="std")) cycle
                  indarrayFlu = indarrayFlu + 1
! Array sizes check and possibile resizing
               if (indarrayFlu>PARTICLEBUFFER) then
                  call diagnostic(arg1=9,arg2=2,arg3=nomsub)
               endif
               Array_Flu(indarrayFlu) = npi
            enddo 
         endif
! Boundary Conditions: end
      endif
! Subroutine for wall BC treatment (DB-SPH)
! Density and pressure updates for wall elements: MUSCL + LPRS  
! + state equation 
      call start_and_stop(2,18)
      if ((Domain%tipo=="bsph").and.(nag>0).and.(DBSPH%n_w>0)) then
         call Gradients_to_MUSCL   
         call BC_wall_elements
      endif
      call start_and_stop(3,18)
! Update of the time stage
      if ((Domain%RKscheme>1).and.(Domain%time_split==0)) then
         Domain%time_stage = MODULO(Domain%time_stage,Domain%RKscheme)
         Domain%time_stage = Domain%time_stage + 1
            else
               done_flag = .true.         
      endif
   enddo
! Post-processing
   if (allocated(Z_fluid_max)) then
      if ((int(simulation_time/Domain%depth_dt_out)>Domain%depth_it_out_last)  &
         .or.(on_going_time_step==1)) then
         Domain%depth_it_out_last = int(simulation_time / Domain%depth_dt_out)
         call Update_Zmax_at_grid_vert_columns(1) 
         else
            call Update_Zmax_at_grid_vert_columns(0)
      endif
   endif
   if (Domain%time_split==0) dtvel = dt
   if (nout>0) then
      call Print_Results(it_eff,it_print,'loop__')
   endif
   if (nres>0) then
      call Memo_Results(it_eff,it_memo,it_rest,dtvel,'loop__')
   endif
! Computing variables at the monitoring points
   call CalcVarp
   if (Domain%icpoi_fr>0) then 
      if ((mod(it,Domain%icpoi_fr)==0).and.npointst>0) then
         call memo_ctl
      endif
      if (n_bodies>0) then
         if (mod(it,Domain%icpoi_fr)==0)  call Body_dynamics_output      
      endif
      elseif (Domain%cpoi_fr>zero) then
         if ((mod(simulation_time,Domain%cpoi_fr)<=dtvel).and.npointst>0) then
            call memo_ctl
         endif
         if (n_bodies>0) then
            if (mod(simulation_time,Domain%cpoi_fr)<=dtvel) then
               call Body_dynamics_output
            endif         
         endif
! DB-SPH post-processing and update of the variable "wet" in the array "pg"
         if ((Domain%tipo=="bsph").and.(mod(simulation_time,Domain%cpoi_fr)<=  &
            dtvel)) then
            if ((DBSPH%n_monitor_points>0).or.(DBSPH%n_monitor_regions==1))    &
               call wall_elements_pp
!$omp parallel do default(none) shared(DBSPH,pg_w) private(npi)
            do npi=1,DBSPH%n_w
               pg_w(npi)%wet = 0 
            enddo    
!$omp end parallel do
         endif
   endif
   if (Domain%ipllb_fr>0) then
      if ((mod(it,Domain%ipllb_fr)==0).and.nlines>0) then
         call calc_pelo
      endif
      elseif (Domain%pllb_fr>zero) then
         if ((mod(simulation_time,Domain%pllb_fr) <= dtvel).and.nlines>0) then
            call calc_pelo
         endif
   endif
   if (Granular_flows_options%monitoring_lines>0) then
      if ((int(simulation_time/Granular_flows_options%dt_out)>                 &
         Granular_flows_options%it_out_last).or.(on_going_time_step==1)) then
         Granular_flows_options%it_out_last = int(simulation_time /            &
            Granular_flows_options%dt_out)
         call interface_post_processing
      endif
   endif
   if (Q_sections%n_sect>0) then
      if ((int(simulation_time/Q_sections%dt_out)>Q_sections%it_out_last).or.  &
         (on_going_time_step==1)) then
         Q_sections%it_out_last = int(simulation_time / Q_sections%dt_out)
         call sub_Q_sections
      endif
   endif
! Post-processing for the water front
   if (Domain%imemo_fr>0) then
      if (mod(it,Domain%imemo_fr)==0) then
         xmax = -1.0d30
         ymax = -1.0d30
         zmax = -1.0d30
         do npi=1,nag
            if ((pg(npi)%vel_type/="std").or.(pg(npi)%cella==0)) cycle
            xmax = max(xmax,pg(npi)%coord(1))
            ymax = max(ymax,pg(npi)%coord(2))
            zmax = max(zmax,pg(npi)%coord(3))
         enddo
         write (nfro,'(4g14.7)') simulation_time,xmax,ymax,zmax
      endif
      elseif (Domain%memo_fr>zero) then
         if (it>1.and.mod(simulation_time,Domain%memo_fr) <= dtvel) then
            xmax = - 1.0d30
            ymax = - 1.0d30
            zmax = - 1.0d30
            do npi=1,nag
               if (pg(npi)%vel_type/="std".or.pg(npi)%cella==0) cycle
               xmax = max(xmax,pg(npi)%coord(1))
               ymax = max(ymax,pg(npi)%coord(2))
               zmax = max(zmax,pg(npi)%coord(3))
            enddo
            write (nfro,'(4g14.7)') simulation_time, xmax ,ymax, zmax
         endif
   endif
! Paraview output and .txt file concatenation
   if (vtkconv) then
      call result_converter ('loop__')
   endif
! If the "kill file" exists, then the run is stopped and last results are saved.
   inquire (file=nomefilekill,EXIST=kill_flag)
   if (kill_flag) exit ITERATION_LOOP
   if (simulation_time>=Domain%tmax) exit ITERATION_LOOP
enddo ITERATION_LOOP 
! Post-processing: log file
if (it_eff/=it_print.and.nout>0) then
   it_print = it_eff
   call Print_Results(it_eff,it_print,'fine__')
endif
! Post-processing: restart file
if (it_eff/=it_memo.and.nres>0) then
   call Memo_Results(it_eff,it_memo,it_rest,dtvel,'fine__')
endif
if (vtkconv) then
   call result_converter ('fine__')
endif
if (nout>0) then
   write (nout,*) " "
   write (nout,'(a)')                                                          &
"----------------------------------------------------------------------------------------"
   write (nout,*) " "
   SpCountot = 0
   do i=1,NMedium
      SpCountot = SpCountot + SpCount(i)
      write (nout,'(a,i15,a,a)')                                               &
         "Number of source particles        :  SpCount = ",SpCount(i),         &
         " medium ",Med(i)%tipo
   enddo
   write (nout,'(a,i15)') "Number of total source particles  :  SpCountot = ", &
      SpCountot
   write (nout,*) " "
   OpCountot = 0
   do i=1,NMedium
      OpCountot = OpCountot + OpCount(i)
      write (nout,'(a,i15,a,a)')                                               &
         "Number of outgone particles       :  OpCount = ",OpCount(i),         &
         " medium ",Med(i)%tipo
   enddo
   write (nout,'(a,i15)')                                                      &
      "Number of total outgone particles :  OpCountot = ",OpCountot
   write (nout,*) " "
   EpCountot = 0
   do i=1,NMedium
      EpCountot = EpCountot + EpCount(i)
      write (nout,'(a,i15,a,a)')                                               &
         "Number of escaped particles       :  EpCount = ",EpCount(i),         &
         " medium ",Med(i)%tipo
   enddo
   write (nout,'(a,i15)')                                                      &
      "Number of total escaped particles :  EpCountot = ",EpCountot
   write (nout,*) " "
   EpOrdGridtot = 0
   do i=1,NMedium
      EpOrdGridtot = EpOrdGridtot + EpOrdGrid(i)
      write (nout,'(a,i15,a,a)')                                               &
         "Number of escaped particles (OrdGrid1)       :  EpOrdGrid = ",       & 
         EpOrdGrid(i)," medium ",Med(i)%tipo
   enddo
   write (nout,'(a,i15)')                                                      &
      "Number of total escaped particles (OrdGrid1) :  EpOrdGridtot = ",       &

      EpOrdGridtot
   write (nout,*) " "
   write (nout,*) " "
   write (nout,*) "Final number of particles:       NAG = ",nag
   if (Domain%tipo=="bsph") write (nout,*)                                     &
      "Final number of wall particles:       DBSPH%n_w = ",DBSPH%n_w
   write (nout,*) " "
   write (nout,'(a)')                                                          &
"----------------------------------------------------------------------------------------"
   write (nout,*) " "
endif
if ((Domain%tipo=="bsph").and.(DBSPH%n_w>0)) deallocate(pg_w)
!------------------------
! Deallocations
!------------------------
return
end subroutine Loop_Irre_3D

