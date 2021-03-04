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
! Program unit: Memo_Results               
! Description: To write detailed results for restart.
!-------------------------------------------------------------------------------
subroutine Memo_Results(it,it_memo,it_rest,dtvel,str)
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
integer(4),intent(in) :: it
double precision,intent(in) :: dtvel
character(6),intent(in) :: str
integer(4),intent(inout) :: it_memo,it_rest
integer(4) :: nrecords,restartcode,i,i_zone,size_aux
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
restartcode = 0
!------------------------
! Statements
!------------------------
! Post-processing for first (and last) step
if (index(str,'inizio')/=0) then
   nrecords = 5
   if (NumVertici>0) nrecords = nrecords + 1
#ifdef SPACE_3D
      if (NumFacce>0) nrecords = nrecords + 2
#endif
   if (NumTratti>0) nrecords = nrecords + 1
   if (NPartZone>0) nrecords = nrecords + 1
   if (NumBVertices>0) nrecords = nrecords + 1
#ifdef SPACE_2D
      if (NumBSides>0) nrecords = nrecords + 1
#endif
#ifdef SPACE_3D
      if (Domain%tipo=="semi") then
         if ((allocated(GCBFVector)).and.(GCBFVecDim>0)) nrecords = nrecords + 1
         if ((allocated(GCBFPointers)).and.(Grid%nmax>1)) nrecords = nrecords  &
                                                                     + 1
      endif
#endif
   write(nres) version,nrecords
#ifdef SPACE_3D
   write(nres) nag,NPartZone,NumVertici,NumFacce,NumTratti,                    &   
      NumBVertices,GCBFVecDim,Grid%nmax,npointst,NPoints,NPointsl,             &
      NPointse,NLines,doubleh
#elif defined SPACE_2D
   write(nres) nag,NPartZone,NumVertici,NumTratti,                             &   
      NumBVertices,NumBSides,Grid%nmax,npointst,NPoints,NPointsl,              &
      NPointse,NLines,doubleh
#endif
   write(nres) domain
   write(nres) Grid
   if (NumVertici>0) write(nres) Vertice(1:SPACEDIM,1:NumVertici)
#ifdef SPACE_3D
   if (NumFacce>0) then
      write(nres) BoundaryFace(1:NumFacce)
      write(nres) BFaceList(1:NumFacce)
   endif
#endif
   if (NumTratti>0) write(nres) Tratto(1:NumTratti)
   do i_zone=1,NPartZone
      size_aux = 0
#ifdef SPACE_3D
      if (allocated(Partz(i_zone)%BC_zmax_vertices)) then
         size_aux = size(Partz(i_zone)%BC_zmax_vertices,1)
      endif
#endif
      write(nres) size_aux
      write(nres) Partz(i_zone)%DBSPH_fictitious_reservoir_flag,               &
         Partz(i_zone)%ipool,Partz(i_zone)%icol,                               &
         Partz(i_zone)%Medium,Partz(i_zone)%npointv,                           &
         Partz(i_zone)%IC_source_type,Partz(i_zone)%Car_top_zone,              &
         Partz(i_zone)%slip_coefficient_mode,                                  &
#ifdef SPACE_3D
         Partz(i_zone)%plan_reservoir_points,Partz(i_zone)%ID_first_vertex_sel,&
         Partz(i_zone)%ID_last_vertex_sel,Partz(i_zone)%dam_zone_ID,           &
         Partz(i_zone)%dam_zone_n_vertices,Partz(i_zone)%dx_CartTopog,         &
         Partz(i_zone)%H_res,                                                  &
#endif
         Partz(i_zone)%BC_shear_stress_input,                                  &
         Partz(i_zone)%avg_comp_slip_coeff,Partz(i_zone)%avg_ni_T_SASPH,       &
         Partz(i_zone)%avg_tau_wall_f,Partz(i_zone)%pool,Partz(i_zone)%valp,   &
         Partz(i_zone)%limit(1:2),Partz(i_zone)%vel(1:3),                      &
         Partz(i_zone)%coordMM(1:3,1:2),                                       &
         Partz(i_zone)%plan_reservoir_pos(1:4,1:2)
#ifdef SPACE_3D
      write(nres) Partz(i_zone)%dam_zone_vertices(1:4,1:2)
      if (allocated(Partz(i_zone)%BC_zmax_vertices)) then
         write(nres) Partz(i_zone)%BC_zmax_vertices(1:size_aux,1:3)
      endif
#endif
      write(nres) Partz(i_zone)%vlaw(0:3,MAXPOINTSVLAW),Partz(i_zone)%shape,   &
         Partz(i_zone)%bend,Partz(i_zone)%pressure,Partz(i_zone)%move,         &
         Partz(i_zone)%tipo,Partz(i_zone)%label
   enddo
   if (NumBVertices>0) write(nres) BoundaryVertex(1:NumBVertices)
#ifdef SPACE_2D
   if (NumBSides>0) write(nres) BoundarySide(1:NumBSides)
#endif
#ifdef SPACE_3D
      if (Domain%tipo=="semi") then
         if ((allocated(GCBFVector)).and.(GCBFVecDim>0)) write(nres)           &
            GCBFVector(1:GCBFVecDim)
         if ((allocated(GCBFPointers)).and.(Grid%nmax>1)) write(nres)          &
            GCBFPointers(1:Grid%nmax,1:2)
      endif
#endif
   flush(nres)
   write(ulog,'(a)')                                                           &
" ----------------------------------------------------------------------------"
   write(ulog,'(a,i10,a,g15.6)') " Results and restart heading saved step: ",  &
      it," time: ",simulation_time
   write(ulog,'(a)')                                                           &
" ----------------------------------------------------------------------------"
endif
if (Domain%irest_fr>0) then
   if (mod(it,Domain%irest_fr)==0) then
      it_rest = it
   endif
! Case with restart  
   elseif (input_any_t%rest_fr>zero) then
      if ((it>1).and.(mod(simulation_time,input_any_t%rest_fr)<=dtvel)) then
         it_rest = it
     endif
endif
if (Domain%imemo_fr>0) then
   if (mod(it,Domain%imemo_fr)==0) then
        it_memo = it
   endif
   elseif (input_any_t%memo_fr>zero) then
      if ((it>1).and.(mod(simulation_time,input_any_t%memo_fr)<=dtvel)) then
         it_memo = it
      endif
endif
if ((it_rest==it).or.(index(str,'inizio')/=0).or.(index(str,'fine')/=0)) then
! If restartcode=1, then to save the whole arrays "pg","pg_w"
   restartcode = 1
   write(nres) it,simulation_time,dt,nag,restartcode
   write(nres) pg(1:nag)
   if (allocated(pg_w)) write(nres)                                            &
      pg_w(1:DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet)
   do i=1,n_bodies
      write(nres) body_arr(i)%npart,body_arr(i)%Ic_imposed,                    &
         body_arr(i)%imposed_kinematics,body_arr(i)%n_records,body_arr%mass,   &
         body_arr(i)%umax,body_arr(i)%pmax,body_arr(i)%x_CM(1:3),              &
         body_arr(i)%alfa(1:3),body_arr(i)%u_CM(1:3),body_arr(i)%omega(1:3),   &
         body_arr(i)%Force(1:3),body_arr(i)%Moment(1:3),                       &
         body_arr(i)%Ic(1:3,1:3),body_arr(i)%Ic_inv(1:3,1:3)
      if ((allocated(body_arr(i)%body_kinematics)).and.                        &
         (body_arr(i)%n_records>0)) then
         write(nres) body_arr(i)%body_kinematics(1:body_arr(i)%n_records,1:7)
      endif
   enddo
   if (allocated(bp_arr)) write(nres) bp_arr(1:n_body_part)
   if (allocated(surf_body_part)) write(nres) surf_body_part(1:n_surf_body_part)
#ifdef SPACE_3D
      if (allocated(Z_fluid_max)) write(nres)                                  &
         Z_fluid_max(1:Grid%ncd(1)*Grid%ncd(2),1:2)
      if (allocated(q_max)) write(nres) q_max(1:size(q_max))
   if (allocated(substations%sub)) then
      write(nres) substations%sub(1:substations%n_sub)%POS_fsum(1),            &
         substations%sub(1:substations%n_sub)%POS_fsum(2)
      write(nres) substations%sub(1:substations%n_sub)%Ymax(1),                &
         substations%sub(1:substations%n_sub)%Ymax(2)
      write(nres) substations%sub(1:substations%n_sub)%EOT(1),                 &
         substations%sub(1:substations%n_sub)%EOT(2)
   endif
#endif
   if (allocated(Granular_flows_options%minimum_saturation_flag)) write(nres)  &
      Granular_flows_options%minimum_saturation_flag(1:Grid%ncd(1),            &
      1:Grid%ncd(2))
   if (allocated(Granular_flows_options%maximum_saturation_flag)) write(nres)  &
      Granular_flows_options%maximum_saturation_flag(1:Grid%ncd(1),            &
      1:Grid%ncd(2))
   flush(nres)
   if (index(str,'inizio')==0) then
      write(ulog,'(a)')                                                        &
      " --------------------------------------------------------------------"
      write(ulog,'(a,i10,a,g15.6)') " Results and restart saved step: ",it,    &
         " time: ",simulation_time
      write(ulog,'(a)')                                                        &
      " --------------------------------------------------------------------"
      call system("cp -p *.ris restart_tjob_0.ris0")
   endif
   elseif (it_memo==it) then
! If restartcode=0, then to save "pg" only for visualizations
      restartcode = 0
      write(nres) it,simulation_time,dt,nag,restartcode
      write(nres) pg(1:nag)%coord(1),pg(1:nag)%coord(2),pg(1:nag)%coord(3),    &
         pg(1:nag)%vel(1),pg(1:nag)%vel(2),pg(1:nag)%vel(3),pg(1:nag)%pres,    &
         pg(1:nag)%dens,pg(1:nag)%mass,pg(1:nag)%kin_visc,pg(1:nag)%imed,      &
         pg(1:nag)%icol
      flush(nres)
      if (index(str,'inizio')==0) then
         write(ulog,'(a)')                                                     &
            " --------------------------------------------------------"
         write(ulog,'(a,i10,a,g15.6)') " Results saved step: ",it," time: ",   &
            simulation_time
         write(ulog,'(a)')                                                     &
            " --------------------------------------------------------"
         call system("cp -p *.ris restart_tjob_0.ris0")
      endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine Memo_Results
