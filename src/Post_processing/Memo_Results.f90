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
integer(4),intent(IN) :: it
double precision,intent(IN) :: dtvel
character(6),intent(IN) :: str
integer(4),intent(INOUT) :: it_memo
integer(4),intent(INOUT) :: it_rest
integer(4) :: nrecords,restartcode,i
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
   if (NumFacce>0) nrecords = nrecords + 1
   if (NumFacce>0) nrecords = nrecords + 1
   if (NumTratti>0) nrecords = nrecords + 1
   if (NPartZone>0) nrecords = nrecords + 1
   if (NumBVertices>0) nrecords = nrecords + 1
   if (NumBSides>0) nrecords = nrecords + 1
   if (Domain%tipo=="semi") then
      if ((allocated(GCBFVector)).and.(GCBFVecDim>0)) nrecords = nrecords + 1
      if ((allocated(GCBFPointers)).and.(Grid%nmax>1)) nrecords = nrecords + 1
   endif
   write(nres) version,nrecords
   write(nres) ncord,Nag,NMedium,NPartZone,NumVertici,NumFacce,NumTratti,      &   
      NumBVertices,NumBSides,GCBFVecDim,Grid%nmax,NPointst,NPoints,NPointsl,   &
      NPointse,NLines,NSections,doubleh
   write(nres) domain
   write(nres) grid
   write(nres) Med(1:NMedium)
   if (NumVertici>0) write(nres) Vertice(1:SPACEDIM,1:NumVertici)
   if (NumFacce>0) write(nres) BoundaryFace(1:NumFacce)
   if (NumFacce>0) write(nres) BFaceList(1:NumFacce)
   if (NumTratti>0) write(nres) Tratto(1:NumTratti)
   if (NPartZone>0) write(nres) Partz(1:NPartZone)
   if (NumBVertices>0) write(nres) BoundaryVertex(1:NumBVertices)
   if (NumBSides>0) write(nres) BoundarySide(1:NumBSides)
   if (Domain%tipo=="semi") then
      if ((allocated(GCBFVector)).and.(GCBFVecDim>0)) write(nres)              &
         GCBFVector(1:GCBFVecDim)
      if ((allocated(GCBFPointers)).and.(Grid%nmax>1)) write(nres)             &
         GCBFPointers(1:Grid%nmax,1:2)
   endif
   flush(nres)
   write(nout,'(a,i10,a,f15.5)')                                               &
" ----------------------------------------------------------------------------"
   write(nout,'(a,i10,a,f15.5)') " Results and restart heading saved   step: ",&
      it,"   time: ",simulation_time
   write(nout,'(a,i10,a,f15.5)')                                               &
" ----------------------------------------------------------------------------"
endif
if (Domain%irest_fr>0) then
   if (mod(it,Domain%irest_fr)==0) then
      it_rest = it
   endif
! Case with restart  
   elseif (Domain%rest_fr>zero) then
      if ((it>1).and.(mod(simulation_time,Domain%rest_fr)<=dtvel)) then
         it_rest = it
     endif
endif
if (Domain%imemo_fr>0) then
   if (mod(it,Domain%imemo_fr)==0) then
        it_memo = it
   endif
   elseif (Domain%memo_fr>zero) then
      if ((it>1).and.(mod(simulation_time,Domain%memo_fr)<=dtvel)) then
         it_memo = it
      endif
endif
if ((it_rest==it).or.(index(str,'inizio')/=0).or.(index(str,'fine')/=0)) then
! If restartcode=1, then to save the whole arrays "pg","pg_w"
   restartcode = 1
   write(nres) it,simulation_time,dt,nag,ncord,restartcode
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
   if (allocated(Z_fluid_max)) write(nres)                                     &
      Z_fluid_max(1:Grid%ncd(1)*Grid%ncd(2))
   if (allocated(q_max)) write(nres) q_max(1:size(q_max))
   if (allocated(Granular_flows_options%minimum_saturation_flag)) write(nres)  &
      Granular_flows_options%minimum_saturation_flag(1:Grid%ncd(1),            &
      1:Grid%ncd(2))
   if (allocated(Granular_flows_options%maximum_saturation_flag)) write(nres)  &
      Granular_flows_options%maximum_saturation_flag(1:Grid%ncd(1),            &
      1:Grid%ncd(2))
   flush(nres)
   if (index(str,'inizio')==0) then
      write(nout,'(a,i10,a,f15.5)')                                            &
      " --------------------------------------------------------------------"
      write(nout,'(a,i10,a,f15.5)') " Results and restart saved   step: ",it,  &
         "   time: ",simulation_time
      write(nout,'(a,i10,a,f15.5)')                                            &
      " --------------------------------------------------------------------"
   endif
   elseif (it_memo==it) then
! If restartcode=0, then to save "pg" only for visualizations
      restartcode = 0
      write(nres) it,simulation_time,dt,nag,ncord,restartcode
      write(nres) pg(1:nag)%coord(1),pg(1:nag)%coord(2),pg(1:nag)%coord(3),    &
         pg(1:nag)%vel(1),pg(1:nag)%vel(2),pg(1:nag)%vel(3),pg(1:nag)%pres,    &
         pg(1:nag)%dens,pg(1:nag)%mass,pg(1:nag)%visc,pg(1:nag)%IntEn,         &
         pg(1:nag)%VolFra,pg(1:nag)%imed,pg(1:nag)%icol
      flush(nres)
      if (index(str,'inizio')==0) then
         write(nout,'(a,i10,a,f15.5)')                                         &
            " --------------------------------------------------------"
         write(nout,'(a,i10,a,f15.5)') " Results saved   step: ",it,"   time: "&
            ,simulation_time
         write(nout,'(a,i10,a,f15.5)')                                         &
            " --------------------------------------------------------"
      endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine Memo_Results

