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
! Program unit: deallocation_sequence         
! Description: Sequence of deallocations
!-------------------------------------------------------------------------------
subroutine deallocation_sequence
!------------------------
! Modules
!------------------------
use I_O_file_module
use Static_allocation_module
use Dynamic_allocation_module
use Memory_I_O_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: i
character(100) :: array_name
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
!------------------------
! Deallocations
!------------------------
array_name = "Vertice"
call allocate_de_dp_r2(.false.,Vertice,array_name=array_name)
array_name = "Tratto"
call allocate_de_BouStr_r1(.false.,Tratto,array_name=array_name)
#ifdef SPACE_3D
array_name = "BoundaryFace"
call allocate_de_BouFac_r1(.false.,BoundaryFace,array_name=array_name)
array_name = "BFaceList"
call allocate_de_int4_r1(.false.,BFaceList,array_name=array_name)
#endif
array_name = "BoundaryVertex"
call allocate_de_int4_r1(.false.,BoundaryVertex,array_name=array_name)
#ifdef SPACE_2D
array_name = "BoundarySide"
call allocate_de_BouSid_r1(.false.,BoundarySide,array_name=array_name)
#endif
array_name = "Partz"
call allocate_de_Zon_r1(.false.,Partz,array_name=array_name)
array_name = "Med"
call allocate_de_Med_r1(.false.,Med,array_name=array_name)
array_name = "OpCount"
call allocate_de_int4_r1(.false.,OpCount,array_name=array_name)
array_name = "SpCount"
call allocate_de_int4_r1(.false.,SpCount,array_name=array_name)
array_name = "EpCount"
call allocate_de_int4_r1(.false.,EpCount,array_name=array_name)
array_name = "EpOrdGrid"
call allocate_de_int4_r1(.false.,EpOrdGrid,array_name=array_name)
array_name = "Control_Points"
call allocate_de_CtlPoi_r1(.false.,Control_Points,array_name=array_name)
array_name = "Control_Lines"
call allocate_de_CtlLin_r1(.false.,Control_Lines,array_name=array_name)
array_name = "pg"
call allocate_de_Par_r1(.false.,pg,array_name=array_name)
array_name = "ts0_pg"
call allocate_de_TimSta_r1(.false.,ts0_pg,array_name=array_name)
array_name = "NPartOrd"
call allocate_de_int4_r1(.false.,NPartOrd,array_name=array_name)
array_name = "Icont"
call allocate_de_int4_r1(.false.,Icont,array_name=array_name)
array_name = "nPartIntorno"
call allocate_de_int4_r1(.false.,nPartIntorno,array_name=array_name)
array_name = "PartIntorno"
call allocate_de_int4_r1(.false.,PartIntorno,array_name=array_name)
array_name = "PartKernel"
call allocate_de_dp_r2(.false.,PartKernel,array_name=array_name)
array_name = "rag"
call allocate_de_dp_r2(.false.,rag,array_name=array_name)
array_name = "BoundaryDataTab"
call allocate_de_BouDat_r1(.false.,BoundaryDataTab,array_name=array_name)
array_name = "BoundaryDataPointer"
call allocate_de_int4_r2(.false.,BoundaryDataPointer,array_name=array_name)
array_name = "Array_Flu"
call allocate_de_int4_r1(.false.,Array_Flu,array_name=array_name)
if (Domain%tipo=="bsph") then
   array_name = "NPartOrd_w"
   call allocate_de_int4_r1(.false.,NPartOrd_w,array_name=array_name)
   array_name = "Icont_w"
   call allocate_de_int4_r1(.false.,Icont_w,array_name=array_name)
   array_name = "nPartIntorno_fw"
   call allocate_de_int4_r1(.false.,nPartIntorno_fw,array_name=array_name)
   array_name = "PartIntorno_fw"
   call allocate_de_int4_r1(.false.,PartIntorno_fw,array_name=array_name)
   array_name = "grad_vel_VSL_fw"
   call allocate_de_dp_r2(.false.,grad_vel_VSL_fw,array_name=array_name)
   array_name = "kernel_fw"
   call allocate_de_dp_r2(.false.,kernel_fw,array_name=array_name)
   array_name = "rag_fw"
   call allocate_de_dp_r2(.false.,rag_fw,array_name=array_name)
endif
do i=1,n_bodies
   write(ulog,'(a,i6)') "Deallocation of the body n. ",i
   array_name = "body_arr(i)%body_kinematics"
   call allocate_de_dp_r2(.false.,body_arr(i)%body_kinematics,                 &
      array_name=array_name)
enddo
array_name = "body_arr"
call allocate_de_Bod_r1(.false.,body_arr,array_name=array_name)
array_name = "bp_arr"
call allocate_de_BodPar_r1(.false.,bp_arr,array_name=array_name)
array_name = "Icont_bp"
call allocate_de_int4_r1(.false.,Icont_bp,array_name=array_name)
array_name = "NPartOrd_bp"
call allocate_de_int4_r1(.false.,NPartOrd_bp,array_name=array_name)
array_name = "nPartIntorno_bp_f"
call allocate_de_int4_r1(.false.,nPartIntorno_bp_f,array_name=array_name)
array_name = "PartIntorno_bp_f"
call allocate_de_int4_r1(.false.,PartIntorno_bp_f,array_name=array_name)
array_name = "KerDer_bp_f_cub_spl"
call allocate_de_dp_r1(.false.,KerDer_bp_f_cub_spl,array_name=array_name)
array_name = "KerDer_bp_f_Gal"
call allocate_de_dp_r1(.false.,KerDer_bp_f_Gal,array_name=array_name)
array_name = "rag_bp_f"
call allocate_de_dp_r2(.false.,rag_bp_f,array_name=array_name)
array_name = "surf_body_part"
call allocate_de_int4_r1(.false.,surf_body_part,array_name=array_name)
array_name = "nPartIntorno_bp_bp"
call allocate_de_int4_r1(.false.,nPartIntorno_bp_bp,array_name=array_name)
array_name = "PartIntorno_bp_bp"
call allocate_de_int4_r1(.false.,PartIntorno_bp_bp,array_name=array_name)
array_name = "rag_bp_bp"
call allocate_de_dp_r2(.false.,rag_bp_bp,array_name=array_name)
array_name = "impact_vel"
call allocate_de_dp_r2(.false.,impact_vel,array_name=array_name)
#ifdef SPACE_3D
array_name = "BoundaryConvexEdge"
call allocate_de_BouConEdg_r1(.false.,BoundaryConvexEdge,array_name=array_name)
array_name = "GCBFVector"
call allocate_de_int4_r1(.false.,GCBFVector,array_name=array_name)
array_name = "GCBFPointers"
call allocate_de_int4_r2(.false.,GCBFPointers,array_name=array_name)
if(allocated(Q_sections%section)) then
   do i=1,Q_sections%n_sect
      write(ulog,'(a,i6)') "Deallocation of the flow-rate monitoring section ",&
         "n. ",i
      array_name = "Q_sections%section(i)%flow_rate"
      call allocate_de_dp_r1(.false.,Q_sections%section(i)%flow_rate,          &
         array_name=array_name)
   enddo
   array_name = "Q_sections%section"
   call allocate_de_QSec_r1(.false.,Q_sections%section,array_name=array_name)
endif
do i=1,substations%n_sub
   write(ulog,'(a,i6)') "Deallocation of the substation n. ",i
   array_name = "substations%sub(i)%DEMvert"
   call allocate_de_int4_r1(.false.,substations%sub(i)%DEMvert,                &
      array_name=array_name)
enddo
array_name = "substations%sub"
call allocate_de_Sub_r1(.false.,substations%sub,array_name=array_name)
array_name = "CLC%z0"
call allocate_de_dp_r2(.false.,CLC%z0,array_name=array_name)
#endif
array_name = "Granular_flows_options%lines"
call allocate_de_dp_r2(.false.,Granular_flows_options%lines,                   &
   array_name=array_name)
array_name = "Granular_flows_options%minimum_saturation_flag"
call allocate_de_log_r2(.false.,Granular_flows_options%minimum_saturation_flag,&
   array_name=array_name)
array_name = "Granular_flows_options%maximum_saturation_flag"
call allocate_de_log_r2(.false.,Granular_flows_options%maximum_saturation_flag,&
   array_name=array_name)
array_name = "Granular_flows_options%saturation_conditions"
call allocate_de_int4_r2(.false.,Granular_flows_options%saturation_conditions, &
   array_name=array_name)
array_name = "DBSPH%kinematics"
call allocate_de_dp_r3(.false.,DBSPH%kinematics,array_name=array_name)
array_name = "DBSPH%n_kinematics_records"
call allocate_de_int4_r1(.false.,DBSPH%n_kinematics_records,                   &
   array_name=array_name)
array_name = "DBSPH%rotation_centre"
call allocate_de_dp_r2(.false.,DBSPH%rotation_centre,array_name=array_name)
array_name = "DBSPH%inlet_sections"
call allocate_de_dp_r2(.false.,DBSPH%inlet_sections,array_name=array_name)
array_name = "DBSPH%outlet_sections"
call allocate_de_dp_r2(.false.,DBSPH%outlet_sections,array_name=array_name)
return
end subroutine deallocation_sequence
