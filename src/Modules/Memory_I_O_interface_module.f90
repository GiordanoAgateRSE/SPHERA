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
! Program unit: Memory_I_O_interface_module
! Description: Interfaces to the program units of the folder Memory_IO
!-------------------------------------------------------------------------------
module Memory_I_O_interface_module
interface
#ifdef SOLID_BODIES
   subroutine allocate_de_BodPar_r1(allocation_flag,array,extent_1,array_name, &
      ulog_flag)
      use Hybrid_allocation_module
      implicit none
      type (body_particle),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_BodPar_r1
   subroutine allocate_de_Bod_r1(allocation_flag,array,extent_1,array_name,    &
      ulog_flag)
      use Hybrid_allocation_module
      implicit none
      type (body),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_Bod_r1
#endif
#ifdef SPACE_3D
   subroutine allocate_de_BouConEdg_r1(allocation_flag,array,extent_1,         &
      array_name,ulog_flag)
      use Hybrid_allocation_module
      implicit none
      type (TyBoundaryConvexEdge),dimension(:),allocatable,intent(inout) ::    &
         array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_BouConEdg_r1
#endif
   subroutine allocate_de_BouDat_r1(allocation_flag,array,extent_1,array_name, &
      ulog_flag)
      use Hybrid_allocation_module
      implicit none
      type (TyBoundaryData),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_BouDat_r1
#ifdef SPACE_3D
   subroutine allocate_de_BouFac_r1(allocation_flag,array,extent_1,array_name, &
      ulog_flag)
      use Hybrid_allocation_module
      implicit none
      type (TyBoundaryFace),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_BouFac_r1
#endif
#ifdef SPACE_2D
   subroutine allocate_de_BouSid_r1(allocation_flag,array,extent_1,array_name, &
      ulog_flag)
      use Hybrid_allocation_module
      implicit none
      type (TyBoundarySide),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_BouSid_r1
#endif
   subroutine allocate_de_BouStr_r1(allocation_flag,array,extent_1,array_name, &
      ulog_flag)
      use Hybrid_allocation_module
      implicit none
      type (TyBoundaryStretch),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_BouStr_r1
   subroutine allocate_de_ch100_r1(allocation_flag,array,extent_1,array_name,  &
      ulog_flag)
      implicit none
      character(100),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_ch100_r1
#ifdef SPACE_3D
   subroutine allocate_de_CLCp_r1(allocation_flag,array,extent_1,array_name,   &
      ulog_flag)
      use Hybrid_allocation_module
      implicit none
      type (CLC_polygon_der_type),dimension(:),allocatable,intent(inout) ::    &
         array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_CLCp_r1
#endif
   subroutine allocate_de_CtlLin_r1(allocation_flag,array,extent_1,array_name, &
      ulog_flag)
      use Hybrid_allocation_module
      implicit none
      type (TyCtlLine),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_CtlLin_r1
   subroutine allocate_de_CtlPoi_r1(allocation_flag,array,extent_1,array_name, &
      ulog_flag)
      use Hybrid_allocation_module
      implicit none
      type (TyCtlPoint),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_CtlPoi_r1
   subroutine allocate_de_dp_r1(allocation_flag,array,extent_1,array_name,     &
      ulog_flag)
      implicit none
      double precision,dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_dp_r1
   subroutine allocate_de_dp_r2(allocation_flag,array,extent_1,extent_2,       &
      array_name,ulog_flag)
      implicit none
      double precision,dimension(:,:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1,extent_2
      character(100),intent(in) :: array_name
   end subroutine allocate_de_dp_r2
   subroutine allocate_de_dp_r3(allocation_flag,array,extent_1,extent_2,       &
      extent_3,array_name,ulog_flag)
      implicit none
      double precision,dimension(:,:,:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1,extent_2,extent_3
      character(100),intent(in) :: array_name
   end subroutine allocate_de_dp_r3
   subroutine allocate_de_dp_r4(allocation_flag,array,extent_1,extent_2,       &
      extent_3,extent_4,array_name,ulog_flag)
      implicit none
      double precision,dimension(:,:,:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1,extent_2,extent_3,extent_4
      character(100),intent(in) :: array_name
   end subroutine allocate_de_dp_r4
   subroutine allocate_de_int4_r1(allocation_flag,array,extent_1,array_name,   &
      ulog_flag)
      implicit none
      integer(4),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_int4_r1
   subroutine allocate_de_int4_r2(allocation_flag,array,extent_1,extent_2,     &
      array_name,ulog_flag)
      implicit none
      integer(4),dimension(:,:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1,extent_2
      character(100),intent(in) :: array_name
   end subroutine allocate_de_int4_r2
   subroutine allocate_de_log_r2(allocation_flag,array,extent_1,extent_2,      &
      array_name,ulog_flag)
      implicit none
      logical,dimension(:,:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1,extent_2
      character(100),intent(in) :: array_name
   end subroutine allocate_de_log_r2
   subroutine allocate_de_Med_r1(allocation_flag,array,extent_1,array_name,    &
      ulog_flag)
      use Hybrid_allocation_module
      implicit none
      type (TyMedium),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_Med_r1   
   subroutine allocate_de_Par_r1(allocation_flag,array,extent_1,array_name,    &
      ulog_flag)
      use Hybrid_allocation_module
      implicit none
      type (TyParticle),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_Par_r1
#ifdef SPACE_3D
   subroutine allocate_de_QSec_r1(allocation_flag,array,extent_1,array_name,   &
      ulog_flag)
      use Hybrid_allocation_module
      implicit none
      type (tyQ_section_array),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_QSec_r1
#endif
#ifdef SPACE_3D
   subroutine allocate_de_Sub_r1(allocation_flag,array,extent_1,array_name,    &
      ulog_flag)
      use Hybrid_allocation_module
      implicit none
      type (type_substation),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_Sub_r1
#endif      
   subroutine allocate_de_TimSta_r1(allocation_flag,array,extent_1,array_name, &
      ulog_flag)
      use Hybrid_allocation_module
      implicit none
      type (Tytime_stage),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_TimSta_r1
   subroutine allocate_de_Zon_r1(allocation_flag,array,extent_1,array_name,    &
      ulog_flag)
      use Hybrid_allocation_module
      implicit none
      type (TyZone),dimension(:),allocatable,intent(inout) :: array
      logical,intent(in) :: allocation_flag,ulog_flag
      integer(4),intent(in),optional :: extent_1
      character(100),intent(in) :: array_name
   end subroutine allocate_de_Zon_r1
   subroutine check_max_file_unit_ID(max_file_unit_ID,file_unit_requested,     &
      file_name)
      implicit none
      integer(4),intent(inout) :: max_file_unit_ID
      integer(4),intent(in) :: file_unit_requested
      character(100),intent(in) :: file_name
   end subroutine check_max_file_unit_ID
   subroutine diagnostic(arg1,arg2,arg3)
      use Static_allocation_module
      integer(4),intent(in) :: arg1
      integer(4),intent(in),optional :: arg2
      character(len=lencard),intent(in),optional :: arg3
   end subroutine
   subroutine open_close_file(open_flag,I_O_unit,file_name)
      implicit none
      logical,intent(in) :: open_flag
      integer(4),intent(in) :: I_O_unit
      character(100),intent(in) :: file_name
   end subroutine open_close_file
end interface 
end module Memory_I_O_interface_module
