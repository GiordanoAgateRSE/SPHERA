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
! Program unit: fluid_global_quantities
! Description: Record for the output time series of the fluid global quantities:
!              mass, linear momentum, angular momentum (with respect to the 
!              domain origin), volume.
!-------------------------------------------------------------------------------
subroutine fluid_global_quantities
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Dynamic_allocation_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: ip
double precision :: glob_mass,glob_vol
double precision,dimension(3) :: glob_lin_mom,glob_ang_mom,aux_vec
character(100) :: file_name
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine Vector_Product(uu,VV,ww,SPACEDIM)
      implicit none
      integer(4),intent(in) :: SPACEDIM
      double precision,intent(in),dimension(SPACEDIM) :: uu,VV
      double precision,intent(inout),dimension(SPACEDIM) :: ww
   end subroutine
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
glob_mass = 0.d0
glob_lin_mom(1:3) = 0.d0
glob_ang_mom(1:3) = 0.d0
glob_vol = 0.d0
!------------------------
! Statements
!------------------------
do ip=1,nag
   glob_mass = glob_mass + pg(ip)%mass
   glob_lin_mom(1:3) = glob_lin_mom(1:3) + pg(ip)%mass * pg(ip)%vel(1:3)
   call Vector_Product(pg(ip)%coord,pg(ip)%vel,aux_vec,3)
   glob_ang_mom(1:3) = glob_ang_mom(1:3) + pg(ip)%mass * aux_vec(1:3)
   glob_vol = glob_vol + pg(ip)%mass / pg(ip)%dens
enddo
! File creation and heading
write(file_name,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),&
   '_fluid_global_quantities_',on_going_time_step,".txt"
call open_close_file(.true.,ufgl,file_name)
if (on_going_time_step==1) then
   write(ufgl,*) "Fluid global quantities "
   write(ufgl,'((11x,a),(10x,a),3(1x,a),3(a),(8x,a))')                         &
      " Time(s)"," Mass(kg)"," LinMomx(kg*m*s-2)"," LinMomy(kg*m*s-2)",        &
      " LinMomz(kg*m*s-2)"," AngMomx(kg*m2*s-2)"," AngMomy(kg*m2*s-2)",        &
      " AngMomz(kg*m2*s-2)"," Volume(m3)"
endif
! Output record
write(ufgl,'(9(ES19.9E2,1x))') simulation_time,glob_mass,glob_lin_mom(1:3),    &
   glob_ang_mom(1:3),glob_vol
call open_close_file(.false.,ufgl,file_name)
!------------------------
! Deallocations
!------------------------
return
end subroutine fluid_global_quantities
