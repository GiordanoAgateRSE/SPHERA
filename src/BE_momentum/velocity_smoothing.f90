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
! Program unit: velocity_smoothing 
! Description: Partial velocity smoothing (part 1). Full velocity smoothing 
!              within the "zmax" zones (part 1).
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine velocity_smoothing(BC_zmax_flag)
#elif defined SPACE_2D
subroutine velocity_smoothing
#endif
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
#ifdef SPACE_3D
logical,intent(in) :: BC_zmax_flag
#endif
integer(4) :: npi,npj,contj,npartint,ii
double precision :: rhoi,rhoj,amassj,pesoj,moddervel
double precision,dimension(3) :: dervel 
double precision,dimension(:,:),allocatable :: dervel_mat
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
if (n_bodies>0) then  
   allocate(dervel_mat(nag,3))
   dervel_mat(:,:) = 0.d0
endif
!------------------------
! Statements
!------------------------
! Body particle contributions to velocity smoothing
if (n_bodies>0) then
   call start_and_stop(3,7)
   call start_and_stop(2,19)
   call body_to_smoothing_vel(dervel_mat)
   call start_and_stop(3,19)
   call start_and_stop(2,7)
endif
!$omp parallel do default(none)                                                &
!$omp shared(pg,Med,nPartIntorno,NMAXPARTJ,PartIntorno,PartKernel,indarrayFlu) &
#ifdef SPACE_3D
!$omp shared(Array_Flu,Domain,n_bodies,dervel_mat,Partz,BC_zmax_flag)          &
#elif defined SPACE_2D
!$omp shared(Array_Flu,Domain,n_bodies,dervel_mat)                             &
#endif
!$omp private(ii,npi,contj,npartint,npj,rhoi,rhoj,amassj,dervel,moddervel)     &
!$omp private(pesoj)
do ii=1,indarrayFlu
   npi = Array_Flu(ii)
#ifdef SPACE_3D
! Selective partial smoothing only for the new particles in the "zmax" zones 
! each time step (BC condition)
   if ((BC_zmax_flag.eqv..true.).and.(Partz(pg(npi)%izona)%tipo/="zmax")) cycle
#endif
   pg(npi)%var = zero
! The mixture particles, which are temporarily affected by the frictional 
! viscosity threshold are fixed.
   if (pg(npi)%mu==Med(pg(npi)%imed)%mumx) cycle
   do contj=1,nPartIntorno(npi)
      npartint = (npi - 1) * NMAXPARTJ + contj
      npj = PartIntorno(npartint)
#ifdef SPACE_3D
! The newly generated "zmax" particles are not considered as neighbours (in 
! case of partial smoothing no "zmax" type should be active)
      if (Partz(pg(npj)%izona)%tipo=="zmax") cycle
#endif
      rhoi = pg(npi)%dens
      rhoj = pg(npj)%dens
      amassj = pg(npj)%mass
      dervel(:) = pg(npj)%vel(:) - pg(npi)%vel(:)
      if (pg(npj)%vel_type/="std") then
         rhoj = rhoi
         amassj = pg(npi)%mass
         moddervel = - two * (pg(npi)%vel(1) * pg(npj)%zer(1) + pg(npi)%vel(2) &
                     * pg(npj)%zer(2) + pg(npi)%vel(3) * pg(npj)%zer(3))
         dervel(:) = moddervel * pg(npj)%zer(:)    
      endif
      if (Med(pg(npj)%imed)%den0/=Med(pg(npi)%imed)%den0) cycle
      pesoj = amassj * PartKernel(4,npartint) / rhoj
      pg(npi)%var(:) = pg(npi)%var(:) + dervel(:) * pesoj
   enddo
   if (n_bodies>0) then
      pg(npi)%var(:) = pg(npi)%var(:) + dervel_mat(npi,:)
   endif
! Impose boundary conditions at inlet and outlet sections (DB-SPH)
   if (Domain%tipo=="bsph") then
      call DBSPH_inlet_outlet(npi)
      else
#ifdef SPACE_3D
            call velocity_smoothing_SA_SPH_3D(npi)
#elif defined SPACE_2D
               call velocity_smoothing_SA_SPH_2D(npi)
#endif
   endif
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
if (n_bodies>0) deallocate(dervel_mat)
return
end subroutine velocity_smoothing
