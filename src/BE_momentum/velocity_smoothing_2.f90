!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2020 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: velocity_smoothing_2
! Description: Partial smoothing of the fluid velocity field (part 2). Full 
!              velocity smoothing within the "zmax" zones (part 2).
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine velocity_smoothing_2(BC_zmax_flag)
#elif defined SPACE_2D
subroutine velocity_smoothing_2
#endif
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Dynamic_allocation_module
use Hybrid_allocation_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
#ifdef SPACE_3D
logical,intent(in) :: BC_zmax_flag
#endif
integer(4) :: npi,ii
double precision :: TetaV1
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
!$omp parallel do default(none)                                                &
#ifdef SPACE_3D
!$omp shared(pg,Med,Domain,dt,indarrayFlu,Array_Flu,input_any_t,BC_zmax_flag)  &
!$omp shared(Partz,ulog)                                                       &
#elif defined SPACE_2D
!$omp shared(pg,Med,Domain,dt,indarrayFlu,Array_Flu,input_any_t,ulog)          &
#endif
!$omp private(npi,ii,TetaV1)
! Loop over all the active particles
do ii=1,indarrayFlu
   npi = Array_Flu(ii)
#ifdef SPACE_3D
! Selective partial smoothing only for the new particles in the "zmax" zones 
! each time step (BC condition)
   if ((BC_zmax_flag.eqv..true.).and.(Partz(pg(npi)%izona)%tipo/="zmax")) cycle
   if (BC_zmax_flag.eqv..true.) then
! Full smoothing for newly emitted particles in the "zmax" zones
      TetaV1 = 1.d0
! Change the new-particle zone from the current “zmax” zone to the associated 
! boundary zone (“Car_top_zone”) and initialize the "vel_type" so that the 
! newly-initialized particles will be treated as standard computational 
! particles after initialization.
      pg(npi)%izona = Partz(pg(npi)%izona)%Car_top_zone
      pg(npi)%vel_type = "std"
      write(ulog,'(2a,a5,i12)') 'Program unit "velocity_smoothing_2": ',       &
         'pg(npi)%vel_type, npi: ',pg(npi)%vel_type,npi
      else
#endif
         TetaV1 = input_any_t%TetaV
#ifdef SPACE_3D
   endif
#endif
   if (TetaV1>1.d-9) then
! TetaV depending on the time step
      TetaV1 = TetaV1 * Med(pg(npi)%imed)%Celerita * dt / Domain%h
      if (pg(npi)%kodvel==0) then             
! The particle is inside the domain and far from boundaries
         pg(npi)%var(:) = pg(npi)%vel(:) + TetaV1 * pg(npi)%var(:)     
         pg(npi)%vel(:) = pg(npi)%var(:)                                   
         else
! The particle is close to a "source", "level" or "normal velocity boundary
! (kodvel=1 or kodvel=2): the final velocity is kept unmodified
            pg(npi)%var(:) = pg(npi)%vel(:)                                
      endif
      else
         pg(npi)%var(:) = pg(npi)%vel(:)
   endif
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine velocity_smoothing_2
