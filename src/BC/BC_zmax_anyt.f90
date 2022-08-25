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
! Program unit: BC_zmax_anyt
! Description: Dirichlet's BC for the maximum fluid height. The selected zone 
!              is continuously filled with SPH particles from the 
!              non-stationary free surface height to the given input height 
!              "zmax", under hydrostatic or homogeneous pressure conditions. 
!              Only active in 3D with SASPH. No need to impose pressure to the 
!              particles already in the "z_max" zone (ie., the particles not 
!              emitted). Each "zmax" zone needs at least 5 adjacent open 
!              sections in the input files in order to set BCs on the fluid 
!              depth. However, a "zmax" zone can be used for a more generic 
!              purpose. The associated DEM-DTM or bottom has to be locally a 
!              uniform and isotropic Cartesian grid.
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine BC_zmax_anyt
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
integer(4) :: i_zone,aux_factor,i_vertex,n_levels,ii,jj,kk,npi
double precision :: z_FS_max
integer(4),external :: ParticleCellNumber
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine EoS_barotropic_linear(k_bulk,rho_ref,p_ref,rho_in,p_in,rho_out,  &
      p_out)
      implicit none
      double precision,intent(in) :: k_bulk,rho_ref,p_ref
      double precision,intent(in),optional :: rho_in,p_in
      double precision,intent(out),optional :: rho_out,p_out
   end subroutine EoS_barotropic_linear
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
! Loop over the zones
do i_zone=1,NPartZone
! To select any zone of type "zmax"
   if (Partz(i_zone)%tipo/="zmax") cycle
   aux_factor = nint(Partz(i_zone)%dx_CartTopog / Domain%dx)
! Loop over the vertices saved at the beginning of the simulation
   do i_vertex=1,size(Partz(i_zone)%BC_zmax_vertices,1)
! Assess the maximum free surface height in the 9-point stencil around the 
! current vertex
      call z_FS_max_9p_stencil(i_zone,i_vertex,z_FS_max)
! Maybe the function "ceiling" is better than "int" here
      n_levels = int((Partz(i_zone)%H_res - z_FS_max) / Domain%dx)
! Loop over the relative IDs of the coordinates of the new particles which 
! blossom around the same bottom vertex 
!$omp parallel do default(none)                                                &
!$omp shared (Partz,i_zone,Domain,aux_factor,nag,SpCount,PARTICLEBUFFER,uerr)  &
!$omp shared (pg,PgZero,Med,i_vertex,n_levels,input_any_t)                     &
!$omp private (ii,jj,kk,npi)
      do kk=1,n_levels
         do jj=1,aux_factor
            do ii=1,aux_factor
! Generate the new particles: extrusion from the vertices saved at t=0s 
! starting from the free-surface (not from the DEM-DTM), but analogously to the 
! program unit “GeneratePart”. Locally, the extrusion is not precise as the 
! initialization of water bodies from DEM-DTM because the free surface 
! (contrarily to the DEM-DTM) is not discretized by faces (it seems 
! unefficient). However, the discretization voids are negligible and the 
! procedure dynamically adjusts every time step. Plus, no mass penetration is 
! guaranteed.
!$omp critical (omp_zmax_cs)
! Update the total number of particles and the number of particles per medium
               nag = nag + 1
               npi = nag
               SpCount(Partz(i_zone)%Medium) = SpCount(Partz(i_zone)%Medium) + 1
!$omp end critical (omp_zmax_cs)
               if (npi>PARTICLEBUFFER) then
                  write(uerr,*) 'The number of particles is larger than the ', &
                     'maximum number of particles in the program unit ',       &
                     '"BC_zmax_anyt"; the program stops here. '
                  stop
               endif
! Initialization of the quantities of the new particle
               pg(npi) = PgZero
! Coordinates. One avoids particles to be located on the very diagonals 
! of topography (otherwise some particles could not be removed even if below 
! the topography; it is unclear why this rare eventuality would happen, but no 
! problem survives). Here there is no test of visibility. The lines on the 
! particle position are rearranged from "particle_position_extrusion", where 
! particles are extruded from DEM-DTM.
               pg(npi)%coord(1) = Partz(i_zone)%BC_zmax_vertices(i_vertex,1) - &
                                  aux_factor / 2.d0 * Domain%dx + (ii - 1 +    &
                                  0.501d0) * Domain%dx
               pg(npi)%coord(2) = Partz(i_zone)%BC_zmax_vertices(i_vertex,2) - &
                                  aux_factor / 2.d0 * Domain%dx + (jj - 1 +    &
                                  0.5d0) * Domain%dx
               pg(npi)%coord(3) = (Partz(i_zone)%H_res - Domain%dx / 2.d0) -   &
                                  (kk - 1) * Domain%dx
               pg(npi)%izona = i_zone
               pg(npi)%volume = Domain%PVolume
               pg(npi)%mass = pg(npi)%volume * Med(Partz(i_zone)%Medium)%den0
               pg(nag)%dden_ALE12 = 0.d0
               pg(npi)%imed = Partz(i_zone)%Medium  
               pg(npi)%kin_visc = Med(Partz(i_zone)%Medium)%kin_visc
               pg(npi)%mu = Med(Partz(i_zone)%Medium)%kin_visc *               &
                            Med(Partz(i_zone)%Medium)%den0
               if (index(Med(Partz(i_zone)%Medium)%tipo,"liquid")>0) then
                  pg(npi)%state = "flu"
                  elseif (index(Med(Partz(i_zone)%Medium)%tipo,"granular")>0)  &
                     then
                     pg(npi)%state = "sol"
               endif
               pg(npi)%cella = ParticleCellNumber(pg(npi)%coord)
               if (Partz(i_zone)%pressure=="pa") then       
! Uniform pressure 
                  pg(npi)%pres = Partz(i_zone)%valp + Domain%prif
                  elseif (Partz(i_zone)%pressure=="qp") then
! or hydrostatic pressure (as function of "zmax")
                     pg(npi)%pres = Med(Partz(i_zone)%Medium)%den0 *           &
                                    Domain%grav(3) * (pg(npi)%coord(3) -       &
                                    Partz(i_zone)%valp) + Domain%prif
               endif
! EoS inverse
               call EoS_barotropic_linear(Med(Partz(i_zone)%Medium)%eps,       &
                  Med(Partz(i_zone)%Medium)%den0,Domain%prif,p_in=pg(npi)%pres,&
                  rho_out=pg(npi)%dens)
! Mass update
               if (input_any_t%ALE3) then
                  pg(npi)%mass = pg(npi)%dens * pg(npi)%volume
               endif
! Formal null velocity initialization (it will follow a selective partial 
! velocity smoothing)
               pg(npi)%vel(:) = 0.d0
               pg(npi)%vel_type = "std"
               pg(npi)%sect_old_pos(:) = pg(npi)%coord(:)
            enddo
         enddo
      enddo
!$omp end parallel do
   enddo
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine BC_zmax_anyt
#endif
