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
! Program unit: BC_wall_elements
! Description: Wall element density and pressure (Amicarelli et al., 2013, 
!              IJNME). Shepard correction for the kinematic viscosity of the 
!              semi-particles. Shepard correction for the velocity gradient 
!              (times shear viscosity) of the semi-particles in the Viscous  
!              Sub-Layer of the Surface Neutral Boundary Layer.         
!-------------------------------------------------------------------------------
subroutine BC_wall_elements
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
integer(4) :: contj,npartint,npi,npj
double precision :: rhoL,uL,uR,cL,rhostar,pstar,Ww_Shep
double precision :: uCartL(3)
! Array to detect wall elements with fluid neighbours
integer(4), dimension(:), allocatable :: neigh_w
double precision,dimension(:), allocatable :: den
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
allocate(den(DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet))
allocate(neigh_w(DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet))
!------------------------
! Initializations
!------------------------
den = zero
neigh_w = 0
!------------------------
! Statements
!------------------------
! Free-surface condition (a sort of zeroing) for wall pressure and density
!$omp parallel do default(none) private(npi) shared (DBSPH,pg_w,med)
do npi=1,DBSPH%n_w
   pg_w(npi)%pres = zero
   pg_w(npi)%dens = med(1)%den0
   if (pg_w(npi)%sigma>1d-6) then
! Kinematic viscosity of the semi-particles
      pg_w(npi)%kin_visc_semi_part = pg_w(npi)%kin_visc_semi_part /            &
                                     pg_w(npi)%sigma
! Semi-particle velocity gradient (times shear viscosity) in the Viscous  
! Sub-Layer of the Surface Neutral Boundary Layer: normalization
      pg_w(npi)%grad_vel_VSL_times_mu(:) = pg_w(npi)%grad_vel_VSL_times_mu(:) /&
                                        pg_w(npi)%sigma
      else
         pg_w(npi)%kin_visc_semi_part = 0.d0
         pg_w(npi)%grad_vel_VSL_times_mu(:) = 0.d0   
   endif
enddo
! Loop over computational fluid particles for density and pressure 
! contributions to wall elements
!$omp end parallel do
do npi=1,nag
   if (pg(npi)%imed/=pg(1)%imed) cycle
! Loop over the wall neighbnours of the computational particle
   do contj=1,nPartIntorno_fw(npi)
      npartint = (npi - 1) * NMAXPARTJ + contj
      npj = PartIntorno_fw(npartint)
! Zeroing wall density when a wall element has fluid neighbours 
! (just the first time it is encountered)
      if (neigh_w(npj)==0) then
         pg_w(npj)%dens = zero
         neigh_w(npj) = 1
      endif
! MUSCL reconstruction (first order): it provides the initial left side 
! conditions for the Riemann solver (subscripts L).
! The right state (R) is defined by the wall element properties.
! ucartL represents the left state condition in Cartesian coordinates: 
! uL is its projection over the wall normal.
      rhoL = pg(npi)%dens - (pg(npi)%drho(1) * rag_fw(1,npartint) +            &
             pg(npi)%drho(2) * rag_fw(2,npartint) +                            &
             pg(npi)%drho(3) * rag_fw(3,npartint))
      uCartL(1) = pg(npi)%var(1) - (pg(npi)%dvel(1,1) * rag_fw(1,npartint)     &
                  + pg(npi)%dvel(1,2) * rag_fw(2,npartint) +                   &
                  pg(npi)%dvel(1,3)*rag_fw(3,npartint))
      uCartL(2) = pg(npi)%var(2) - (pg(npi)%dvel(2,1) * rag_fw(1,npartint)     &
                  + pg(npi)%dvel(2,2) * rag_fw(2,npartint) +                   &
                  pg(npi)%dvel(2,3) * rag_fw(3,npartint))
      uCartL(3) = pg(npi)%var(3) - (pg(npi)%dvel(3,1) * rag_fw(1,npartint)     &
                  + pg(npi)%dvel(3,2) * rag_fw(2,npartint) +                   &
                  pg(npi)%dvel(3,3) * rag_fw(3,npartint))
      uL = uCartL(1) * pg_w(npj)%normal(1) + uCartL(2) * pg_w(npj)%normal(2)   &
           + uCartL(3) * pg_w(npj)%normal(3)
! Linearized Partial Riemann Solver: start
! The partial Riemann problem is here 1D, oriented along the wall normal. 
! Solution (density) in the central zone (here the wall element): subscript star
      uR = pg_w(npj)%vel(1) * pg_w(npj)%normal(1) + pg_w(npj)%vel(2) *         &
           pg_w(npj)%normal(2) + pg_w(npj)%vel(3) * pg_w(npj)%normal(3)
      cL = dsqrt(Med(pg(npi)%imed)%eps/rhoL)
      rhostar = rhoL + (uL - uR) * rhoL / cL
      if (rhostar<10.d0) then
         rhostar = 10.d0
      endif
! EoS
      call EoS_barotropic_linear(Med(pg(npi)%imed)%eps,Med(pg(npi)%imed)%den0, &
         Domain%prif,rho_in=rhostar,p_out=pstar)
      if (pstar<-99.d3) then
         pstar = -99.d3
      endif
! Linearized Partial Riemann Solver: end
! "den" is an auxiliary vector to add contributions to the density and 
! pressure denominators
      Ww_Shep = pg(npi)%mass / pg(npi)%dens * kernel_fw(1,npartint) 
      den(npj) = den(npj) + Ww_Shep
      pg_w(npj)%dens = pg_w(npj)%dens + rhostar * Ww_Shep
      pg_w(npj)%pres = pg_w(npj)%pres + pstar * Ww_Shep
   enddo
enddo
! Updating wall element pressure and density
!$omp parallel do default(none) private(npi) shared(neigh_w,pg_w,den,DBSPH)
do npi=1,DBSPH%n_w
   if (neigh_w(npi)==1) then
      pg_w(npi)%dens = pg_w(npi)%dens / den(npi)
      pg_w(npi)%pres = pg_w(npi)%pres / den(npi)
! non-negative pressure values
      if ((DBSPH%negative_wall_p_allowed.eqv..false.).and.                     &
         (pg_w(npi)%pres<0.d0)) pg_w(npi)%pres = zero 
   endif
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
deallocate(den)
deallocate(neigh_w) 
return
end subroutine BC_wall_elements
