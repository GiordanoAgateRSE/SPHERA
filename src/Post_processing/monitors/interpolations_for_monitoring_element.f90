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
! Program unit: interpolations_for_monitoring_element
! Description:  SPH interpolation of physical quantities at a monitoring 
!               element point using the neighbouring SPH elements.
!-------------------------------------------------------------------------------
subroutine interpolations_for_monitoring_element(pglocal)
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
integer(4) :: nceli,ncel,igridi,kgridi,jgridi,irang,krang,jrang,fw
integer(4) :: mm,npj,irestocell
#ifdef SOLID_BODIES
integer(4) :: sb
#endif
double precision :: rijlocal,uni, pesoj,plocal,rho
double precision,dimension(3) :: raglocal,vel
type (TyCtlPoint) :: pglocal
integer(4),external :: CellIndices,CellNumber
double precision,external :: w
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
nceli = pglocal%cella
if (nceli==0) return
irestocell = CellIndices(nceli,igridi,jgridi,kgridi)
uni = zero
plocal = zero
rho = zero
vel(:) = zero
!------------------------
! Statements
!------------------------
! Loop over the grid cells
do jrang = jgridi-1,jgridi+1
   do irang = igridi-1,igridi+1 
      do krang = kgridi-1,kgridi+1
         ncel = CellNumber (irang,jrang,krang)
         if (ncel==0) cycle
! Contributions from fluid particles
         if (Icont(ncel+1)>Icont(ncel)) then
! Loop over the cell particles
            do mm=Icont(ncel),Icont(ncel+1)-1
               npj = NPartOrd(mm)
               if (pg(npj)%vel_type=="fix") cycle
               raglocal(:) = pglocal%coord(:) - pg(npj)%coord(:)
               rijlocal = dot_product(raglocal,raglocal)
               if (rijlocal>square_doubleh) cycle
               rijlocal = dsqrt(rijlocal) 
               pesoj = pg(npj)%mass * w(rijlocal,Domain%h,Domain%coefke) /     &
                       pg(npj)%dens
               uni = uni + pesoj
! Update of the monitoring element pressure
               plocal  = plocal  + pg(npj)%pres * pesoj        
! Update of the monitoring element density (mono-phase SPH approximation)
               rho = rho + pg(npj)%dens * pesoj
               vel(:) = vel(:) + pg(npj)%vel(:) * pesoj
            enddo
         endif
! Contributions from wall elements
         if ((Domain%tipo=="bsph").and.(DBSPH%n_w>0)) then
            if (Icont_w(ncel+1)>Icont_w(ncel)) then
! Loop over the neighbouring wall particles in the cell
               do fw=Icont_w(ncel),Icont_w(ncel+1)-1
                  npj = NPartOrd_w(fw)
! Relative positions and distances
                  raglocal(:) = pglocal%coord(:) - pg_w(npj)%coord(:)
                  rijlocal = dot_product(raglocal,raglocal)
! The wall element lies within the "kernel support" of the monitoring element
                  if (rijlocal>square_doubleh) cycle
                  rijlocal = dsqrt(rijlocal)
                  pesoj = pg_w(npj)%mass * w(rijlocal,Domain%h,Domain%coefke) /&
                          pg_w(npj)%dens
                  uni = uni + pesoj
                  plocal = plocal  + pg_w(npj)%pres * pesoj
                  rho = rho + pg_w(npj)%dens * pesoj                  
                  vel(:) = vel(:) + pg_w(npj)%vel(:) * pesoj        
               enddo
            endif
         endif
#ifdef SOLID_BODIES
! Contributions from body particles
            if (Icont_bp(ncel+1)>Icont_bp(ncel)) then
! Loop over the neighbouring solid particles in the cell
               do sb=Icont_bp(ncel),Icont_bp(ncel+1)-1
                  npj = NPartOrd_bp(sb)
! Relative positions and distances
                  raglocal(:) = pglocal%coord(:) - bp_arr(npj)%pos(:)
                  rijlocal = dot_product(raglocal,raglocal)
                  if (rijlocal>square_doubleh) cycle
! The body particle lies within the "kernel support" of the monitoring element
                  rijlocal = dsqrt(rijlocal)
! "Inner" and surface body particles are both useful
                  pesoj = w(rijlocal,Domain%h,Domain%coefke) *                 &
                          bp_arr(npj)%volume
                  uni = uni + pesoj
                  plocal = plocal  + bp_arr(npj)%pres * pesoj
                  vel(:) = vel(:) + bp_arr(npj)%vel_mir(:) * pesoj
! Fluid density cannot be interpolated from solid body particles
               enddo
            endif
#endif
      enddo
   enddo
enddo
pglocal%uni = uni
if (uni>=1.d-1) then
! In case the dicrete Shepard coefficient is smaller than 0.1 the monitor 
! element assumes null pressure, velocity and density.
! SPH interpolations at the position of the monitoring element from the values 
! of the neighbouring fluid particles, wall elements and body particles.
         pglocal%pres = plocal/uni
         pglocal%dens = rho/uni
         pglocal%vel(:) = vel(:)/uni
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine interpolations_for_monitoring_element
