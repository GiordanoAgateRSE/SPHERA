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
! Program unit: Continuity_Equation
! Description: To accumulate contributions for the continuity equation.
!              Computation of velocity gradients and the second invariant of the
!              strain-rate tensor (with 2D renormalization).    
!-------------------------------------------------------------------------------
subroutine Continuity_Equation(npi)
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
integer(4),intent(in) :: npi
integer(4) :: npj,contj,npartint
double precision :: rhoi,rhoj,amassj,moddervel,appo
#ifdef SPACE_3D
double precision :: moddia,modout
#elif defined SPACE_2D
double precision :: det
#endif
double precision,dimension(3) :: pesogradj,dvar
double precision,dimension(9) :: dvdi
#ifdef SPACE_2D
double precision,dimension(9) :: aij
#endif
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
if ((pg(npi)%vel_type/="std").or.(pg(npi)%cella==0)) return
pg(npi)%dden  = zero
dvar(:) = zero
#ifdef SPACE_2D
   aij(:) = zero
#endif
dvdi(:) = zero
! First loop to find interacting particles and saving 
do contj=1,nPartIntorno(npi)
   npartint = (npi - 1) * NMAXPARTJ + contj
   npj = PartIntorno(npartint)
   rhoi = pg(npi)%dens
   rhoj = pg(npj)%dens
   amassj = pg(npj)%mass
   dvar(:) = pg(npj)%var(:) - pg(npi)%var(:)
   if (pg(npj)%vel_type/="std") then          
      rhoj = rhoi
      amassj = pg(npi)%mass
! Particles with "fix" movement, but not still.
      dvar(:) = pg(npj)%vel(:) - pg(npi)%var(:)
      if ((pg(npj)%vel(1)==zero).and.(pg(npj)%vel(2)==zero).and.               &
         (pg(npj)%vel(3)==zero)) then
         if (pg(npj)%slip=="f") then 
            moddervel = - two * (pg(npi)%var(1) * pg(npj)%zer(1) +             &
                        pg(npi)%var(3) * pg(npj)%zer(3))
            dvar(:) = moddervel * pg(npj)%zer(:)
         endif
      endif
   endif
! Continuity equation
   pesogradj(:) = amassj * rag(:,npartint) * PartKernel(1,npartint) / rhoj
   if (Granular_flows_options%KTGF_config.ne.1) then
      appo = amassj * PartKernel(1,npartint) *                                 &
             (dvar(1)*rag(1,npartint) + dvar(2)*rag(2,npartint) +              &
             dvar(3)*rag(3,npartint))
      else
         appo = rhoi * PartKernel(1,npartint) * (amassj / rhoj) *              &
                (dvar(1) * rag(1,npartint) + dvar(2) * rag(2,npartint) +       &
                dvar(3) * rag(3,npartint))
   endif
   pg(npi)%dden = pg(npi)%dden - appo
! Velocity derivatives 
   if (pg(npj)%vel_type/="std") return
   if (Granular_flows_options%KTGF_config==1) then
#ifdef SPACE_3D
!du/dy
      dvdi(2) = dvdi(2) + pesogradj(2) * dvar(1)
!dv/dx
      dvdi(4) = dvdi(4) + pesogradj(1) * dvar(2)
!dv/dy
      dvdi(5) = dvdi(5) + pesogradj(2) * dvar(2)
!dv/dz
      dvdi(6) = dvdi(6) + pesogradj(3) * dvar(2)
!dw/dy
      dvdi(8) = dvdi(8) + pesogradj(2) * dvar(3)
#endif
!du/dx
      dvdi(1) = dvdi(1) + pesogradj(1) * dvar(1)
!du/dz 
      dvdi(3) = dvdi(3) + pesogradj(3) * dvar(1)
!dw/dx
      dvdi(7) = dvdi(7) + pesogradj(1) * dvar(3)
!dw/dz
      dvdi(9) = dvdi(9) + pesogradj(3) * dvar(3)
#ifdef SPACE_2D
         aij(1) = aij(1) + (pg(npj)%coord(1) - pg(npi)%coord(1)) * pesogradj(1)
         aij(2) = aij(2) + (pg(npj)%coord(3) - pg(npi)%coord(3)) * pesogradj(1)
         aij(3) = aij(3) + (pg(npj)%coord(1) - pg(npi)%coord(1)) * pesogradj(3)
         aij(4) = aij(4) + (pg(npj)%coord(3) - pg(npi)%coord(3)) * pesogradj(3)
#endif
   endif
enddo
! Boundary contributions (DB-SPH)
#ifdef SPACE_2D
   if (Domain%tipo=="bsph") then
      do contj=1,nPartIntorno_fw(npi)
         npartint = (npi - 1) * NMAXPARTJ + contj
         npj = PartIntorno_fw(npartint)
         dvar(:) = pg_w(npj)%vel(:) - pg(npi)%var(:)
         appo = - pg_w(npj)%dens * pg_w(npj)%weight * kernel_fw(1,npartint)    &
                * pg(npi)%Gamma * (dvar(1) * pg_w(npj)%normal(1) + dvar(2)     &
                * pg_w(npj)%normal(2) + dvar(3) * pg_w(npj)%normal(3))
         appo = appo - pg_w(npj)%mass * kernel_fw(2,npartint) *                &
                (dvar(1) * rag_fw(1,npartint) + dvar(2) *                      &
                rag_fw(2,npartint) + dvar(3) * rag_fw(3,npartint))
         pg(npi)%dden = pg(npi)%dden + appo
      enddo
   endif
#endif
if (Granular_flows_options%KTGF_config==1) then
#ifdef SPACE_3D
      moddia = (dvdi(1) * dvdi(1) + dvdi(5) * dvdi(5) + dvdi(9) * dvdi(9))
      modout = ((dvdi(2) + dvdi(4)) * (dvdi(2) + dvdi(4)) + (dvdi(3) +         &
                dvdi(7)) * (dvdi(3) + dvdi(7)) + (dvdi(6)+dvdi(8))*(dvdi(6)+   &
                dvdi(8)))
      pg(npi)%secinv = dsqrt( half * moddia + quarter * modout)
#elif defined SPACE_2D
         det = aij(1) * aij(4) - aij(2) * aij(3)
         if (det<1.d-3) det = 1.d-3
         pg(npi)%dudx = (dvdi(1) * aij(4) - dvdi(3) * aij(2)) / det       
         pg(npi)%dudy = ( - dvdi(1) * aij(3) + dvdi(3) * aij(1)) / det       
         pg(npi)%dvdx = ( dvdi(7) * aij(4) - dvdi(9) * aij(2)) / det       
         pg(npi)%dvdy = ( - dvdi(7) * aij(3) + dvdi(9) * aij(1)) / det       
         pg(npi)%secinv = dsqrt(half * pg(npi)%dudx * pg(npi)%dudx + half * &
                          pg(npi)%dvdy * pg(npi)%dvdy + quarter *           &
                          (pg(npi)%dudy + pg(npi)%dvdx) * (pg(npi)%dudy +   &
                          pg(npi)%dvdx)) 
#endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine Continuity_Equation
