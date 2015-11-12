!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-Ricerca di Sistema)



! SPHERA authors and email contact are provided on SPHERA documentation.

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
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: inter_EqCont_2D 
! Description: To accumulate contributions for the 2D continuity equation. Computation of velocity gradients and the second 
!              invariant of the strain-rate tensor.        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine inter_EqCont_2D 
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
integer(4) :: npi,npj,contj,npartint   
double precision :: rhoi,rhoj,amassj,moddervel,det,appo,factdiff, rvw, dervol   
double precision,dimension(3) :: pesogradj,dvar      
double precision,dimension(4) :: derspa,aij
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
!$omp parallel do default(none) &
!$omp private(npi,derspa,dvar,aij,contj,npartint,npj,rhoi,rhoj,amassj)         &         
!$omp private(moddervel,pesogradj,dervol,factdiff,rvw,det,appo)                &
!$omp shared(nag,pg,nPartIntorno,NMAXPARTJ,PartIntorno,PartKernel,rag)         &
!$omp shared(diffusione,pg_w,nPartIntorno_fw,PartIntorno_fw,DBSPH,kernel_fw)   &
!$omp shared(rag_fw,Domain,Granular_flows_options)
! Loop over all the active particles
do npi=1,nag
   if ((pg(npi)%vel_type/="std").or.(pg(npi)%cella==0)) cycle
   pg(npi)%dden  = zero
   pg(npi)%diffu = zero
   derspa(:) = zero
   dvar(:) = zero
   aij(:) = zero
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
! fix = .TRUE.
         if ((pg(npj)%vel(1)==zero).AND.(pg(npj)%vel(2)==zero).AND.            &
            (pg(npj)%vel(3)==zero)) then   
            if (pg(npj)%slip=="f") then 
               moddervel = - two * (pg(npi)%var(1) * pg(npj)%zer(1) +          &
                           pg(npi)%var(3) * pg(npj)%zer(3))
               dvar(:) = moddervel * pg(npj)%zer(:)
            end if
         end if
      end if
! Continuity equation
      pesogradj(1:3) = amassj * rag(1:3,npartint) * PartKernel(1,npartint) /   &
                       rhoj
      if (Granular_flows_options%ID_erosion_criterion.ne.1) then
         appo = amassj * PartKernel(1,npartint) *                              &
               (dvar(1)*rag(1,npartint) + dvar(2)*rag(2,npartint) +            &
               dvar(3)*rag(3,npartint))
         else
            appo = rhoi * PartKernel(1,npartint) * (amassj/rhoj) *             &
                   (dvar(1)*rag(1,npartint) + dvar(2)*rag(2,npartint) +        &
                   dvar(3)*rag(3,npartint))
      endif
      pg(npi)%dden = pg(npi)%dden - appo
      if (diffusione) then
         dervol = pg(npj)%VolFra - pg(npi)%VolFra
         call diffumorris (npi,npj,npartint,dervol,factdiff,rvw)
         pg(npi)%diffu = pg(npi)%diffu + factdiff * rvw
      end if
! Velocity derivatives  
      if (pg(npj)%vel_type/="std") cycle
      if ((Granular_flows_options%ID_erosion_criterion==1).or.                 &
         (Granular_flows_options%ID_erosion_criterion==3)) then     
         derspa(1) = derspa(1) + pesogradj(1) * dvar(1)
         derspa(2) = derspa(2) + pesogradj(3) * dvar(1)
         derspa(3) = derspa(3) + pesogradj(1) * dvar(3)
         derspa(4) = derspa(4) + pesogradj(3) * dvar(3)
         aij(1) = aij(1) + (pg(npj)%coord(1) - pg(npi)%coord(1)) * pesogradj(1)
         aij(2) = aij(2) + (pg(npj)%coord(3) - pg(npi)%coord(3)) * pesogradj(1)
         aij(3) = aij(3) + (pg(npj)%coord(1) - pg(npi)%coord(1)) * pesogradj(3)
         aij(4) = aij(4) + (pg(npj)%coord(3) - pg(npi)%coord(3)) * pesogradj(3)
      endif      
   end do
! Boundary contributions (DB-SPH)
   if (Domain%tipo=="bsph") then
      do contj=1,nPartIntorno_fw(npi)
         npartint = (npi - 1) * NMAXPARTJ + contj
         npj = PartIntorno_fw(npartint)
         dvar(:) = pg_w(npj)%vel(:) - pg(npi)%var(:)
         appo = - pg_w(npj)%dens * pg_w(npj)%weight * kernel_fw(1,npartint) *  &
                pg(npi)%Gamma * (dvar(1) * pg_w(npj)%normal(1) + dvar(2) *     &
                pg_w(npj)%normal(2) + dvar(3) * pg_w(npj)%normal(3))
         appo = appo - pg_w(npj)%mass * kernel_fw(2,npartint) *                &
                (dvar(1) * rag_fw(1,npartint) + dvar(2) * rag_fw(2,npartint) + &
                dvar(3) * rag_fw(3,npartint))
         pg(npi)%dden = pg(npi)%dden + appo
      end do
   endif
   if ((Granular_flows_options%ID_erosion_criterion==1).or.                    &
      (Granular_flows_options%ID_erosion_criterion==3)) then
      det = aij(1) * aij(4) - aij(2) * aij(3)
      if (det<0.001d0) det = 0.001d0
      pg(npi)%dudx = (derspa(1) * aij(4) - derspa(2) * aij(2)) / det       
      pg(npi)%dudy = ( - derspa(1) * aij(3) + derspa(2) * aij(1)) / det       
      pg(npi)%dvdx = ( derspa(3) * aij(4) - derspa(4) * aij(2)) / det       
      pg(npi)%dvdy = ( - derspa(3) * aij(3) + derspa(4) * aij(1)) / det       
      pg(npi)%secinv = Dsqrt(half * pg(npi)%dudx * pg(npi)%dudx + half *       &
                       pg(npi)%dvdy * pg(npi)%dvdy + quarter * (pg(npi)%dudy   & 
                       + pg(npi)%dvdx) * (pg(npi)%dudy + pg(npi)%dvdx)) 
   endif
end do
!$omp end parallel do
call start_and_stop(3,12)
call start_and_stop(2,19)
if (n_bodies>0) call body_particles_to_continuity
call start_and_stop(3,19)
call start_and_stop(2,12)
!------------------------
! Deallocations
!------------------------
return
end subroutine inter_EqCont_2D

