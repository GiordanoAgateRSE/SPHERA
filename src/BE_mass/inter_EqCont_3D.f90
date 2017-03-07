!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
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
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: inter_EqCont_2D 
! Description: To accumulate contributions for the 3D continuity equation.
!              Computation of velocity gradients and the second invariant of the
!              strain-rate tensor.        
!-------------------------------------------------------------------------------
subroutine inter_EqCont_3D
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
double precision :: rhoi,rhoj,amassj,moddervel,det,moddia,modout,appo,factdiff 
double precision :: rvw,dervol
double precision,dimension (3) :: pesogradj,dvar
double precision,dimension (9) :: derspa,aij,bij,dvdi
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
!$omp private(npi,derspa,dvar,aij,bij,dvdi,contj,npartint,npj,rhoi,rhoj)       &
!$omp private(amassj,moddervel,pesogradj,dervol,factdiff,rvw,det,moddia)       &
!$omp private(modout,appo)                                                     &
!$omp shared(nag,pg,nPartIntorno,NMAXPARTJ,PartIntorno,PartKernel,rag)         &
!$omp shared(diffusione,Granular_flows_options)
! Loop over all the active particles
do npi=1,nag
   if ((pg(npi)%vel_type/="std").or.(pg(npi)%cella==0)) cycle
   pg(npi)%dden  = zero
   pg(npi)%diffu = zero
   derspa(:) = zero
   dvar(:)   = zero
   aij(:)    = zero
   bij(:)    = zero
   dvdi(:)   = zero
! First loop to find interacting particles and saving 
   do contj=1,nPartIntorno(npi)
      npartint = (npi - 1)* NMAXPARTJ + contj
      npj = PartIntorno(npartint)
      rhoi   = pg(npi)%dens
      rhoj   = pg(npj)%dens
      amassj = pg(npj)%mass
      dvar(:) = pg(npj)%var(:) - pg(npi)%var(:)
      if (pg(npj)%vel_type/="std") then          
         rhoj = rhoi
         amassj = pg(npi)%mass
! Particles with "fix" movement, but not still.
         dvar(:) = pg(npj)%vel(:) - pg(npi)%var(:)        
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
         appo = amassj * PartKernel(1,npartint) * &
                (dvar(1) * rag(1,npartint) + dvar(2) * rag(2,npartint) +       &
                dvar(3) * rag(3,npartint))
         else
            appo = rhoi * PartKernel(1,npartint) * (amassj/rhoj) * &
                   (dvar(1) * rag(1,npartint) + dvar(2) * rag(2,npartint) +    &
                   dvar(3) * rag(3,npartint))
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
         dvdi(1) = dvdi(1) + pesogradj(1) * dvar(1)      !du/dx  
         dvdi(2) = dvdi(2) + pesogradj(2) * dvar(1)      !du/dy  
         dvdi(3) = dvdi(3) + pesogradj(3) * dvar(1)      !du/dz  
         dvdi(4) = dvdi(4) + pesogradj(1) * dvar(2)      !dv/dx
         dvdi(5) = dvdi(5) + pesogradj(2) * dvar(2)      !dv/dy
         dvdi(6) = dvdi(6) + pesogradj(3) * dvar(2)      !dv/dz
         dvdi(7) = dvdi(7) + pesogradj(1) * dvar(3)      !dw/dx
         dvdi(8) = dvdi(8) + pesogradj(2) * dvar(3)      !dw/dy
         dvdi(9) = dvdi(9) + pesogradj(3) * dvar(3)      !dw/dz
      endif 
! End continuity equation
   enddo
   if ((Granular_flows_options%ID_erosion_criterion==1).or.                    &
      (Granular_flows_options%ID_erosion_criterion==3)) then
! Wrong renormalization: start
! This renormalization presents several errors; it should be re-written)
!   bij(1) = (aij(5) * aij(9) - aij(6) * aij(8))
!   bij(2) =-(aij(4) * aij(9) - aij(7) * aij(6))
!   bij(3) = (aij(4) * aij(8) - aij(5) * aij(7))
!   bij(4) =-(aij(2) * aij(9) - aij(3) * aij(8))
!   bij(5) = (aij(1) * aij(9) - aij(7) * aij(7))
!   bij(6) =-(aij(1) * aij(8) - aij(2) * aij(2))
!   bij(7) = (aij(2) * aij(6) - aij(3) * aij(5))
!   bij(8) =-(aij(1) * aij(6) - aij(3) * aij(4))
!   bij(9) = (aij(1) * aij(5) - aij(2) * aij(4))
!   det = aij(1)*aij(5)*aij(9)+aij(2)*aij(6)*aij(7)+aij(3)*aij(4)*aij(8) - &
!         aij(1)*aij(6)*aij(8)-aij(4)*aij(2)*aij(9)-aij(7)*aij(5)*aij(3)
!   if (det < 0.001d0) det = 0.001d0
!dudx
!   dvdi(1) = (bij(1)*derspa(1) - bij(2)*derspa(2) + bij(3)*derspa(3)) / det     
!dudy
!   dvdi(2) = (bij(1)*derspa(4) - bij(2)*derspa(5) + bij(3)*derspa(6)) / det     
!dudz
!   dvdi(3) = (bij(1)*derspa(7) - bij(2)*derspa(8) + bij(3)*derspa(9)) / det     
!dvdx
!   dvdi(4) = (bij(4)*derspa(1) - bij(5)*derspa(2) + bij(6)*derspa(3)) / det     
!dvdy
!   dvdi(5) = (bij(4)*derspa(4) - bij(5)*derspa(5) + bij(6)*derspa(6)) / det     
!dvdz
!   dvdi(6) = (bij(4)*derspa(7) - bij(5)*derspa(8) + bij(6)*derspa(9)) / det     
!dwdx
!   dvdi(7) = (bij(7)*derspa(1) - bij(8)*derspa(2) + bij(9)*derspa(3)) / det     
!dwdy
!   dvdi(8) = (bij(7)*derspa(4) - bij(8)*derspa(5) + bij(9)*derspa(6)) / det     
!dwdz
!   dvdi(9) = (bij(7)*derspa(7) - bij(8)*derspa(8) + bij(9)*derspa(9)) / det     
! modifica solo dei termini utili per tensore vel def
!   dvdi(2) = half * (dvdi(2) + dvdi(4))
!   dvdi(3) = half * (dvdi(3) + dvdi(7))
!   dvdi(6) = half * (dvdi(6) + dvdi(8))
!   moddia  = (dvdi(1)*dvdi(1) + dvdi(5)*dvdi(5) + dvdi(9)*dvdi(9))
!   modout  = (dvdi(2)*dvdi(2) + dvdi(3)*dvdi(3) + dvdi(6)*dvdi(6))
!   pg(npi)%secinv = Dsqrt( half*moddia + modout)
! Wrong renormalization: end
      moddia  = (dvdi(1) * dvdi(1) + dvdi(5) * dvdi(5) + dvdi(9) * dvdi(9))
      modout  = ((dvdi(2) + dvdi(4)) * (dvdi(2) + dvdi(4)) + (dvdi(3) +        &
                dvdi(7)) * (dvdi(3) + dvdi(7)) + (dvdi(6)+dvdi(8))*(dvdi(6)+   &
                dvdi(8)) )
      pg(npi)%secinv = Dsqrt( half * moddia + quarter * modout)
   endif   
enddo
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
end subroutine inter_EqCont_3D

