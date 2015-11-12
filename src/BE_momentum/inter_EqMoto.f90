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
! Program unit: inter_EqMoto
! Description: Computation of the momentum equation RHS (with DB-SPH boundary treatment scheme, Shepard's coefficient and gravity 
!              are added at a later stage) and the energy equation RHS (this last equation is not validated).   
!----------------------------------------------------------------------------------------------------------------------------------

subroutine inter_EqMoto (npi,tpres,tdiss,tvisc)
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
integer(4),intent(IN) :: npi
double precision,intent(INOUT),dimension(1:SPACEDIM) :: tpres,tdiss,tvisc
integer(4) :: npj,contj,npartint,index_rij_su_h
double precision :: rhoi,rhoj,amassj,pi,pj,alpha,veln,velti,veltj,deltan,pre   
double precision :: coeff,secinv,nupa,nu,modderveln,moddervelt,moddervel
double precision :: dvtdn,denorm,rij_su_h,ke_coef,kacl_coef,rij_su_h_quad
double precision :: vol_Shep,Ww_Shep,rijtemp,rijtemp2
double precision :: gradmod,gradmodwacl,wu,denom
double precision,dimension(3) :: dervel,dervelmorr,appopres,appodiss,rvw
double precision,dimension(3) :: rvwalfa,rvwbeta,ragtemp
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
pg(npi)%dEdT = zero
tpres(:) = zero
tdiss(:) = zero
tvisc(:) = zero
dervel(:) = zero
dervelmorr(:) = zero
deltan = 1.d+07
ke_coef = Domain%coefke / Domain%h
kacl_coef = Domain%coefkacl / Domain%h
!------------------------
! Statements
!------------------------
! Loop to find the distance from the fix wall (if present; SA-SPH). 
! This "fix boundary" is composed by particles whose movement is imposed.
if (Domain%Slip) then
   do contj=1,nPartIntorno(npi)
      npartint = (npi - 1)* NMAXPARTJ + contj
      npj = PartIntorno(npartint)
      if (pg(npj)%vel_type=="std") cycle
      denorm = max(zero,(rag(1,npartint) * pg(npj)%zer(1) + rag(2,npartint)    &
               * pg(npj)%zer(2) + rag(3,npartint) * pg(npj)%zer(3)))
      deltan = min(deltan,denorm)
   end do
end if 
do contj=1,nPartIntorno(npi)
   npartint = (npi - 1)* NMAXPARTJ + contj
   npj = PartIntorno(npartint)
! To check if this term is needed in the previous loop 
   if (npi==npj) cycle  
! Kernel computations only for intermediate time stages
   if (Domain%time_stage>1) then
! To compute inter-particle distance vector 
      ragtemp(1:3) = pg(npi)%coord(1:3) - pg(npj)%coord(1:3)
      if ((abs(ragtemp(1))>doubleh).or.(abs(ragtemp(2))>doubleh).or.           &
         (abs(ragtemp(3))>doubleh)) cycle
! Square distance is preferred to improve accuracy 
      rijtemp = ragtemp(1) * ragtemp(1) + ragtemp(2) * ragtemp(2) + ragtemp(3) &
                *ragtemp(3)
      if (rijtemp>doublesquareh) cycle
! Saving inter-particle distance 
      rijtemp2 = rijtemp
      rijtemp = Dsqrt(rijtemp)
      rij_su_h = rijtemp / Domain%h
      rij_su_h_quad = rijtemp2 / squareh
      index_rij_su_h = int(rij_su_h)
      denom = one / (rijtemp + eta)
      rag(1:3,npartint) = ragtemp(1:3)
! Kernel computations
      gradmod = zero
      gradmodwacl = zero
      wu = zero
      PartKernel(1:4,npartint) = zero
      if (index_rij_su_h>=2) cycle
      gradmod = -two * rij_su_h + 1.5d0 * rij_su_h_quad
      gradmodwacl = -12.0d0 - 3.0d0 * rij_su_h_quad + 12.0d0 * rij_su_h 
      wu = 0.666666667d0 + rij_su_h_quad * (rij_su_h * half - one)
      if (index_rij_su_h>0) then
         gradmod = -gradmod + rij_su_h_quad - two
         wu = (two - rij_su_h) * (two - rij_su_h) * (two - rij_su_h) *         &
         0.166666666667d0
      end if 
      gradmod = gradmod * ke_coef
      gradmodwacl = gradmodwacl * kacl_coef
      PartKernel(1,npartint) = gradmod * denom 
      PartKernel(2,npartint) = PartKernel(1,npartint) / (rijtemp2 + eta2)
      PartKernel(3,npartint) = gradmodwacl * denom
      PartKernel(4,npartint) = wu * Domain%coefke
   end if
   rhoi = pg(npi)%dens
   rhoj = pg(npj)%dens
   amassj = pg(npj)%mass
   pi = pg(npi)%pres
   pj = pg(npj)%pres
   dervel(:) = pg(npj)%vel(:) - pg(npi)%vel(:)
   dervelmorr(:) = pg(npj)%vel(:) - pg(npi)%vel(:)
   if (pg(npj)%vel_type/="std") then  
      pj = pi
      rhoj = rhoi
      amassj = pg(npi)%mass
      moddervel = - two* (pg(npi)%vel(1) * pg(npj)%zer(1) + pg(npi)%vel(2) *   &
                  pg(npj)%zer(2) + pg(npi)%vel(3) * pg(npj)%zer(3))
      dervel(:) = moddervel * pg(npj)%zer(:)   
      if (pg(npj)%slip=="f") then
         dervelmorr(:) = dervel(:)
         else if (pg(npj)%slip=="c") then
            denorm = max(zero,abs(rag(1,npartint) * pg(npj)%zer(1) +           &
                     rag(3,npartint) * pg(npj)%zer(3)))
            veln = (pg(npi)%vel(1) * pg(npj)%zer(1) + pg(npi)%vel(3) *         &
                   pg(npj)%zer(3))
            velti = (pg(npi)%vel(1) * pg(npj)%zer(3) - pg(npi)%vel(3) *        &
                    pg(npj)%zer(1))
            secinv = abs(velti / (deltan + 0.0001d0))
            nu = Med(pg(npi)%imed)%visc
            if (index(Med(pg(npi)%imed)%tipo,"liquid")>0) then
               nu = Med(pg(npi)%imed)%visc
               else if (index(Med(pg(npi)%imed)%tipo,"gas")>0) then
                  nu = Med(pg(npi)%imed)%visc
                  else if (index(Med(pg(npi)%imed)%tipo,"general")>0) then
                     nupa = Med(pg(npi)%imed)%taucri / (secinv + 0.0001d0) +   &
                            Med(pg(npi)%imed)%visc * ((secinv + 0.0001d0) **   &
                            (Med(pg(npi)%imed)%cuin-one))
                     nu = min(Med(pg(npi)%imed)%numx,nupa)
                     else if (index(Med(pg(npi)%imed)%tipo,"granular")>0) then
                        pre = (max(zero,pg(npi)%pres)) / pg(npi)%dens
                        coeff = sin (Med(pg(npi)%imed)%phi)
                        nupa = (pre*coeff) / (secinv + 0.0001d0) +             &
                               Med(pg(npi)%imed)%visc
                        nu = min(nupa,Med(pg(npi)%imed)%numx)
            end if
            dvtdn = (sin(pg(npj)%ang)) * (pg(npi)%dudy + pg(npi)%dvdx) +       &
                    (cos(pg(npj)%ang)) * (pg(npi)%dudx - pg(npi)%dvdy)
            veltj = ( - two * (deltan / (denorm + 0.0001d0)) * nu /            &
                    (pg(npi)%visc + 0.0001d0) + one) * velti + two * dvtdn     &
                    * deltan
            moddervelt = veltj - velti
            modderveln = - two * veln               
            dervelmorr(1) = modderveln * pg(npj)%zer(1) + moddervelt *         &
                            pg(npj)%zer(3)
            dervelmorr(3) = modderveln * pg(npj)%zer(3) - moddervelt *         &
                            pg(npj)%zer(1)
            else if (pg(npj)%slip=="n") then
               dervelmorr(:) = - pg(npi)%vel(:)
      end if
   end if
! Momentum equation: start 
   if (Granular_flows_options%ID_erosion_criterion/=1) then
      alpha = pi / (rhoi * rhoi) + pj / (rhoj * rhoj) 
      else
         alpha = (pi + pj) / rhoi 
   endif
   if (Domain%tipo=="semi") then
      if (Granular_flows_options%ID_erosion_criterion/=1) then
         appopres(:) = ( - amassj * alpha * rag(:,npartint) *                  &
                       PartKernel(3,npartint) )
         else
            appopres(:) = ( - amassj / rhoj * alpha * rag(:,npartint) *        &
                          PartKernel(3,npartint) )
      endif
      else
! DB-SPH: the equation has to use the same kernel type when using the kernel 
! function (cubic spline) and its derivative
         if (Domain%tipo=="bsph") then
            appopres(:) = ( - amassj * alpha * rag(:,npartint) *               &
                          PartKernel(1,npartint) ) 
         endif
   endif
   tpres(:) = tpres(:) + appopres(:)
   call viscomon (npi,npj,npartint,dervel,rvwalfa,rvwbeta)
   appodiss(:) = rvwalfa(:) + rvwbeta(:)
! Monaghan term (artificial viscosity)
   tdiss(:) = tdiss(:) + appodiss(:)   
   call viscomorris (npi,npj,npartint,dervel,rvw)
   tvisc(:) = tvisc(:) + rvw(:)
! Momentum equation: end
   if (esplosione) &
      pg(npi)%dEdT = pg(npi)%dEdT + half * (dervel(1) * (appopres(1) +         &
                     appodiss(1)) + dervel(2) * (appopres(2) + appodiss(2)) +  &
                     dervel(3) * (appopres(3) + appodiss(3)))
end do
! Boundary contributions (DB-SPH)
if ((Domain%tipo=="bsph").and.(DBSPH%n_w > 0)) then
   do contj=1,nPartIntorno_fw(npi)
      npartint = (npi - 1)* NMAXPARTJ + contj
      npj = PartIntorno_fw(npartint)
      dervel(:) = pg_w(npj)%vel(:) - pg(npi)%vel(:)
      appopres(:) = - pg_w(npj)%dens * pg_w(npj)%weight * kernel_fw(1,npartint)&
                    * (pg(npi)%pres / (pg(npi)%dens * pg(npi)%dens) +          &
                    pg_w(npj)%pres / (pg_w(npj)%dens * pg_w(npj)%dens))        &
                    * pg_w(npj)%normal(:)
      appopres(:) = appopres(:) - pg_w(npj)%mass * (pg(npi)%pres /             &
                    (pg(npi)%dens * pg(npi)%dens) + pg_w(npj)%pres /           &
                    (pg_w(npj)%dens * pg_w(npj)%dens)) * rag_fw(:,npartint) *  &
                    kernel_fw(2,npartint)  
      tpres(:) = tpres(:) + appopres(:)
      call viscomon_wall_elements(npi,npj,npartint,dervel,rvwalfa,rvwbeta)
      call viscomorris_wall_elements(npi,npj,npartint,dervel,rvw)
      tvisc(:) = tvisc(:) + rvw(:)
   end do
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine inter_EqMoto

