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
! Program unit: inter_EqMoto
! Description: Computation of the momentum equation RHS (with DB-SPH boundary 
!              treatment scheme, Shepard's coefficient and gravity are added at 
!              a later stage) merged with the control-volume velocity equation, 
!              included contributions to the ALE1 velocity increment.
!-------------------------------------------------------------------------------
subroutine inter_EqMoto(npi,tpres,tdiss,tvisc)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: npi
double precision,intent(inout),dimension(1:SPACEDIM) :: tpres,tdiss,tvisc
integer(4) :: npj,contj,npartint,index_rij_su_h
double precision :: rhoi,rhoj,amassj,pi,pj,alpha,veln,velti,veltj,deltan,pre   
double precision :: coeff,secinv,nupa,nu,modderveln,moddervelt,moddervel
double precision :: dvtdn,denorm,rij_su_h,ke_coef,kacl_coef,rij_su_h_quad
double precision :: rijtemp,rijtemp2,numx
double precision :: gradmod,gradmodwacl,wu,denom,absv_pres_grav_inner
double precision :: absv_Morris_inner,Morris_inner_weigth,kernel_der
double precision :: dervel(3),dervelmorr(3),appopres(3),appodiss(3),rvw(3)
double precision :: rvwalfa(3),rvwbeta(3),ragtemp(3),rvw_sum(3),rvw_semi_part(3)
double precision :: DBSPH_wall_she_vis_term(3),t_visc_semi_part(3),aux_vec(3)
double precision :: t_pres_aux(3),ALE1_term_sum(3)
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine viscomorris(mass_comput_part,dens_comput_part,                   &
      kin_visc_comput_part,mass_neighbour,dens_neighbour,kin_visc_neighbour,   &
      kernel_der,vel_type,rel_dis,dervel,rvw)
   implicit none
   double precision,intent(in) :: mass_comput_part,dens_comput_part
   double precision,intent(in) :: kin_visc_comput_part,mass_neighbour
   double precision,intent(in) :: dens_neighbour,kin_visc_neighbour,kernel_der
   double precision,intent(in) :: rel_dis(3),dervel(3)
   character(3),intent(in) :: vel_type
   double precision,intent(out) :: rvw(3)
   end subroutine viscomorris
   subroutine MatrixProduct(AA,BB,CC,nr,nrc,nc)
      implicit none
      integer(4),intent(in) :: nr,nrc,nc
      double precision,intent(in),dimension(nr,nrc) :: AA
      double precision,intent(in),dimension(nrc,nc) :: BB
      double precision,intent(inout),dimension(nr,nc) :: CC
   end subroutine MatrixProduct
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
dervel(:) = zero
dervelmorr(:) = zero
deltan = 1.d7
ke_coef = Domain%coefke / Domain%h
kacl_coef = Domain%coefkacl / Domain%h
rvw_sum(:) = zero
DBSPH_wall_she_vis_term(:) = 0.d0
t_visc_semi_part(:) = 0.d0
Morris_inner_weigth = 0.d0
t_pres_aux(1:3) = 0.d0
ALE1_term_sum(1:3) = 0.d0
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
   enddo
endif
do contj=1,nPartIntorno(npi)
   npartint = (npi - 1) * NMAXPARTJ + contj
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
      if (rijtemp>square_doubleh) cycle
! Saving inter-particle distance 
      rijtemp2 = rijtemp
      rijtemp = dsqrt(rijtemp)
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
      gradmodwacl = -12.d0 - 3.d0 * rij_su_h_quad + 12.d0 * rij_su_h 
      wu = 0.666666667d0 + rij_su_h_quad * (rij_su_h * half - one)
      if (index_rij_su_h>0) then
         gradmod = -gradmod + rij_su_h_quad - two
         wu = (two - rij_su_h) * (two - rij_su_h) * (two - rij_su_h) *         &
         0.166666666667d0
      endif 
      gradmod = gradmod * ke_coef
      gradmodwacl = gradmodwacl * kacl_coef
      PartKernel(1,npartint) = gradmod * denom 
      PartKernel(2,npartint) = PartKernel(1,npartint) / (rijtemp2 + eta2)
      PartKernel(3,npartint) = gradmodwacl * denom
      PartKernel(4,npartint) = wu * Domain%coefke
   endif
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
         elseif (pg(npj)%slip=="c") then
            denorm = max(zero,abs(rag(1,npartint) * pg(npj)%zer(1) +           &
                     rag(3,npartint) * pg(npj)%zer(3)))
            veln = (pg(npi)%vel(1) * pg(npj)%zer(1) + pg(npi)%vel(3) *         &
                   pg(npj)%zer(3))
            velti = (pg(npi)%vel(1) * pg(npj)%zer(3) - pg(npi)%vel(3) *        &
                    pg(npj)%zer(1))
            secinv = abs(velti / (deltan + 1.d-4))
            if (index(Med(pg(npi)%imed)%tipo,"liquid")>0) then
               nu = Med(pg(npi)%imed)%kin_visc
               elseif (index(Med(pg(npi)%imed)%tipo,"granular")>0) then
                  pre = (max(zero,pg(npi)%pres)) / pg(npi)%dens
                  coeff = sin (Med(pg(npi)%imed)%phi)
                  nupa = (pre*coeff) / (secinv + 1.d-4) +                      &
                         Med(pg(npi)%imed)%kin_visc
                  numx = Med(pg(npi)%imed)%mumx / Med(pg(npi)%imed)%den0
                  nu = min(nupa,numx)
            endif
            dvtdn = (sin(pg(npj)%ang)) * (pg(npi)%dudy + pg(npi)%dvdx) +       &
                    (cos(pg(npj)%ang)) * (pg(npi)%dudx - pg(npi)%dvdy)
            veltj = ( - two * (deltan / (denorm + 1.d-4)) * nu /               &
                    (pg(npi)%kin_visc + 1.d-4) + one) * velti + two * dvtdn    &
                    * deltan
            moddervelt = veltj - velti
            modderveln = - two * veln               
            dervelmorr(1) = modderveln * pg(npj)%zer(1) + moddervelt *         &
                            pg(npj)%zer(3)
            dervelmorr(3) = modderveln * pg(npj)%zer(3) - moddervelt *         &
                            pg(npj)%zer(1)
            elseif (pg(npj)%slip=="n") then
               dervelmorr(:) = - pg(npi)%vel(:)
      endif
   endif
! Momentum equation: start
! "grad_p + ALE1" terms: start
   if (Domain%tipo=="semi") then
      if (Granular_flows_options%KTGF_config/=1) then
! Liquids
! For the pressure-gradient term
         alpha = (pj - pi) / (rhoi * rhoj)
         appopres(1:3) = -amassj * alpha * rag(1:3,npartint) *                 &
                         PartKernel(3,npartint)
! Contribution to the ALE1 term in the ME-VC
         alpha = (pi * (rhoj / rhoi + 1.d0) + pj * (rhoi / rhoj - 1.d0)) /     &
                 (rhoi * rhoj)
         ALE1_term_sum(1:3) = ALE1_term_sum(1:3) - amassj * alpha *            &
                              rag(1:3,npartint) * PartKernel(3,npartint)
         else
! Dense granular flows
            alpha = (pi + pj) / (rhoi * rhoj)
            appopres(1:3) = -amassj * alpha * rag(1:3,npartint) *              &
                            PartKernel(3,npartint)            
      endif
      else
! DBSPH: the equation has to use the same kernel type when using the kernel 
! function (cubic spline) and its derivative
         if (Domain%tipo=="bsph") then
            alpha = pi / (rhoi * rhoi) + pj / (rhoj * rhoj)
            appopres(1:3) = -amassj * alpha * rag(1:3,npartint) *              &
                            PartKernel(1,npartint)            
         endif
   endif
   t_pres_aux(1:3) = t_pres_aux(1:3) + appopres(1:3)
! "grad_p + ALE1" terms: end
! To compute Monaghan term (artificial viscosity)
   call viscomon(npi,npj,npartint,dervel,rvwalfa,rvwbeta)
   appodiss(:) = rvwalfa(:) + rvwbeta(:)
! To add Monaghan term (artificial viscosity)
   tdiss(:) = tdiss(:) + appodiss(:)
! To compute Morris term (interaction with neighbouring fluid particle)
   if ((pg(npi)%kin_visc>1.d-12).and.(pg(npj)%kin_visc>1.d-12)) then
      call viscomorris(pg(npi)%mass,pg(npi)%dens,pg(npi)%kin_visc,             &
         pg(npj)%mass,pg(npj)%dens,pg(npj)%kin_visc,PartKernel(2,npartint),    &
         pg(npj)%vel_type,rag(1:3,npartint),dervel,rvw)
! To add  Morris term (interaction with neighbouring fluid particle)   
      tvisc(:) = tvisc(:) + rvw(:)
      rvw_sum(:) = rvw_sum(:) + rvw(:)
   endif
! Momentum equation: end
enddo
if (input_any_t%C1_BE) then
! 1st-order consistency
   call MatrixProduct(pg(npi)%B_ren_gradp,BB=t_pres_aux,CC=aux_vec,nr=3,nrc=3, &
      nc=1)
   t_pres_aux(1:3) = -aux_vec(1:3)
endif
! grad_p term and ALE1 term are summed to the acceleration (ALE1 term is here 
! formally pre-added to the grad_p term)
tpres(1:3) = tpres(1:3) + t_pres_aux(1:3) + ALE1_term_sum(1:3)
if (input_any_t%ALE3) then
! ALE1 term (always used) is explicitly saved, printed and used 
! also for the CE only if ALE3==.true. 
   pg(npi)%dvel_ALE1(1:3) = pg(npi)%dvel_ALE1(1:3) + ALE1_term_sum(1:3)
endif
pg(npi)%laminar_flag = 0
if (pg(npi)%kin_visc>0.d0) then
   absv_pres_grav_inner = dsqrt(dot_product(tpres,tpres)) + GI
   absv_Morris_inner = dsqrt(dot_product(rvw_sum,rvw_sum))
   if (absv_pres_grav_inner/=0.) Morris_inner_weigth = absv_Morris_inner /     &
                                                       absv_pres_grav_inner *  &
                                                       100.d0
! Taking into account the condition (absv_pres_grav_inner==0.d0) is only due to 
! graphical reasons.
   if ((Morris_inner_weigth>5.d0).or.(absv_pres_grav_inner==0.d0)) then
      pg(npi)%laminar_flag = 1
      else
         pg(npi)%laminar_flag = 0
   endif
endif
! Boundary contributions (DB-SPH), only in case of a simulated local laminar
! regime or imposed no-slip conditions
if (DBSPH%n_w>0) then
   do contj=1,nPartIntorno_fw(npi)
      npartint = (npi - 1) * NMAXPARTJ + contj
      npj = PartIntorno_fw(npartint)
      dervel(:) = pg_w(npj)%vel(:) - pg(npi)%vel(:)
      appopres(:) = - pg_w(npj)%dens * pg_w(npj)%weight *                      &
                    kernel_fw(1,npartint) * (pg(npi)%pres / (pg(npi)%dens *    &
                    pg(npi)%dens) + pg_w(npj)%pres / (pg_w(npj)%dens *         &
                    pg_w(npj)%dens)) * pg_w(npj)%normal(:)
      appopres(:) = appopres(:) - pg_w(npj)%mass * (pg(npi)%pres /             &
                    (pg(npi)%dens * pg(npi)%dens) + pg_w(npj)%pres /           &
                    (pg_w(npj)%dens * pg_w(npj)%dens)) * rag_fw(:,npartint)*   &
                    kernel_fw(2,npartint)  
      tpres(:) = tpres(:) + appopres(:)
      if ((pg(npi)%laminar_flag==1).or.(DBSPH%slip_ID==1)) then
         call DBSPH_BC_shear_viscosity_term(npi,npj,npartint,                  &
            DBSPH_wall_she_vis_term)
         if ((pg(npi)%kin_visc>1.d-12).and.                                    &
            (pg_w(npj)%kin_visc_semi_part>1.d-12)) then
! To compute Morris term (interaction with neighbouring semi-particle)
            kernel_der = kernel_fw(2,npartint)/(dot_product(rag_fw(:,npartint),&
                         rag_fw(:,npartint)) + eta2)
            call viscomorris(pg(npi)%mass,pg(npi)%dens,pg(npi)%kin_visc,       &
               pg_w(npj)%mass,pg_w(npj)%dens,pg_w(npj)%kin_visc_semi_part,     &
               kernel_der,"std",rag_fw(1:3,npartint),dervel,rvw_semi_part)
            t_visc_semi_part(:) = t_visc_semi_part(:) + rvw_semi_part(:)
         endif
      endif
   enddo
   if ((pg(npi)%laminar_flag==1).or.(DBSPH%slip_ID==1)) then
! Computation of the boundary shear viscosity term in DB-SPH-NS
      DBSPH_wall_she_vis_term(:) = DBSPH_wall_she_vis_term(:) / pg(npi)%dens   &
                                   / pg(npi)%Gamma
   endif
! Update the overall (inner+BC) shear viscosity term in DB-SPH-NS
   tvisc(:) = tvisc(:) + DBSPH_wall_she_vis_term(:) + t_visc_semi_part(:)   
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine inter_EqMoto
