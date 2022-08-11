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
! Program unit: AddBoundaryContributions_to_ME2D                                
! Description: To compute boundary terms for the 2D momentum equation 
!              (gradPsuro,ViscoF). Equations refer to particle npi. In case of 
!              a neighbouring inlet section, the particle velocity is assigned 
!              (Di Monaco et al., 2011, EACFM).
!              Inversion of the renormalization matrix in 2D, even in the 
!              absence of SASPH neighbours (to fasten the algorithm under 
!              general conditions).                    
!-------------------------------------------------------------------------------
#ifdef SPACE_2D
subroutine AddBoundaryContributions_to_ME2D(npi,IntNcbs,tpres,tdiss,tvisc,     &
   slip_coeff_counter)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use SA_SPH_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
double precision,parameter :: pibmink = 10.d0
double precision,parameter :: eps = 0.1d0
integer(4),intent(in) :: npi,IntNcbs
double precision,intent(inout),dimension(1:SPACEDIM) :: tpres,tdiss,tvisc
double precision,intent(inout),dimension(1:size(Partz)) :: slip_coeff_counter
integer(4) :: imed,i,j,icbs,ibdt,ibdp,iside,sidestr
double precision :: IntWds,IntWdV,cinvisci,Monvisc,roi,ro0,celeri,pressi,u_t_0     
double precision :: alfaMon,xpi,ypi,SVforce,slip_coefficient,TN,ypimin,ypigradN 
double precision :: GradNpsuro,IntWd1s0,IntWd1s2,IntWd3s0,DvelN,viscN  
double precision :: GravN,PressB,distpi,distpimin,pressib,pressibmin
double precision :: QiiIntWdS,level,pressj,Qsi,Qsj,velix,veliz,veliq,hcrit
double precision :: hcritmin,zbottom,FlowRate1,Lb,L,minquotanode,maxquotanode
double precision :: SomQsiQsj,DiffQsiQsj,d_50
integer(4),dimension(1:PLANEDIM) :: acix
double precision,dimension(1:PLANEDIM) :: RHS,RG,ss,nnlocal,gradbPsuro
double precision,dimension(1:PLANEDIM) :: ViscoMon,ViscoShear,sidevel,TT,Dvel
double precision,dimension(1:PLANEDIM) :: RG_like,gradbPsuro_like,RG_sum
double precision,dimension(1:SPACEDIM) :: u_t_0_vector,aux_vec_2,aux_vec
! Unit vector of the unity vector: direction of (1,1)
double precision,dimension(1:SPACEDIM) :: one_vec_dir
type (TyBoundarySide) :: RifBoundarySide
character(4):: strtype
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine wall_function_for_SASPH(u_t_0,d_50,r_0w,slip_coefficient_0w,     &
      ni_T_0w)
      implicit none
         double precision,intent(in) :: u_t_0,d_50,r_0w
         double precision,intent(out) :: slip_coefficient_0w,ni_T_0w
   end subroutine wall_function_for_SASPH
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
! active coordinate indexes
acix(1) = 1  
acix(2) = 3
if ((Domain%time_stage==1).or.(Domain%time_split==1)) then 
   pg(npi)%kodvel = 0
   pg(npi)%velass(:) = zero
endif
imed = pg(npi)%imed
ro0 = Med(imed)%den0
roi = pg(npi)%dens
pressi = pg(npi)%pres
distpimin = Domain%dx
RHS(1) = -tpres(1)
RHS(2) = -tpres(3)
ViscoMon(1) = zero
ViscoMon(2) = zero
ViscoShear(1) = zero
ViscoShear(2) = zero
ibdt = BoundaryDataPointer(3,npi)
one_vec_dir(1:3) = 0.d0
do i=1,PLANEDIM
   one_vec_dir(acix(i)) = 1.d0 / dsqrt(2.d0)
enddo
RG_sum(1:2) = 0.d0
!------------------------
! Statements
!------------------------
do icbs=1,IntNcbs
   ibdp = ibdt + icbs - 1
   iside = BoundaryDataTab(ibdp)%CloBoNum
   RifBoundarySide = BoundarySide(iside)
   sidestr = RifBoundarySide%stretch
   strtype = Tratto(sidestr)%tipo
   if (strtype=="open") cycle
   xpi = BoundaryDataTab(ibdp)%LocXYZ(1)
   ypi = BoundaryDataTab(ibdp)%LocXYZ(2)
   IntWdS = BoundaryDataTab(ibdp)%BoundaryIntegral(1)
   IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(3)
   IntWd1s0 = BoundaryDataTab(ibdp)%BoundaryIntegral(6)
   IntWd3s0 = BoundaryDataTab(ibdp)%BoundaryIntegral(7)
   IntWd1s2 = BoundaryDataTab(ibdp)%BoundaryIntegral(8)
   GravN = zero
   do i=1,PLANEDIM
      ss(i) = RifBoundarySide%T(acix(i),acix(1))
      nnlocal(i) = RifBoundarySide%T(acix(i),acix(2))
      gradbPsuro(i) = Domain%grav(acix(i))
! For the renormalization at SASPH frontiers
      if (input_any_t%ME_gradp_cons==3) then
         gradbPsuro_like(i) = one_vec_dir(acix(i))
      endif
      GravN = GravN + Domain%grav(acix(i)) * nnlocal(i)
   enddo
   select case (strtype)
      case ("fixe")
         do i=1,PLANEDIM
            sidevel(i) = zero
         enddo
      case ("tapi")
         do i=1,PLANEDIM
            sidevel(i) = RifBoundarySide%velocity(acix(i))
         enddo
   endselect
   if ((strtype=="fixe").or.(strtype=="tapi"))  then
      do i=1,PLANEDIM
         Dvel(i) = two * (pg(npi)%var(acix(i)) - sidevel(i))   
      enddo
      DvelN = Dvel(1) * nnlocal(1) + Dvel(2) * nnlocal(2)
! SASPH "PPST term": start
      distpi = max(ypi,distpimin)
      pressib = pressi
      pressibmin = pibmink * ro0 * distpimin
      if (DvelN<=zero) then
         if (pressib<pressibmin) pressib = pressibmin
      endif
      PressB = pressib - ro0 * GravN * distpi
      QiiIntWdS = IntWdS * (pressib + PressB) / roi
! SASPH "PPST term": end
! Sub-critical flow
      elseif (strtype=="leve") then   
         Qsi = pressi / roi
! The BC_shear_stress_input (formerly Shear_coeff) was used before SPHERA v.7.0 
! for other purposes in case of boundaries of type "leve"
         level = Partz(Tratto(sidestr)%zone)%BC_shear_stress_input
         if (pg(npi)%coord(3)<level) then    
            pressj = Med(imed)%den0 * Domain%grav(3) * (pg(npi)%coord(3) -     &
                     level)
            Qsj = pressj / Med(imed)%den0  
            SomQsiQsj = Qsi + Qsj
            DiffQsiQsj = Qsi - Qsj       
            else            
               SomQsiQsj = zero
               DiffQsiQsj = zero    
         endif         
         QiiIntWdS = SomQsiQsj * IntWds       
         ypimin = eps * Domain%h               
         ypigradN = ypi          
         if (ypigradN<ypimin) ypigradN = ypimin     
         GradNpsuro = (Qsi - Qsj) / ypigradN       
         gradbPsuro(1) = gradbPsuro(1) + GradNpsuro * nnlocal(1)    
         gradbPsuro(2) = gradbPsuro(2) + GradNpsuro * nnlocal(2)     
         pg(npi)%velass(:) = zero        
! Imposed normal velocity
         elseif (strtype=="velo") then     
            if ((Domain%time_stage==1).or.(Domain%time_split==1)) then 
               pg(npi)%kodvel = 2
               pg(npi)%velass(1) = Tratto(sidestr)%NormVelocity * nnlocal(1)
               pg(npi)%velass(2) = zero            
               pg(npi)%velass(3) = Tratto(sidestr)%NormVelocity * nnlocal(2)
            endif
            return             
            elseif (strtype=="flow") then     
! Imposing the flow rate 
               if ((Domain%time_stage==1).or.(Domain%time_split==1)) then 
                  pg(npi)%kodvel = 2
! Normal velocity 
                  if (BoundarySide(iside)%CloseParticles>0) then
                     minquotanode = min(vertice(3,                             &
                        BoundarySide(iside)%vertex(1)),                        &
                        vertice(3,BoundarySide(iside)%vertex(2)))
                     maxquotanode = max                                        &
                        (vertice(3,BoundarySide(iside)%vertex(1)),             &
                        vertice(3,BoundarySide(iside)%vertex(2)))
                     Lb = BoundarySide(iside)%CloseParticles_maxQuota -        &
                        minquotanode
                     L = maxquotanode - minquotanode
                     FlowRate1 = Tratto(sidestr)%FlowRate * Lb / L
                     Tratto(sidestr)%NormVelocity = doubleh * FlowRate1 /      &
                        (BoundarySide(iside)%CloseParticles * Domain%PVolume)
                     else
                        Tratto(sidestr)%NormVelocity = zero
                  endif
                  pg(npi)%velass(1) = Tratto(sidestr)%NormVelocity * nnlocal(1)
                  pg(npi)%velass(2) = zero            
                  pg(npi)%velass(3) = Tratto(sidestr)%NormVelocity * nnlocal(2)
               endif
               return
               elseif (strtype=="sour") then
                  if (xpi>=zero.and.xpi<=RifBoundarySide%length) then
                     if ((Domain%time_stage==1).or.(Domain%time_split==1)) then 
                        pg(npi)%kodvel = 2
                        pg(npi)%velass(1) = Tratto(sidestr)%NormVelocity *     &
                                            nnlocal(1)
                        pg(npi)%velass(2) = zero            
                        pg(npi)%velass(3) = Tratto(sidestr)%NormVelocity *     &
                                            nnlocal(2)
                     endif
                     return
                  endif
! Critical flow
                  elseif (strtype=="crit") then     
! Non-stationary critical flow
                     Qsi = pressi / roi                          
                     velix = pg(npi)%vel(1)
                     veliz = pg(npi)%vel(3)         
                     veliq = velix * velix + veliz * veliz  
                     hcrit = veliq / dabs(Domain%grav(3))
! The BC_shear_stress_input (formerly Shear_coeff) was used before SPHERA v.7.0 
! for other purposes in case of boundaries of type "leve"
                     hcritmin =                                                &
                        Partz(Tratto(sidestr)%zone)%BC_shear_stress_input       
                     if (hcritmin<Domain%h) hcritmin = Domain%h
                     if (hcrit<hcritmin) hcrit = hcritmin
                     zbottom = Vertice(3,RifBoundarySide%vertex(1))  
                     if (zbottom>Vertice(3,RifBoundarySide%vertex(2)))         &
                        zbottom = Vertice(3,RifBoundarySide%vertex(2))
                     level = zbottom + hcrit          
                     if (pg(npi)%coord(3)<level) then          
                        pressj = Med(imed)%den0 * Domain%grav(3) *             &
                           (pg(npi)%coord(3) - level)
                        Qsj = pressj / Med(imed)%den0
                        SomQsiQsj = Qsi + Qsj              
                        DiffQsiQsj = Qsi - Qsj             
                        else                                       
                           SomQsiQsj = zero            
                           DiffQsiQsj = zero
                     endif
                     QiiIntWdS = SomQsiQsj * IntWds              
                     ypimin = eps * Domain%h                      
                     ypigradN = ypi                                    
                     if (ypigradN<ypimin) ypigradN = ypimin           
                     GradNpsuro = DiffQsiQsj / ypigradN                    
                     gradbPsuro(1) = gradbPsuro(1) + GradNpsuro * nnlocal(1)           
                     gradbPsuro(2) = gradbPsuro(2) + GradNpsuro * nnlocal(2)              
! Zeroing the vertical component of velocity 
                     if ((Domain%time_stage==1).or.(Domain%time_split==1)) then 
                        pg(npi)%kodvel = 1                                                 
                        pg(npi)%velass(:) = zero
                     endif             
                     elseif (strtype=="open") then       
                        cycle
   endif
! SASPH contribution to "grad_p" and renormalization matrix: start
   RG = zero
   RG_like(1:2) = 0.d0
   do i=1,PLANEDIM
      do j=1,PLANEDIM
         RG(i) = RG(i) + RifBoundarySide%RN(acix(i),acix(j)) * gradbPsuro(j)
         if (input_any_t%ME_gradp_cons==3) then          
! Renormalization at SASPH frontiers
            RG_like(i) = RG_like(i) + RifBoundarySide%RN(acix(i),acix(j)) *    &
                         gradbPsuro_like(j)
         endif
      enddo
      RG(i) = RG(i) * IntWdV
      if (input_any_t%ME_gradp_cons==3) then
! Renormalization at SASPH frontiers
         RG_like(i) = RG_like(i) * IntWdV
      endif
   enddo
! Collecting the temporary SASPH contributions to the grad_p term
   RG_sum(1:2) = RG_sum(1:2) + RG(1:2)
   if (input_any_t%ME_gradp_cons==3) then
! Renormalization at SASPH frontiers
! Notice that the contributions to RHS and B_ren_fp have different signs as RHS 
! will be subtracted from the acceleration
      do i=1,PLANEDIM
            pg(npi)%B_ren_fp(acix(i),1:3) = pg(npi)%B_ren_fp(acix(i),1:3) -    &
                                            RG_like(i)
      enddo
   endif
! SASPH contribution to "grad_p" and renormalization matrix: end
! Collection of the SASPH PPST terms directly in the RHS vector
   RHS(1:2) = RHS(1:2) - nnlocal(1:2) * QiiIntWdS
! Volume viscosity force (with changed sign) and artificial viscosity term
   if (strtype=="fixe".or.strtype=="tapi") then
      if (xpi>=zero.and.xpi<=RifBoundarySide%length) then
         if (Partz(Tratto(sidestr)%zone)%slip_coefficient_mode==1) then
! Slip coefficient and molecular viscosity from input
            slip_coefficient = Partz(Tratto(sidestr)%zone)%BC_shear_stress_input
            cinvisci = pg(npi)%kin_visc
            elseif (Partz(Tratto(sidestr)%zone)%slip_coefficient_mode==2) then
! Slip coefficient computed
! Particle tangential (relative) velocity (vector)
! Both "DvelN" and "T" are defined with an opposite direction
               u_t_0_vector(:) = (pg(npi)%var(:) - RifBoundarySide%velocity(:))&
                                 - (0.5d0 * DvelN * RifBoundarySide%T(:,3))
! Particle tangential (relative) velocity (absolute value)
               u_t_0 = dsqrt(dot_product(u_t_0_vector(:),u_t_0_vector(:)))
! To assess the slip coefficient and the turbulent viscosity
               d_50 = Partz(Tratto(sidestr)%zone)%BC_shear_stress_input * 10.d0
               call wall_function_for_SASPH(u_t_0,d_50,                        &
                  BoundaryDataTab(ibdp)%LocXYZ(2),slip_coefficient,cinvisci)
               if (slip_coefficient>1.d-12) then
!$omp critical (avg_slip_coefficient_2D)         
! Update of the incremental sum for the slip coefficient
                  Partz(Tratto(sidestr)%zone)%avg_comp_slip_coeff =            &
                     Partz(Tratto(sidestr)%zone)%avg_comp_slip_coeff +         &
                     slip_coefficient
! Update of the incremental sum for the turbulent viscosity
                  Partz(Tratto(sidestr)%zone)%avg_ni_T_SASPH =                 &
                     Partz(Tratto(sidestr)%zone)%avg_ni_T_SASPH + cinvisci
! Update of the incremental sum for the wall-function shear stress
                  Partz(Tratto(sidestr)%zone)%avg_tau_wall_f =                 &
                     Partz(Tratto(sidestr)%zone)%avg_tau_wall_f +              &
                     slip_coefficient * pg(npi)%dens * cinvisci * u_t_0 /      &
                     BoundaryDataTab(ibdp)%LocXYZ(2)
! Update the counter for both the slip coefficient, the turbulent viscosity and 
! the wall-function shear stress
                  slip_coeff_counter(Tratto(sidestr)%zone) =                   &
                     slip_coeff_counter(Tratto(sidestr)%zone) + 1
!$omp end critical (avg_slip_coefficient_2D)
               endif
         endif
         celeri = Med(imed)%celerita
         alfaMon = Med(imed)%alfaMon
         Monvisc = alfaMon * celeri * Domain%h
         viscN = Monvisc
         TN = viscN * DvelN * IntWd3s0
         do i=1,PLANEDIM
            TT(i) = TN * nnlocal(i)
         enddo
         if ((pg(npi)%laminar_flag==1).or.                                     &
            (Tratto(sidestr)%laminar_no_slip_check.eqv..false.)) then
            SVforce = 2.d0 * cinvisci * slip_coefficient * IntWd1s0
            else
               SVforce = 0.d0
         endif
         do i=1,PLANEDIM
            ViscoMon(i) = ViscoMon(i) + TT(i)
            ViscoShear(i) = ViscoShear(i) + Dvel(i) * SVforce
         enddo
      endif
   endif
enddo
! The renormalization matrix for the pressure-gradient term is inverted 
! just after all its components are collected and just before the 
! 1st-order consistency scheme applies to the summation of all the 
! particle-boundary contributions
if (input_any_t%ME_gradp_cons>0) then
! Inversion of the renormalization matrix
   call B_ren_gradp_inversion(npi)
endif
if (IntNcbs==0) return
! grad_p (renormalization at boundaries): start
if (input_any_t%ME_gradp_cons==3) then
   aux_vec_2(1:3) = 0.d0
   do i=1,PLANEDIM
      aux_vec_2(acix(i)) = RG_sum(i)
   enddo
   call MatrixProduct(pg(npi)%B_ren_fp,BB=aux_vec_2,CC=aux_vec,nr=3,nrc=3,     &
      nc=1)
   do i=1,PLANEDIM
      RG_sum(i) = -aux_vec(acix(i))
   enddo
endif
! grad_p (renormalization at boundaries): end
! SASPH contributions to the grad_p term are collected in RHS (both with and 
! without renormalization)
RHS(1:2) = RHS(1:2) + RG_sum(1:2)
do i=1,PLANEDIM
   tpres(acix(i)) = - RHS(i)
   tdiss(acix(i)) = tdiss(acix(i)) - ViscoMon(i)
   tvisc(acix(i)) = tvisc(acix(i)) - ViscoShear(i)
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine AddBoundaryContributions_to_ME2D
#endif
