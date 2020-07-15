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
! Program unit: AddBoundaryContributions_to_ME2D                                
! Description: To compute boundary terms for the 2D momentum equation 
!              (gradPsuro,ViscoF). Equations refer to particle npi. (Di Monaco 
!              et al., 2011, EACFM).                        
!-------------------------------------------------------------------------------
#ifdef SPACE_2D
subroutine AddBoundaryContributions_to_ME2D(npi,IntNcbs,tpres,tdiss,tvisc) 
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
integer(4) :: imed,i,j,icbs,ibdt,ibdp,iside,sidestr
double precision :: IntWds,IntWdV,cinvisci,Monvisc,roi,ro0,celeri,pressi,u_t_0     
double precision :: alfaMon,xpi,ypi,SVforce,slip_coefficient,TN,ypimin,ypigradN 
double precision :: GradNpsuro,IntWd1s0,IntWd1s2,IntWd3s0,DvelN,viscN  
double precision :: GravN,PressB,distpi,distpimin,pressib,pressibmin
double precision :: QiiIntWdS,level,pressj,Qsi,Qsj,velix,veliz,veliq,hcrit
double precision :: hcritmin,zbottom,FlowRate1,Lb,L,minquotanode,maxquotanode
double precision :: SomQsiQsj,DiffQsiQsj
integer(4),dimension(1:PLANEDIM) :: acix
double precision,dimension(1:PLANEDIM) :: IntLocXY,RHS,RG,ss,nnlocal,gradbPsuro
double precision,dimension(1:PLANEDIM) :: ViscoMon,ViscoShear,sidevel,TT,Dvel
double precision,dimension(1:SPACEDIM) :: u_t_0_vector
type (TyBoundarySide) :: RifBoundarySide
character(4):: strtype
double precision,dimension(1:size(Partz)) :: slip_coeff_counter
!------------------------
! Explicit interfaces
!------------------------
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
slip_coeff_counter(:) = 0
!------------------------
! Statements
!------------------------
do icbs=1,IntNcbs
   ibdp = ibdt + icbs - 1
   IntLocXY(1:PLANEDIM) = BoundaryDataTab(ibdp)%LocXYZ(1:PLANEDIM)
   iside = BoundaryDataTab(ibdp)%CloBoNum
   RifBoundarySide = BoundarySide(iside)
   sidestr = RifBoundarySide%stretch
   strtype = Tratto(sidestr)%tipo
   if (strtype=="open") cycle
   xpi = IntLocXY(1)
   ypi = IntLocXY(2)
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
! Boundary contribution to "gradP" 
! (pressure gradient term in the momentum equation) 
   if ((strtype=="fixe").or.(strtype=="tapi"))  then
      do i=1,PLANEDIM
         Dvel(i) = two * (pg(npi)%var(acix(i)) - sidevel(i))   
      enddo
      DvelN = Dvel(1) * nnlocal(1) + Dvel(2) * nnlocal(2)
      distpi = max(ypi,distpimin)
      pressib = pressi
      pressibmin = pibmink * ro0 * distpimin
      if (DvelN<=zero) then
         if (pressib<pressibmin) pressib = pressibmin
      endif
      PressB = pressib - ro0 * GravN * distpi 
      QiiIntWdS = IntWdS * (pressib + PressB) / roi
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
                     hcrit = veliq / abs(Domain%grav(3))
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
   RG = zero
   do i=1,PLANEDIM
      do j=1,PLANEDIM
         RG(i) = RG(i) + RifBoundarySide%RN(acix(i),acix(j)) * gradbPsuro(j) 
      enddo
      RG(i) = RG(i) * IntWdV
   enddo
   do i=1,PLANEDIM
      RHS(i) = RHS(i) - nnlocal(i) * QiiIntWdS + RG(i)
   enddo
! Volume viscosity force (with changed sign) and artificial viscosity term
   if (strtype=="fixe".or.strtype=="tapi") then
      if (xpi>=zero.and.xpi<=RifBoundarySide%length) then
         if (Partz(Tratto(sidestr)%zone)%slip_coefficient_mode==1) then
! Slip coefficient from input
            cinvisci = pg(npi)%kin_visc
            elseif (Partz(Tratto(sidestr)%zone)%slip_coefficient_mode==2) then
! Slip coefficient computed
! Particle tangential (relative) velocity (vector)
               u_t_0_vector(:) = (pg(npi)%var(:) - RifBoundarySide%velocity(:))&
                                 - (0.5d0 * DvelN * RifBoundarySide%T(:,3))
! Particle tangential (relative) velocity (absolute value)
               u_t_0 = dsqrt(dot_product(u_t_0_vector(:),u_t_0_vector(:)))
! To assess the slip coefficient and the turbulent viscosity
               call wall_function_for_SASPH(u_t_0,                             &
                  Partz(Tratto(sidestr)%zone)%BC_shear_stress_input,           &
                  pg(npi)%dens,BoundaryDataTab(ibdp)%LocXYZ(3),                &
                  slip_coefficient,cinvisci)
! Update of the incremental sum of the slip coefficient values
                  Partz(Tratto(sidestr)%zone)%avg_comp_slip_coeff =            &
                     Partz(Tratto(sidestr)%zone)%avg_comp_slip_coeff +         &
                     slip_coefficient
                  slip_coeff_counter(Tratto(sidestr)%zone) =                   &
                     slip_coeff_counter(Tratto(sidestr)%zone) + 1
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
do i=1,PLANEDIM
   tpres(acix(i)) = - RHS(i)
   tdiss(acix(i)) = tdiss(acix(i)) - ViscoMon(i)
   tvisc(acix(i)) = tvisc(acix(i)) - ViscoShear(i)
enddo
! Update of the average slip coefficient for each boundary zone
do i=1,size(Partz)
   if (slip_coeff_counter(i)>0) then
      Partz(i)%avg_comp_slip_coeff = Partz(i)%avg_comp_slip_coeff /            &
                                     slip_coeff_counter(i)
   endif
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine AddBoundaryContributions_to_ME2D
#endif
