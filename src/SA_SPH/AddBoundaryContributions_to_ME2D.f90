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
! Program unit: AddBoundaryContributions_to_ME2D                                
! Description: To compute boundary terms for the 2D momentum equation (gradPsuro, 
!              ViscoF). Equations refer to particle npi. (Di Monaco et al., 2011, EACFM).                        
!----------------------------------------------------------------------------------------------------------------------------------

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
double precision,parameter :: pibmink = 10.0d0
double precision,parameter :: eps = 0.1d0
integer(4),intent(IN) :: npi,IntNcbs
double precision,intent(INOUT),dimension(1:SPACEDIM) :: tpres,tdiss,tvisc
integer(4) :: imed,i,j,icbs,ibdt,ibdp,iside,sidestr
double precision :: IntWds,IntWdV,cinvisci,Monvisc,roi,ro0,celeri,pressi      
double precision :: alfaMon,xpi,ypi,SVforce,SVcoeff,TN,ypimin,ypigradN 
double precision :: GradNpsuro,IntWd1s0,IntWd1s2,IntWd3s0,DvelN,viscN  
double precision :: GravN,PressB,distpi,distpimin,pressib,pressibmin
double precision :: QiiIntWdS,level,pressj,Qsi,Qsj,velix,veliz,veliq,hcrit
double precision :: hcritmin,zbottom,FlowRate1,Lb,L,minquotanode,maxquotanode
double precision :: SomQsiQsj,DiffQsiQsj,tpres_save1,ViscoMon_save1
integer(4),dimension(1:PLANEDIM) :: acix
double precision,dimension(1:PLANEDIM) :: IntLocXY,RHS,RG,ss,nnlocal,gradbPsuro
double precision,dimension(1:PLANEDIM) :: ViscoMon,ViscoShear,sidevel,TT,Dvel
type (TyBoundarySide) :: RifBoundarySide
character(4):: strtype
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
cinvisci = pg(npi)%visc
roi = pg(npi)%dens
pressi = pg(npi)%pres
distpimin = Domain%dd
tpres_save1 = zero
ViscoMon_save1 = zero
RHS(1) = -tpres(1)
RHS(2) = -tpres(3)
ViscoMon(1) = zero
ViscoMon(2) = zero
ViscoShear(1) = zero
ViscoShear(2) = zero
ibdt = BoundaryDataPointer(3,npi)
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
   do i = 1, PLANEDIM
      ss(i) = RifBoundarySide%T(acix(i), acix(1))
      nnlocal(i) = RifBoundarySide%T(acix(i), acix(2))
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
   end select
! Boundary contribution to "gradP" 
! (pressure gradient term in the momentum equation) 
   if ((strtype=="fixe").OR.(strtype=="tapi"))  then
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
         level = Tratto(sidestr)%ShearCoeff  
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
               if (tempo<Tratto(sidestr)%trampa) then
                  pg(npi)%velass(1) = Tratto(sidestr)%NormVelocity*nnlocal(1) *&
                     tempo / Tratto(sidestr)%trampa
                  pg(npi)%velass(2) = zero           
                  pg(npi)%velass(3) = Tratto(sidestr)%NormVelocity*nnlocal(2) *&
                     tempo / Tratto(sidestr)%trampa
                  else
                     pg(npi)%velass(1) = Tratto(sidestr)%NormVelocity*nnlocal(1)
                     pg(npi)%velass(2) = zero            
                     pg(npi)%velass(3) = Tratto(sidestr)%NormVelocity*nnlocal(2)
               endif
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
                  if (tempo<Tratto(sidestr)%trampa) then
                     pg(npi)%velass(1) = Tratto(sidestr)%NormVelocity *        &
                        nnlocal(1) * tempo / Tratto(sidestr)%trampa
                     pg(npi)%velass(2) = zero           
                     pg(npi)%velass(3) = Tratto(sidestr)%NormVelocity *        &
                        nnlocal(2) * tempo / Tratto(sidestr)%trampa
                     else
                        pg(npi)%velass(1) = Tratto(sidestr)%NormVelocity *     &
                           nnlocal(1)
                        pg(npi)%velass(2) = zero            
                        pg(npi)%velass(3) = Tratto(sidestr)%NormVelocity *     &
                           nnlocal(2)
                  endif
               endif        
               return             
               elseif (strtype=="sour") then
                  if (xpi>=zero.AND.xpi<=RifBoundarySide%length) then
                     if ((Domain%time_stage==1).or.(Domain%time_split==1)) then 
                        pg(npi)%kodvel = 2
                        if (tempo<Tratto(sidestr)%trampa) then
                           pg(npi)%velass(1) = Tratto(sidestr)%NormVelocity *  &
                                               nnlocal(1) * tempo /            &
                                               Tratto(sidestr)%trampa
                           pg(npi)%velass(2) = zero           
                           pg(npi)%velass(3) = Tratto(sidestr)%NormVelocity *  &
                                               nnlocal(2) * tempo /            &
                                               Tratto(sidestr)%trampa
                           else
                              pg(npi)%velass(1) = Tratto(sidestr)%NormVelocity &
                                 * nnlocal(1)
                              pg(npi)%velass(2) = zero            
                              pg(npi)%velass(3) = Tratto(sidestr)%NormVelocity &
                                 * nnlocal(2)
                        endif
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
                     hcrit = veliq / Abs(Domain%grav(3))        
                     hcritmin = Tratto(sidestr)%ShearCoeff       
                     if (hcritmin<Domain%h) hcritmin = Domain%h       
                     if (hcrit<hcritmin) hcrit = hcritmin          
                     zbottom = Vertice(3, RifBoundarySide%vertex(1))    
                     if (zbottom>Vertice(3, RifBoundarySide%vertex(2)))        &
                        zbottom = Vertice(3, RifBoundarySide%vertex(2))
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
! explosion
      tpres_save1 = tpres_save1 - (nnlocal(i) * QiiIntWdS + RG(i)) * DvelN *   &
         nnlocal(i)
! explosion
      RHS(i) = RHS(i) - nnlocal(i) * QiiIntWdS + RG(i)
   enddo
! Volume viscosity force (with changed sign)
   if (strtype=="fixe".OR.strtype=="tapi") then
      if (xpi>=zero.AND.xpi<=RifBoundarySide%length) then
         SVcoeff = Tratto(sidestr)%ShearCoeff
         do i=1,PLANEDIM
            Dvel(i) = two * (pg(npi)%var(acix(i)) - sidevel(i))
         enddo
         DvelN = Dvel(1) * nnlocal(1) + Dvel(2) * nnlocal(2)
         viscN = 0.3333333d0 * cinvisci
         celeri = Med(imed)%celerita
         alfaMon = Med(imed)%alfaMon
         Monvisc = alfaMon * celeri * Domain%h
         viscN = Monvisc
         TN = viscN * DvelN * IntWd3s0
         do i=1,PLANEDIM
            TT(i) = TN * nnlocal(i)
         enddo
! Shear viscosity force (with changed sign)
         SVforce = SVcoeff * (cinvisci + cinvisci) * IntWd1s0
         do i=1,PLANEDIM
            ViscoMon(i) = ViscoMon(i) + TT(i)
! explosion
            ViscoMon_save1 = ViscoMon_save1 + ViscoMon(i) * DvelN * nnlocal(i)
! explosion
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
! contribution for specific internal energy
if (esplosione) then
   pg(npi)%dEdT = - half * (tpres_save1 - ViscoMon_save1)
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine AddBoundaryContributions_to_ME2D

