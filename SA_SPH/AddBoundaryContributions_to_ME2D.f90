!cfile AddBoundaryContributions_to_ME2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : AddBoundaryContributions_to_ME2D
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
!AA504
! 02  Amicarelli        08/04/2014     (v5.04) Monaghan's term always active (also for separating particles)   
!
!************************************************************************************
! Module purpose : Module to compute boundary contributions to rodivV,
!                  gradPsuro and ViscoF relative to particle npi
!
! Calling routine: Loop_Irre_2D
!
! Called routines: 
!
!************************************************************************************
!
subroutine AddBoundaryContributions_to_ME2D (npi, IntNcbs, tpres, tdiss, tvisc) !, Ncbs
!
!.. Computes boundary contributions to rodivV, gradPsuro and ViscoF relative to particle npi
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use BoundIntegralTab_Module
use FILES_ENTITIES
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
double precision,parameter :: pibmink = 10.0d0
double precision,parameter :: eps = 0.1d0
!
!.. Formal Arguments ..
integer(4),      intent(IN)    :: npi
!integer(4),      intent(IN)    :: Ncbs
integer(4),      intent(IN)    :: IntNcbs
double precision,intent(INOUT),dimension(1:SPACEDIM) :: tpres
double precision,intent(INOUT),dimension(1:SPACEDIM) :: tdiss
double precision,intent(INOUT),dimension(1:SPACEDIM) :: tvisc
!
!.. Local Scalars ..
integer(4)       :: imed  !, izonelocal
integer(4)       :: i, j, icbs, ibdt, ibdp, iside, sidestr
double precision :: IntWds, IntWdV, cinvisci, Monvisc, &
                    roi, ro0, celeri, pressi, alfaMon, xpi, ypi, &   !Qii, stiffi, 
                    SVforce, SVcoeff, TN, ypimin, ypigradN, GradNpsuro, &
                    IntWd1s0, IntWd1s2, IntWd3s0, DvelN, viscN, &             !viscS, 
                    GravN, PressB, distpi, distpimin, pressib, pressibmin
double precision :: QiiIntWdS, level, pressj, Qsi, Qsj, velix, veliz, veliq, hcrit, hcritmin, zbottom
double precision :: FlowRate1, Lb , L, minquotanode, maxquotanode, &
                    SomQsiQsj, DiffQsiQsj
double precision :: tpres_save1, ViscoMon_save1
character(4):: strtype
!
!.. Local Arrays ..
integer(4),      dimension(1:PLANEDIM) :: acix
double precision,dimension(1:PLANEDIM) :: IntLocXY
double precision,dimension(1:PLANEDIM) :: RHS
double precision,dimension(1:PLANEDIM) :: RG
double precision,dimension(1:PLANEDIM) :: ss
double precision,dimension(1:PLANEDIM) :: nnlocal
double precision,dimension(1:PLANEDIM) :: gradbPsuro, ViscoMon, ViscoShear, sidevel, TT, Dvel !, tpres_save, ViscoMon_save

type (TyBoundarySide) :: RifBoundarySide
!
!.. Executable Statements ..
!
 acix(1) = 1  !active coordinate indexes
 acix(2) = 3
!
!AA404 sub
          if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
   pg(npi)%kodvel = 0
   pg(npi)%velass(:) = zero
 end if
 imed = pg(npi)%imed
 ro0 = Med(imed)%den0
 cinvisci = pg(npi)%visc
 roi = pg(npi)%dens
 pressi = pg(npi)%pres
! Qii = (pressi + pressi) / roi !AdM 15-10-08
 distpimin = Domain%dd
!
 tpres_save1 = zero
! tpres_save = zero
! ViscoMon_save = zero
 ViscoMon_save1 = zero
 RHS(1) = -tpres(1)
 RHS(2) = -tpres(3)
! Amatr(1, 1) = one  !AdM 15-10-08
! Amatr(1, 2) = zero !AdM 15-10-08
! Amatr(2, 1) = zero !AdM 15-10-08
! Amatr(2, 2) = one  !AdM 15-10-08
 ViscoMon(1) = zero
 ViscoMon(2) = zero
 ViscoShear(1) = zero
 ViscoShear(2) = zero
!
 ibdt = BoundaryDataPointer(3,npi)
!
 do icbs = 1, IntNcbs
!
   ibdp = ibdt + icbs - 1
   IntLocXY(1:PLANEDIM) = BoundaryDataTab(ibdp)%LocXYZ(1:PLANEDIM)
   iside = BoundaryDataTab(ibdp)%CloBoNum
!
   RifBoundarySide = BoundarySide(iside)
   sidestr = RifBoundarySide%stretch
   strtype = Tratto(sidestr)%tipo
!
   if (strtype == "open") cycle
!
   xpi = IntLocXY(1)
   ypi = IntLocXY(2)
!
   IntWdS = BoundaryDataTab(ibdp)%BoundaryIntegral(1)
!   IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(2)
   IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(3)
   IntWd1s0 = BoundaryDataTab(ibdp)%BoundaryIntegral(6)
   IntWd3s0 = BoundaryDataTab(ibdp)%BoundaryIntegral(7)
   IntWd1s2 = BoundaryDataTab(ibdp)%BoundaryIntegral(8)
!
   GravN = zero !AdM 15-10-08
   do i = 1, PLANEDIM
     ss(i) = RifBoundarySide%T(acix(i), acix(1))
     nnlocal(i) = RifBoundarySide%T(acix(i), acix(2))
     gradbPsuro(i) = Domain%grav(acix(i))
     GravN = GravN + Domain%grav(acix(i)) * nnlocal(i) !AdM 15-10-08
   end do
!
!   sidevel = zero
   Select Case (strtype)
   Case ("fixe")
     do i = 1, PLANEDIM
       sidevel(i) = zero
     end do
   Case ("tapi")
     do i = 1, PLANEDIM
       sidevel(i) = RifBoundarySide%velocity(acix(i))
     end do
   End Select
!
!*******  Boundary contribution to gradP  ******************************************
! 
!   QiiIntWdS = zero
   if (strtype == "fixe" .OR. strtype == "tapi")  then
!
!     QiiIntWdS = Qii * IntWds !AdM 15-10-08
     do i = 1, PLANEDIM
!       Dvel(i) = pg(npi)%vel(acix(i)) - sidevel(i)   ! vel velocita' mediata
       Dvel(i) = two * (pg(npi)%var(acix(i)) - sidevel(i))   ! var velocita' mediata
     end do
     DvelN = Dvel(1) * nnlocal(1) + Dvel(2) * nnlocal(2)
!     distpi = max(ypi,Domain%h) !AdM 15-10-08
     distpi = max(ypi,distpimin)
     pressib = pressi
     pressibmin = pibmink * ro0 * distpimin
     if (DvelN <= zero) then
       if (pressib < pressibmin) pressib = pressibmin
     end if
!    
!     PressB = pressi - ro0 * GravN * distpi !sostituisce quella sopra  !AdM 15-10-08
!     QiiIntWdS = IntWds * (pressi + PressB) / roi  !sostituisce quella sopra  !AdM 15-10-08
     PressB = pressib - ro0 * GravN * distpi 
     QiiIntWdS = IntWdS * (pressib + PressB) / roi  !sostituisce quella sopra  !AdM 15-10-08
!
   else if (strtype == "leve") then   !flusso subcritico
!
     Qsi = pressi / roi
!     Qsj = zero
     level = Tratto(sidestr)%ShearCoeff  
     if (pg(npi)%coord(3) < level) then    
       pressj = Med(imed)%den0 * Domain%grav(3) * (pg(npi)%coord(3) - level)
       Qsj = pressj / Med(imed)%den0  
       SomQsiQsj = Qsi + Qsj        
       DiffQsiQsj = Qsi - Qsj       
     else            
       SomQsiQsj = zero   
       DiffQsiQsj = zero    
     end if         
!     pressj = Med(imed)%den0 * Domain%grav(3) * (pg(npi)%coord(3) - level)  
!     if (pressj < -pressi)  pressj = -pressi   
!     Qsj = pressj / roi        
!     QiiIntWdS = (Qsi + Qsj) * IntWds       
     QiiIntWdS = SomQsiQsj * IntWds       
     ypimin = eps * Domain%h               
     ypigradN = ypi          
     if (ypigradN < ypimin) ypigradN = ypimin     
     GradNpsuro = (Qsi - Qsj) / ypigradN       
     !if (GradNpsuro < zero) GradNpsuro = zero  
!     gradbPsuro(1) = GradNpsuro * nnlocal(1)     
!     gradbPsuro(2) = GradNpsuro * nnlocal(2)       
     !gradbPsuro(1) = zero         
     !gradbPsuro(2) = zero            
     gradbPsuro(1) = gradbPsuro(1) + GradNpsuro * nnlocal(1)    
     gradbPsuro(2) = gradbPsuro(2) + GradNpsuro * nnlocal(2)     
!  
! Imposizione a zero della componente verticale di velocita'
     !pg(npi)%kodvel = 1    !Prova per confronto con Fresurf_SPH_2D   
!     pg(npi)%kodvel = 1               
     pg(npi)%velass(:) = zero             
!
   else if (strtype == "velo") then   !velocita' normale assegnata  
!  
! Imposizione delle componenti di velocita'     
!
!AA404 sub
    if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
       pg(npi)%kodvel = 2
       if (tempo < Tratto(sidestr)%trampa) then
         pg(npi)%velass(1) = Tratto(sidestr)%NormVelocity*nnlocal(1) * tempo / Tratto(sidestr)%trampa
         pg(npi)%velass(2) = zero           
         pg(npi)%velass(3) = Tratto(sidestr)%NormVelocity*nnlocal(2) * tempo / Tratto(sidestr)%trampa
       else
         pg(npi)%velass(1) = Tratto(sidestr)%NormVelocity*nnlocal(1)
         pg(npi)%velass(2) = zero            
         pg(npi)%velass(3) = Tratto(sidestr)%NormVelocity*nnlocal(2)
       end if
     end if
     return             
!
   else if (strtype == "flow") then   ! portata normale assegnata  
!
! Imposizione delle componenti di velocita'     
!
!AA404 sub
       if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
       pg(npi)%kodvel = 2
! calcolo velocita' normale
       if (BoundarySide(iside)%CloseParticles > 0) then
         minquotanode = min (vertice(3,BoundarySide(iside)%vertex(1)),vertice(3,BoundarySide(iside)%vertex(2)))
         maxquotanode = max (vertice(3,BoundarySide(iside)%vertex(1)),vertice(3,BoundarySide(iside)%vertex(2)))
!
         Lb = BoundarySide(iside)%CloseParticles_maxQuota - minquotanode
         L = maxquotanode - minquotanode
!         L =  BoundarySide(iside)%length
         FlowRate1 = Tratto(sidestr)%FlowRate * Lb / L
         Tratto(sidestr)%NormVelocity = doubleh * FlowRate1 / (BoundarySide(iside)%CloseParticles * Domain%PVolume)
       else
         Tratto(sidestr)%NormVelocity = zero
       end if
!
       if (tempo < Tratto(sidestr)%trampa) then
         pg(npi)%velass(1) = Tratto(sidestr)%NormVelocity*nnlocal(1) * tempo / Tratto(sidestr)%trampa
         pg(npi)%velass(2) = zero           
         pg(npi)%velass(3) = Tratto(sidestr)%NormVelocity*nnlocal(2) * tempo / Tratto(sidestr)%trampa
       else
         pg(npi)%velass(1) = Tratto(sidestr)%NormVelocity*nnlocal(1)
         pg(npi)%velass(2) = zero            
         pg(npi)%velass(3) = Tratto(sidestr)%NormVelocity*nnlocal(2)
       end if
     end if
     return             
!
   else if (strtype == "sour") then
!
     if (xpi >= zero .AND. xpi <= RifBoundarySide%length) then
!
!AA404 sub
       if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
         pg(npi)%kodvel = 2
         if (tempo < Tratto(sidestr)%trampa) then
           pg(npi)%velass(1) = Tratto(sidestr)%NormVelocity*nnlocal(1) * tempo / Tratto(sidestr)%trampa
           pg(npi)%velass(2) = zero           
           pg(npi)%velass(3) = Tratto(sidestr)%NormVelocity*nnlocal(2) * tempo / Tratto(sidestr)%trampa
         else
           pg(npi)%velass(1) = Tratto(sidestr)%NormVelocity*nnlocal(1)
           pg(npi)%velass(2) = zero            
           pg(npi)%velass(3) = Tratto(sidestr)%NormVelocity*nnlocal(2)
         end if
       end if
       return
     end if
!
   else if (strtype == "crit") then  !flusso critico   
!
! Condizione di flusso critico non stazionario
     Qsi = pressi / roi                          
     velix = pg(npi)%vel(1)         
     veliz = pg(npi)%vel(3)         
     veliq = velix * velix + veliz * veliz  
     hcrit = veliq / Abs(Domain%grav(3))        
     hcritmin = Tratto(sidestr)%ShearCoeff       
     if (hcritmin < Domain%h) hcritmin = Domain%h       
     if (hcrit < hcritmin) hcrit = hcritmin          
     zbottom = Vertice(3, RifBoundarySide%vertex(1))    
     if (zbottom > Vertice(3, RifBoundarySide%vertex(2))) &
     zbottom = Vertice(3, RifBoundarySide%vertex(2))
     level = zbottom + hcrit          
     if (pg(npi)%coord(3)< level) then          
       pressj = Med(imed)%den0 * Domain%grav(3) * (pg(npi)%coord(3) - level)
       Qsj = pressj / Med(imed)%den0       
       SomQsiQsj = Qsi + Qsj              
       DiffQsiQsj = Qsi - Qsj             
     else                                       
       SomQsiQsj = zero            
       DiffQsiQsj = zero                   
     end if                                          
!     if (pressj < -pressi) pressj = -pressi           
!     Qsj = pressj / roi                       
!     QiiIntWdS = (Qsi + Qsj) * IntWds          
     QiiIntWdS = SomQsiQsj * IntWds              
     ypimin = eps * Domain%h                      
     ypigradN = ypi                                    
     if (ypigradN < ypimin) ypigradN = ypimin           
     GradNpsuro = DiffQsiQsj / ypigradN                    
     !if (GradNpsuro < zero) GradNpsuro = zero                  
!     gradbPsuro(1) = GradNpsuro * nnlocal(1)          
!     gradbPsuro(2) = GradNpsuro * nnlocal(2)             
     !gradbPsuro(1) = zero                                      
     !gradbPsuro(2) = zero                                    
     gradbPsuro(1) = gradbPsuro(1) + GradNpsuro * nnlocal(1)           
     gradbPsuro(2) = gradbPsuro(2) + GradNpsuro * nnlocal(2)              
!
! Imposizione a zero della componente verticale di velocita'
!
!AA404 sub
     if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
       pg(npi)%kodvel = 1                                                 
       pg(npi)%velass(:) = zero
     end if             
!
   else if (strtype == "open") then       
!
!! Imposizione a zero della componente verticale di velocita'
!     pg(npi)%kodvel = 1                        
!     pg(npi)%velass(:) = zero                    
!     cycle                               
     cycle                                 
!
   end if
!
   RG = zero
   do i = 1, PLANEDIM
     do j = 1, PLANEDIM
       RG(i) = RG(i) + RifBoundarySide%RN(acix(i), acix(j)) * gradbPsuro(j) 
     end do
     RG(i) = RG(i) * IntWdV
   end do
!
   do i = 1, PLANEDIM
!!     tpres_save(i) = tpres_save(i) - nnlocal(i) * QiiIntWdS + RG(i)
!     tpres_save(i) = tpres_save(i) - (nnlocal(i) * QiiIntWdS + RG(i)) * DvelN * nnlocal(i)
!.. explosion
     tpres_save1 = tpres_save1 - (nnlocal(i) * QiiIntWdS + RG(i)) * DvelN * nnlocal(i)
!.. explosion
     RHS(i) = RHS(i) - nnlocal(i) * QiiIntWdS + RG(i)
!------------------- !AdM 15-10-08 ----------------------------------------------
!     do j = 1, PLANEDIM
!       Amatr(i, j) = Amatr(i, j) - RifBoundarySide%R(acix(i), acix(j)) * IntWdV
!     end do
!------------------- !AdM 15-10-08 ----------------------------------------------
   end do
!
!
!*******  Boundary contribution to ViscoF  ******************************************
!
! Volume viscosity force (with changed sign)
   if (strtype == "fixe" .OR. strtype == "tapi") then
     if (xpi >= zero .AND. xpi <= RifBoundarySide%length) then
       SVcoeff = Tratto(sidestr)%ShearCoeff
       do i = 1, PLANEDIM
!         Dvel(i) = pg(npi)%vel(acix(i)) - sidevel(i)
         Dvel(i) = two * (pg(npi)%var(acix(i)) - sidevel(i))
       end do
!       DvelS = Dvel(1) * ss(1) + Dvel(2) * ss(2)
       DvelN = Dvel(1) * nnlocal(1) + Dvel(2) * nnlocal(2)
       viscN = 0.3333333d0 * cinvisci
!       viscN = zero
!       viscS = viscN * SVcoeff
!       viscS = zero
!AA504 removed line in order to keep Monaghan's term, even when a particle leaves a boundary
!         stiffi = Med(imed)%eps
         celeri = Med(imed)%celerita
         alfaMon = Med(imed)%alfaMon
         Monvisc = alfaMon * celeri * Domain%h
         viscN = Monvisc
!       TS = viscS * DvelS * IntWd1s2
       TN = viscN * DvelN * IntWd3s0
       do i = 1, PLANEDIM
!         TT(i) = TS * ss(i) + TN * nnlocal(i)
         TT(i) = TN * nnlocal(i)
       end do
! Shear viscosity force (with changed sign)
       SVforce = SVcoeff * (cinvisci + cinvisci) * IntWd1s0
       do i = 1, PLANEDIM
         ViscoMon(i) = ViscoMon(i) + TT(i)
!         ViscoMon_save(i) = ViscoMon_save(i) + ViscoMon(i) * DvelN * nnlocal(i)
!.. explosion
         ViscoMon_save1 = ViscoMon_save1 + ViscoMon(i) * DvelN * nnlocal(i)
!.. explosion
         ViscoShear(i) = ViscoShear(i) + Dvel(i) * SVforce
       end do
     end if
   end if
!
 end do
!
!.. inlining della subroutine RisolviSistema_2x2. Risolve il sistema con la regola di Cramer
!
!------------------- !AdM 15-10-08 ----------------------------------------------
!  detAmatr = amatr(1,1) * amatr(2,2) - amatr(2,1) * amatr(1,2)
!  if (abs(detAmatr) < epsilon(detAmatr)) then
!
!.. la matrice è singolare
!
!    write(nout,*) "la matrice e'' singolare. detA=",detAmatr
!    write(nscr,*) "la matrice e'' singolare. detA=",detAmatr
!  else
!    gradPsuro(1) = (RHS(1) * amatr(2,2) - RHS(2) * amatr(1,2))/detAmatr
!    gradPsuro(2) = (amatr(1,1) * RHS(2) - amatr(2,1) * RHS(1))/detAmatr
!  end if
!------------------- !AdM 15-10-08 ----------------------------------------------
!
  do i = 1, PLANEDIM
!    tpres(acix(i)) = -gradPsuro(i) !AdM 15-10-08
    tpres(acix(i)) = - RHS(i) !AdM 15-10-08
    tdiss(acix(i)) = tdiss(acix(i)) - ViscoMon(i)
    tvisc(acix(i)) = tvisc(acix(i)) - ViscoShear(i)
  end do
!
!............................................ 2011 mar 08
!.. contribution for Specific Internal Energy
!
  if (esplosione) then
!    do i = 1, PLANEDIM
!!!      pg(npi)%dEdT = pg(npi)%dEdT - half * ( dvel(i)*(tpres_save(i) - ViscoMon(i)) )
!!      pg(npi)%dEdT = pg(npi)%dEdT - half * ( tpres_save(i) - ViscoMon_save(i) )
!!      pg(npi)%dEdT = pg(npi)%dEdT - half * ( tpres_save(i) )
!      pg(npi)%dEdT = pg(npi)%dEdT - half * ( tpres_save(i) - dvel(i)*ViscoMon(i) )
!    end do
      pg(npi)%dEdT = - half * ( tpres_save1 - ViscoMon_save1 )
  end if
!
!.. end contribution for Specific Internal Energy
!................................................
!
!prova!
!  tpres = floor(tpres * azzeramento) / azzeramento
!  where (dabs(tpres) < arrotondamento) tpres = zero
!  tdiss = floor(tdiss * azzeramento) / azzeramento
!  where (dabs(tdiss) < arrotondamento) tdiss = zero
!  tvisc = floor(tvisc * azzeramento) / azzeramento
!  where (dabs(tvisc) < arrotondamento) tvisc = zero
!prova!
!
  return
  end subroutine AddBoundaryContributions_to_ME2D
!---split

