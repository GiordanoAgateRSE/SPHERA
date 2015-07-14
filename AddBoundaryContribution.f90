!cfile AddBoundaryContribution_to_CE2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : AddBoundaryContribution_to_CE2D
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!
!************************************************************************************
! Module purpose : Module to compute boundary contributions to rodivV
!
! Calling routine: Loop_Irre_2D
!
! Called routines: /
!
!************************************************************************************
!
  subroutine AddBoundaryContribution_to_CE2D (npi, IntNcbs, BCrodivV)  !, Ncbs
!
!Computes boundary contributions to rodivV
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
  use Diagnostic_MODULE
  use BoundIntegralTab_Module
!
!.. Implicit Declarations ..
  implicit none
!
!  double precision,parameter :: eps = 0.005d0 !AdM 15-10-08
!
!.. Formal Arguments ..
  integer(4),      intent(IN)      :: npi
!  integer(4),      intent(IN)      :: Ncbs
  integer(4),      intent(IN)      :: IntNcbs
  double precision,intent(INOUT)   :: BCrodivV
!
!.. Local Scalars ..
  integer(4)       :: pd, icbs, iside, sidestr, ibdt, ibdp
  double precision :: IntWds, roi, vin   !,IntWdV,  xpi, ypi,xpmin, xpmax, interlen,  etalocal !AdM 15-10-08
  character(4)     :: strtype
  type (TyBoundarySide) :: RifBoundarySide
!
!.. Local Arrays ..
  double precision,dimension(1:PLANEDIM)    :: IntLocXY
  integer(4),      dimension(1:PLANEDIM)    :: acix
  double precision,dimension(1:PLANEDIM)    :: nnlocal
  double precision,dimension(1:PLANEDIM)    :: Dvel
!
!.. Executable Statements ..
!
  acix(1) = 1        !active coordinate indexes
  acix(2) = 3
!
  BCrodivV = zero
!
  if (IntNcbs <= 0) return
!
  roi = pg(npi)%dens
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
    if (strtype == "fixe" .OR. strtype == "tapi" .OR. strtype == "velo" .OR. strtype == "flow" .OR. strtype == "sour") then 
!
      IntWdS = BoundaryDataTab(ibdp)%BoundaryIntegral(1)
!      IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(3) !AdM 15-10-08
      vin = zero
!
      do pd = 1, PLANEDIM
        nnlocal(pd) = RifBoundarySide%T(acix(pd), acix(2))
      end do
!
      select case (strtype)
!
        case ("fixe")
          do pd = 1, PLANEDIM
            Dvel(pd) = pg(npi)%var(acix(pd))
          end do
          vin = Dvel(1) * nnlocal(1) + Dvel(2) * nnlocal(2)
!          etalocal = eps * Domain%h !AdM 15-10-08
!          SIntWds2 = IntWdV / (ypi + etalocal) !AdM 15-10-08
!
        case ("tapi")
          do pd = 1, PLANEDIM
            Dvel(pd) = pg(npi)%var(acix(pd))-RifBoundarySide%velocity(acix(pd))
          end do
          vin = Dvel(1) * nnlocal(1) + Dvel(2) * nnlocal(2)
!          etalocal = eps * Domain%h !AdM 15-10-08
!          SIntWds2 = IntWdV / (ypi + etalocal) !AdM 15-10-08
!
        case ("velo", "flow", "sour")
!AA404 sub
          if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
            pg(npi)%koddens = 2
            pg(npi)%densass = roi
            BCrodivV = zero
          end if
          return
!
      end select
!
        !****  Boundary contribution to rodivV  ****
!
!      BCrodivV = BCrodivV + vin * roi * ( IntWds + SIntWds2 ) !AdM 15-10-08
      BCrodivV = BCrodivV + two * vin * roi * IntWdS
!
    end if 
!
  end do
!
  return
  end subroutine AddBoundaryContribution_to_CE2D
!---split

!cfile AddBoundaryContribution_to_CE3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : AddBoundaryContribution_to_CE3D
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!
!************************************************************************************
! Module purpose : Module to compute boundary contributions to rodivV, in 3D geometry,
!                  relative to particle npi
!                  Additional term due to a linear distribution of normal velocity
!                  has been introduced
!
! Calling routine: Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
!
  subroutine AddBoundaryContribution_to_CE3D (npi, Ncbf, BCtorodivV)
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
  integer(4),      intent(IN)    :: npi
  integer(4),      intent(IN)    :: Ncbf
  double precision,intent(INOUT) :: BCtorodivV
!
!.. Local Scalars ..
  integer(4)       :: sd, sdj, icbf, iface, ibdt, ibdp
  integer(4)       :: stretch        
  double precision :: roi, scaprod
  character(4)     :: boundtype            
!
!.. Local Arrays ..
  double precision,dimension(1:SPACEDIM) :: vb
  double precision,dimension(1:SPACEDIM) :: vi
  double precision,dimension(1:SPACEDIM) :: dvij
  double precision,dimension(1:SPACEDIM) :: LocPi, LocDvij
!
!.. Executable Statements ..
!
  roi = pg(npi)%dens
  vi(:) = pg(npi)%var(:)
  BCtorodivV = zero
!
  if (Ncbf <= 0) return
!
  ibdt = BoundaryDataPointer(3,npi)
!
  do icbf = 1, Ncbf
!
    ibdp = ibdt + icbf - 1
    LocPi(1:SPACEDIM) = BoundaryDataTab(ibdp)%LocXYZ(1:SPACEDIM)
    iface = BoundaryDataTab(ibdp)%CloBoNum
!
    stretch = BoundaryFace(iface)%stretch
    boundtype = Tratto(stretch)%tipo
!
    if (boundtype == "fixe" .OR. boundtype == "tapi") then
!
      if (LocPi(3) > zero) then          !La faccia "iface" interagisce con la particella Pi
!
        vb(:) = BoundaryFace(iface)%velocity(:)
        dvij(:) = two * (vi(:) - vb(:))
!
!!!        scaprod = zero
!!!        do sd = 1, SPACEDIM
!!!          Locdvij(sd) = zero  !componenti locali del vettore 2*(vi-vb)
!!!          do sdj = 1, SPACEDIM
!!!            Locdvij(sd) = Locdvij(sd) + dvij(sdj) * BoundaryFace(iface)%T(sdj,sd)
!!!          end do
!!!          scaprod = scaprod + Locdvij(sd) * BoundaryDataTab(ibdp)%BoundaryIntegral(3+sd)
!!!        end do
        scaprod = zero
        sd = 3
        Locdvij(sd) = zero  !componenti locali del vettore 2*(vi-vb)
        do sdj = 1, SPACEDIM
          Locdvij(sd) = Locdvij(sd) + dvij(sdj) * BoundaryFace(iface)%T(sdj,sd)
        end do
        scaprod = scaprod + Locdvij(sd) * BoundaryDataTab(ibdp)%BoundaryIntegral(3+sd)        
!
        !*******  Boundary contribution to divV  ***********************************
        BCtorodivV = BCtorodivV + roi * scaprod
!
      end if
!
    else if (boundtype == "velo" .or. boundtype == "flow" .or. boundtype == "sour") then
!
!AA404 sub
          if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
        pg(npi)%koddens = 2
        pg(npi)%densass = roi
        BCtorodivV = zero
      end if
      return
!
    end if
!
  end do
!
  return
  end subroutine AddBoundaryContribution_to_CE3D
!---split

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

!cfile AddBoundaryContributions_to_ME3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name    : AddBoundaryContributions_to_ME3D
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07      Graphic windows calls removed
! 01  Agate/Flamini    08/10/07      Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!AA504
! 03 Amicarelli         08/04/2014     (v5.04) Viscosity term depending on velocity divergence is removed (negligible); boundary term depending 
!                                      on molecular viscosity is disabled because of errors for granular flows and general lack of validation.
!
!************************************************************************************
! Module purpose : Module to compute boundary contributions to rodivV,
!              gradPsuro and ViscoF relative to particle npi
!              Performs implicit computation of gradPsuro
!
! Calling routine: Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
!
  subroutine AddBoundaryContributions_to_ME3D (npi, Ncbf, tpres, tdiss, tvisc)
!
!Computes boundary contributions to rodivV, gradPsuro and ViscoF relative to particle npi
!Performs implicit computation of gradPsuro
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
  use FILES_ENTITIES
!
!.. Implicit Declarations ..
  implicit none
!
!.. Local Parameters ..
!  double precision,parameter :: etaratio  = 0.1d0
!  double precision,parameter :: minIntWdV = 0.001d0  
!  double precision,parameter :: pibmink = 10.0d0
!
!.. Formal Arguments ..
  integer(4),      intent(IN)    :: npi
  integer(4),      intent(IN)    :: Ncbf
  double precision,intent(INOUT),dimension(1:SPACEDIM) :: tpres, tdiss, tvisc
!
!.. Local Scalars ..
  integer(4)       :: sd, icbf, iface, ibdp
  integer(4)       :: sdj, i, j, mati, stretch
  double precision :: IntdWrm1dV
  double precision :: cinvisci, Monvisc, cinviscmult, pressi, dvn
  double precision :: FlowRate1, Lb, L, minquotanode, maxquotanode
  double precision :: Qii, roi, celeri, alfaMon, Mmult, IntGWZrm1dV
  double precision :: tpres_save1, ViscoMon_save1
!
  character(4)     :: stretchtype
!
!.. local Arrays ..
  double precision,dimension(1:SPACEDIM) :: vb, vi, dvij, RHS, nnlocal, Grav_Loc, Gpsurob_Loc, Gpsurob_Glo
  double precision,dimension(1:SPACEDIM) :: ViscoMon, ViscoShear, LocPi
!
!.. Executable Statements ..
!
!.. initializations
!
  mati = pg(npi)%imed
  cinvisci = pg(npi)%visc
  roi = pg(npi)%dens
  pressi = pg(npi)%pres
  Qii = (pressi + pressi) / roi
  tpres_save1 = zero
  ViscoMon_save1 = zero
!
  vi(:) = pg(npi)%var(:)
!
  RHS(:) = zero
  ViscoMon(:) = zero
  ViscoShear(:) = zero
!
!AA404 sub
  if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
    pg(npi)%kodvel = 0
    pg(npi)%velass(:) = zero
  end if
!
  face_loop: do icbf = 1,Ncbf
!
    ibdp = BoundaryDataPointer(3,npi) + icbf - 1 
    iface = BoundaryDataTab(ibdp)%CloBoNum
    stretch = BoundaryFace(iface)%stretch
    stretchtype = Tratto(stretch)%tipo
    nnlocal(:) = BoundaryFace(iface)%T(:,3)
!
!.. skips for the open boundary condition (no constraint must be applied in this case)
    if (stretchtype == "open") cycle face_loop
!
!!!    LocPi(:) = BoundaryDataTab(ibdp)%LocXYZ(:) !Coordinate locali della particella reale Pi
!!!    if (LocPi(3) > zero) then                  !La faccia "iface" interagisce con la particella Pi
!
    if (stretchtype == "fixe" .Or. stretchtype == "tapi") then
!
      LocPi(:) = BoundaryDataTab(ibdp)%LocXYZ(:) !Coordinate locali della particella reale Pi
      if (LocPi(3) > zero) then                  !La faccia "iface" interagisce con la particella Pi
        dvn = zero
        do SD = 1,SPACEDIM
          vb(SD) = BoundaryFace(iface)%velocity(SD)
          dvij(SD) = two * (vi(SD) - vb(SD))
          dvn = dvn + BoundaryFace(iface)%T(SD, 3) * dvij(SD)
          Grav_Loc(SD) = zero                           !Componenti locali della gravita'
          do sdj = 1,SPACEDIM
            Grav_Loc(SD) = Grav_Loc(SD) + BoundaryFace(iface)%T(sdj, SD) * Domain%grav(sdj)
          end do
        end do
           
        !*******  Boundary contribution to gradPsuro  ***********************************
            
        !Componenti locali
        do i = 1,SPACEDIM
          Gpsurob_Loc(i) = -Qii * BoundaryDataTab(ibdp)%BoundaryIntegral(3+i)
          do j = 1,SPACEDIM
            Gpsurob_Loc(i) = Gpsurob_Loc(i) - BoundaryDataTab(ibdp)%IntGiWrRdV(i, j) * Grav_Loc(j)
          end do
        end do
        !Componenti globali
        do i = 1,SPACEDIM
          Gpsurob_Glo(i) = zero
          do j = 1,SPACEDIM
            Gpsurob_Glo(i) = Gpsurob_Glo(i) + BoundaryFace(iface)%T(i, j) * Gpsurob_Loc(j)
!.. explosion
          if (esplosione) then
            tpres_save1 = tpres_save1 - (nnlocal(i) * BoundaryDataTab(ibdp)%IntGiWrRdV(i, j) &
                        + Gpsurob_Glo(i)) * dvn * nnlocal(i)
          end if
!.. explosion
          end do
          RHS(i) = RHS(i) + Gpsurob_Glo(i)
        end do
            
        !*******  Boundary contribution to ViscoF  ***********************************
        IntGWZrm1dV = BoundaryDataTab(ibdp)%BoundaryIntegral(7)
        IntdWrm1dV = BoundaryDataTab(ibdp)%BoundaryIntegral(3)
        alfaMon = Med(mati)%alfaMon
        if (alfaMon > zero .or. cinvisci > zero) then
!AA504 sub start: the molecular viscosity term, depending on velocity divergence, can be neglected (otherwise it may cause several problems); 
! further Monaghan's term is now activated even for detaching particles
            celeri = Med(mati)%celerita
            Monvisc = alfaMon * celeri * Domain%h 
!AA504 sub end
          Mmult = -Monvisc * dvn * IntGWZrm1dV
          ViscoMon(:) = ViscoMon(:) + Mmult * nnlocal(:)
!.. explosion
          do i = 1, SPACEDIM
            ViscoMon_save1 = ViscoMon_save1 + ViscoMon(i) * dvn * nnlocal(i)
          end do
!.. explosion
        end if
        if (cinvisci > zero) then
          cinviscmult = two * cinvisci * IntdWrm1dV * Tratto(stretch)%ShearCoeff
          ViscoShear(:) = ViscoShear(:) + cinviscmult * dvij(:)
        end if
      end if
!
    else if (stretchtype == "velo") then 
!
!AA404 sub
      if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
        pg(npi)%kodvel = 2
        if (tempo < Tratto(stretch)%trampa) then
          pg(npi)%velass(:) = Tratto(stretch)%NormVelocity*nnlocal(:) * tempo / Tratto(stretch)%trampa
        else
          pg(npi)%velass(:) = Tratto(stretch)%NormVelocity*nnlocal(:) 
        end if
      end if
      return   
! 
    else if (stretchtype == "flow") then 
!
!AA404 sub
       if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
        pg(npi)%kodvel = 2
! calcolo velocita' normale
        if (BoundaryFace(iface)%CloseParticles > 0) then
          minquotanode = 9999.0d0
          maxquotanode = const_m_9999
          do i = 1,BoundaryFace(iface)%nodes
            if (minquotanode > vertice(3,BoundaryFace(iface)%node(i)%name)) minquotanode = vertice(3,BoundaryFace(iface)%node(i)%name)
            if (maxquotanode < vertice(3,BoundaryFace(iface)%node(i)%name)) maxquotanode = vertice(3,BoundaryFace(iface)%node(i)%name)
          end do
          Lb = BoundaryFace(iface)%CloseParticles_maxQuota - minquotanode
          L = maxquotanode - minquotanode
          FlowRate1 = Tratto(stretch)%FlowRate * Lb / L
          Tratto(stretch)%NormVelocity = doubleh * FlowRate1 / (BoundaryFace(iface)%CloseParticles * Domain%PVolume)
        else
          Tratto(stretch)%NormVelocity = zero
        end if
        if (tempo < Tratto(stretch)%trampa) then
          pg(npi)%velass(:) = Tratto(stretch)%NormVelocity*nnlocal(:) * tempo / Tratto(stretch)%trampa
        else
          pg(npi)%velass(:) = Tratto(stretch)%NormVelocity*nnlocal(:) 
        end if
      end if
      return   
! 
!.. for the source boundary that assume three velocity components
!                                                        
    else if (stretchtype == "sour") then         
!
!AA404 sub
      if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
        pg(npi)%kodvel = 2
        if (tempo < Tratto(stretch)%trampa) then
          pg(npi)%velass(:) = Tratto(stretch)%NormVelocity*nnlocal(:) * tempo / Tratto(stretch)%trampa
        else
          pg(npi)%velass(:) = Tratto(stretch)%NormVelocity*nnlocal(:) 
        end if
      end if
      return                                                           
    end if
!
!!    end if
!
  end do face_loop

!*******  Addition of the boundary contribution to gradP  ***********************************
!.. Boundary contributions to momentum equation  
!
!!!!if (it_corrente ==1 .and. (npi== 10447 .or. npi==10425)) then
!!!!write (99,*) 'me3d',tpres,tdiss,tvisc
!!!!write (99,*) ' '
!!!!write (99,*) 'me3d',RHS,ViscoMon,ViscoShear
!!!!end if
!
  tpres(:) = tpres(:) - RHS(:)
  tdiss(:) = tdiss(:) - ViscoMon(:)
!AA504 rm and comm: this 3D boundary term has not been tested and seems not to work properly (it is erased at the moment)
! It seems useless at this stage to comment all the other lines involved as they are sparse and do not cause relevant computational time.
!  tvisc(:) = tvisc(:) - ViscoShear(:)
!............................................ 2011 mar 08
!.. contribution for Specific Internal Energy
!
  if (esplosione) then
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
  end subroutine AddBoundaryContributions_to_ME3D
!---split

!cfile AddElasticBoundaryReaction_2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : AddElasticBoundaryReaction_2D
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!
!************************************************************************************
! Module purpose : Module to compute the boundary integral IntWdS
!
! Calling routine: Loop_Irre_2D
!
! Called routines: 
!
!************************************************************************************
!
subroutine AddElasticBoundaryReaction_2D (npi, Ncbs, BoundReaction)

! Adds supplementariìy normal boundary reaction to reinforce insufficient
! pressure gradient in case of few neighbouring particles and presence of
! normal component of mass force (gravity)
! The normal reaction is computed with the formula R=(c0^2/d) ln(zi/d) [for zi<d],
! stemming from the compressible reaction of the fluid, where:
! c0^2 = E/ro0 is the square celerity of the fluid;
! zi is the distance of the particle Pi from the boundary face
! d is a reference distance from which the reaction is added

!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
integer(4),      intent(IN)    :: npi, Ncbs
double precision,intent(INOUT),dimension(1:SPACEDIM) :: BoundReaction
!
!.. Local Parameters ..
double precision,parameter :: ymincoeff = 0.25d0
double precision,parameter :: reafactor = 1.0d0
!
!.. Local Scalars ..
integer(4) :: sd, icbs, iside, nt, ibdt, ibdp, mate
double precision :: xpi, ypi, ypimin, celer02, vin, normreact
!
!.. Executable Statements ..
!
  mate = pg(npi)%imed
  ypimin = ymincoeff * Domain%dd
  celer02 = Med(mate)%eps / Med(mate)%den0
!
  ibdt = BoundaryDataPointer(3,npi)
  do icbs = 1, Ncbs
!
    ibdp = ibdt + icbs - 1
    iside = BoundaryDataTab(ibdp)%CloBoNum
    nt = BoundarySide(iside)%stretch
!
    if (Tratto(nt)%tipo == "fixe" .or. Tratto(nt)%tipo == "tapi") then
!
      xpi = BoundaryDataTab(ibdp)%LocXYZ(1)
      ypi = BoundaryDataTab(ibdp)%LocXYZ(2)
!
      if (ypi < ypimin) then
!
        if (xpi > zero .and. xpi < BoundarySide(iside)%Length) then  !La proiezione normale della particella npi sul piano
                                                                     ! del lato "iside" è interna al lato "iside"
          vin = zero
          do sd = 1, SPACEDIM
            vin = vin + pg(npi)%var(sd) * BoundarySide(iside)%T(sd, 3)
          end do
          if (vin < zero) then
!            normreact = -reafactor * celer02 * DLog(ypi / ypimin) / ypimin
            normreact = -reafactor * celer02 * DLog((Domain%h + ypi - ypimin) / Domain%h) / Domain%h
            do sd = 1, SPACEDIM
              BoundReaction(sd) = BoundReaction(sd) + normreact * BoundarySide(iside)%T(sd, 3)
            end do
          end if
        end if
      end if
    end if
  end do
!
!prova!
!  BoundReaction = floor(BoundReaction * azzeramento) / azzeramento
!  where (dabs(BoundReaction) < arrotondamento) BoundReaction = zero
!prova!
!
return 
end subroutine AddElasticBoundaryReaction_2D
!---split

!cfile AddElasticBoundaryReaction_3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : AddElasticBoundaryReaction_3D
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!
!************************************************************************************
! Module purpose : Module to compute the boundary integral IntWdS
!
! Calling routine: AddBoundaryContribution_to_CE3D
!                  AddBoundaryContributions_to_ME3D
!
! Called routines: 
!
!************************************************************************************
!
subroutine AddElasticBoundaryReaction_3D (npi, Ncbf, BoundReaction)

! Adds supplementari normal boundary reaction to reinforce insufficient
! pressure gradient in case of few neighbouring particles and presence of
! normal component of mass force (gravity)
! The normal reaction is computed with the formula R=(c0^2/d) ln(zi/d) [for zi<d],
! stemming from the compressible reaction of the fluid, where:
! c0^2 = E/ro0 is the square celerity of the fluid;
! zi is the distance of the particle Pi from the boundary face
! d is a reference distance from which the reaction is added

!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
integer(4),      intent(IN)    :: npi, Ncbf
double precision,intent(INOUT),dimension(1:SPACEDIM) :: BoundReaction
!
!.. Local Parameters ..
double precision,parameter :: zmincoeff = 0.25d0
double precision,parameter :: reafactor = 1.0d0
!
!.. Local Scalars ..
integer(4) :: sd, icbf, iface, nt, ibdt, ibdp, mate, fkod
integer(4) :: ne, NCloseEdgeF
double precision :: zi, zimin, celer02, vin, normreact
double precision :: scaprod, edgelen2, tau, edgedist2
!
!.. Local Arrays ..
double precision,dimension(1:SPACEDIM)   :: PXLoc, csi
double precision,dimension(1:SPACEDIM)   :: XQ, QP, QPcosdir
logical, dimension(1:Domain%MAXCLOSEBOUNDFACES) :: ReaFace
!
! External functions and subrotuines
logical, external    :: IsPointInternal
!
!.. Executable Statements ..
!
  BoundReaction = zero
  mate = pg(npi)%imed
  zimin = zmincoeff * Domain%dd
  celer02 = Med(mate)%eps / Med(mate)%den0
!
  ibdt = BoundaryDataPointer(3,npi)
  do icbf = 1, Ncbf
!
    ibdp = ibdt + icbf - 1
    iface = BoundaryDataTab(ibdp)%CloBoNum
    nt = BoundaryFace(iface)%stretch
!
    ReaFace(icbf) = .false.
!
    if (Tratto(nt)%tipo == "fixe" .or. Tratto(nt)%tipo == "tapi") then
!
      PXLoc(:) = BoundaryDataTab(ibdp)%LocXYZ(:)
      zi = PXLoc(3)
!
      if (zi < zimin) then
!
        call LocalNormalCoordinates (PXLoc, csi, iface)
        fkod = BoundaryFace(iface)%nodes - 2
!
        if (IsPointInternal(fkod, csi)) then    !La proiezione normale della particella npi sul piano
                                                ! della faccia "iface" è interna alla faccia "iface"
          vin = zero
          do sd = 1,SPACEDIM
            vin = vin + pg(npi)%var(sd) * BoundaryFace(iface)%T(sd, 3)
          end do
          if (vin < zero) then
!!!            normreact = -reafactor * celer02 * DLog(zi / zimin) / zimin
            normreact = -reafactor * celer02 * DLog((Domain%h + zi - zimin) / Domain%h) / Domain%h
            BoundReaction(:) = BoundReaction(:) + normreact * BoundaryFace(iface)%T(:, 3)
          end if
          ReaFace(icbf) = .true.
        else
          ReaFace(icbf) = .false.
        end if
      end if
    end if
  end do
!
!.. Reazione da eventuali spigoli (edges) vicini
  do ne = 1, NumBEdges
    !.. Verifica se la particella è vicina ad almeno una delle facce comuni allo spigolo ne
    NCloseEdgeF = 0
    
    ibdt = BoundaryDataPointer(3,npi)
    do icbf = 1, Ncbf
!
      ibdp = ibdt + icbf - 1
      iface = BoundaryDataTab(ibdp)%CloBoNum
      if (iface == BoundaryConvexEdge(ne)%face(1) .and. .Not. ReaFace(icbf)) then
        NCloseEdgeF = NCloseEdgeF + 1
      else if (iface == BoundaryConvexEdge(ne)%face(2) .and. .Not. ReaFace(icbf)) then
        NCloseEdgeF = NCloseEdgeF + 1
      end if
    end do
    if (NCloseEdgeF /= 2) return
    
    !..Distanza zi della particella npi dallo spigolo ne
    scaprod = zero
    do sd = 1,SPACEDIM
      scaprod = scaprod + (pg(npi)%Coord(sd) - BoundaryConvexEdge(ne)%node(1)%GX(sd)) * BoundaryConvexEdge(ne)%component(sd)
    end do
    edgelen2 = BoundaryConvexEdge(ne)%length * BoundaryConvexEdge(ne)%length
    tau = scaprod / edgelen2
    if (tau >= zero .and. tau <= one) then
      edgedist2 = zero
      XQ(:) = BoundaryConvexEdge(ne)%node(1)%GX(:) + BoundaryConvexEdge(ne)%component(:) * tau
      QP(:) = pg(npi)%Coord(:) - XQ(:)
      do sd = 1,SPACEDIM
        edgedist2 = edgedist2 + QP(sd) * QP(sd)
      end do
      zi = Dsqrt(edgedist2)
      if (zi < zimin) then
        vin = zero
        QPcosdir(:) = QP(:) / zi
        do sd = 1,SPACEDIM
          vin = vin + pg(npi)%var(sd) * QPcosdir(sd)
        end do
        if (vin < zero) then
          normreact = -reafactor * celer02 * DLog(zi / zimin) / zimin
          BoundReaction(:) = BoundReaction(:) +  normreact * QPcosdir(:)
        end if
      end if
    end if
  end do
!
!prova!
!  BoundReaction = floor(BoundReaction * azzeramento) / azzeramento
!  where (dabs(BoundReaction) < arrotondamento) BoundReaction = zero
!prova!
!
return 
end subroutine AddElasticBoundaryReaction_3D
!---split

