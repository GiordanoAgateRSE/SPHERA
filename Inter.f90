!cfile inter_EqMoto.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : inter_EqMoto
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Amicarelli/Agate  05/07/2011     Time stage parameters
! 04  Amicarelli/Agate  30Nov11        BSPH: wall element contributions,
!                                            density and velocity gradients
!AA501b
! 05  Amicarelli-Agate  13nov12        Body dynamics
!AA504
! 06Amicarelli          08Apr14        (v5.04) Modifications for granular flows
!
!************************************************************************************
!AA406 sub
! Module purpose : Computation of the terms of the momentum equation 
!                  (BSPH: Shepard's coefficient and gravity are applied later on)
!                  and the energy equation. 
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: viscomon
!                  viscomorris
!AA406
!                  viscomon_wall_elements (BSPH)
!                  viscomorris_wall_elements (BSPH)
!
!************************************************************************************
!
subroutine inter_EqMoto (npi,tpres,tdiss,tvisc)
! ex inter2
!* implementa il meccanismo di ricerca delle particelle che agiscono
!* su quella i-esima.
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
!integer(4),parameter :: local_d = 500  ! num max part entro 2h
!
!.. Formal Arguments ..
integer(4),intent(IN) :: npi
double precision,intent(INOUT),dimension(1:SPACEDIM) :: tpres, tdiss, tvisc
!
!.. Local Scalars ..
!
integer(4)       :: npj,contj,npartint
!
!!!double precision :: unity,pesoj
double precision :: rhoi,rhoj,amassj
! quantita' step2 eq del moto
double precision :: pi,pj,alpha,rvw(3),rvwalfa(3),rvwbeta(3)  !,factvis
double precision :: veln,velti,veltj,deltan,pre,coeff,secinv,nupa,nu,modderveln,moddervelt,moddervel
double precision :: dvtdn,denorm
!
!.. Local Arrays ..
! quantita' step eq continuita'
double precision,dimension(3) :: dervel, dervelmorr, appopres, appodiss  !, appovisc      ! dvel per eq standard, dvar usa vel smoot di monaghan
!
!AA402 start
integer(4)       :: index_rij_su_h
double precision :: rij_su_h,ke_coef,kacl_coef,rij_su_h_quad,vol_Shep,Ww_Shep
double precision :: rijtemp,rijtemp2,ragtemp(3)
double precision :: gradmod,gradmodwacl,wu,denom
!AA402 end
!
!.. Executable Statements ..
!
!* azzeramento quantita generali
 pg(npi)%dEdT  = zero
 tpres(:)      = zero
 tdiss(:)      = zero
 tvisc(:)      = zero
 dervel(:)     = zero
 dervelmorr(:) = zero
 deltan        = 1.d+07
!
 ke_coef = Domain%coefke / Domain%h
 kacl_coef = Domain%coefkacl / Domain%h
!
!!! unity         = zero
!
!*_______________________________________________________________
!*prima passata per trovare la distanza della parete fissa se c'è
!AA406 comm
! This fix boundary is intended to be active when the semi-analytic approach works
! It is composed by particles whose movement is imposed
 if ( Domain%Slip ) then
!
   do contj = 1, nPartIntorno(npi)
!
     npartint = (npi-1)* NMAXPARTJ + contj
     npj = PartIntorno(npartint)
!
     if (pg(npj)%vel_type == "std") cycle
!
     denorm = max(zero,(rag(1,npartint)*pg(npj)%zer(1)+rag(2,npartint)*pg(npj)%zer(2)+rag(3,npartint)*pg(npj)%zer(3)))
     deltan = min(deltan,denorm)
!
   end do
!
 end if 
!
!.. prima passata calcolo delle interazioni
!
 do contj = 1, nPartIntorno(npi)
!
   npartint = (npi-1)* NMAXPARTJ + contj
   npj = PartIntorno(npartint)
!
   if (npi == npj) cycle   ! da verificare se e' da aggiungere anche nel ciclo precedente.
!
!AA402 start
!
!AA406 sub
!  kernel computations (parts from the subroutine CalcVarLength): only for intermediate time stages
!
   if (Domain%time_stage > 1) then
!
!.. calcola le componenti deltaX, deltaY e deltaZ della distanza tra le due particelle (corrente ed adiacente)
     ragtemp(1:3) = pg(npi)%coord(1:3) - pg(npj)%coord(1:3)
!
!.. se una sola delle componenti supera la sfera di influenza di una particella, scarta la coppia (non interagiscono)
     if (abs(ragtemp(1)) > doubleh .or. abs(ragtemp(2)) > doubleh .or. abs(ragtemp(3)) > doubleh ) cycle
!
!.. calcola la distanza effettiva (considera direttamente il quadrato per aumentare l'accuratezza)
     rijtemp = ragtemp(1)*ragtemp(1) + ragtemp(2)*ragtemp(2) + ragtemp(3)*ragtemp(3)
!
!..se il valore della distanza tra le particelle e' maggiore del raggio di influenza, scarta la coppia
      if (rijtemp > doublesquareh) cycle
!
!.. memorizzo le distanze deltax,deltay e deltaz e la distanza totale rij
     rijtemp2 = rijtemp
     rijtemp = Dsqrt(rijtemp)
     rij_su_h = rijtemp / Domain%h
     rij_su_h_quad = rijtemp2 / squareh
     index_rij_su_h = int(rij_su_h)
     denom = one / (rijtemp + eta)
     rag(1:3,npartint) = ragtemp(1:3)
!
! kernel computations
     gradmod = zero
     gradmodwacl = zero
     wu = zero
     PartKernel(1:4,npartint) = zero
     if (index_rij_su_h >= 2) cycle
     gradmod = -two * rij_su_h + 1.5d0 * rij_su_h_quad
     gradmodwacl = -12.0d0 - 3.0d0 * rij_su_h_quad + 12.0d0 * rij_su_h 
     wu = 0.666666667d0 + rij_su_h_quad * (rij_su_h * half - one)
     if (index_rij_su_h > 0) then
       gradmod = -gradmod + rij_su_h_quad - two
       wu = (two - rij_su_h) * (two - rij_su_h) * (two - rij_su_h) * 0.166666666667d0
     end if 
     gradmod = gradmod * ke_coef
     gradmodwacl = gradmodwacl * kacl_coef
     PartKernel(1,npartint) = gradmod * denom 
     PartKernel(2,npartint) = PartKernel(1,npartint) / (rijtemp2 + eta2)
     PartKernel(3,npartint) = gradmodwacl * denom
     PartKernel(4,npartint) = wu * Domain%coefke
   end if
!AA402 end
!
   rhoi   = pg(npi)%dens
   rhoj   = pg(npj)%dens
   amassj = pg(npj)%mass
   pi = pg(npi)%pres
   pj = pg(npj)%pres
!
   dervel(:)     = pg(npj)%vel(:) - pg(npi)%vel(:)
   dervelmorr(:) = pg(npj)%vel(:) - pg(npi)%vel(:)
!
   if ( pg(npj)%vel_type /= "std" ) then          !non part fix o altro
     pj = pi
     rhoj   = rhoi
     amassj = pg(npi)%mass
     moddervel = -two* (pg(npi)%vel(1)*pg(npj)%zer(1)+pg(npi)%vel(2)*pg(npj)%zer(2)+pg(npi)%vel(3)*pg(npj)%zer(3))
     dervel(:) = moddervel * pg(npj)%zer(:)   !+pg(npj)%vstart(:)
     if (pg(npj)%slip == "f") then
       dervelmorr(:) = dervel(:)
     else if (pg(npj)%slip == "c") then
       denorm = max(zero,abs(rag(1,npartint)*pg(npj)%zer(1)+rag(3,npartint)*pg(npj)%zer(3)))
!       denorm = max(zero,abs(rag(1,contj,npi)*pg(npj)%zer(1)+rag(3,contj,npi)*pg(npj)%zer(3)))
       veln   = (pg(npi)%vel(1)*pg(npj)%zer(1)+pg(npi)%vel(3)*pg(npj)%zer(3))
       velti  = (pg(npi)%vel(1)*pg(npj)%zer(3)-pg(npi)%vel(3)*pg(npj)%zer(1))
       secinv = abs(velti/(deltan+0.0001d0))
       nu     = Med(pg(npi)%imed)%visc
       if (index(Med(pg(npi)%imed)%tipo,"liquid") > 0) then
         nu    = Med(pg(npi)%imed)%visc
       else if (index(Med(pg(npi)%imed)%tipo,"gas") > 0) then
         nu    = Med(pg(npi)%imed)%visc
       else if (index(Med(pg(npi)%imed)%tipo,"general") > 0) then
         nupa  = Med(pg(npi)%imed)%taucri/(secinv+0.0001d0)+Med(pg(npi)%imed)%visc*((secinv+0.0001d0)**(Med(pg(npi)%imed)%cuin-one))
         nu    = min(Med(pg(npi)%imed)%numx,nupa)
       else if (index(Med(pg(npi)%imed)%tipo,"granular") > 0) then
         pre   = (max(zero,pg(npi)%pres))/pg(npi)%dens
         coeff = sin (Med(pg(npi)%imed)%phi)
         nupa  = (pre*coeff)/(secinv+0.0001d0)+Med(pg(npi)%imed)%visc
         nu    = min(nupa,Med(pg(npi)%imed)%numx)
       end if
       dvtdn = (sin (pg(npj)%ang))*(pg(npi)%dudy+pg(npi)%dvdx)+(cos (pg(npj)%ang))*(pg(npi)%dudx-pg(npi)%dvdy)
       veltj = (-two*(deltan/(denorm+0.0001d0))*nu/(pg(npi)%visc+0.0001d0)+one)*velti +two*dvtdn*deltan
       moddervelt = veltj - velti
       modderveln = -two* veln               
       dervelmorr(1) = modderveln*pg(npj)%zer(1)+moddervelt*pg(npj)%zer(3)
       dervelmorr(3) = modderveln*pg(npj)%zer(3)-moddervelt*pg(npj)%zer(1)
     else if (pg(npj)%slip == "n") then
       dervelmorr(:) = - pg(npi)%vel(:)
     end if
   end if

!=============  EQUAZIONE DEL MOTO  ================

!* calcolo unita'
!
!   if ( pg(npj)%vel_type /= "std" .OR. pg(npj)%imed /= pg(npi)%imed ) then
!     pj = 1.2d0 * pj
!   end if
!   alpha = (pi+pj)/(rhoi*rhoj) 
!

!AA504 sub start
     if (Granular_flows_options%ID_erosion_criterion.ne.1) then
        alpha = pi / (rhoi*rhoi) + pj / (rhoj*rhoj) 
        else
           alpha = (pi+pj)/rhoi 
     endif
!AA504 sub end
                      
!
!AA406 sub start
   if (Domain%tipo == "semi") then
!AA504 sub start
     if (Granular_flows_options%ID_erosion_criterion.ne.1) then
        appopres(:) = ( -amassj * alpha * rag(:,npartint) * PartKernel(3,npartint) )
        else
         appopres(:) = ( -amassj/rhoj * alpha * rag(:,npartint) * PartKernel(3,npartint) )
     endif
!AA504 sub end
      else
! BSPH: the equation has to use the same kernel type when using the kernel function (cubic spline) and the derivative
         if (Domain%tipo == "bsph") then
            appopres(:) = ( -amassj * alpha * rag(:,npartint) * PartKernel(1,npartint) ) 
         endif
   endif
!AA406 sub end
!
   tpres(:) = tpres(:) + appopres(:)
!
!!!   pesoj = amassj * PartKernel(4,npartint) / rhoj
!
!!!   if (pg(npj)%vel_type /= "std") unity = unity + pesoj  !??????????????
!
   call viscomon (npi,npj,npartint,dervel,rvwalfa,rvwbeta)
!
   appodiss(:) = rvwalfa(:) + rvwbeta(:)
   tdiss(:) = tdiss(:) + appodiss(:)   ! termine dissipativo di Monaghan
!
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!if (it_corrente == 2 .and. npi > 16900 .and. npi < 17200) then
!  write (99,*) ' it = ',it_corrente,' n. particella = ',npi,npj
!  write (99,'(a,9d23.15)') ' inter_eqMoto   pg(npi)%  vel,var,dens,mass,pres ',pg(npi)%vel,pg(npi)%var,pg(npi)%dens,pg(npi)%mass,pg(npi)%pres
!  write (99,'(a,9d23.15)') ' inter_eqMoto   dervel,appodiss,tdiss ',dervel,appodiss,tdiss
!end if
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!
!???????????????????? j > nag   dervelmorr non utilizzato
!   if ( npj > nag .AND. pg(npj)%velmorr(1) /= zero) dervelmorr(1) = pg(npj)%velmorr(1) - pg(npi)%vel(1)
!
   call viscomorris (npi,npj,npartint,dervel,rvw)
!
   tvisc(:) = tvisc(:) + rvw(:)
!
!============= FINE EQUAZIONE DEL MOTO  ================
!
!........................................... 2011 mar 08
!.. calculation for Specific Internal Energy
   if (esplosione) &
     pg(npi)%dEdT = pg(npi)%dEdT + half * ( dervel(1) * (appopres(1)+appodiss(1)) + &
                  dervel(2)*(appopres(2)+appodiss(2)) + dervel(3)*(appopres(3)+appodiss(3)) )
!.. end calculation for Specific Internal Energy
!...............................................
!
 end do
!
!AA406 start
!Boundary contributions (BSPH)
if ((Domain%tipo == "bsph").and.(DBSPH%n_w > 0)) then
    do contj=1,nPartIntorno_fw(npi)
       npartint = (npi-1)* NMAXPARTJ + contj
       npj = PartIntorno_fw(npartint)
!
!AA501test
!       if (pg_w(npj)%wet == 1) then
!
       dervel(:) = pg_w(npj)%vel(:) - pg(npi)%vel(:)
!
!AA406test
!      if (nPartIntorno(npi) > 0) then
         appopres(:) = -pg_w(npj)%dens * pg_w(npj)%weight * kernel_fw(1,npartint) &
                       *(pg(npi)%pres/(pg(npi)%dens*pg(npi)%dens)+pg_w(npj)%pres/(pg_w(npj)%dens*pg_w(npj)%dens)) &
                       *pg_w(npj)%normal(:)
!         else 
!            appopres(:) = -pg_w(npj)%dens * pg_w(npj)%weight * kernel_fw(3,npartint) &
!                          *(pg(npi)%pres/(pg(npi)%dens*pg(npi)%dens)+pg_w(npj)%pres/(pg_w(npj)%dens*pg_w(npj)%dens)) &
!                          *pg_w(npj)%normal(:)
!      endif
!       appopres(:) = -pg(npi)%dens * pg_w(npj)%weight * kernel_fw(1,npartint) &
!                     *(pg(npi)%pres/(pg(npi)%dens*pg(npi)%dens)+pg(npi)%pres/(pg(npi)%dens*pg(npi)%dens)) &
!                     *pg_w(npj)%normal(:)
!       appopres(:) = -kernel_fw(2,npartint)* pg_w(npj)%weight * kernel_fw(1,npartint) &
!                     *(pg(npi)%pres/(pg(npi)%dens*pg(npi)%dens)+kernel_fw(3,npartint)/(kernel_fw(2,npartint)*kernel_fw(2,npartint))) &
!                     *pg_w(npj)%normal(:)
!
!AA406 contributions from semi-particles
!AA406test
 !     if (nPartIntorno(npi) > 0) then
          appopres(:) = appopres(:) - pg_w(npj)%mass * (pg(npi)%pres/(pg(npi)%dens*pg(npi)%dens)+ &
                        pg_w(npj)%pres/(pg_w(npj)%dens*pg_w(npj)%dens)) * rag_fw(:,npartint) * kernel_fw(2,npartint)  
 !     else
 !         appopres(:) = appopres(:) - pg_w(npj)%mass * (pg(npi)%pres/(pg(npi)%dens*pg(npi)%dens)+ &
 !                       pg_w(npj)%pres/(pg_w(npj)%dens*pg_w(npj)%dens)) * rag_fw(:,npartint) * kernel_fw(4,npartint)
 !     endif
!
       tpres(:) = tpres(:) + appopres(:)
       call viscomon_wall_elements(npi,npj,npartint,dervel,rvwalfa,rvwbeta)
!AA406test
!       tdiss(:) = tdiss(:) + rvwalfa(:) + rvwbeta(:)
!       
       call viscomorris_wall_elements(npi,npj,npartint,dervel,rvw)
       tvisc(:) = tvisc(:) + rvw(:)
!
!AA501test
!       endif
!
    end do
 endif
!AA406 end
!
!!! pg(npi)%uni = unity
!
!prova!
!tpres = floor(tpres * azzeramento) / azzeramento
!where (dabs(tpres) < arrotondamento) tpres = zero
!tdiss = floor(tdiss * azzeramento) / azzeramento
!where (dabs(tdiss) < arrotondamento) tdiss = zero
!tvisc = floor(tvisc * azzeramento) / azzeramento
!where (dabs(tvisc) < arrotondamento) tvisc = zero
!prova!
!
return
end subroutine inter_EqMoto
!---split

!cfile inter_EqCont_2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : inter_EqCont_2D
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Amicarelli/Agate  30Nov11        BSPH: wall element parameters and subroutines
!AA501b
! 04  Amicarelli-Agate  13nov12        Body dynamics
!AA504
! 06Amicarelli          08Apr14        (v5.04) Modifications for granular flows; second invariant of the strain rate tensor is computed only when necessary.
!
!************************************************************************************
! Module purpose : Module to accumulate the contributions of density variation in the
!                  continuity equation of the particles that are in the sphere of
!                  influence of the particle considered for the 2D case. The velocity 
!                  spatial derivatives and the second invariant of the deformation   
!                  velocity tensor are calculated.
! 
! Calling routine: Loop_Irre_2D
!
! Called routines: diffumorris
!
!************************************************************************************
!
subroutine inter_EqCont_2D 
! ex inter32d
!* implementa il meccanismo di ricerca delle particelle che agiscono
!* su quella i-esima.
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Local Parameters ..
!integer(4), parameter :: local_d = 500  ! num max part entro 2h
!
!.. Formal Arguments ..
!
!.. Local Scalars ..
integer(4) :: npi,npj,contj,npartint    !!!,npar2h
double precision :: rhoi,rhoj,amassj,moddervel,det,appo       !!! unity,pesoj,  
!modello bifluido
double precision :: factdiff, rvw, dervol   !, tdiff, cuei
!modello bifluido
!!!logical :: fix
!
!.. Local Arrays ..
double precision,dimension(3) :: pesogradj
! quantita' step eq continuita'
double precision,dimension(3) :: dvar      ! dervel, dvel per eq standard, dvar usa vel smoot di monaghan
double precision,dimension(4) :: derspa, aij
!
!.. Executable Statements ..
!
!!! !$omp private(i,unity,npar2h,derspa,dvar,aij,contj,npartint,j,fix,rhoi,rhoj,amassj,moddervel,pesoj,pesogradj) &
!
!AA504 sub omp directives
!AA406 sub
!$omp parallel do default(none) &
!$omp private(npi,derspa,dvar,aij,contj,npartint,npj,rhoi,rhoj,amassj,moddervel,pesogradj) &
!$omp private(dervol,factdiff,rvw,det,appo) &
!$omp shared(nag,pg,nPartIntorno,NMAXPARTJ,PartIntorno,PartKernel,rag,diffusione) &
!$omp shared(pg_w,nPartIntorno_fw,PartIntorno_fw,DBSPH,kernel_fw,rag_fw,Domain,Granular_flows_options)
!
!.. loops on all the active particles
!
 do npi = 1,nag
!
   if ( pg(npi)%vel_type /= "std" .or. pg(npi)%cella == 0) cycle
!
!* azzeramento quantita generali
!!!   npar2h    = 0
!!!   unity     = zero
   pg(npi)%dden  = zero
   pg(npi)%diffu = zero
   derspa(:) = zero
   dvar(:)   = zero
   aij(:)    = zero
!*_______________________________________________________________
!*prima passata per trovare celle interagenti e memorizzazione
!
   do contj = 1, nPartIntorno(npi)
!
     npartint = (npi-1)* NMAXPARTJ + contj
     npj = PartIntorno(npartint)
!
!!!     fix = .FALSE.
!
     rhoi    = pg(npi)%dens
     rhoj    = pg(npj)%dens
     amassj  = pg(npj)%mass
     dvar(:) = pg(npj)%var(:) - pg(npi)%var(:)
!
     if ( pg(npj)%vel_type /= "std" ) then          !non part fix o altro
       rhoj    = rhoi
       amassj  = pg(npi)%mass
       dvar(:) = pg(npj)%vel(:) - pg(npi)%var(:)        !!! resta definiti cosi' per particelle con moto fix ma non ferme
       if ( pg(npj)%vel(1) == zero .AND. pg(npj)%vel(2) == zero .AND. pg(npj)%vel(3) == zero ) then   !!! fix = .TRUE.
!!!       if ( pg(npj)%slip == "f" .AND. fix ) then 
         if ( pg(npj)%slip == "f" ) then 
           moddervel = -two * (pg(npi)%var(1) * pg(npj)%zer(1) + pg(npi)%var(3) * pg(npj)%zer(3))
           dvar(:) = moddervel * pg(npj)%zer(:)
         end if
       end if
     end if

!=============  EQUAZIONE DI CONTINUITA' ================
!
!* calcolo unita'
!!!     pesoj = amassj * PartKernel(4,npartint) / rhoj
     pesogradj(1:3) = amassj * rag(1:3,npartint) * PartKernel(1,npartint) / rhoj
!!!     unity = unity + pesoj
!!!     if ( pg(npj)%vel_type /= "fix" ) npar2h = npar2h + 1
!
!*eq. di continuita' calcolo drho ++++++++++++++++++++++++++++++++++++++++++++++++
!     pg(npi)%dden = pg(npi)%dden - amassj * PartKernel(1,npartint) * &
!                  (dervel(1)*rag(1,npartint) + dervel(2)*rag(2,npartint) + dervel(3)*rag(3,npartint)) 
!AA504 sub start
     if (Granular_flows_options%ID_erosion_criterion.ne.1) then
        appo = amassj * PartKernel(1,npartint) * &
               (dvar(1)*rag(1,npartint) + dvar(2)*rag(2,npartint) + dvar(3)*rag(3,npartint))
        else
             appo = rhoi * PartKernel(1,npartint) * (amassj/rhoj) * &
                    (dvar(1)*rag(1,npartint) + dvar(2)*rag(2,npartint) + dvar(3)*rag(3,npartint))
     endif
!AA504 sub end

!AA501btest
     pg(npi)%dden = pg(npi)%dden - appo

!*calcolo termine diffusivo
! modello bifluido
     if (diffusione) then
       dervol = pg(npj)%VolFra - pg(npi)%VolFra
       call diffumorris (npi,npj,npartint,dervol,factdiff,rvw)
       pg(npi)%diffu = pg(npi)%diffu + factdiff * rvw
     end if
! modello bifluido
!
!*calcolo derivate spaziali 
     if ( pg(npj)%vel_type /= "std" ) cycle
     
!AA504
     if ((Granular_flows_options%ID_erosion_criterion==1).or.(Granular_flows_options%ID_erosion_criterion==3)) then     
        derspa(1) = derspa(1) + pesogradj(1) * dvar(1)
        derspa(2) = derspa(2) + pesogradj(3) * dvar(1)
        derspa(3) = derspa(3) + pesogradj(1) * dvar(3)
        derspa(4) = derspa(4) + pesogradj(3) * dvar(3)
        aij(1) = aij(1) + (pg(npj)%coord(1) - pg(npi)%coord(1)) * pesogradj(1)
        aij(2) = aij(2) + (pg(npj)%coord(3) - pg(npi)%coord(3)) * pesogradj(1)
        aij(3) = aij(3) + (pg(npj)%coord(1) - pg(npi)%coord(1)) * pesogradj(3)
        aij(4) = aij(4) + (pg(npj)%coord(3) - pg(npi)%coord(3)) * pesogradj(3)
!AA504
     endif      
   end do
!
!AA406 start
! boundary contributions (BSPH)
    if (Domain%tipo == "bsph") then
       do contj = 1,nPartIntorno_fw(npi)
          npartint = (npi-1)* NMAXPARTJ + contj
          npj = PartIntorno_fw(npartint)
          dvar(:) = pg_w(npj)%vel(:) - pg(npi)%var(:)
!AA406test occhio nche al subtest Gamma
          appo = -pg_w(npj)%dens * pg_w(npj)%weight * kernel_fw(1,npartint) * pg(npi)%Gamma * &
                (dvar(1) * pg_w(npj)%normal(1) + dvar(2) * pg_w(npj)%normal(2) + dvar(3) * pg_w(npj)%normal(3))
          appo = appo - pg_w(npj)%mass * kernel_fw(2,npartint) * &
                (dvar(1) * rag_fw(1,npartint) + dvar(2) * rag_fw(2,npartint) + dvar(3) * rag_fw(3,npartint))
          pg(npi)%dden = pg(npi)%dden + appo
       end do
    endif
!AA406 end
! 
!prova!
!   pg(npi)%dden = floor(pg(npi)%dden * azzeramento) / azzeramento
!   if (dabs(pg(npi)%dden) < arrotondamento) pg(npi)%dden = zero 
!prova!
!
!AA504
   if ((Granular_flows_options%ID_erosion_criterion==1).or.(Granular_flows_options%ID_erosion_criterion==3)) then
      det = aij(1) * aij(4) - aij(2) * aij(3)
!
      if (det < 0.001d0) det = 0.001d0
!
      pg(npi)%dudx = ( derspa(1)*aij(4) - derspa(2)*aij(2)) / det       ! derivate norm accoppiate
      pg(npi)%dudy = (-derspa(1)*aij(3) + derspa(2)*aij(1)) / det       ! derivate norm accoppiate
      pg(npi)%dvdx = ( derspa(3)*aij(4) - derspa(4)*aij(2)) / det       ! derivate norm accoppiate
      pg(npi)%dvdy = (-derspa(3)*aij(3) + derspa(4)*aij(1)) / det       ! derivate norm accoppiate
!
!ç   pg(npi)%dudy = half * (pg(npi)%dudy + pg(npi)%dvdx)
!ç   pg(npi)%dvdx = pg(npi)%dudy
!
      pg(npi)%secinv = Dsqrt(half * pg(npi)%dudx*pg(npi)%dudx + half * pg(npi)%dvdy*pg(npi)%dvdy + quarter * &
                       (pg(npi)%dudy + pg(npi)%dvdx) * (pg(npi)%dudy + pg(npi)%dvdx)) !SaMa: forse perche' dudy=dvdx
!ç   pg(npi)%secinv = Dsqrt(half * pg(npi)%dudx*pg(npi)%dudx + half * pg(npi)%dvdy*pg(npi)%dvdy + pg(npi)%dudy*pg(npi)%dudy)
!AA504
   endif
!
!!!   pg(npi)%npar2h = npar2h
!
!!!   pg(npi)%uni = unity
!
 end do  ! nag
!
!$omp end parallel do

!AA501b
 call start_and_stop(3,12)
 call start_and_stop(2,19)
 if (n_bodies > 0) call body_particles_to_continuity
 call start_and_stop(3,19)
 call start_and_stop(2,12)

!
return
end subroutine inter_EqCont_2D
!---split

!cfile inter_EqCont_3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : inter_EqCont_3D
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!AA501b
! 03  Amicarelli-Agate  13nov12        Body dynamics
!AA504
! 06Amicarelli          08Apr14        (v5.04) Modifications for granular flows; second invariant of the strain rate tensor is computed only when necessary.
!                                      Renormalization deactivated (because of errors). Correction of the formula for second invariant of the strain rate tensor. 
!
!************************************************************************************
!AA406 sub
! Module purpose : Module to accumulate the contributions of density variation in the
!
!                  continuity equation of the particles that are in the sphere of
!                  influence of the particle considered for the 3D case. The velocity 
!                  spatial derivatives and the second invariant of the deformation   
!                  velocity tensor are calculated.
!
! Calling routine: Loop_Irre_3D
!
! Called routines: diffumorris
!
!************************************************************************************
!
subroutine inter_EqCont_3D
! ex inter33d
!* implementa il meccanismo di ricerca delle particelle che agiscono
!* su quella i-esima.
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Local Parameters ..
!integer(4), parameter :: local_d = 500  ! num max part entro 2h
!
!.. Local Scalars ..
integer(4) :: npi,npj,contj,npartint    !!!,npar2h
double precision :: rhoi,rhoj,amassj,moddervel,det,moddia,modout,appo    !!! unity,pesoj,   
! quantita' step2 eq del moto
!double precision pi,pj
! modello bifluido
double precision :: factdiff, rvw, dervol   !, tdiff, cuei
! modello bifluido
!!!logical :: fix
!
!.. Local Arrays ..
double precision,dimension (3) :: pesogradj
! quantita' step eq continuita'
double precision,dimension (3) :: dvar     ! dervel, dvel per eq standard, dvar usa vel smoot di monaghan
double precision,dimension (9) :: derspa,aij,bij,dvdi
!
!.. Executable Statements ..
!
!!! !$omp private(i,unity,npar2h,derspa,dvar,aij,bij,dvdi,contj,npartint,j,fix,rhoi,rhoj,amassj,pi,pj,moddervel) &
!!! !$omp private(pesoj,pesogradj,dervol,factdiff,rvw,det,moddia,modout) &.
!pi,pj,
!AA504 sub omp directives
!$omp parallel do default(none) &
!$omp private(npi,derspa,dvar,aij,bij,dvdi,contj,npartint,npj,rhoi,rhoj,amassj,moddervel) &
!$omp private(pesogradj,dervol,factdiff,rvw,det,moddia,modout,appo) &
!$omp shared(nag,pg,nPartIntorno,NMAXPARTJ,PartIntorno,PartKernel,rag,diffusione,Granular_flows_options)
!
!.. loops on all the active particles
!
 do npi = 1,nag
!
   if ( pg(npi)%vel_type /= "std" .or. pg(npi)%cella == 0) cycle
!
!* azzeramento quantita generali
!!!   npar2h    = 0
!!!   unity     = zero
   pg(npi)%dden  = zero
   pg(npi)%diffu = zero
   derspa(:) = zero
   dvar(:)   = zero
   aij(:)    = zero
   bij(:)    = zero
   dvdi(:)   = zero
!
!*_______________________________________________________________
!*prima passata per trovare celle interagenti e memorizzazione
!
   do contj = 1, nPartIntorno(npi)
!
     npartint = (npi-1)* NMAXPARTJ + contj
     npj = PartIntorno(npartint)
!
!!!     fix = .FALSE.
!
     rhoi   = pg(npi)%dens
     rhoj   = pg(npj)%dens
     amassj = pg(npj)%mass
!     pj     = pg(npj)%pres
!     pi     = pg(npi)%pres
!
     dvar(:) = pg(npj)%var(:) - pg(npi)%var(:)
!
     if ( pg(npj)%vel_type /= "std" ) then          !non part fix o altro
!       pj    = pi
       rhoj    = rhoi
       amassj  = pg(npi)%mass
       dvar(:) = pg(npj)%vel(:) - pg(npi)%var(:)        !!! resta definiti così per particelle con moto fix ma non ferme
       if ( pg(npj)%vel(1) == zero .AND. pg(npj)%vel(2) == zero .AND. pg(npj)%vel(3) == zero ) then  !!! fix = .TRUE.
!!!       if ( pg(npj)%slip == "f" .AND. fix ) then 
         if ( pg(npj)%slip == "f") then 
           moddervel = -two * (pg(npi)%var(1) * pg(npj)%zer(1) + pg(npi)%var(3) * pg(npj)%zer(3))
           dvar(:) = moddervel * pg(npj)%zer(:)
         end if
       end if
     end if
!
!=============  EQUAZIONE DI CONTINUITA' ================
!      
!* calcolo unita'
!!!     pesoj = amassj * PartKernel(4,npartint) / rhoj
     pesogradj(1:3) = amassj * rag(1:3,npartint) * PartKernel(1,npartint) / rhoj
!!!     unity = unity + pesoj  
!!!     if ( pg(npj)%vel_type /= "fix" ) npar2h = npar2h + 1

!*eq. di continuita' calcolo drho ++++++++++++++++++++++++++++++++++++++++++++++++
!     pg(npi)%dden = pg(npi)%dden - amassj * PartKernel(1,npartint) * &
!                  (dervel(1)*rag(1,npartint) + dervel(2)*rag(2,npartint) + dervel(3)*rag(3,npartint))
!AA504 sub start
     if (Granular_flows_options%ID_erosion_criterion.ne.1) then
        appo = amassj * PartKernel(1,npartint) * &
           (dvar(1)*rag(1,npartint) + dvar(2)*rag(2,npartint) + dvar(3)*rag(3,npartint))
     else
              appo = rhoi * PartKernel(1,npartint) * (amassj/rhoj) * &
                     (dvar(1)*rag(1,npartint) + dvar(2)*rag(2,npartint) + dvar(3)*rag(3,npartint))
     endif
!AA504 sub end
     pg(npi)%dden = pg(npi)%dden - appo
!
!*calcolo termine diffusivo
! modello bifluido
     if (diffusione) then
       dervol = pg(npj)%VolFra - pg(npi)%VolFra
       call diffumorris (npi,npj,npartint,dervol,factdiff,rvw)
       pg(npi)%diffu = pg(npi)%diffu + factdiff * rvw
     end if
! modello bifluido

!*calcolo derivate spaziali 
     if ( pg(npj)%vel_type /= "std" ) cycle 
!AA504
     if ((Granular_flows_options%ID_erosion_criterion==1).or.(Granular_flows_options%ID_erosion_criterion==3)) then
!AA504 rm start     
!     derspa(1) = derspa(1) + pesogradj(1) * dvar(1)      !du/dx  !! LA MATRICE DERSPA
!     derspa(2) = derspa(2) + pesogradj(2) * dvar(1)      !du/dy  !! è ORDINATA 
!     derspa(3) = derspa(3) + pesogradj(3) * dvar(1)      !du/dz  !! PER COLONNE
!     derspa(4) = derspa(4) + pesogradj(1) * dvar(2)      !dv/dx
!     derspa(5) = derspa(5) + pesogradj(2) * dvar(2)      !dv/dy
!     derspa(6) = derspa(6) + pesogradj(3) * dvar(2)      !dv/dz
!     derspa(7) = derspa(7) + pesogradj(1) * dvar(3)      !dw/dx
!     derspa(8) = derspa(8) + pesogradj(2) * dvar(3)      !dw/dy
!     derspa(9) = derspa(9) + pesogradj(3) * dvar(3)      !dw/dz
!!.. normalizzazione del gradiente di velocita'
!     aij(1) = aij(1) + (pg(npj)%coord(1) - pg(npi)%coord(1)) * pesogradj(1)   !a11
!     aij(2) = aij(2) + (pg(npj)%coord(2) - pg(npi)%coord(2)) * pesogradj(1)   !a12
!     aij(3) = aij(3) + (pg(npj)%coord(3) - pg(npi)%coord(3)) * pesogradj(1)   !a13
!     aij(4) = aij(4) + (pg(npj)%coord(1) - pg(npi)%coord(1)) * pesogradj(2)   !a21
!     aij(5) = aij(5) + (pg(npj)%coord(2) - pg(npi)%coord(2)) * pesogradj(2)   !a22
!     aij(6) = aij(6) + (pg(npj)%coord(3) - pg(npi)%coord(3)) * pesogradj(2)   !a23
!     aij(7) = aij(7) + (pg(npj)%coord(1) - pg(npi)%coord(1)) * pesogradj(3)   !a31
!     aij(8) = aij(8) + (pg(npj)%coord(2) - pg(npi)%coord(2)) * pesogradj(3)   !a32
!     aij(9) = aij(9) + (pg(npj)%coord(3) - pg(npi)%coord(3)) * pesogradj(3)   !a33
!AA504 rm end
!AA504 start
     dvdi(1) = dvdi(1) + pesogradj(1) * dvar(1)      !du/dx  
     dvdi(2) = dvdi(2) + pesogradj(2) * dvar(1)      !du/dy  
     dvdi(3) = dvdi(3) + pesogradj(3) * dvar(1)      !du/dz  
     dvdi(4) = dvdi(4) + pesogradj(1) * dvar(2)      !dv/dx
     dvdi(5) = dvdi(5) + pesogradj(2) * dvar(2)      !dv/dy
     dvdi(6) = dvdi(6) + pesogradj(3) * dvar(2)      !dv/dz
     dvdi(7) = dvdi(7) + pesogradj(1) * dvar(3)      !dw/dx
     dvdi(8) = dvdi(8) + pesogradj(2) * dvar(3)      !dw/dy
     dvdi(9) = dvdi(9) + pesogradj(3) * dvar(3)      !dw/dz
!AA504 end
!AA504
     endif 
!============= FINE EQUAZIONE DI CONTINUITA' ================
!
   end do
!
!prova!
!   pg(npi)%dden = floor(pg(npi)%dden * azzeramento) / azzeramento
!   if (dabs(pg(npi)%dden) < arrotondamento) pg(npi)%dden = zero 
!prova!
!AA504
    if ((Granular_flows_options%ID_erosion_criterion==1).or.(Granular_flows_options%ID_erosion_criterion==3)) then
!AA504 rm start (this renormalization presents several errors; it should be re-written)
!   bij(1) = (aij(5) * aij(9) - aij(6) * aij(8))
!   bij(2) =-(aij(4) * aij(9) - aij(7) * aij(6))
!   bij(3) = (aij(4) * aij(8) - aij(5) * aij(7))
!   bij(4) =-(aij(2) * aij(9) - aij(3) * aij(8))
!   bij(5) = (aij(1) * aij(9) - aij(7) * aij(7))
!   bij(6) =-(aij(1) * aij(8) - aij(2) * aij(2))
!   bij(7) = (aij(2) * aij(6) - aij(3) * aij(5))
!   bij(8) =-(aij(1) * aij(6) - aij(3) * aij(4))
!   bij(9) = (aij(1) * aij(5) - aij(2) * aij(4))
!!
!   det = aij(1)*aij(5)*aij(9)+aij(2)*aij(6)*aij(7)+aij(3)*aij(4)*aij(8) - &
!         aij(1)*aij(6)*aij(8)-aij(4)*aij(2)*aij(9)-aij(7)*aij(5)*aij(3)
!!
!   if (det < 0.001d0) det = 0.001d0
!!
!   dvdi(1) = (bij(1)*derspa(1) - bij(2)*derspa(2) + bij(3)*derspa(3)) / det     !dudx
!   dvdi(2) = (bij(1)*derspa(4) - bij(2)*derspa(5) + bij(3)*derspa(6)) / det     !dudy
!   dvdi(3) = (bij(1)*derspa(7) - bij(2)*derspa(8) + bij(3)*derspa(9)) / det     !dudz
!   dvdi(4) = (bij(4)*derspa(1) - bij(5)*derspa(2) + bij(6)*derspa(3)) / det     !dvdx
!   dvdi(5) = (bij(4)*derspa(4) - bij(5)*derspa(5) + bij(6)*derspa(6)) / det     !dvdy
!   dvdi(6) = (bij(4)*derspa(7) - bij(5)*derspa(8) + bij(6)*derspa(9)) / det     !dvdz
!   dvdi(7) = (bij(7)*derspa(1) - bij(8)*derspa(2) + bij(9)*derspa(3)) / det     !dwdx
!   dvdi(8) = (bij(7)*derspa(4) - bij(8)*derspa(5) + bij(9)*derspa(6)) / det     !dwdy
!   dvdi(9) = (bij(7)*derspa(7) - bij(8)*derspa(8) + bij(9)*derspa(9)) / det     !dwdz
!AA504 rm end
!
! modifica solo dei termini utili per tensore vel def
!   dvdi(2) = half * (dvdi(2) + dvdi(4))
!   dvdi(3) = half * (dvdi(3) + dvdi(7))
!   dvdi(6) = half * (dvdi(6) + dvdi(8))
! 
!   moddia  = (dvdi(1)*dvdi(1) + dvdi(5)*dvdi(5) + dvdi(9)*dvdi(9))
!   modout  = (dvdi(2)*dvdi(2) + dvdi(3)*dvdi(3) + dvdi(6)*dvdi(6))
!   pg(npi)%secinv = Dsqrt( half*moddia + modout)
! 
   moddia  = (dvdi(1) * dvdi(1) + dvdi(5) * dvdi(5) + dvdi(9) * dvdi(9))
!AA504 sub 
   modout  = ( (dvdi(2)+dvdi(4))*(dvdi(2)+dvdi(4)) + (dvdi(3)+dvdi(7))*(dvdi(3)+dvdi(7)) + (dvdi(6)+dvdi(8))*(dvdi(6)+dvdi(8)) )
! old formula
! modout  = ( (dvdi(2)*dvdi(4))*(dvdi(2)*dvdi(4)) + (dvdi(3)*dvdi(7))*(dvdi(3)*dvdi(7)) + (dvdi(6)*dvdi(8))*(dvdi(6)*dvdi(8)) )
   
   pg(npi)%secinv = Dsqrt( half * moddia + quarter * modout)
!AA504
   endif   
! 
!!!   pg(npi)%npar2h = npar2h
!
!!!   pg(npi)%uni = unity
!
 end do  ! nag
!
!$omp end parallel do
!

!AA501b
 call start_and_stop(3,12)
 call start_and_stop(2,19)
 if (n_bodies > 0) call body_particles_to_continuity
 call start_and_stop(3,19)
 call start_and_stop(2,12)

return
end subroutine inter_EqCont_3D
!---split

!cfile inter_SmoothVelo_2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : inter_SmoothVelo_2D
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code

!AA501b
! 03  Amicarelli-Agate  13nov12        Body particle contributions 

!
!************************************************************************************
! Module purpose : Module to calculate the corrective term of the velocity
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D, Euler, Heun
!

!AA501b modified
!AA601 sub
! Called routines: body_to_smoothing_vel,DBSPH_inlet_outlet
 
!
!************************************************************************************
!
subroutine inter_SmoothVelo_2D 
! ex inter42d
!* implementa il meccanismo di ricerca delle particelle che agiscono su quella npi-esima.
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Local Parameters ..
!integer(4), parameter :: local_d = 500  ! num max part entro 2h
!
!.. Local Scalars ..
  integer(4)       :: npi,i,ii,npj,contj,npartint
  integer(4)       :: icbs, Ncbs, IntNcbs, ibdt, ibdp, iside, sidestr   !Ncols, 
!AA406 sub
  double precision :: rhoi,rhoj,amassj,pesoj,  moddervel,unity
!
  double precision :: IntWdV  !xpmin,xpmax,interlen,xpi,ypi, deltai,sidelen, 
  character(4)     :: strtype
!
!.. Local Arrays ..
  integer(4),      dimension(1:PLANEDIM) :: acix
  double precision,dimension(1:PLANEDIM) :: sss, nnn, DVLoc, DVGlo, BCLoc, BCGlo
  double precision,dimension(1:PLANEDIM) :: IntLocXY
! quantita' step eq continuita'
  double precision,dimension(3) :: dervel     ! dvel per eq standard, dvar usa vel smoot di monaghan
  
!AA501b start
  double precision,dimension(:,:),allocatable :: dervel_mat
  double precision,dimension(:),allocatable :: unity_vec

  if (n_bodies > 0) then  
!Allocations
     allocate(dervel_mat(nag,3))
     allocate(unity_vec(nag))
!Initializations
     dervel_mat = 0.
     unity_vec = 0.
  endif
!AA501b end  
  
!.. Executable Statements ..
!
  acix(1) = 1
  acix(2) = 3
!

!AA406
  unity = zero
  
!AA501b
! Body particle contributions to pressure smoothing
  if (n_bodies > 0) then
     call start_and_stop(3,14)
     call start_and_stop(2,19)
     call body_to_smoothing_vel(dervel_mat,unity_vec)
     call start_and_stop(3,19)
     call start_and_stop(2,14)
  endif
!AA501b end   
  
!
!!! !$omp private(ii,npi,unity,contj,npartint,npj,rhoi,rhoj,amassj,dervel,moddervel,pesoj) &   !sidelen,
!AA501b modified
!$omp parallel do default(none) &
!$omp private(ii,npi,contj,npartint,npj,rhoi,rhoj,amassj,dervel,moddervel,pesoj) &
!$omp private(Ncbs,IntNcbs,ibdt,ibdp,icbs,IntLocXY,iside,sidestr,strtype,i) &
!$omp private(sss,nnn,DVGlo,DVLoc,IntWdV,BCLoc,BCGlo) &
!$omp shared(nag,Pg,Domain,Med,Tratto,acix,nPartIntorno,NMAXPARTJ,PartIntorno,PartKernel) &
!$omp shared(BoundaryDataPointer,BoundaryDataTab,BoundarySide,indarrayFlu,Array_Flu,esplosione,kernel_fw,unity,dervel_mat,unity_vec,n_bodies)
!
!.. loops on all the active particles in order to calculate the average velocity depending on 
!.. the smoothing coefficient of velocity ThetaV
!
!!!!!!  do npi = 1,nag
!!!!!!!
!!!!!!    if ( pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" .or. pg(npi)%state == "sol") cycle
!$$$$$
  do ii = 1,indarrayFlu
    npi = Array_Flu(ii)
!$$$$$
!
!* azzeramento quantita generali
    pg(npi)%var = zero
    pg(npi)%Envar = zero
!!!    unity   = zero
!
!.. contributo delle particelle interne
!
    do contj = 1, nPartIntorno(npi)
!
      npartint = (npi-1)* NMAXPARTJ + contj
      npj = PartIntorno(npartint)
!
      rhoi   = pg(npi)%dens
      rhoj   = pg(npj)%dens
      amassj = pg(npj)%mass
!
      dervel(:) = pg(npj)%vel(:) - pg(npi)%vel(:)
!
      if ( pg(npj)%vel_type /= "std" ) then          !non part fix o altro
        rhoj   = rhoi
        amassj = pg(npi)%mass
        moddervel = -two * (pg(npi)%vel(1)*pg(npj)%zer(1)+pg(npi)%vel(2)*pg(npj)%zer(2)+pg(npi)%vel(3)*pg(npj)%zer(3))
        dervel(:) = moddervel * pg(npj)%zer(:)    !+pg(npj)%vstart(:)
      end if
!!============= CORREZIONE VELOCITA' ===================
      if ( Med(pg(npj)%imed)%den0 /= Med(pg(npi)%imed)%den0 ) cycle
!
      pesoj = amassj * PartKernel(4,npartint) / rhoj
!
!AA406
      unity = unity + pesoj
!
!!!      unity = unity + pesoj  
      pg(npi)%var(:) = pg(npi)%var(:) + dervel(:) * pesoj   
! calcolo velocita' per modello diffusivo
!   pg(npj)%veldif(:) = pg(npj)%veldif(:) + pg(npj)%vel(:) * pesoj
!
!...................................... 2011 mar 08
!.. update of Specific Internal Energy
      if (esplosione) pg(npi)%Envar = pg(npi)%Envar + (pg(npj)%IntEn - pg(npi)%IntEn) * pesoj
!.. end update of Specific Internal Energy
!.........................................
!
    end do
!
!

!AA501b start
      if (n_bodies > 0) then
         pg(npi)%var(:) = pg(npi)%var(:) + dervel_mat(npi,:)
         unity = unity + unity_vec(npi)
      endif

!AA406 start
      if (Domain%tipo == "bsph") then
!AA406test
          pg(npi)%var(:) = pg(npi)%var(:)
!          pg(npi)%var(:) = pg(npi)%var(:) / unity
!AA601 sub start
! Impose boundary conditions at inlet and outlet sections (DB-SPH)
          call DBSPH_inlet_outlet(npi)
!AA601 sub end
     else
!AA406 end
!
!.. contributo dei contorni
!
    Ncbs = BoundaryDataPointer(1,npi)
    IntNcbs = BoundaryDataPointer(2,npi)
    ibdt = BoundaryDataPointer(3,npi)
!
    if (IntNcbs > 0) then
!
      do icbs = 1, IntNcbs
!
        ibdp = ibdt + icbs - 1
        IntLocXY(1:PLANEDIM) = BoundaryDataTab(ibdp)%LocXYZ(1:PLANEDIM)
        iside = BoundaryDataTab(ibdp)%CloBoNum
        sidestr = BoundarySide(iside)%stretch
        strtype = Tratto(sidestr)%tipo
!
        if (strtype == 'sour' .or. strtype == 'velo' .or. strtype == 'flow') then
          pg(npi)%var(:) = zero   
!!!          unity = one
          exit  
        end if
!
!        sidelen = BoundarySide(iside)%Length
        do i = 1, PLANEDIM
          sss(i) = BoundarySide(iside)%T(acix(i),1)
          nnn(i) = BoundarySide(iside)%T(acix(i),3)
          DVGlo(i) = two * (Tratto(sidestr)%velocity(acix(i)) - pg(npi)%vel(acix(i)))
        end do
        DVLoc(1) = sss(1) * DVGlo(1) + sss(2) * DVGlo(2)
        DVLoc(2) = nnn(1) * DVGlo(1) + nnn(2) * DVGlo(2)
!
!        IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(2)
        IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(3)
!
        if (strtype == 'fixe' .or. strtype == 'tapi') then
!
          BCLoc(1) = DVLoc(1) * IntWdV * Tratto(sidestr)%ShearCoeff
          BCLoc(2) = DVLoc(2) * IntWdV
          BCGlo(1) = sss(1) * BCLoc(1) + nnn(1) * BCLoc(2)
          BCGlo(2) = sss(2) * BCLoc(1) + nnn(2) * BCLoc(2)
!
          pg(npi)%var(1) = pg(npi)%var(1) + BCGlo(1)   
          pg(npi)%var(3) = pg(npi)%var(3) + BCGlo(2)   
!
!!!          unity = unity + IntWdV
!
        end if
      end do
    end if
!
!AA406
    endif
!
!
!!!    pg(npi)%uni = unity
!
!prova!
!    pg(npi)%var = floor(pg(npi)%var * azzeramento) / azzeramento
!    where (dabs(pg(npi)%var) < arrotondamento) pg(npi)%var = zero 
!prova!
!
  end do
!
!$omp end parallel do
!

!AA501b start
!Deallocations
  if (n_bodies > 0) then
     deallocate(dervel_mat)
     deallocate(unity_vec)
  endif
!AA501b end

return
end subroutine inter_SmoothVelo_2D
!---split

!cfile inter_SmoothVelo_3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : inter_SmoothVelo_3D
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code

!AA501b
! 03  Amicarelli-Agate  13nov12        Body dynamics 

!
!************************************************************************************
! Module purpose : Module to calculate the corrective term of the velocity
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D, Euler, Heun
!
!AA501b modified
!AA601 sub
! Called routines: body_to_smoothing_vel,DBSPH_inlet_outlet
!
!************************************************************************************
!
subroutine inter_SmoothVelo_3D
! ex inter43d
!* implementa il meccanismo di ricerca delle particelle che agiscono su quella npi-esima.
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
!integer(4), parameter :: local_d = 500  ! num max part entro 2h
!
!.. Local Scalars ..
integer(4)       :: npi,i,j,npj,contj,npartint, ibdt, ibdp, ii
integer(4)       :: Ncbf, icbf, iface, facestr
double precision :: rhoi,rhoj,amassj,pesoj,  moddervel !unity
double precision :: IntWdV
character(4)     :: strtype
!
!.. Local Arrays ..
double precision,dimension(1:SPACEDIM) :: DVLoc, DVGlo, BCLoc, BCGlo
double precision,dimension(1:SPACEDIM) :: LocX
! quantita' step eq continuita'
double precision,dimension(3) :: dervel     ! dvel per eq standard, dvar usa vel smoot di monaghan

!AA501btest start
  double precision,dimension(:,:),allocatable :: dervel_mat
  double precision,dimension(:),allocatable :: unity_vec
  double precision :: unity

  if (n_bodies > 0) then  
!Allocations
     allocate(dervel_mat(nag,3))
     allocate(unity_vec(nag))
!Initializations
     dervel_mat = 0.
     unity_vec = 0.
     unity = 0.
  endif
!AA501btest end  

!
!.. Executable Statements ..
!!! !$omp private(ii,npi,unity,contj,npartint,npj,rhoi,rhoj,amassj,dervel,moddervel,pesoj) &
!

!AA501b
! Body particle contributions to pressure smoothing
  if (n_bodies > 0) then
     call start_and_stop(3,14)
     call start_and_stop(2,19)
     call body_to_smoothing_vel(dervel_mat,unity_vec)
     call start_and_stop(3,19)
     call start_and_stop(2,14)
  endif
!AA501b end  

!AA501b modified
!AA501 sub
!$omp parallel do default(none) &
!$omp private(ii,npi,contj,npartint,npj,rhoi,rhoj,amassj,dervel,moddervel,pesoj) &
!$omp private(Ncbf,icbf,ibdt,ibdp,LocX,iface,facestr,strtype,DVGlo,DVLoc,i,j,IntWdV,BCLoc,BCGlo,unity) &
!$omp shared(nag,pg,Med,nPartIntorno,NMAXPARTJ,PartIntorno,PartKernel) &
!$omp shared(BoundaryDataPointer,BoundaryDataTab,BoundaryFace,Tratto,indarrayFlu,Array_Flu,esplosione,Domain,n_bodies,unity_vec,dervel_mat)
!
!.. loops on all the active particles in order to calculate the average velocity depending on 
!.. the smoothing coefficient of velocity ThetaV
!
!!!!!!  do npi = 1,nag
!!!!!!!
!!!!!!    if ( pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" .or. pg(npi)%state == "sol") cycle
!
!$$$$$
  do ii = 1,indarrayFlu
    npi = Array_Flu(ii)
!$$$$$
!
!* azzeramento quantita generali
    pg(npi)%var = zero
!!!    unity   = zero
!
!.. contributo delle particelle interne
!
    do contj = 1, nPartIntorno(npi)
!
      npartint = (npi-1)* NMAXPARTJ + contj
      npj = PartIntorno(npartint)
!
      rhoi   = pg(npi)%dens
      rhoj   = pg(npj)%dens
      amassj = pg(npj)%mass
!
      dervel(:) = pg(npj)%vel(:) - pg(npi)%vel(:)
!
      if ( pg(npj)%vel_type /= "std" ) then          !non part fix o altro
        rhoj   = rhoi
        amassj = pg(npi)%mass
        moddervel = -two * (pg(npi)%vel(1)*pg(npj)%zer(1)+pg(npi)%vel(2)*pg(npj)%zer(2)+pg(npi)%vel(3)*pg(npj)%zer(3))
        dervel(:) = moddervel * pg(npj)%zer(:)    !+pg(npj)%vstart(:)
      end if
!!============= CORREZIONE VELOCITA' ===================
      if ( Med(pg(npj)%imed)%den0 /= Med(pg(npi)%imed)%den0 ) cycle
!
      pesoj = amassj * PartKernel(4,npartint) / rhoj
!
!!!      unity = unity + pesoj  
      pg(npi)%var(:) = pg(npi)%var(:) + dervel(:) * pesoj   
!
! calcolo velocita' per modello diffusivo
!      pg(npj)%veldif(:) = pg(npj)%veldif(:) + pg(npj)%vel(:) * pesoj
!
    end do
!
!...................................... 2011 mar 08
!.. update of Specific Internal Energy
      if (esplosione) pg(npi)%Envar = pg(npi)%Envar + (pg(npj)%IntEn - pg(npi)%IntEn) * pesoj
!.. end update of Specific Internal Energy
!.........................................
!
!AA406test
!

!AA501b start
      if (n_bodies > 0) then
         pg(npi)%var(:) = pg(npi)%var(:) + dervel_mat(npi,:)
         unity = unity + unity_vec(npi)
      endif
 
!AA601 sub start
! Impose boundary conditions at inlet and outlet sections (DB-SPH)
      if (Domain%tipo == "bsph") then
         call DBSPH_inlet_outlet(npi)
         else
!AA601 sub end
      
!.. contributo dei contorni
!
    ncbf = BoundaryDataPointer(1,npi)
    ibdt = BoundaryDataPointer(3,npi)
!
    if (Ncbf > 0) then        !Ci sono facce di contorno vicine
!
      do icbf = 1, Ncbf
!
        ibdp = ibdt + icbf - 1
        LocX(1:SPACEDIM) = BoundaryDataTab(ibdp)%LocXYZ(1:SPACEDIM)
        iface = BoundaryDataTab(ibdp)%CloBoNum
!
        facestr = BoundaryFace(iface)%stretch
        strtype = Tratto(facestr)%tipo
!
        if (strtype == 'sour' .or. strtype == 'velo' .or. strtype == 'flow') then
          pg(npi)%var(:) = zero   
!!!          unity = one
          exit  
        end if
!       
!       do i = 1, SPACEDIM
!         rrr(i) = BoundaryFace(iface)%T(i,1)
!         sss(i) = BoundaryFace(iface)%T(i,2)
!         nnn(i) = BoundaryFace(iface)%T(i,3)
!         DVGlo(i) = two * (Tratto(facestr)%velocity(i) - pg(npi)%vel(i))
!       end do
!       DVLoc(1) = rrr(1) * DVGlo(1) + rrr(2) * DVGlo(2) + rrr(3) * DVGlo(3)
!       DVLoc(2) = sss(1) * DVGlo(1) + sss(2) * DVGlo(2) + sss(3) * DVGlo(3)
!       DVLoc(3) = nnn(1) * DVGlo(1) + nnn(2) * DVGlo(2) + nnn(3) * DVGlo(3)
!       
        DVGlo(:) = two * (Tratto(facestr)%velocity(:) - pg(npi)%vel(:))
        do i = 1,SPACEDIM
          DVLoc(i) = zero
          do j = 1,SPACEDIM
            DVLoc(i) = DVLoc(i) + BoundaryFace(iface)%T(j,i) * DVGlo(j)
          end do
        end do
!
        IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(2)
!
        if (strtype == 'fixe' .or. strtype == 'tapi') then
!
          BCLoc(1) = DVLoc(1) * IntWdV * Tratto(facestr)%ShearCoeff
          BCLoc(2) = DVLoc(2) * IntWdV * Tratto(facestr)%ShearCoeff
          BCLoc(3) = DVLoc(3) * IntWdV
!          BCGlo(:) = rrr(:) * BCLoc(1) + sss(:) * BCLoc(2) + nnn(:) * BCLoc(3)
          do i = 1, SPACEDIM
            BCGlo(i) = zero
            do j = 1,SPACEDIM
              BCGlo(i) = BCGlo(i) + BoundaryFace(iface)%T(i,j) * BCLoc(j)
            end do
          end do
!
          pg(npi)%var(:) = pg(npi)%var(:) + BCGlo(:)   
!
        end if
      end do
    end if

!AA601
     endif   
    
!!!    pg(npi)%uni = unity
!
!prova!
!    pg(npi)%var = floor(pg(npi)%var * azzeramento) / azzeramento
!    where (dabs(pg(npi)%var) < arrotondamento) pg(npi)%var = zero 
!prova!
!
  end do
!
!$omp end parallel do
!

!AA501b start
!Deallocations
  if (n_bodies > 0) then
     deallocate(dervel_mat)
     deallocate(unity_vec)
  endif
!AA501b end

return
end subroutine inter_SmoothVelo_3D
!---split

!cfile inter_CoefDif.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : inter_CoefDif
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
! Module purpose : Module to calculate the corrective term of the velocity
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D, time_integration
!
! Called routines: 
!
!************************************************************************************
!
subroutine inter_CoefDif (npi)
! ex inter6
!* implementa il meccanismo di ricerca delle particelle che agiscono
!* su quella i-esima.
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Local Parameters ..
!integer(4), parameter :: local_d = 500  ! num max part entro 2h
!
!.. Formal Arguments ..
integer(4),intent(IN) :: npi
!
!.. Local Scalars ..
integer(4)       :: npj,contj,npartint
double precision :: unity,rhoj,amassj,pesoj   !moddervel,  rhoi,
!
!.. Executable Statements ..
!
!* azzeramento quantita generali
 unity = zero
 pg(npi)%veldif(:) = zero
!
!*_______________________________________________________________
!*prima passata per trovare celle interagenti e memorizzazione
!
 do contj = 1, nPartIntorno(npi)
!
   npartint = (npi-1)* NMAXPARTJ + contj
   npj = PartIntorno(npartint)
!
!   rhoi   = pg(npi)%dens
   rhoj   = pg(npj)%dens
   amassj = pg(npj)%mass
!
!============= CALCOLO VELOCITA' PER COEFFICIENTE DIFFUSIVO ===================

   pesoj = amassj * PartKernel(4,npartint) / rhoj

   unity = unity + pesoj  

! calcolo velocita' per modello diffusivo
   pg(npi)%veldif(:) = pg(npi)%veldif(:) + pg(npj)%vel(:) * pesoj  

!============= FINE CALCOLO VELOCITA' PER COEFFICIENTE DIFFUSIVO ================
!
!prova!
!   pg(npi)%veldif = floor(pg(npi)%veldif * azzeramento) / azzeramento
!   where (dabs(pg(npi)%veldif) < arrotondamento) pg(npi)%veldif = zero 
!prova!
!
 end do
!
 pg(npi)%uni = unity
!
return
end subroutine inter_CoefDif
!---split

!cfile inter_SmoothPres.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : inter_SmoothPres
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!
!************************************************************************************
! Module purpose : Module to calculate the corrective term of the pressure
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
!
subroutine inter_SmoothPres
! ex inter7
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Local Scalars ..
integer(4)       :: npi,npj,contj,npartint,ii
double precision :: unity,presi,presj,rhoj,amassj,pesoj,appo1,appo2,TetaP1
!
!.. Executable Statements ..
!
!$omp parallel do default(none) &
!$omp private(npi,ii,unity,appo1,appo2,contj,npartint,npj,presi,rhoj,presj,amassj,pesoj) &
!$omp shared(nag,pg,Med,Domain,nPartIntorno,NMAXPARTJ,PartIntorno,PartKernel,indarrayFlu,Array_Flu)
!
!.. loops on all the particles
!
!!!!!!  do npi = 1,nag
!!!!!!!
!!!!!!    if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" .or. pg(npi)%state == "sol") cycle
!$$$$$
  do ii = 1,indarrayFlu
    npi = Array_Flu(ii)
!$$$$$
!
!* azzeramento quantita generali
    unity = zero
    appo1 = zero
    appo2 = zero
!
    do contj = 1, nPartIntorno(npi)
!
      npartint = (npi-1)* NMAXPARTJ + contj
      npj = PartIntorno(npartint)
!
      if ( pg(npj)%vel_type /= "std" ) cycle          !non part fix o altro
!
      presi  = pg(npi)%pres
      rhoj   = pg(npj)%dens    
      presj  = pg(npj)%pres    
      amassj = pg(npj)%mass
!
!.. calcola i pesi
      pesoj = amassj * PartKernel(4,npartint) / rhoj
!
      unity = unity + pesoj  
!
      appo1 = appo1 + (presj - presi) * pesoj  
      appo2 = appo2 - Domain%grav(3)*Med(pg(npi)%imed)%den0 * (pg(npj)%coord(3) - pg(npi)%coord(3)) * pesoj     !aggiunto SaMa
!
    end do
!
!--------------- SaMa -----------------
    if (Domain%Psurf /= 's') then
      pg(npi)%vpres = appo1
    else if (unity > 0.8d0) then
      pg(npi)%vpres = appo1
    else
      pg(npi)%vpres = appo1 + appo2
    end if
!
    pg(npi)%uni = unity
!--------------- SaMa -----------------
!
  end do
!
!$omp end parallel do
!
!  call cpu_time(cpu_loop1a)
!  call cpu_time(cpu_loop2a)
!
!$omp parallel do default(none) private(npi,ii,TetaP1) shared(nag,Pg,Med,Domain,dt,indarrayFlu,Array_Flu,esplosione)
!
!.. applies the correction to all the particles
!
!!!!!!         do npi = 1,nag
!!!!!!!
!!!!!!           if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" .or. pg(npi)%state == "sol") cycle
!$$$$$
  do ii = 1,indarrayFlu
    npi = Array_Flu(ii)
!$$$$$
!
    if (esplosione) then
!.. con Csound al posto di Celerita e' circa uguale
      TetaP1 = Domain%TetaP * pg(npi)%Csound * dt / Domain%h
    else
! calcolo TetaP adeguato al passo temporale
      TetaP1 = Domain%TetaP * Med(pg(npi)%imed)%Celerita * dt / Domain%h
    end if
!
!    if (pg(npi)%densass == 0) then
!.. updates the pressure 
    pg(npi)%pres = pg(npi)%pres + TetaP1 * pg(npi)%vpres / pg(npi)%uni
!
!.. updates the density depending on the local pressure, reference medium density and the comprimibility eps
    pg(npi)%dens = Med(pg(npi)%imed)%den0 * (one + pg(npi)%pres / Med(pg(npi)%imed)%eps)
!    end if
  end do

!$omp end parallel do
!
return
end subroutine inter_SmoothPres
!---split

!cfile InterFix.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : InterFix
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
! Module purpose : Module to accumulate the contributions of the particles that are
!                  in the sphere of influence of the particle considered
!
! Calling routine: NormFix
!
! Called routines: 
!
!************************************************************************************
!
subroutine InterFix (npi,appo,unity)
!* implementa il meccanismo di ricerca delle particelle che agiscono
!* su quella i-esima.
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
integer(4), parameter :: local_d = 500  ! num max part entro 2h
!
!.. Formal Arguments ..
integer(4),      intent(IN)    :: npi
double precision,intent(INOUT) :: unity
double precision,intent(INOUT),dimension(3) :: appo
!
!.. Local Scalars ..
integer(4) :: npj,contj,npartint   !!!,nfix
double precision :: rhoj,amassj,pesoj    !rhoi,
!
!.. Local Arrays ..
double precision,dimension(3) :: pesogradj
!
!.. Executable Statements ..
!
!* azzeramento quantita generali
 unity   = zero
 appo(:) = zero
!!! nfix  = zero
!
!*_______________________________________________________________
!*prima passata per trovare celle interagenti e memorizzazione
!
 do contj = 1, nPartIntorno(npi)
!
   npartint = (npi-1)* NMAXPARTJ + contj
   npj = PartIntorno(npartint)
!
   if ( pg(npj)%vel_type == "std" ) cycle      !non part fix o altro 
!!!   nfix = nfix + 1  
!
!   rhoi   = pg(npi)%dens
   rhoj   = pg(npj)%dens
   amassj = pg(npj)%mass
!
!* calcolo unita'
   pesoj = amassj * Partkernel(4,npartint) / rhoj
   pesogradj(1:3) = amassj * rag(1:3,npartint) * PartKernel(1,npartint) / rhoj
!
   unity = unity + pesoj  
   appo(:) = appo(:) + pesogradj(:)  
!
 end do
!
 appo(:) = -appo(:)
!
return
end subroutine InterFix 
!---split

!cfile inter_SmoothVF.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : inter_SmoothVF
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!
!************************************************************************************
! Module purpose : Module to calculate the corrective term of the volume fraction
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
!
subroutine inter_SmoothVF (npi,appo1,unity)
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
double precision,intent(INOUT) :: appo1, unity
!
!.. Local Scalars ..
integer(4)       :: npj,contj,npartint
double precision :: voli,volj,rhoj,amassj,pesoj
!
!.. Executable Statements ..
!
!* azzeramento quantita generali
 unity = zero
 appo1 = zero
!
 do contj = 1, nPartIntorno(npi)
!
   npartint = (npi-1)* NMAXPARTJ + contj
   npj = PartIntorno(npartint)
!
   voli   = pg(npi)%VolFra
   rhoj   = pg(npj)%dens    
   volj   = pg(npj)%VolFra    
   amassj = pg(npj)%mass
!
   if ( pg(npj)%vel_type /= "std" ) cycle          !non part fix o altro
!
!.. calcola i pesi
!
   pesoj = amassj * PartKernel(4,npartint) / rhoj
!
   unity = unity + pesoj  
!
   appo1 = appo1 + (volj - voli ) * pesoj  
!
!============= FINE CALCOLO VELOCITA' PER COEFFICIENTE DIFFUSIVO ================
!
 end do
!
 pg(npi)%uni = unity
!
return
end subroutine inter_SmoothVF
!---split

