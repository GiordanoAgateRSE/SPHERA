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

