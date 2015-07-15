!cfile Shields.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : Shields
!
! Last updating : April 08, 2014
!
! Improvement traceback:
!
! ..  G. Agate, S. Manenti           Initial development of the code
! 01  G. Agate, S. Manenti  15/05/12 Check on max number of iterations to compute Ustar
!                                    Check on beta angle for tau correction (Shields)
! 02  G. Agate, S. Manenti  22/05/12 DistZmin and Kbeta calculation modified
!AA504
! 03  Amicarelli            08/04/14 (v5.04) Main modifications for making the granular mixture to erode the fixed bed and calling the 3D erosion criterion.
!                                    Minor modifications: k=3d_90, Shields threshold for low Re*, convergence criterion and iterations from input data, ... 
!
!************************************************************************************
! Module purpose : Shields Erosion Model.
! 
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
!AA504 sub
! Called routines: diagnostic,fixed_bed_slope_limited,
!                  compute_k_BetaGamma,initialization_fixed_granular_particle
!
!************************************************************************************
!
!AA504 
subroutine Shields (npi)
!
!.. assign modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
!!double precision, parameter :: MinTotUnit = 0.95d0
!
!.. Formal Arguments ..
!
!AA504 sub start
integer(4),intent(in) :: npi 
integer(4) :: imed, intliq_id,intsol_id,indpeloloc,aux_ID,ncelcorr,igridi,jgridi,kgridi,iappo,contatore,ind_neigh,ind_ref_interface 
!AA504 sub end
double precision :: interf_liq, interf_sol, partLiq_pres 
double precision :: Ks, Ustar, Ustarold, Z0, Velocity, Velocity2
double precision :: Taub, Restar, RestarC, Tetacr, Taubcror, Taubcr, phi, beta, gamma, Kbeta, Kgamma, nu_ustar, deltanu, DistZmin
double precision :: pretot, DifRel
!AA504 start
double precision :: abs_vel_rij,aux_scal
double precision :: rel_dis(3),vel_rij(3),vel_s(3),aux_vec(3),aux_vec_2(3)
logical          :: flagz0,bed_slope_test 
!AA504 end
character(len=lencard)  :: nomsub = "Shields"
!
!.. External Routines ..
integer(4),external :: ParticleCellNumber, CellIndices, CellNumber
!
!.. Executable Statements ..
!

!AA504 rm part
!AA504 sub
 if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") return !cycle
!.. controllo movimento granulare e calcolo viscosita'
 imed = pg(npi)%imed
!AA504 sub start    
! Criterion computed only for granular particles (otherwise return)
 if (index(Med(imed)%tipo,"granular") <= 0) return !cycle      
 ncelcorr = ParticleCellNumber(pg(npi)%coord)
 iappo = CellIndices(ncelcorr,igridi,jgridi,kgridi) 
 ! In case of no erosion: simple management of SPH granualr particles (with or without a velocity threshold)
 if (Granular_flows_options%erosion_flag==1) then
    if (ind_interfaces(igridi,jgridi,4).ne.0) then 
       if (pg(npi)%coord(3)>=(pg(ind_interfaces(igridi,jgridi,4))%coord(3)-2.d0*Domain%h)) then 
          if (pg(npi)%state=="flu") return   
          pg(npi)%state = "flu"
          else
             pg(npi)%state="sol"
             pg(npi)%var = 0.d0
             pg(npi)%vel = 0.d0
             pg(npi)%sigma_prime = 0.0d0
       endif
       else
          if (pg(npi)%state=="flu") return   
          pg(npi)%state = "flu" 
    endif
    aux_ID = ind_interfaces(igridi,jgridi,4)
    if (aux_ID==0) aux_ID = ind_interfaces(igridi,jgridi,3)
    aux_vec(:) = pg(aux_ID)%vel_old(:) - pg(npi)%vel_old(:)
    Velocity2 = dot_product(aux_vec,aux_vec) 
    pretot = pg(aux_ID)%pres  + (pg(aux_ID)%coord(3) - pg(npi)%coord(3)) * (-Domain%grav(3)) * med(pg(aux_ID)%imed)%den0 + Velocity2 * half * Med(pg(aux_ID)%imed)%den0
    pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita*Med(imed)%celerita))   
    return
 endif
 if (Granular_flows_options%ID_erosion_criterion==1) then
! SPH granular particles below the fixed bed (more than 2h) are fixed (this involve Beta_max=45°)
    bed_slope_test = .true.
    call fixed_bed_slope_limited(npi,igridi,jgridi,bed_slope_test) 
    if (bed_slope_test==.false.) return
! Deposition for mobile particles without fixed neighbours, with no fixed bed in the column and close to frontiers (option in input: deposition_at_frontiers);
! no deposition for mobile particles far from the fixed bed.    
    if ((pg(npi)%state=="flu").and.(pg(npi)%ind_neigh_mix_bed==0)) then
       if (((BoundaryDataPointer(2,npi)>0).and.(Granular_flows_options%deposition_at_frontiers==1)) .and. (ind_interfaces(igridi,jgridi,4)==0)) then
          pg(npi)%state="sol"
          pg(npi)%var = 0.d0
          pg(npi)%vel = 0.d0
          pg(npi)%sigma_prime = 0.0d0
          if (pg(npi)%indneighliqsol.ne.0) then
             aux_ID = pg(npi)%indneighliqsol
             else 
                aux_ID = pg(npi)%ind_neigh_mob_for_granmob    
          endif
          if (aux_ID.ne.0) then
             aux_vec(:) = pg(aux_ID)%vel_old(:) - pg(npi)%vel_old(:)
             Velocity2 = dot_product(aux_vec,aux_vec) 
             pretot = pg(aux_ID)%pres  + (pg(aux_ID)%coord(3) - pg(npi)%coord(3)) * (-Domain%grav(3)) * med(pg(aux_ID)%imed)%den0 + Velocity2 * half * Med(pg(aux_ID)%imed)%den0
             pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita*Med(imed)%celerita))   
             else
                indpeloloc = ind_interfaces(igridi,jgridi,3)
                if (indpeloloc.ne.0) then
                   pretot = pg(indpeloloc)%pres  + (pg(indpeloloc)%coord(3) - pg(npi)%coord(3)) * (-Domain%grav(3)) * med(imed)%den0
                   pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita*Med(imed)%celerita))   
                   else
                      pg(npi)%dens = med(imed)%den0
                endif 
          endif 
       endif
       return 
    endif
!AA504 sub end
!AA504 start 
! Particles with no free surface in the grid column: with fluid neighbours are set free to move, with no fluid neighbours are set fixed.
    if (Granular_flows_options%erosion_flag==2) then
       if (ind_interfaces(igridi,jgridi,1)==0) then
          if (pg(npi)%indneighliqsol==1)  then
             pg(npi)%state = "flu" 
             aux_ID = pg(npi)%indneighliqsol
             aux_vec(:) = pg(aux_ID)%vel_old(:) - pg(npi)%vel_old(:)
             Velocity2 = dot_product(aux_vec,aux_vec) 
             pretot = pg(aux_ID)%pres  + (pg(aux_ID)%coord(3) - pg(npi)%coord(3)) * (-Domain%grav(3)) * med(pg(aux_ID)%imed)%den0 + Velocity2 * half * Med(pg(aux_ID)%imed)%den0
             pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita*Med(imed)%celerita))   
             else
                pg(npi)%state="sol"
                pg(npi)%var = 0.d0
                pg(npi)%vel = 0.d0
                pg(npi)%sigma_prime = 0.0d0
                indpeloloc = ind_interfaces(igridi,jgridi,3)
                if (indpeloloc.ne.0) then
                   pretot = pg(indpeloloc)%pres  + (pg(indpeloloc)%coord(3) - pg(npi)%coord(3)) * (-Domain%grav(3)) * med(imed)%den0
                   pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita*Med(imed)%celerita))   
                   else
                      pg(npi)%dens = med(imed)%den0
                endif 
          endif
          return 
       endif
    endif
!AA504 end
 endif
!AA504 end 
!AA504 rm part
!.. controllo se la particella si muove in caso di esplosione non applico il modello di erosione  
 if (esplosione) then
    Velocity2 = pg(npi)%vel(1)*pg(npi)%vel(1) + pg(npi)%vel(2)*pg(npi)%vel(2) + pg(npi)%vel(3)*pg(npi)%vel(3)
!AA504 sub        
    if (Velocity2 > 1.0e-3) return 
 end if
!AA504 part moved above
!.. controllo esistenza particella liquid vicina
!AA504 sub start
! Choice of the representative neighbour (priority to pure fluid neighbours)
 if (Granular_flows_options%ID_erosion_criterion==1) then
    if (pg(npi)%indneighliqsol.ne.0) then
       ind_neigh = pg(npi)%indneighliqsol
       elseif (pg(npi)%state=="sol") then
          ind_neigh = pg(npi)%ind_neigh_mix_bed
          else
             ind_neigh = pg(npi)%ind_neigh_mob_for_granmob
    endif
    else
       ind_neigh = pg(npi)%indneighliqsol
 endif 
 ind_ref_interface = ind_interfaces(igridi,jgridi,3)
 if (ind_neigh == 0) then
!AA504 sub end
    pg(npi)%state = "sol"
    pg(npi)%vel = zero
    pg(npi)%var = zero
!AA504 sub start
! Deposition or no erosion for particles with no neighbours 
    if (Granular_flows_options%ID_erosion_criterion==1) then
       pg(npi)%sigma_prime = 0.0d0
       indpeloloc = ind_interfaces(igridi,jgridi,3)
       else
          indpeloloc = ind_interfaces(igridi,jgridi,1)
    endif
    if (indpeloloc /= 0) then
       intsol_id = ind_ref_interface
       if (intsol_id == 0) intsol_id = npi
       if (Granular_flows_options%ID_erosion_criterion==1) then
          Velocity2 = pg(indpeloloc)%vel_old(1)*pg(indpeloloc)%vel_old(1) + pg(indpeloloc)%vel_old(2)*pg(indpeloloc)%vel_old(2) + &
                      pg(indpeloloc)%vel_old(3)*pg(indpeloloc)%vel_old(3)
          pretot = pg(indpeloloc)%pres + (pg(indpeloloc)%coord(3) - pg(intsol_id)%coord(3)) * (-Domain%grav(3)) * &
                   med(pg(indpeloloc)%imed)%den0 + (pg(intsol_id)%coord(3) - pg(npi)%coord(3)) * (-Domain%grav(3)) * med(imed)%den0
          else
             Velocity2 = pg(indpeloloc)%vel(1)*pg(indpeloloc)%vel(1) + pg(indpeloloc)%vel(2)*pg(indpeloloc)%vel(2) + &
                         pg(indpeloloc)%vel(3)*pg(indpeloloc)%vel(3)
             pretot = pg(indpeloloc)%pres  + Velocity2 * half * med(pg(indpeloloc)%imed)%den0 + &
                      (pg(indpeloloc)%coord(3) - pg(intsol_id)%coord(3)) * (-Domain%grav(3)) * med(pg(indpeloloc)%imed)%den0 + &
                      (pg(intsol_id)%coord(3) - pg(npi)%coord(3)) * (-Domain%grav(3)) * med(imed)%den0     
       endif
       pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita*Med(imed)%celerita))
    end if
       return 
 end if
!AA504 sub end
!.. calcolo interfaccia liquida e solida
!AA504 sub start
if (Granular_flows_options%ID_erosion_criterion==1) then
   intliq_id = ind_neigh
   intsol_id = npi
   else
       intliq_id = pg(npi)%indneighliqsol
       if (pg(npi)%indneighliqsol.ne.0) then
          intsol_id = pg(pg(npi)%indneighliqsol)%indneighliqsol
          else
             if (ind_interfaces(igridi,jgridi,3).ne.0) then
                intsol_id = ind_interfaces(igridi,jgridi,3)
                else
                   intsol_id = npi
             endif
       endif
 endif 
!AA504 sub end
 interf_liq = pg(intliq_id)%coord(3)
 interf_sol = pg(intsol_id)%coord(3)
!AA504 sub start 
! Compute pretot with kinetic term (depeding on the relative velocity) 
 if (Granular_flows_options%ID_erosion_criterion==1) then      
    Velocity2 = dot_product(pg(intliq_id)%vel_old,pg(intliq_id)%vel_old)
    Velocity = dsqrt(Velocity2)
    PartLiq_pres = (interf_liq - interf_sol) * med(pg(intliq_id)%imed)%den0 * (-Domain%grav(3)) + pg(intliq_id)%pres 
    pretot = PartLiq_pres + Velocity2 * half * med(pg(intliq_id)%imed)%den0    
!    pretot = PartLiq_pres + Med(pg(intliq_id)%imed)%celerita * Velocity * med(pg(intliq_id)%imed)%den0 
    else
       Velocity2 = pg(intliq_id)%vel(1)*pg(intliq_id)%vel(1) + pg(intliq_id)%vel(2)*pg(intliq_id)%vel(2) + &
                   pg(intliq_id)%vel(3)*pg(intliq_id)%vel(3)
       PartLiq_pres = (interf_liq - interf_sol) * med(pg(intliq_id)%imed)%den0 * (-Domain%grav(3)) + &
                      pg(intliq_id)%pres + Velocity2 * half * med(pg(intliq_id)%imed)%den0
       pretot  = PartLiq_pres + med(pg(intsol_id)%imed)%den0 * (-Domain%grav(3)) * (interf_sol - pg(npi)%Coord(3)) 
 endif
!AA504 sub end  
!.. calcolo velocita'
 if (Velocity2 == zero) then
    pg(npi)%state = "sol"
!AA504 rm line       
    pg(npi)%vel = zero
    pg(npi)%var = zero
    pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita*Med(imed)%celerita))
!AA504
    if (Granular_flows_options%ID_erosion_criterion==1) pg(npi)%sigma_prime = 0.0d0
       return 
 end if
!.. ipotizzo RE*<2 (laminare) e Z<deltanu (particella entro lo strato laminare)
!AA504 sub start
! Compute roughness (k): c_4=k/d_50 (from input, Shields 2D) or c_4=3 (Shields 3D)
 if (Granular_flows_options%ID_erosion_criterion == 1) then
    Ks = 3.d0 * med(imed)%d_90
    else
       Ks = med(imed)%RoughCoef * med(imed)%D50
 endif
 if (Granular_flows_options%ID_erosion_criterion > 1) then     
    Velocity = Dsqrt(Velocity2)
    else
!       if (ind_interfaces(igridi,jgridi,4).ne.0) pg(npi)%normal_int(:) = pg(ind_interfaces(igridi,jgridi,4))%normal_int_old(:)  !this line needs normal_int_old as a particle can be non-representative of its column, but representative for a close column  
       if (Granular_flows_options%erosion_flag==2) then
! default value for interface normal (in case of no free surface)           
          if (ind_interfaces(igridi,jgridi,1)==0) then ! at this point the neighbour can be only a granular particle
             pg(npi)%normal_int(1) = 0.d0
             pg(npi)%normal_int(2) = 0.d0
             pg(npi)%normal_int(3) = 1.d0
          endif
       endif
! Compute velocity component perpendicular to the interface normal       
       vel_s(:) = pg(intliq_id)%vel_old(:) - dot_product(pg(intliq_id)%vel_old,pg(npi)%normal_int) * pg(npi)%normal_int(:)  
       Velocity = dsqrt(dot_product(vel_s,vel_s))
 endif             
!AA504 sub end      
!.. utilizzare rijtempmin e non la distanza in z (non funziona per piano inclinato vedi Flushing)
!AA504 sub comm
!..  DistZmin è la minima distanza ortogonale all'interfaccia di riferimento
!      DistZmin = Dabs( pg(npi)%coord(3) - pg(pg(npi)%indneighliqsol)%coord(3) )
!AA504 sub start
! Compute distance between the compuattional particle and the interacting particle (the height difference does not work for inclined beds)
 if (Granular_flows_options%ID_erosion_criterion==1) then
    if (pg(npi)%indneighliqsol.ne.0) then
       DistZmin = pg(npi)%rijtempmin(1)
       elseif (pg(npi)%state=="sol") then
          DistZmin = pg(npi)%rijtempmin(2)
          else
             DistZmin = pg(npi)%rijtempmin(3)
    endif
    else
       DistZmin = pg(npi)%rijtempmin(1)
 endif 
 Ustar = dsqrt(Velocity * pg(intliq_id)%mu / (Med(pg(intliq_id)%imed)%den0 * DistZmin))
!AA504 sub end
!.. ciclo determinazione Z0 e Ustar da Zhou Liu 'Sediment transport' pag.9
 flagz0 = .true.
 DifRel = one
 contatore = 0
!      iter_ustar: do while (flagz0 .and. (DifRel > 0.01d0))
!AA504 sub
 iter_ustar: do while (flagz0 .and. (DifRel > Granular_flows_options%conv_crit_erosion))
!.. controllo su numero iterazioni
    contatore = contatore + 1
!AA504 sub        
    if (contatore > Granular_flows_options%n_max_iterations) then
       inquire (file=nomefileerr, EXIST=err_flag)
       if (.not. err_flag) open(unit=uniterr,file=nomefileerr,form='formatted',status='unknown')
       write (nout,*) ' WARNING! Superato il numero massimo di iterazioni per il calcolo di Ustar. ' 
       write (nscr,*) ' WARNING! Superato il numero massimo di iterazioni per il calcolo di Ustar. ' 
       write (uniterr,*) '   ' 
       write (uniterr,*) ' WARNING! Superato il numero massimo di iterazioni per il calcolo di Ustar. ' 
       write (uniterr,*) '  Provare a ridurre il valore di Ks e rilanciare (for Shields 2D or Mohr erosion criterion). ' 
       write (uniterr,*) '  tempo = ',tempo,'  Velocity',Velocity
       write (uniterr,*) '  Ustar',Ustar,'  Ustarold',Ustarold
!AA504 sub
       write (uniterr,'(a,i6,a,3f10.7,a,e10.3)') '  npi = ',npi,'      pg(npi)%coord = ',pg(npi)%coord(1),pg(npi)%coord(2),pg(npi)%coord(3)," residual = ",DifRel
       write (uniterr,*) '  DistZmin',DistZmin,'  Z0',Z0
!.. Hydraulically rough flow
       Z0 = 0.033d0 * Ks
       flagz0 = .false.
       if (DistZmin > Z0) Ustar = vKconst * Velocity / Dlog(DistZmin / Z0)
      exit iter_ustar
    end if
!AA504 sub
    nu_ustar = pg(intliq_id)%mu / (Med(pg(intliq_id)%imed)%den0 * Ustar)
    RestarC = Ks / nu_ustar
    deltanu = 11.6d0 * nu_ustar
    if (deltanu > DistZmin .and. RestarC < 5.0d0) then
       flagz0 = .false.
       else
          if (RestarC >= 70.0d0) then
!.. Hydraulically rough flow
             Z0 = 0.033d0 * Ks
             else if (RestarC <= 5.0d0) then
!.. Hydraulically smooth flow
                Z0 = 0.11d0 * nu_ustar
                else
!.. Hydraulically transition flow
                   Z0 = 0.11d0 * nu_ustar + 0.033d0 * Ks
          end if
          Ustarold = Ustar
!.. inizio prova 06giu2011
!          if (Z0 <= Ks) then
!.. fine prova 06giu2011
!.. inizio prova 07giu2011
!          if (Z0 <= med(imed)%D50) then
!            Ustar = dsqrt(Velocity * pg(intliq_id)%mu / (pg(intliq_id)%dens * DistZmin))
!            flagz0 = .false.
!          else if (Z0 < DistZmin) then
!.. fine prova 07giu2011
          if (DistZmin > Z0) then
             Ustar = vKconst * Velocity / Dlog(DistZmin / Z0)
!AA504 sub              
             DifRel = dabs(Ustar-Ustarold) / Ustar
!            DifRel = dabs(Ustar-Ustarold) / Ustarold
             else
                flagz0 = .false.
                inquire (file=nomefileerr, EXIST=err_flag)
                if (.not. err_flag) open(unit=uniterr,file=nomefileerr,form='formatted',status='unknown')
                write (nout,*) ' WARNING! Z0>rijtempmin per il calcolo di Ustar. ' 
                write (nscr,*) ' WARNING! Z0>rijtempmin per il calcolo di Ustar. ' 
                write (uniterr,*) '   ' 
                write (uniterr,*) ' WARNING! Z0>rijtempmin per il calcolo di Ustar. ' 
                write (uniterr,*) '  tempo = ',tempo,'  Velocity',Velocity
                write (uniterr,*) '  Ustar',Ustar,'  Ustarold',Ustarold
                write (uniterr,'(a,i6,a,3f10.7)') '  npi = ',npi,'      pg(npi)%coord = ',pg(npi)%coord(1),pg(npi)%coord(2),pg(npi)%coord(3)
                write (uniterr,*) '  DistZmin',DistZmin,'  Z0',Z0
          end if
    endif
 enddo iter_ustar 
 !!!!if (isnan(Ustar)==.true.) then
 if (isnan(Ustar)) then
    call diagnostic (arg1=11,arg2=1,arg3=nomsub)
 end if
!.. calcolo teta critico secondo Browlie 1981
!!!!SM      Rp = med(imed)%D50 / pg(intflu_id)%visc * Dsqrt((pg(intsol_id)%dens/pg(intflu_id)%dens - one) * GI * med(imed)%D50)
!!!!AG      Rp = med(imed)%D50 / pg(intflu_id)%visc * Dsqrt((pg(npi)%dens/pg(intflu_id)%dens - one) * GI * med(imed)%D50)
!!!      Rp = med(imed)%D50 / pg(intflu_id)%visc * Dsqrt((pg(intsol_id)%dens / pg(intflu_id)%dens - one) * GI * med(imed)%D50)
!!!      Rp = Rp ** 0.6D0
!!!      Tetacr = 0.22D0 / Rp + 0.06D0 * Dexp(-17.77D0 / Rp)
!.. calcolo taub
!AA504 sub
 Taub = Ustar * Ustar * Med(pg(intliq_id)%imed)%den0
!.. calcolo Restar con nuovo ustar e teta critico secondo Shield 1936
!AA504 sub
 Restar = Med(pg(intliq_id)%imed)%den0 * Ustar * med(imed)%D50 / pg(intliq_id)%mu
!AA504 sub start      
 if (Restar >= 500.0D0) then
    Tetacr = 0.068D0
    elseif (Restar<=1.d0) then
       Tetacr = 0.1d0    
       else
          Tetacr = 0.010595D0 * Dlog(Restar) + 0.110476D0 / Restar + 0.0027197D0
 end if
!AA504 sub end           
!.. calcolo taub critico orizzontale
!AA504 sub start
 Taubcror = Tetacr * GI * med(imed)%D50 * (Med(pg(intsol_id)%imed)%den0_s - Med(pg(intliq_id)%imed)%den0)
 if (Granular_flows_options%ID_erosion_criterion == 1) then
    pg(npi)%u_star = Ustar
    call compute_k_BetaGamma(npi,intliq_id,DistZmin)
    Taubcr = Taubcror * pg(npi)%k_BetaGamma
    else
!AA504 sub end          
!.. calcolo angolo inclinazione
!.. Attenzione! vel(1) e' la velocita' nella direzione prevalente
!!      beta = Datan2(pg(intliq_id)%vel(3) , pg(intliq_id)%vel(1))
!      beta = Datan(pg(intliq_id)%vel(3) / Dabs(pg(intliq_id)%vel(1)))   ! funziona anche con X invertite ??!!
       beta = Dabs(Datan(pg(intliq_id)%vel(3) / pg(intliq_id)%vel(1)))   ! funziona anche con X invertite
       phi = med(imed)%phi
!      
!      if (Dabs(beta) > phi) then
       if (beta > phi) then
          Kbeta = zero 
!      else if (Dabs(beta) < PIGRECO/18.0d0) then
!        kbeta = 0.4
!!!      else if (pg(intliq_id)%vel(3) <= zero .and. beta < PIGRECO/18.0d0) then
!!!        kbeta = 0.4 * Dsin(beta)/Dsin(PIGRECO/18.0d0)
          else
             if (pg(intliq_id)%vel(3) <= zero) then 
                Kbeta = Dsin(phi - beta) / Dsin(phi)  ! downsloping
                else
                   Kbeta = Dsin(phi + beta) / Dsin(phi)  ! upsloping 
             end if
       end if
       if (ncord == 2) then
          Kgamma = one
          else
!!        gamma = Datan2(pg(intliq_id)%vel(3) , pg(intliq_id)%vel(2))
!        gamma = Datan(pg(intliq_id)%vel(3) / Dabs(pg(intliq_id)%vel(2)))   ! funziona anche con X invertite ??!!
             gamma = Dabs(Datan(pg(intliq_id)%vel(3) / pg(intliq_id)%vel(2)))   ! funziona anche con X invertite
             Kgamma = one - Dtan(gamma)*Dtan(gamma) / (Dtan(phi)*Dtan(phi))
             if (kgamma > zero) then
                Kgamma = Dcos(gamma) * Dsqrt(kgamma)
                else
                   kgamma = zero
             end if
       end if
!$$$$
!.. ATTENZIONE al calcolo di gamma usare i coseni direttori 
!..    il calcolo precedente e' errato???!!!!!!!!!!!!!!!!!!!!!!!!
       kgamma = one
!$$$$
!$$$$
!.. ATTENZIONE al calcolo di kbeta e kgamma per caso traversa 
!!! caso 7 e 8      kbeta = 0.4
!!! caso 9      kgamma = 0.4
!      kbeta = 0.4
!$$$$
!
!.. calcolo taub critico da Zhang 2006
       Taubcr = Taubcror * Kbeta * Kgamma
!
!      if (Taub > Taubcr ) then                                              !10
!      if (Taub > Taubcr .or. pg(npi)%coord(3) <= pg(intliq_id)%coord(3)) then !10.2
!      if (Taub > Taubcr .or. pg(npi)%coord(3) > pg(intliq_id)%coord(3)) then  !10.3
!      if (Taub > Taubcr .or. pg(npi)%punta) then                               !10.41 e 10.42 25gen2011

!AA504 start
 endif
 if (Taubcr.ne.0.d0) then
    pg(npi)%tau_tauc = Taub/Taubcr
    if (pg(npi)%tau_tauc>99999.d0) pg(npi)%tau_tauc=99999.d0
    else
       pg(npi)%tau_tauc = 99999.d0  
 endif
!AA504 end    
!AA504 rm line
!AA504 start 
 if ((Taub > Taubcr).and.(it_corrente > Med(imed)%NIterSol)) then 
    if (pg(npi)%state.ne."flu") then
       pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita*Med(imed)%celerita))
       pg(npi)%state = "flu"
    endif   
!AA504end           
    else
!.. controllo se le particelle sono vicine all'uscita
       if (pg(npi)%CloseBcOut == 0) then
          pg(npi)%state = "sol"
          pg(npi)%vel = zero
          pg(npi)%var = zero
!AA504 mv part 
!AA504 sub start
          pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita*Med(imed)%celerita))
          if (diffusione) pg(npi)%dens = med(imed)%den0 
!AA504 sub end          
!AA504 start  
! Initializing some SPH granular parameters for fixed particles
          if (Granular_flows_options%ID_erosion_criterion==1) pg(npi)%sigma_prime = 0.0d0
!AA504 end 
       end if
 end if
!AA504 rm part
return
end subroutine Shields
!---split

