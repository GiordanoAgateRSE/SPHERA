!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
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
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: Shields 
! Description: 3D erosion criterion based on the formulation of both Shields-van
!              Rijn 2D criterion and Seminara et al.(2002) 3D criterion.
!              2D Shields erosion criterion based on pure fluid - fixed bed 
!              interactions  (Manenti et al., 2012, JHE). 
!              Extension for bed load transport layer - fixed bed interactions 
!              (Amicarelli et al., CAF, submitted).
!              Extension to the third dimension (Amicarelli et al., CAF, 
!              submitted).
!              k=3d_90 (Manenti et al., 2012, JHE; Amicarelli et al., CAF, 
!              submitted).
!              Shields threshold for low Re* (Amicarelli et al., CAF,submitted).  
!-------------------------------------------------------------------------------
subroutine Shields(npi)
!------------------------
! Modules
!------------------------ 
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: npi 
logical :: flagz0,bed_slope_test
integer(4) :: imed,intliq_id,intsol_id,indpeloloc,aux_ID,ncelcorr,igridi 
integer(4) :: jgridi,kgridi,iappo,contatore,ind_neigh,ind_ref_interface 
double precision :: interf_liq,interf_sol,partLiq_pres,Ks,Ustar
double precision :: Ustarold,Z0,Velocity,Velocity2,Taub,Restar,RestarC
double precision :: Tetacr,Taubcror,Taubcr,phi,beta,gamma,Kbeta,Kgamma
double precision :: nu_ustar,deltanu,DistZmin,pretot,DifRel,aux_scal
double precision :: rel_dis(3),vel_rij(3),vel_s(3),aux_vec(3),aux_vec_2(3)
character(len=lencard)  :: nomsub = "Shields"
integer(4),external :: ParticleCellNumber,CellIndices,CellNumber
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
if ((pg(npi)%cella==0).or.(pg(npi)%vel_type/="std")) return 
! Check motion status and viscosity computation
imed = pg(npi)%imed
! Criterion computed only for granular particles (otherwise "return")
if (index(Med(imed)%tipo,"granular")<=0) return       
ncelcorr = ParticleCellNumber(pg(npi)%coord)
iappo = CellIndices(ncelcorr,igridi,jgridi,kgridi) 
! In case of no erosion: simple management of SPH mixture particles
! (eventual influence of both a positive velocity threshold and/or the 
! frictional viscosity threshold).
if (Granular_flows_options%erosion_flag==1) then
   if (ind_interfaces(igridi,jgridi,4).ne.0) then
! The mixture particles which define the fixed bed have null velocities and 
! the variable "state" equal to '"flu"', even though they are fixed.
      if (pg(npi)%coord(3)>=(pg(ind_interfaces(igridi,jgridi,4))%coord(3)      &
         -2.d0 * Domain%h)) then 
         if (pg(npi)%state=="flu") return
         pg(npi)%state = "flu"
         else
            pg(npi)%state = "sol"
            pg(npi)%var = 0.d0
            pg(npi)%vel = 0.d0
            pg(npi)%sigma_prime_m = 0.0d0
            pg(npi)%pres_fluid = 0.0d0
      endif
      else
         if (pg(npi)%state=="flu") return   
         pg(npi)%state = "flu" 
   endif
   aux_ID = ind_interfaces(igridi,jgridi,4)
   if (aux_ID==0) aux_ID = ind_interfaces(igridi,jgridi,3)
   aux_vec(:) = pg(aux_ID)%vel_old(:) - pg(npi)%vel_old(:)
   Velocity2 = dot_product(aux_vec,aux_vec) 
   pretot = pg(aux_ID)%pres  + (pg(aux_ID)%coord(3) - pg(npi)%coord(3)) *      &
            (-Domain%grav(3)) * med(pg(aux_ID)%imed)%den0 + Velocity2 * half * &
            Med(pg(aux_ID)%imed)%den0
   pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita *             &
                  Med(imed)%celerita))  
   return
endif
if (Granular_flows_options%ID_erosion_criterion==1) then
! SPH granular particles below the fixed bed (more than 2h) are fixed 
! (this involves Beta_max=45°)
   bed_slope_test = .true.
   call fixed_bed_slope_limited(npi,igridi,jgridi,bed_slope_test) 
   if (bed_slope_test.eqv..false.) return
! Deposition for mobile particles without fixed neighbours, with no fixed bed 
! in the column and close to frontiers (option in input: 
! deposition_at_frontiers);
! no deposition for mobile particles far from the fixed bed.    
   if ((pg(npi)%state=="flu").and.(pg(npi)%ind_neigh_mix_bed==0)) then
      if (((BoundaryDataPointer(2,npi)>0).and.                                 &
         (Granular_flows_options%deposition_at_frontiers==1)).and.             &
         (ind_interfaces(igridi,jgridi,4)==0)) then
         pg(npi)%state="sol"
         pg(npi)%var = 0.d0
         pg(npi)%vel = 0.d0
         pg(npi)%sigma_prime_m = 0.0d0
         pg(npi)%pres_fluid = 0.0d0
         if (pg(npi)%indneighliqsol.ne.0) then
            aux_ID = pg(npi)%indneighliqsol
            else 
               aux_ID = pg(npi)%ind_neigh_mob_for_granmob    
         endif
         if (aux_ID.ne.0) then
            aux_vec(:) = pg(aux_ID)%vel_old(:) - pg(npi)%vel_old(:)
            Velocity2 = dot_product(aux_vec,aux_vec) 
            pretot = pg(aux_ID)%pres + (pg(aux_ID)%coord(3) - pg(npi)%coord(3))&
                     * ( - Domain%grav(3)) * med(pg(aux_ID)%imed)%den0 +       &
                     Velocity2 * half * Med(pg(aux_ID)%imed)%den0
            pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita *    &
                           Med(imed)%celerita))   
            else
               indpeloloc = ind_interfaces(igridi,jgridi,3)
               if (indpeloloc.ne.0) then
                  pretot = pg(indpeloloc)%pres  + (pg(indpeloloc)%coord(3) -   &
                           pg(npi)%coord(3)) * (-Domain%grav(3)) *             &
                           med(imed)%den0
                  pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita&
                                 * Med(imed)%celerita))   
                  else
                     pg(npi)%dens = med(imed)%den0
               endif 
         endif 
      endif
      return 
   endif
! Particles with no free surface in the grid column: they are set free to move 
! in the presence of fluid neighbours; they are set fixed otherwise.
   if (Granular_flows_options%erosion_flag==2) then
      if (ind_interfaces(igridi,jgridi,1)==0) then
         if (pg(npi)%indneighliqsol==1)  then
            pg(npi)%state = "flu" 
            aux_ID = pg(npi)%indneighliqsol
            aux_vec(:) = pg(aux_ID)%vel_old(:) - pg(npi)%vel_old(:)
            Velocity2 = dot_product(aux_vec,aux_vec) 
            pretot = pg(aux_ID)%pres  + (pg(aux_ID)%coord(3) -                 &
                     pg(npi)%coord(3)) * ( - Domain%grav(3)) *                 &
                     med(pg(aux_ID)%imed)%den0 + Velocity2 * half *            &
                     Med(pg(aux_ID)%imed)%den0
            pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita *    &
                           Med(imed)%celerita))   
            else
               pg(npi)%state="sol"
               pg(npi)%var = 0.d0
               pg(npi)%vel = 0.d0
               pg(npi)%sigma_prime_m = 0.0d0
               pg(npi)%pres_fluid = 0.0d0
               indpeloloc = ind_interfaces(igridi,jgridi,3)
               if (indpeloloc.ne.0) then
                  pretot = pg(indpeloloc)%pres  + (pg(indpeloloc)%coord(3) -   &
                           pg(npi)%coord(3)) * ( - Domain%grav(3)) *           &
                           med(imed)%den0
                  pg(npi)%dens = med(imed)%den0 + (pretot /                    &
                                 (Med(imed)%celerita*Med(imed)%celerita))   
                  else
                     pg(npi)%dens = med(imed)%den0
               endif 
         endif
         return 
      endif
   endif
endif
! In case of explosion, erosion criterion is not active.  
if (esplosione) then
   Velocity2 = pg(npi)%vel(1) * pg(npi)%vel(1) + pg(npi)%vel(2) *              &
               pg(npi)%vel(2) + pg(npi)%vel(3) * pg(npi)%vel(3)
   if (Velocity2>1.0e-3) return 
end if
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
if (ind_neigh==0) then
   pg(npi)%state = "sol"
   pg(npi)%vel = zero
   pg(npi)%var = zero
! Deposition or no erosion for particles with no neighbours 
   if (Granular_flows_options%ID_erosion_criterion==1) then
      pg(npi)%sigma_prime_m = 0.0d0
      pg(npi)%pres_fluid = 0.0d0
      indpeloloc = ind_interfaces(igridi,jgridi,3)
      else
         indpeloloc = ind_interfaces(igridi,jgridi,1)
   endif
   if (indpeloloc/=0) then
      intsol_id = ind_ref_interface
      if (intsol_id==0) intsol_id = npi
      if (Granular_flows_options%ID_erosion_criterion==1) then
         Velocity2 = pg(indpeloloc)%vel_old(1) * pg(indpeloloc)%vel_old(1) +   &
                     pg(indpeloloc)%vel_old(2) * pg(indpeloloc)%vel_old(2) +   &
                     pg(indpeloloc)%vel_old(3) * pg(indpeloloc)%vel_old(3)
         pretot = pg(indpeloloc)%pres + (pg(indpeloloc)%coord(3) -             &
                  pg(intsol_id)%coord(3)) * ( - Domain%grav(3)) *              &
                  med(pg(indpeloloc)%imed)%den0 + (pg(intsol_id)%coord(3) -    &
                  pg(npi)%coord(3)) * ( - Domain%grav(3)) * med(imed)%den0
         else
            Velocity2 = pg(indpeloloc)%vel(1) * pg(indpeloloc)%vel(1) +        &
                        pg(indpeloloc)%vel(2) * pg(indpeloloc)%vel(2) +        &
                        pg(indpeloloc)%vel(3) * pg(indpeloloc)%vel(3)
            pretot = pg(indpeloloc)%pres + Velocity2 * half *                  &
                     med(pg(indpeloloc)%imed)%den0 +                           &
                     (pg(indpeloloc)%coord(3) - pg(intsol_id)%coord(3)) *      &
                     ( - Domain%grav(3)) * med(pg(indpeloloc)%imed)%den0 +     &
                     (pg(intsol_id)%coord(3) - pg(npi)%coord(3)) *             &
                     ( - Domain%grav(3)) * med(imed)%den0     
      endif
      pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita *          &
                     Med(imed)%celerita))
   end if
   return 
end if
! Pure fluid - bed load transport layer interface 
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
interf_liq = pg(intliq_id)%coord(3)
interf_sol = pg(intsol_id)%coord(3)
! To compute "pretot" with kinetic term (depeding on the relative velocity) 
if (Granular_flows_options%ID_erosion_criterion==1) then      
   Velocity2 = dot_product(pg(intliq_id)%vel_old,pg(intliq_id)%vel_old)
   Velocity = dsqrt(Velocity2)
   PartLiq_pres = (interf_liq - interf_sol) * med(pg(intliq_id)%imed)%den0 *   &
                  ( - Domain%grav(3)) + pg(intliq_id)%pres 
   pretot = PartLiq_pres + Velocity2 * half * med(pg(intliq_id)%imed)%den0    
   else
      Velocity2 = pg(intliq_id)%vel(1) * pg(intliq_id)%vel(1) +                &
                  pg(intliq_id)%vel(2) * pg(intliq_id)%vel(2) +                &
                  pg(intliq_id)%vel(3) * pg(intliq_id)%vel(3)
      PartLiq_pres = (interf_liq - interf_sol) * med(pg(intliq_id)%imed)%den0  &
                     * ( - Domain%grav(3)) + pg(intliq_id)%pres + Velocity2 *  &
                     half * med(pg(intliq_id)%imed)%den0
      pretot  = PartLiq_pres + med(pg(intsol_id)%imed)%den0 *                  &
                ( - Domain%grav(3)) * (interf_sol - pg(npi)%Coord(3)) 
endif
! To compute velocity 
if (Velocity2==zero) then
   pg(npi)%state = "sol"
   pg(npi)%vel = zero
   pg(npi)%var = zero
   pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita *             &
                  Med(imed)%celerita))
   if (Granular_flows_options%ID_erosion_criterion==1) then
      pg(npi)%sigma_prime_m = 0.0d0
      pg(npi)%pres_fluid = 0.0d0
   endif
   return 
end if
! To compute roughness (k): c_4=k/d_50 (from input, Shields 2D) or c_4=3 
! (Shields 3D)
if (Granular_flows_options%ID_erosion_criterion==1) then
   Ks = 3.d0 * med(imed)%d_90
   else
      Ks = med(imed)%RoughCoef * med(imed)%d50
endif
if (Granular_flows_options%ID_erosion_criterion>1) then     
   Velocity = Dsqrt(Velocity2)
   else
      if (Granular_flows_options%erosion_flag==2) then
! Default value for interface normal (in case of no free surface)           
         if (ind_interfaces(igridi,jgridi,1)==0) then 
! At this point the neighbour can be only a mixture particle
            pg(npi)%normal_int(1) = 0.d0
            pg(npi)%normal_int(2) = 0.d0
            pg(npi)%normal_int(3) = 1.d0
         endif
      endif
! To compute the velocity component perpendicular to the interface normal       
      vel_s(:) = pg(intliq_id)%vel_old(:) -                                    &
                 dot_product(pg(intliq_id)%vel_old,pg(npi)%normal_int) *       &
                 pg(npi)%normal_int(:)  
      Velocity = dsqrt(dot_product(vel_s,vel_s))
endif             
! To compute the distance between the computational particle and the interacting 
! particle (the height difference does not work for inclined beds)
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
Ustar = dsqrt(Velocity * pg(intliq_id)%mu / (Med(pg(intliq_id)%imed)%den0 *    &
        DistZmin))
! Cycle to determine "Z0" and "Ustar", acording to Zhou Liu's "Sediment 
! transport" (page 9).
flagz0 = .true.
DifRel = one
contatore = 0
iter_ustar: do while ((flagz0).and.                                            &
   (DifRel>Granular_flows_options%conv_crit_erosion))
! Check the iteration number 
   contatore = contatore + 1
   if (contatore>Granular_flows_options%n_max_iterations) then
      inquire(file=nomefileerr,EXIST=err_flag)
      if (.not.err_flag) open(unit=uniterr,file=nomefileerr,form='formatted',  &
         status='unknown')
      write (nout,*) ' WARNING! Superato il numero massimo di iterazioni per', &
         ' il calcolo di Ustar. ' 


      write (nscr,*) ' WARNING! Superato il numero massimo di iterazioni per', &
         ' il calcolo di Ustar. ' 
      write (uniterr,*) '   ' 
      write (uniterr,*) ' WARNING! Superato il numero massimo di iterazioni ', &
         ' per il calcolo di Ustar. ' 
      write (uniterr,*) '  Provare a ridurre il valore di Ks e rilanciare ',   &
         '(for Shields 2D or Mohr erosion criterion). ' 
      write (uniterr,*) '  simulation_time = ',simulation_time,'  Velocity',   &
         Velocity
      write (uniterr,*) '  Ustar',Ustar,'  Ustarold',Ustarold
      write (uniterr,'(a,i6,a,3f10.7,a,e10.3)') '  npi = ',npi,                &
         '      pg(npi)%coord = ',pg(npi)%coord(1),pg(npi)%coord(2),           &
         pg(npi)%coord(3)," residual = ",DifRel
      write (uniterr,*) '  DistZmin',DistZmin,'  Z0',Z0
! Rough boundary regime
      Z0 = 0.033d0 * Ks
      flagz0 = .false.
      if (DistZmin > Z0) Ustar = vKconst * Velocity / Dlog(DistZmin / Z0)
      exit iter_ustar
   end if
   nu_ustar = pg(intliq_id)%mu / (Med(pg(intliq_id)%imed)%den0 * Ustar)
   RestarC = Ks / nu_ustar
   deltanu = 11.6d0 * nu_ustar
   if ((deltanu>DistZmin).and.(RestarC<5.0d0)) then
      flagz0 = .false.
      else
         if (RestarC>=70.0d0) then
! Rough boundary regime (Nikuradse's formulas)
            Z0 = 0.033d0 * Ks
            else if (RestarC<=5.0d0) then
! Smooth boundary regime (Prandtl's law)
               Z0 = 0.11d0 * nu_ustar
               else
! Internediate boundary regime (Manenti et al., 2012, JHE)
                  Z0 = 0.11d0 * nu_ustar + 0.033d0 * Ks
         end if
         Ustarold = Ustar
         if (DistZmin>Z0) then
            Ustar = vKconst * Velocity / Dlog(DistZmin / Z0)
            DifRel = dabs(Ustar-Ustarold) / Ustar
            else
               flagz0 = .false.
               inquire (file=nomefileerr,EXIST=err_flag)
               if (.not. err_flag) open(unit=uniterr,file=nomefileerr,         &
                  form='formatted',status='unknown')
               write (nout,*) ' WARNING! Z0>rijtempmin per il calcolo di Ustar.' 
               write (nscr,*) ' WARNING! Z0>rijtempmin per il calcolo di Ustar.' 
               write (uniterr,*) '   ' 
               write (uniterr,*) ' WARNING! Z0>rijtempmin per il calcolo di ', &
                  ' Ustar. ' 
               write (uniterr,*) '  simulation_time = ',simulation_time,       &
                  '  Velocity',Velocity
               write (uniterr,*) '  Ustar',Ustar,'  Ustarold',Ustarold
               write (uniterr,'(a,i6,a,3f10.7)') '  npi = ',npi,               &
                  '      pg(npi)%coord = ',pg(npi)%coord(1),pg(npi)%coord(2),  &
                  pg(npi)%coord(3)
               write (uniterr,*) '  DistZmin',DistZmin,'  Z0',Z0
         end if
   endif
enddo iter_ustar 
if (isnan(Ustar)) then
   call diagnostic(arg1=11,arg2=1,arg3=nomsub)
end if
Taub = Ustar * Ustar * Med(pg(intliq_id)%imed)%den0
Restar = Med(pg(intliq_id)%imed)%den0 * Ustar * med(imed)%d50 / pg(intliq_id)%mu
if (Restar>=500.0D0) then
   Tetacr = 0.068D0
   elseif (Restar<=1.d0) then
      Tetacr = 0.1d0    
      else
         Tetacr = 0.010595D0 * Dlog(Restar) + 0.110476D0 / Restar + 0.0027197D0
end if
! To compute Shields critical parameter for flat beds 
Taubcror = Tetacr * GI * med(imed)%d50 * (Med(pg(intsol_id)%imed)%den0_s -     &
           Med(pg(intliq_id)%imed)%den0)
if (Granular_flows_options%ID_erosion_criterion==1) then
   pg(npi)%u_star = Ustar
   call compute_k_BetaGamma(npi,intliq_id,DistZmin)
   Taubcr = Taubcror * pg(npi)%k_BetaGamma
   else
! Slope angle
! "vel(1)" is the velocity component aligned with the main flow.
      beta = Dabs(Datan(pg(intliq_id)%vel(3) / pg(intliq_id)%vel(1)))  
      phi = med(imed)%phi
      if (beta>phi) then
         Kbeta = zero 
         else
            if (pg(intliq_id)%vel(3)<=zero) then 
               Kbeta = Dsin(phi - beta) / Dsin(phi)  
               else
                  Kbeta = Dsin(phi + beta) / Dsin(phi)   
            end if
      end if
      if (ncord==2) then
         Kgamma = one
         else
            gamma = Dabs(Datan(pg(intliq_id)%vel(3) / pg(intliq_id)%vel(2)))
            Kgamma = one - Dtan(gamma) * Dtan(gamma) / (Dtan(phi) *            &
                     Dtan(phi))
            if (kgamma>zero) then
               Kgamma = Dcos(gamma) * Dsqrt(kgamma)
               else
                  kgamma = zero
            end if
      end if
      kgamma = one
! 2D critical Shields parameter, according to Zhang (2006).
      Taubcr = Taubcror * Kbeta * Kgamma
endif
if (Taubcr>1.d-9) then
    pg(npi)%tau_tauc = Taub / Taubcr
    if (pg(npi)%tau_tauc>99999.d0) pg(npi)%tau_tauc = 99999.d0
    else
       pg(npi)%tau_tauc = 99999.d0  
endif
if ((Taub>Taubcr).and.(on_going_time_step>Med(imed)%NIterSol)) then 
   if (pg(npi)%state.ne."flu") then
      pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita *          &
                     Med(imed)%celerita))
      pg(npi)%state = "flu"
   endif   
   else
! To check if particles are close to outlet sections.
      if (pg(npi)%CloseBcOut==0) then
         pg(npi)%state = "sol"
         pg(npi)%vel = zero
         pg(npi)%var = zero
         pg(npi)%dens = med(imed)%den0 + (pretot / (Med(imed)%celerita *       &
                        Med(imed)%celerita))
         if (diffusione) pg(npi)%dens = med(imed)%den0 
! Initializing some SPH parameters for fixed mixture particles.
         if (Granular_flows_options%ID_erosion_criterion==1) then
            pg(npi)%sigma_prime_m = 0.0d0
            pg(npi)%pres_fluid = 0.0d0
         endif
      endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine Shields

