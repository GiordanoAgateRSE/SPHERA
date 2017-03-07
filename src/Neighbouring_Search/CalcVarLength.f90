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
! Program unit: CalcVarLength            
! Description: Neighbouring search (pre-conditioned dynamic vector), relative 
!              positions, kernel functions/derivatives, Shepard's coefficient, 
!              position of the fluid-sediment interfaces along each background 
!              grid column.                                             
!-------------------------------------------------------------------------------
subroutine CalcVarLength
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
integer(4) :: nceli,igridi,kgridi,jgridi,irang,krang,jrang,ncelj,jgrid1
integer(4) :: jgrid2,contliq,mm,npi,npj,npartint,index_rij_su_h
integer(4) :: irestocell,celleloop,fw,i_grid,j_grid,bp,bp_f,aux2,i
double precision :: rij_su_h,ke_coef,kacl_coef,rij_su_h_quad,rijtemp,rijtemp2
double precision :: gradmod,gradmodwacl,wu,denom,normal_int_abs,abs_vel
double precision :: min_sigma_Gamma,dis_fp_dbsph_inoutlet
double precision :: dbsph_inoutlet_threshold,normal_int_mixture_top_abs
double precision :: normal_int_sat_top_abs
double precision :: ragtemp(3)
integer(4),dimension(:),allocatable :: bounded
double precision,dimension(:),allocatable :: dShep_old
character(len=lencard)  :: nomsub = "CalcVarLength"
integer(4),external :: ParticleCellNumber, CellIndices, CellNumber
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
celleloop = 4 + (ncord - 2) * 10 + (3 - ncord)
ke_coef = Domain%coefke / Domain%h
kacl_coef = Domain%coefkacl / Domain%h
nPartIntorno = 0
ind_interfaces = 0
if ((Domain%tipo=="bsph").and.(nag>0)) then
   nPartIntorno_fw = 0
   grad_vel_VSL_fw = 0.d0
   allocate(bounded(nag))
   bounded = 0
   allocate(dShep_old(nag))
   pg_w(:)%sigma = 0.d0
   pg_w(:)%kin_visc_semi_part = 0.d0
   pg_w(:)%grad_vel_VSL_times_mu(1) = 0.d0
   pg_w(:)%grad_vel_VSL_times_mu(2) = 0.d0
   pg_w(:)%grad_vel_VSL_times_mu(3) = 0.d0
endif
if (n_bodies>0) then 
   nPartIntorno_bp_f = 0
   nPartIntorno_bp_bp = 0
   aux2 = 0
endif
!------------------------
! Statements
!------------------------
!$omp parallel do default(none)                                                &
!$omp private(npi,nceli,irestocell,igridi,jgridi,kgridi,jgrid1,jgrid2,irang)   &
!$omp private(jrang,krang,mm,npj,npartint,ncelj,ragtemp,rijtemp,rijtemp2)      &
!$omp private(rij_su_h,rij_su_h_quad,denom,index_rij_su_h,gradmod,gradmodwacl) &
!$omp private(wu,contliq,fw,normal_int_abs,abs_vel,dis_fp_dbsph_inoutlet)      &
!$omp private(dbsph_inoutlet_threshold,normal_int_mixture_top_abs)             &
!$omp private(normal_int_sat_top_abs)                                          &                 
!$omp shared(nag,pg,Domain,Med,Icont,Npartord,NMAXPARTJ,rag,nPartIntorno)      &
!$omp shared(Partintorno,PartKernel,ke_coef,kacl_coef,Doubleh,DoubleSquareh)   &
!$omp shared(squareh,nomsub,ncord,eta,eta2,nout,nscr,ind_interfaces,DBSPH,pg_w)&
!$omp shared(Icont_w,Npartord_w,rag_fw,nPartIntorno_fw,Partintorno_fw)         &
!$omp shared(kernel_fw,dShep_old,Granular_flows_options,NMedium)
loop_nag: do npi=1,nag
   if (Domain%tipo=="bsph") then
      pg(npi)%rhoSPH_old = pg(npi)%rhoSPH_new
      pg(npi)%rhoSPH_new = zero
      dShep_old(npi) = pg(npi)%dShep
      pg(npi)%dShep = zero 
      pg(npi)%sigma = zero 
      pg(npi)%sigma_same_fluid = zero
      pg(npi)%FS = 0
      pg(npi)%DBSPH_inlet_ID = 0
      pg(npi)%DBSPH_outlet_ID = 0
      do npj=DBSPH%n_w+1,DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet  
         dis_fp_dbsph_inoutlet = dsqrt(dot_product(pg(npi)%coord -             &
                                 pg_w(npj)%coord,pg(npi)%coord -               &
                                 pg_w(npj)%coord)) 
         if (npj<=(DBSPH%n_w+DBSPH%n_inlet)) then
            dbsph_inoutlet_threshold =                                         &
               DBSPH%inlet_sections(npj-DBSPH%n_w,10) / dsqrt(2.d0) 
            else
               dbsph_inoutlet_threshold =                                      &
                  DBSPH%outlet_sections(npj-DBSPH%n_w-DBSPH%n_inlet,7) /       &
                  dsqrt(2.d0) 
         endif
         if (dis_fp_dbsph_inoutlet<=dbsph_inoutlet_threshold) then
            if (npj<=(DBSPH%n_w+DBSPH%n_inlet)) then
               pg(npi)%DBSPH_inlet_ID = npj
               else
                  pg(npi)%DBSPH_outlet_ID = npj    
            endif
         endif      
      enddo
   endif
   if (Granular_flows_options%ID_erosion_criterion>0) then
      pg(npi)%normal_int(:) = 0.d0
      pg(npi)%normal_int_mixture_top(:) = 0.d0
      pg(npi)%normal_int_sat_top(:) = 0.d0
   endif
   nceli = pg(npi)%cella
   if (nceli==0) cycle
   irestocell = CellIndices (nceli,igridi,jgridi,kgridi)
   contliq = 0
   pg(npi)%indneighliqsol = 0
   pg(npi)%ind_neigh_mix_bed = 0
   pg(npi)%ind_neigh_mob_for_granmob = 0
   pg(npi)%blt_flag = 0
   pg(npi)%rijtempmin = 99999.d0
   jgrid1 = jgridi - (ncord - 2)
   jgrid2 = jgridi + (ncord - 2)
! Loop over the 9 cells in x-direction 
   loop_jrang: do jrang=jgrid1,jgrid2     
! Loop over the 9/1 cells in y-direction in 3/2D
      loop_irang: do irang=igridi-1,igridi+1    
! Loop over the 9 cells in z-direction  
         loop_krang: do krang=kgridi-1,kgridi+1   
            ncelj = CellNumber (irang,jrang,krang)
! Cell out of domain
            if (ncelj==0) cycle    
! Loop over the particles in the cell 
            loop_mm: do mm=Icont(ncelj),Icont(ncelj+1)-1 
! Warning. A computational particle has reached the maximum number of   
! fluid particle neighbours allowed.                
               if (nPartIntorno(npi)>=NMAXPARTJ) then
                  write (nscr,'(1x,a,i12,a,i12,a)')                            &
                     ' The computational particle ',npi,' has reached ',       &
                     NMAXPARTJ,' fluid particle neighbours.'
                  write (nscr,'(1x,a,3f25.10)') '        Coordinate: ',        &
                     pg(npi)%coord(:)
                  cycle
               endif
! Warning. A computational particle has reached the maximum number of   
! wall element neighbours allowed.
               if ((Domain%tipo=="bsph").and.(DBSPH%n_w>0)) then
                  if (nPartIntorno_fw(npi)>=NMAXPARTJ) then
                     write (nscr,'(1x,a,i12,a,i12,a)')                         &
                        ' The computational particle ',npi,' has reached ',    &
                        NMAXPARTJ,' wall element neighbours.'
                     write (nscr,'(1x,a,3f15.10)') '        Coordinate: ',     &
                        pg(npi)%coord(:)
                     cycle
                  endif
               endif
               if (Icont(ncelj+1)<=Icont(ncelj)) cycle
               npj = NPartOrd(mm)
! Relative distance 
               ragtemp(1:3) = pg(npi)%coord(1:3) - pg(npj)%coord(1:3)
! Square of the relative distance (to improve accuracy later on)
               rijtemp = ragtemp(1)*ragtemp(1) + ragtemp(2)*ragtemp(2) +       &
                         ragtemp(3)*ragtemp(3)
! Check if the second particle is a neighbour of the first one 
               if (rijtemp>doublesquareh) cycle
! Increase of the number of neighburs  
               nPartIntorno(npi) = nPartIntorno(npi) + 1
! Increase an array pointer 
               npartint = (npi - 1) * NMAXPARTJ + nPartIntorno(npi)
! Error. A computational particle has exceeded the maximum number of fluid  
! particle neighbours allowed. Simulation is killed.               
               if (nPartIntorno(npi)>NMAXPARTJ) then
                  write (nscr,'(1x,a,i12,a)')   ' ERROR! The particle ',npi,   &
                     ' has too many fluid particle neighbours.'
                  write (nscr,'(1x,a,3f15.10)') '        Coordinate: ',        &
                     pg(npi)%coord(:)
                  call diagnostic(arg1=10,arg2=1,arg3=nomsub)
               endif
               PartIntorno(npartint) = npj
               rijtemp2 = rijtemp
               rijtemp = Dsqrt(rijtemp)
               rij_su_h = rijtemp / Domain%h
               rij_su_h_quad = rijtemp2 / squareh
               index_rij_su_h = int(rij_su_h)
               denom = one / (rijtemp + eta)
               rag(1:3,npartint) = ragtemp(1:3)
               gradmod = zero
               gradmodwacl = zero
               wu = zero
               PartKernel(1:4,npartint) = zero
               if (index_rij_su_h>=2) cycle
               gradmod = -two * rij_su_h + 1.5d0 * rij_su_h_quad
               gradmodwacl = -12.0d0 - 3.0d0 * rij_su_h_quad + 12.0d0 * rij_su_h 
               wu = 0.666666667d0 + rij_su_h_quad * (rij_su_h * half - one)
               if (index_rij_su_h>0) then
                  gradmod = -gradmod + rij_su_h_quad - two
                  wu = (two - rij_su_h) * (two - rij_su_h) * (two - rij_su_h) *& 
                       0.166666666667d0
               endif 
               gradmod = gradmod * ke_coef
               gradmodwacl = gradmodwacl * kacl_coef
               PartKernel(1,npartint) = gradmod * denom 
               PartKernel(2,npartint) = PartKernel(1,npartint) / (rijtemp2 +   &
                                        eta2)
               PartKernel(3,npartint) = gradmodwacl * denom
               PartKernel(4,npartint) = wu * Domain%coefke
! WendlandC4 kernel, interesting test: start
! PartKernel(1,npartint) = (3.d0/(4.d0*PIGRECO*(2.d0*Domain%h**3))) *          &
!    ((1-rij_su_h/2.d0)**5) * (-280.d0 * (rij_su_h/2.d0)**2 - 56.d0 *          &
!    (rij_su_h/2.d0))
! if (rij_su_h/=0.d0) PartKernel(1,npartint) = PartKernel(1,npartint) * denom
! PartKernel(4,npartint) = (3.d0/(4.d0*PIGRECO*(Domain%h**2))) *               &
!    ((1-rij_su_h/2.d0)**6) * (35.d0 * ((rij_su_h/2.d0)**2) + 18.d0 *          &
!    (rij_su_h/2.d0) + 3.d0)
! WendlandC4 kernel, interesting test: end
! Gallati anti-cluster kernel, interesting test: start
!    kernel_fw(2,npartint) = (5.d0/(16.d0*PIGRECO*Domain%h**2)) *              &
!    ((2.d0 - rij_su_h)**3) 
! Gallati anti-cluster kernel, interesting test: end
               if (Domain%tipo=="bsph") then
!AA!!! test start
                  pg(npi)%sigma = pg(npi)%sigma + pg(npj)%mass *               &
                                  PartKernel(4,npartint) / pg(npj)%dens 
                  if (pg(npi)%imed==pg(npj)%imed) then
                     pg(npi)%rhoSPH_new = pg(npi)%rhoSPH_new + pg(npj)%mass *  &
                                          PartKernel(4,npartint)
                     else
                        pg(npi)%rhoSPH_new = pg(npi)%rhoSPH_new + pg(npj)%mass &
                                             / pg(npj)%dens * pg(npi)%dens *   &
                                             PartKernel(4,npartint)
                  endif
!                  pg(npi)%sigma = pg(npi)%sigma + pg(npj)%mass *               &
!                                  PartKernel(4,npartint) / pg(npj)%dens 
!                  if (pg(npi)%imed==pg(npj)%imed) then
!                     pg(npi)%rhoSPH_new = pg(npi)%rhoSPH_new + pg(npj)%mass *  &
!                                          PartKernel(4,npartint)
!                     if (NMedium>1) then
!                        pg(npi)%sigma_same_fluid = pg(npi)%sigma_same_fluid +  &
!                                                   pg(npj)%mass *              &
!                                                   PartKernel(4,npartint) /    &
!                                                   pg(npj)%dens
!                     endif
!                  endif
!AA!!! test end
               endif
! In case of bed-load transport with/without any erosion criterion
               if (Granular_flows_options%ID_erosion_criterion>0) then
! Searching for the nearest fluid/mixture SPH particle 
                  if (Med(pg(npi)%imed)%tipo/=Med(pg(npj)%imed)%tipo) then
                     if ((rijtemp<pg(npi)%rijtempmin(1)).or.                   &
                        ((rijtemp==pg(npi)%rijtempmin(1)).and.                 &
                        (npj>pg(npi)%indneighliqsol))) then   
                        if ((index(Med(pg(npi)%imed)%tipo,"liquid")>0).and.    &
                           (pg(npi)%coord(3)>pg(npj)%coord(3))) then
                           pg(npi)%indneighliqsol = npj
                           pg(npi)%rijtempmin(1) = rijtemp
                           elseif ((index(Med(pg(npi)%imed)%tipo,"granular")>0)&
                              .and.(pg(npi)%coord(3)<pg(npj)%coord(3))) then
                              pg(npi)%indneighliqsol = npj
                              pg(npi)%rijtempmin(1) = rijtemp
                        endif
                     endif
                  endif
! Searching for the nearest mobile/fixed particle (if cycle apparently less 
! efficient than the previous one because we have to compute the interface 
! normal)              
                  if (pg(npi)%state/=pg(npj)%state) then
                     if (pg(npi)%state=="flu") then
                        if (index(Med(pg(npi)%imed)%tipo,"granular")>0)  &
                           pg(npi)%normal_int(:) = pg(npi)%normal_int(:)       &
                                                - (pg(npj)%coord(:) -          &
                                                pg(npi)%coord(:)) *            &
                                                PartKernel(4,npartint) 
                        if ((rijtemp<pg(npi)%rijtempmin(2)).or.                &
                           ((rijtemp==pg(npi)%rijtempmin(2)).and.              &
                           (npj>pg(npi)%ind_neigh_mix_bed))) then  
                           pg(npi)%ind_neigh_mix_bed = npj
                           pg(npi)%rijtempmin(2) = rijtemp
                        endif
                        elseif (pg(npi)%state=="sol") then
                           pg(npi)%normal_int(:) = pg(npi)%normal_int(:) +     &
                                                   (pg(npj)%coord(:) -         &
                                                   pg(npi)%coord(:)) *         &
                                                   PartKernel(4,npartint)  
                           if ((rijtemp<pg(npi)%rijtempmin(2)).or.             &
                              ((rijtemp==pg(npi)%rijtempmin(2)).and.           &
                              (npj>pg(npi)%ind_neigh_mix_bed))) then 
                              pg(npi)%ind_neigh_mix_bed = npj
                              pg(npi)%rijtempmin(2) = rijtemp
                           endif
                     endif
                  endif
! Searching for the nearest mobile granular/mobile particle
                  if ((index(Med(pg(npi)%imed)%tipo,"granular")>0).and.        &
                     (pg(npi)%state=="flu").and.(pg(npj)%state=="flu")) then
                     if ((rijtemp<pg(npi)%rijtempmin(3)).or.                   &
                        ((rijtemp==pg(npi)%rijtempmin(3)).and.                 &
                        (npj>pg(npi)%ind_neigh_mob_for_granmob))) then 
                        if (pg(npi)%coord(3)<pg(npj)%coord(3)) then
                           pg(npi)%ind_neigh_mob_for_granmob = npj
                           pg(npi)%rijtempmin(3) = rijtemp
                        endif
                     endif
                  endif              
! To estimate the normal to the mixture top
                  if (index(Med(pg(npi)%imed)%tipo,"granular")>0) then
                     if (Med(pg(npi)%imed)%tipo==Med(pg(npj)%imed)%tipo) then
                        pg(npi)%normal_int_mixture_top(:) =                    &
                           pg(npi)%normal_int_mixture_top(:) -                 &
                           (pg(npj)%coord(:) - pg(npi)%coord(:)) *             &
                           PartKernel(4,npartint)
                     endif
                  endif
! To estimate the normal to the top of the saturated zone
                  if ((Med(pg(npi)%imed)%saturated_medium_flag.eqv..true.).and.&
                     (Med(pg(npi)%imed)%saturated_medium_flag.eqv.             &
                     Med(pg(npj)%imed)%saturated_medium_flag)) then
                     pg(npi)%normal_int_sat_top(:) =                           &
                        pg(npi)%normal_int_sat_top(:) - (pg(npj)%coord(:) -    &
                        pg(npi)%coord(:)) * PartKernel(4,npartint)
                  endif
               endif
            enddo loop_mm
! Loop over the neighbouring wall particles in the cell
            if ((Domain%tipo=="bsph").and.(DBSPH%n_w>0)) then
               loop_fw: do fw=Icont_w(ncelj),Icont_w(ncelj+1)-1
                  npj = NPartOrd_w(fw)
! Relative positions and distances
                  ragtemp(1:3) = pg(npi)%coord(1:3) - pg_w(npj)%coord(1:3)
                  rijtemp = ragtemp(1) * ragtemp(1) + ragtemp(2) * ragtemp(2) +&
                            ragtemp(3)*ragtemp(3)
! Distance check
                  if (rijtemp>doublesquareh) cycle
! Counter increases 
                  nPartIntorno_fw(npi) = nPartIntorno_fw(npi) + 1
! Error. A computational particle has exceeded the maximum number of wall  
! element neighbours allowed. Simulation is killed.               
                  if (nPartIntorno_fw(npi)>NMAXPARTJ) then
                     write (nscr,'(1x,a,i12,a)')   ' ERROR! The particle ',npi,&
                        ' has too many wall element neighbours.'
                     write (nscr,'(1x,a,3f15.10)') '        Coordinate: ',     &
                        pg(npi)%coord(:)
                     call diagnostic(arg1=10,arg2=9,arg3=nomsub)
                  endif                  
                  npartint = (npi - 1) * NMAXPARTJ + nPartIntorno_fw(npi)
! Savings
                  PartIntorno_fw(npartint) = npj
                  rijtemp2 = rijtemp
                  rij_su_h = Dsqrt(rijtemp) / Domain%h
                  rij_su_h_quad = rijtemp2 / squareh
                  index_rij_su_h = int(rij_su_h)
                  rag_fw(1:3,npartint) = ragtemp(1:3)
! Kernel computation
                  if (index_rij_su_h>=2) cycle
                  wu = 0.666666667d0 + rij_su_h_quad * (rij_su_h * half - one)
                  if (index_rij_su_h>0) then
                     wu = (two - rij_su_h) * (two - rij_su_h) *                &
                          (two - rij_su_h) * 0.166666666667d0
                  endif 
                  kernel_fw(1,npartint) = wu * Domain%coefke
                  pg(npi)%dShep = pg(npi)%dShep + kernel_fw(1,npartint) *      &
                                  pg_w(npj)%weight * (pg_w(npj)%normal(1) *    &
                                  (pg_w(npj)%vel(1) - pg(npi)%var(1)) +        &
                                  pg_w(npj)%normal(2) * (pg_w(npj)%vel(2) -    &
                                  pg(npi)%var(2)) + pg_w(npj)%normal(3) *      &
                                  (pg_w(npj)%vel(3) - pg(npi)%var(3))) 
                  if ((rij_su_h*Domain%h)<=                                    &
                     (1.3d0*Domain%dx+pg_w(npj)%weight/2.d0)) then
                     pg_w(npj)%wet = 1
                  endif
                  denom = one / (Dsqrt(rijtemp) + eta)
                  gradmod = -two * rij_su_h + 1.5d0 * rij_su_h_quad
                  if (index_rij_su_h>0) then
                     gradmod = - gradmod + rij_su_h_quad - two
                  endif 
                  kernel_fw(2,npartint) = gradmod * ke_coef * denom
                  pg(npi)%sigma = pg(npi)%sigma +  pg_w(npj)%volume *          &
                                  kernel_fw(1,npartint)
!AA!!!test start  
                  if (NMedium==1) then
                     pg(npi)%rhoSPH_new = pg(npi)%rhoSPH_new + pg_w(npj)%mass *&
                                          kernel_fw(1,npartint)
                     else
                        pg(npi)%rhoSPH_new = pg(npi)%rhoSPH_new + pg(npi)%dens &
                                             * pg_w(npj)%volume *              &
                                             kernel_fw(1,npartint)                        
                  endif                                
!                  if (NMedium==1) then
!                     pg(npi)%rhoSPH_new = pg(npi)%rhoSPH_new + pg_w(npj)%mass *&
!                                          kernel_fw(1,npartint)
!                     else
!                        pg(npi)%rhoSPH_new = pg(npi)%rhoSPH_new + pg(npi)%dens &
!                                             * pg_w(npj)%volume *              &
!                                             kernel_fw(1,npartint)                        
!                        pg(npi)%sigma_same_fluid = pg(npi)%sigma_same_fluid +  &
!                                                   pg_w(npj)%volume *          &
!                                                   kernel_fw(1,npartint)
!                  endif
!AA!!!test end                  
! Gallati anti-cluster kernel, interesting test: start
! kernel_fw(3,npartint) =(5.d0/(16.d0*PIGRECO*Domain%h**2))*((2.d0 - rij_su_h)**3) 
! kernel_fw(4,npartint) = (-12.0d0 - 3.0d0 * rij_su_h_quad + 12.0d0 * rij_su_h)&
!    * kacl_coef * denom
! Gallati anti-cluster kernel, interesting test: end
! WendlandC4 kernel, interesting test: start
! kernel_fw(1,npartint) = (3.d0/(4.d0*PIGRECO*(Domain%h**2))) *                &
!    ((1.d0 - rij_su_h/2.d0)**6) * (35.d0 * ((rij_su_h/2.d0)**2) + 18.d0*        &
!    (rij_su_h/2.d0) + 3.d0)
! kernel_fw(2,npartint) = (3.d0/(4.d0*PIGRECO*(2.d0*Domain%h**3))) *           &
!    ((1.d0 - rij_su_h/2.d0)**5) * (-280.d0 * (rij_su_h/2.d0)**2 - 56.d0 *       &
!    (rij_su_h/2.d0))
! if (rij_su_h/=0.d0) kernel_fw(2,npartint) = kernel_fw(2,npartint) * denom
! WendlandC4 kernel, interesting test: end
                  if (DBSPH%slip_ID>0)                                         &
                     call DBSPH_velocity_gradients_VSL_SNBL(npi,npj,npartint)
               enddo loop_fw
            endif
         enddo loop_krang      
      enddo loop_irang       
   enddo loop_jrang   
   if (Granular_flows_options%ID_erosion_criterion>0) then
! Free surface detection along the grid column
      if (index(Med(pg(npi)%imed)%tipo,"liquid")>0) then
!$omp critical (free_surface_detection)
         if (ind_interfaces(igridi,jgridi,1)==0) then
            ind_interfaces(igridi,jgridi,1) = npi   
            else
               if (pg(npi)%coord(3)>                                           &
                  pg(ind_interfaces(igridi,jgridi,1))%coord(3)) then
                  ind_interfaces(igridi,jgridi,1) = npi 
                  elseif (pg(npi)%coord(3)==                                   &
                     pg(ind_interfaces(igridi,jgridi,1))%coord(3)) then
                     if (npi>ind_interfaces(igridi,jgridi,1))                  &
                        ind_interfaces(igridi,jgridi,1) = npi
               endif
         endif
!$omp end critical (free_surface_detection)        
         else
! In case of no free surface in the column and no erosion criterion, updating 
! the mixture - fixed bed interface       
            if (Granular_flows_options%erosion_flag<2) then
!$omp critical (fixed_bed_detection)
               abs_vel = dsqrt(dot_product(pg(npi)%vel,pg(npi)%vel))
               if (abs_vel<=Granular_flows_options%velocity_fixed_bed) then
                  if (ind_interfaces(igridi,jgridi,4)==0) then
                     ind_interfaces(igridi,jgridi,4) = npi   
                     else
                        if (pg(npi)%coord(3)>                                  &
                           pg(ind_interfaces(igridi,jgridi,4))%coord(3)) then
                              ind_interfaces(igridi,jgridi,4) = npi 
                              elseif (pg(npi)%coord(3)==                       &
                                 pg(ind_interfaces(igridi,jgridi,4))%coord(3)) &
                                 then
                                 if (npi>ind_interfaces(igridi,jgridi,4))      &
                                    ind_interfaces(igridi,jgridi,4) = npi
                        endif
                  endif
               endif
!$omp end critical (fixed_bed_detection) 
            endif
!$omp critical (soil_bottom_detection)
            if (ind_interfaces(igridi,jgridi,5)==0) then
               ind_interfaces(igridi,jgridi,5) = npi   
               else
                  if (pg(npi)%coord(3)<                                        &
                     pg(ind_interfaces(igridi,jgridi,5))%coord(3)) then
                        ind_interfaces(igridi,jgridi,5) = npi 
                        elseif (pg(npi)%coord(3)==                             &
                           pg(ind_interfaces(igridi,jgridi,5))%coord(3))       &
                           then
                           if (npi>ind_interfaces(igridi,jgridi,5))            &
                              ind_interfaces(igridi,jgridi,5) = npi
                  endif
            endif
!$omp end critical (soil_bottom_detection)
            if (Med(pg(npi)%imed)%saturated_medium_flag.eqv..true.) then
!$omp critical (saturated_zone_top_detection)
               if (ind_interfaces(igridi,jgridi,6)==0) then
                  ind_interfaces(igridi,jgridi,6) = npi   
                  else
                     if (pg(npi)%coord(3)>                                     &
                        pg(ind_interfaces(igridi,jgridi,6))%coord(3)) then
                           ind_interfaces(igridi,jgridi,6) = npi 
                           elseif (pg(npi)%coord(3)==                          &
                              pg(ind_interfaces(igridi,jgridi,6))%coord(3))    &
                              then
                              if (npi>ind_interfaces(igridi,jgridi,6))         &
                                 ind_interfaces(igridi,jgridi,6) = npi
                     endif
               endif
!$omp end critical (saturated_zone_top_detection)
            endif
      endif                                     
   endif
! In case of bed-load transport with an erosion criterion
   if ((Granular_flows_options%ID_erosion_criterion>0).and.                    &
      (Granular_flows_options%erosion_flag/=1)) then
! Normalization of the interface normal between the granular mixture and the 
! fixed bed (bed-load transport)
      if (index(Med(pg(npi)%imed)%tipo,"granular")>0) then 
         normal_int_abs = dsqrt(dot_product(pg(npi)%normal_int,                &
                          pg(npi)%normal_int)) 
         if (normal_int_abs>zero) then
            pg(npi)%normal_int(:) = pg(npi)%normal_int(:)/normal_int_abs
         endif
      endif    
   endif
! In case of bed-load transport
   if (Granular_flows_options%ID_erosion_criterion>0) then
      if (index(Med(pg(npi)%imed)%tipo,"granular")>0) then 
! Normalization of the mixture top interface normal 
         normal_int_mixture_top_abs = dsqrt(dot_product(                       &
            pg(npi)%normal_int_mixture_top,                                    &
            pg(npi)%normal_int_mixture_top)) 
         if (normal_int_mixture_top_abs>zero) then
            pg(npi)%normal_int_mixture_top(:) =                                &
               pg(npi)%normal_int_mixture_top(:) / normal_int_mixture_top_abs
            else
               pg(npi)%normal_int_mixture_top(1:2) = 0.d0
               pg(npi)%normal_int_mixture_top(3) = 1.d0
         endif    
! Normalization of the interface normal of the saturated zone top 
         normal_int_sat_top_abs = dsqrt(dot_product(                           &
            pg(npi)%normal_int_sat_top,pg(npi)%normal_int_sat_top)) 
         if (normal_int_sat_top_abs>zero) then
            pg(npi)%normal_int_sat_top(:) = pg(npi)%normal_int_sat_top(:) /    &
                                            normal_int_sat_top_abs
            else
               pg(npi)%normal_int_sat_top(1:2) = 0.d0
               pg(npi)%normal_int_sat_top(3) = 1.d0
         endif
      endif
!$omp critical (interface_definition)  
! Update the local position of the upper interface of the bed-load transport 
! layer (mixture side)
      if ((index(Med(pg(npi)%imed)%tipo,"liquid")>0).and.                      &
         (pg(npi)%indneighliqsol.ne.0)) then 
         if (ind_interfaces(igridi,jgridi,3)==0) then
            ind_interfaces(igridi,jgridi,3) = pg(npi)%indneighliqsol   
            elseif (pg(pg(npi)%indneighliqsol)%coord(3)>                       &
               pg(ind_interfaces(igridi,jgridi,3))%coord(3)) then
               ind_interfaces(igridi,jgridi,3) = pg(npi)%indneighliqsol 
               elseif (pg(pg(npi)%indneighliqsol)%coord(3)==                   &
                  pg(ind_interfaces(igridi,jgridi,3))%coord(3)) then
                  if (pg(npi)%indneighliqsol>ind_interfaces(igridi,jgridi,3))  &
                     ind_interfaces(igridi,jgridi,3) = pg(npi)%indneighliqsol
         endif
      endif        
! To update the local position of the upper interface of the bed-load transport 
! layer (liquid side)   
      if ((index(Med(pg(npi)%imed)%tipo,"granular")>0).and.                    &
         (pg(npi)%indneighliqsol.ne.0)) then   
         if (ind_interfaces(igridi,jgridi,2)==0) then
            ind_interfaces(igridi,jgridi,2) = pg(npi)%indneighliqsol   
            elseif (pg(pg(npi)%indneighliqsol)%coord(3)>                       &
               pg(ind_interfaces(igridi,jgridi,2))%coord(3)) then
               ind_interfaces(igridi,jgridi,2) = pg(npi)%indneighliqsol 
               elseif (pg(pg(npi)%indneighliqsol)%coord(3)==                   &
                  pg(ind_interfaces(igridi,jgridi,2))%coord(3)) then
                  if (pg(npi)%indneighliqsol>ind_interfaces(igridi,jgridi,2))  &
                     ind_interfaces(igridi,jgridi,2) = pg(npi)%indneighliqsol
         endif
      endif
! To update the local position of the fixed bed    
      if ((pg(npi)%state=="flu").and.(pg(npi)%ind_neigh_mix_bed.ne.0)) then  
         if (ind_interfaces(igridi,jgridi,4)==0) then
            ind_interfaces(igridi,jgridi,4) = pg(npi)%ind_neigh_mix_bed    
            elseif (pg(pg(npi)%ind_neigh_mix_bed)%coord(3)>                    &
               pg(ind_interfaces(igridi,jgridi,4))%coord(3)) then
               ind_interfaces(igridi,jgridi,4) = pg(npi)%ind_neigh_mix_bed
               elseif (pg(pg(npi)%ind_neigh_mix_bed)%coord(3)==                &
                  pg(ind_interfaces(igridi,jgridi,4))%coord(3)) then
                  if (pg(npi)%ind_neigh_mix_bed>                               &
                     ind_interfaces(igridi,jgridi,4))                          &
                     ind_interfaces(igridi,jgridi,4) = pg(npi)%ind_neigh_mix_bed
         endif
      endif
!$omp end critical (interface_definition) 
   endif
enddo loop_nag
!$omp end parallel do
! In case of bed-load transport 
if (Granular_flows_options%ID_erosion_criterion>0) then
! To compute the interface flags
   do npi=1,nag
      nceli = ParticleCellNumber(pg(npi)%coord)
      irestocell = CellIndices(nceli,igridi,jgridi,kgridi)
      if (index(Med(pg(npi)%imed)%tipo,"granular")>0) then
         if (ind_interfaces(igridi,jgridi,3)==0) then
            ind_interfaces(igridi,jgridi,3) = -npi   
            else
               if (ind_interfaces(igridi,jgridi,3)<0) then
                  if (pg(npi)%coord(3)>                                        &
                     pg(-ind_interfaces(igridi,jgridi,3))%coord(3)) then
                     ind_interfaces(igridi,jgridi,3) = -npi 
                     elseif (pg(npi)%coord(3)==                                &
                        pg(-ind_interfaces(igridi,jgridi,3))%coord(3)) then
                        if (npi>(-ind_interfaces(igridi,jgridi,3)))            &
                           ind_interfaces(igridi,jgridi,3) = -npi
                  endif
               endif
         endif
      endif
   enddo
! Initialization of the saturation flags
   if (Granular_flows_options%saturation_scheme==2) then
      if (Granular_flows_options%time_minimum_saturation>=simulation_time) then
         Granular_flows_options%minimum_saturation_flag = .false.
      endif
      if (Granular_flows_options%time_maximum_saturation>=simulation_time) then
         Granular_flows_options%maximum_saturation_flag = .false.
      endif
   endif
!$omp parallel do default(none)                                                &
!$omp shared(Grid,pg,ind_interfaces,nout,Granular_flows_options,simulation_time)&
!$omp private(i_grid,j_grid)
   do i_grid=1,Grid%ncd(1)
      do j_grid=1,Grid%ncd(2)
         if (ind_interfaces(i_grid,j_grid,3)<0) ind_interfaces(i_grid,j_grid,3)&
            = -ind_interfaces(i_grid,j_grid,3)
         if ((ind_interfaces(i_grid,j_grid,3).ne.0).and.                       &
            (ind_interfaces(i_grid,j_grid,1).ne.0)) then
            if (pg(ind_interfaces(i_grid,j_grid,3))%coord(3)>                  &
               pg(ind_interfaces(i_grid,j_grid,1))%coord(3)) then
               ind_interfaces(i_grid,j_grid,3) = ind_interfaces(i_grid,j_grid,1)
            endif
         endif
         if (ind_interfaces(i_grid,j_grid,4).ne.0)                             &
            pg(ind_interfaces(i_grid,j_grid,4))%blt_flag = 3
         if (ind_interfaces(i_grid,j_grid,3).ne.0)                             &
            pg(ind_interfaces(i_grid,j_grid,3))%blt_flag = 2
         if (ind_interfaces(i_grid,j_grid,1).ne.0)                             &
            pg(ind_interfaces(i_grid,j_grid,1))%blt_flag = 1
         if (ind_interfaces(i_grid,j_grid,5).ne.0)                             &
            pg(ind_interfaces(i_grid,j_grid,5))%blt_flag = 4
         if (ind_interfaces(i_grid,j_grid,6).ne.0)                             &
            pg(ind_interfaces(i_grid,j_grid,6))%blt_flag = 5
! Saturation flags
         if (Granular_flows_options%saturation_scheme==2) then
            if ((Granular_flows_options%time_minimum_saturation>=              &
               simulation_time).and.(ind_interfaces(i_grid,j_grid,1)/=0)) then
               Granular_flows_options%minimum_saturation_flag(i_grid,j_grid) = &
                  .true.
            endif
            if ((Granular_flows_options%time_maximum_saturation>=              &
               simulation_time).and.(ind_interfaces(i_grid,j_grid,1)/=0)) then
               Granular_flows_options%maximum_saturation_flag(i_grid,j_grid) = &
                  .true.
            endif
! Saturation conditions
            if (simulation_time<=Granular_flows_options%time_minimum_saturation&
               )then
! Time lower than / equal to time at minimum saturation
               if (Granular_flows_options%minimum_saturation_flag(i_grid,j_grid&
                  ).eqv..true.) then
! Phreatic zone
                  Granular_flows_options%saturation_conditions(i_grid,j_grid) =&
                     1
                  else
! Dry soil
                     Granular_flows_options%saturation_conditions(i_grid,      &
                        j_grid) = 3
               endif
               elseif (simulation_time<                                        &
                  Granular_flows_options%time_maximum_saturation) then
! Time lower than time at minimum saturation and higher than time at maximum 
! saturation 
                  if (Granular_flows_options%minimum_saturation_flag(i_grid,   &
                     j_grid).eqv..true.) then
! Phreatic zone
                     Granular_flows_options%saturation_conditions(i_grid,      &
                        j_grid) = 1
                     elseif (Granular_flows_options%maximum_saturation_flag    &
                        (i_grid,j_grid).eqv..true.) then
! Unsaturated zone
                        Granular_flows_options%saturation_conditions(i_grid,   &
                           j_grid) = 2
                        else
! Dry soil
                           Granular_flows_options%saturation_conditions(i_grid,&
                              j_grid) = 3
                  endif
                  else
! Time higher than time at maximum saturation         
                     if (Granular_flows_options%maximum_saturation_flag(i_grid,&
                        j_grid).eqv..true.) then
! Phreatic zone
                        Granular_flows_options%saturation_conditions(i_grid,   &
                           j_grid) = 1
                        else
! Dry soil
                           Granular_flows_options%saturation_conditions(i_grid,&
                              j_grid) = 3
                     endif
            endif
         endif       
      enddo
   enddo
!$omp end parallel do
endif
if (Domain%tipo=="bsph") then
   do npi=1,nag
! Gamma (integral Shepard coefficient) intialization (but for inlet sections)
      if (on_going_time_step==-2) then
         if (nPartIntorno_fw(npi)==0) then 
            pg(npi)%Gamma = one
            else
!AA!!! test start            
               pg(npi)%Gamma = pg(npi)%sigma
!               if (NMedium==1) then
!                  pg(npi)%Gamma = pg(npi)%sigma
!                  else
!                     pg(npi)%Gamma = pg(npi)%sigma_same_fluid
!               endif               
!AA!!! test end                  
               if (DBSPH%Gamma_limiter_flag.eqv..true.) pg(npi)%Gamma =         &
                  min(pg(npi)%Gamma,one)    
         endif
      endif
      if (on_going_time_step>-2) then
!AA!!! test start      
         min_sigma_Gamma = min((pg(npi)%sigma+0.05),pg(npi)%Gamma)
!         if (NMedium==1) then
!            min_sigma_Gamma = min((pg(npi)%sigma+0.05),pg(npi)%Gamma)
!            else
!               min_sigma_Gamma = min((pg(npi)%sigma_same_fluid+0.05),          &
!                                 pg(npi)%Gamma)
!         endif
!AA!!! test end
         if ((DBSPH%FS_allowed.eqv..true.).and.                                &
            (min_sigma_Gamma/=pg(npi)%Gamma)) then
            pg(npi)%Gamma_last_active = zero
            pg(npi)%FS = 1
!AA!!! test start
            pg(npi)%uni = pg(npi)%sigma            
!            if (NMedium==1) then
!               pg(npi)%uni = pg(npi)%sigma
!               else
!                  pg(npi)%uni = pg(npi)%sigma_same_fluid
!            endif
!AA!!! test end            
            else
               pg(npi)%uni = pg(npi)%Gamma
               pg(npi)%Gamma_last_active = zero
         endif
         if (pg(npi)%rhoSPH_old==zero) then
            pg(npi)%DensShep = pg(npi)%rhoSPH_new * pg(npi)%Gamma
!AA!!! test start
!            if (NMedium==1) then
!               if (pg(npi)%FS==1) then
!                  pg(npi)%dens = pg(npi)%rhoSPH_new / pg(npi)%sigma 
!                  else
!                     pg(npi)%dens = pg(npi)%rhoSPH_new / pg(npi)%Gamma
!               endif
!               elseif(NMedium>1) then
!                  pg(npi)%dens = pg(npi)%rhoSPH_new / pg(npi)%sigma_same_fluid                  
!            endif
            if (pg(npi)%FS==0) then
               pg(npi)%dens = pg(npi)%rhoSPH_new / pg(npi)%Gamma
               else
!AA!!! test2 start
                  pg(npi)%dens = pg(npi)%rhoSPH_new / pg(npi)%sigma  
!                  if (NMedium==1) then
!                     pg(npi)%dens = pg(npi)%rhoSPH_new / pg(npi)%sigma 
!                     else
!                        pg(npi)%dens = pg(npi)%rhoSPH_new /                    &
!                                       pg(npi)%sigma_same_fluid
!                  endif
!AA!!! test2 end
            endif
!AA!!! test end
         endif
      endif
   enddo
   if (allocated(bounded)) deallocate(bounded)
   if (allocated(dShep_old)) deallocate(dShep_old)
endif 
! SPH parameters for body transport in fluid flows
if (n_bodies>0) then
! Loop over the body particles
!$omp parallel do default(none)                                                &
!$omp private(npi,nceli,irestocell,igridi,jgridi,kgridi,jgrid1,jgrid2,irang)   &
!$omp private(jrang,krang,bp_f,npj,npartint,ncelj,ragtemp,rijtemp,rijtemp2)    &
!$omp private(rij_su_h,rij_su_h_quad,denom,index_rij_su_h,gradmod,bp)          &
!$omp private(gradmodwacl,i,aux2)                                              &
!$omp shared(NMAXPARTJ,ke_coef,kacl_coef,DoubleSquareh,squareh,ncord,Domain)   &
!$omp shared(Icont_bp,NPartOrd_bp,bp_arr,nPartIntorno_bp_bp,PartIntorno_bp_bp) &
!$omp shared(rag_bp_bp,n_body_part,Icont,n_bodies,nPartIntorno_bp_f)           &
!$omp shared(PartIntorno_bp_f,rag_bp_f,KerDer_bp_f_cub_spl,KerDer_bp_f_Gal)    &
!$omp shared(NPartOrd,eta,pg)
   do npi=1,n_body_part
! Computation of the ID of the surface body particles
      i = 0
      aux2 = 0
      do while (i<npi) 
         i = i+1 
         if (bp_arr(i)%area>0.d0) aux2 = aux2+1
      enddo
      nceli = bp_arr(npi)%cell
      if (nceli==0) cycle
      irestocell = CellIndices (nceli,igridi,jgridi,kgridi)
      jgrid1 = jgridi - (ncord-2)
      jgrid2 = jgridi + (ncord-2)
! Loop over the adjacent cells
      do jrang = jgrid1,jgrid2         
         do irang = igridi-1,igridi+1      
            do krang = kgridi-1,kgridi+1      
               ncelj = CellNumber (irang,jrang,krang)
               if (ncelj==0) cycle    
! Parameters for body particle - fluid particle  interactions 
! (for body and fluid dynamics)
! Loop over the neighbouring body particles in the cell
               loop_bp_f: do bp_f=Icont(ncelj),Icont(ncelj+1)-1
                  npj = NPartOrd(bp_f)
! Relative positions and distances
! Sign inversion because body particle acts here as a computational particle
                  ragtemp(1:3) = - pg(npj)%coord(1:3) + bp_arr(npi)%pos(1:3)  
                  rijtemp = ragtemp(1) * ragtemp(1) + ragtemp(2) * ragtemp(2)  &
                            + ragtemp(3) * ragtemp(3)
! Distance check
                  if (rijtemp>doublesquareh) cycle
! neighbouring lists 
                  nPartIntorno_bp_f(npi) = nPartIntorno_bp_f(npi) + 1
                  npartint = (npi - 1) * NMAXPARTJ + nPartIntorno_bp_f(npi)
                  PartIntorno_bp_f(npartint) = npj
! Relative distance
                  rijtemp2 = rijtemp
                  rijtemp = Dsqrt(rijtemp)
                  rij_su_h = rijtemp / Domain%h
                  rij_su_h_quad = rijtemp2 / squareh
                  index_rij_su_h = int(rij_su_h)
                  denom = one / (rijtemp + eta)
                  rag_bp_f(1:3,npartint) = ragtemp(1:3)  
                  if (ncord==2) rag_bp_f(2,npartint) = 0.d0
! Kernel gradients (cubic spline)
                  gradmod = zero
                  KerDer_bp_f_cub_spl(npartint) = zero
                  if (index_rij_su_h>=2) cycle
                  gradmod = -two * rij_su_h + 1.5d0 * rij_su_h_quad
                  if (index_rij_su_h>0) then
                     gradmod = -gradmod + rij_su_h_quad - two
                  endif 
                  gradmod = gradmod * ke_coef
                  KerDer_bp_f_cub_spl(npartint) = gradmod * denom 
! Kernel gradients (Gallati's kernel derivative)
                  gradmodwacl = - 12.0d0 - 3.0d0 * rij_su_h_quad + 12.0d0 *    &
                                rij_su_h 
                  gradmodwacl = gradmodwacl * kacl_coef
                  KerDer_bp_f_Gal(npartint) = gradmodwacl * denom 
               enddo loop_bp_f
! End Loop over the neighbouring body particles in the cell
! Loop over the neighbouring body particles in the cell
               loop_bp: do bp=Icont_bp(ncelj),Icont_bp(ncelj+1)-1
                  npj = NPartOrd_bp(bp)
! Only neighbours belonging to a surface of another body
                  if ((bp_arr(npi)%area>0.d0).and.(bp_arr(npj)%area>0.d0)      &
                     .and.(bp_arr(npi)%body/=bp_arr(npj)%body)) then
! Relative positions and distances
                     ragtemp(1:3) = bp_arr(npi)%pos(1:3) - bp_arr(npj)%pos(1:3)  
                     rijtemp = ragtemp(1) * ragtemp(1) + ragtemp(2) *          &
                               ragtemp(2) + ragtemp(3) * ragtemp(3)
! Distance check
                     if (rijtemp>doublesquareh) cycle
! Neighbouring lists 
                     nPartIntorno_bp_bp(aux2) = nPartIntorno_bp_bp(aux2) + 1
                     npartint = (aux2 - 1) * NMAXPARTJ +                       &
                                nPartIntorno_bp_bp(aux2)
                     PartIntorno_bp_bp(npartint) = npj
! Relative distance
                     rag_bp_bp(1:3,npartint) = ragtemp(1:3)  
                     if (ncord==2) rag_bp_bp(2,npartint) = 0.d0   
                  endif
               enddo loop_bp
            enddo 
         enddo 
      enddo 
! End Loop over the neighbouring body particles in the cell       
   enddo
!$omp end parallel do    
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine CalcVarLength

