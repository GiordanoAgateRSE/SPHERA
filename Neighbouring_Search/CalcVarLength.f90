!cfile CalcVarLength.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : CalcVarLength
!
! Last updating : November 14, 2012
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Amicarelli/Agate  30Nov11        BSPH: Shepard's coefficient, kernel function 
!                                            at wall, relative distances particle-elements
! 04  Amicarelli-Agate  13nov12        (AA501b) Body dynamics     
!AA504
! 05  Amicarelli        08Apr14        (v5.04) Modifications for granular flows: corrections of parallelization errors; mixture - fixed bed interactions; interface flags.
!AA601
! 06  Amicarelli        26Jan15        Inlet/outlet treatment for DBSPH and minor modifications to DBSPH lines.
!
!************************************************************************************
! Module purpose : Module to calculate array to store length variables for the loop
!                  Computation and storage of the interacting particle index,
!AA406 sub
!                  the relative positions, the neighbouring array of the kernel parameter vectors, 
!                  the neighbouring list vectors,Shepard's coefficient,
!                  the position of the fluid-sediment interfaces along each columns. 
!AA501b
!                  Computation of the neighbouring arrays for body particle - fluid particle 
!                  and body particle - body particle contributions (body dynamics)
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
!
subroutine CalcVarLength
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
!.. Local Scalars ..
integer(4)       :: nceli,igridi,kgridi,jgridi,irang,krang,jrang,ncelj,jgrid1,jgrid2,contliq
!
!AA504 sub 
integer(4)       :: mm,npi,npj,npartint,index_rij_su_h,irestocell,celleloop,fw,i_grid,j_grid
double precision :: rij_su_h,ke_coef,kacl_coef,rij_su_h_quad
double precision :: rijtemp,rijtemp2,ragtemp(3)
double precision :: gradmod,gradmodwacl,wu,denom
!AA504 
double precision :: normal_int_abs,abs_vel
character(len=lencard)  :: nomsub = "CalcVarLength"
!
!AA406 start
double precision :: min_sigma_Gamma
!
!AA501b sub
integer(4) :: bp,bp_f,aux2,i !,contj
!
integer(4), dimension(:), allocatable :: bounded
double precision, dimension(:), allocatable :: dShep_old
!AA406 end
!AA601
double precision :: dis_fp_dbsph_inoutlet,dbsph_inoutlet_threshold
!
!.. Local Arrays ..
!integer(4)                          :: AppoPartintorno,npartintj,ier
!integer(4),dimension(:),allocatable :: PartIntornoAppo
!double precision,dimension(3)       :: Apporag
!double precision,dimension(4)       :: AppoPartKernel
!
!.. External routines ..
integer(4),external :: ParticleCellNumber, CellIndices, CellNumber
!
!.. Executable Statements ..
!
  celleloop = 4 + (ncord-2)*10 + (3-ncord)
!
!.. azzeramento quantita' generali e preparazione costanti
!
  ke_coef = Domain%coefke / Domain%h
  kacl_coef = Domain%coefkacl / Domain%h
!
!----------------------------------------------------------------------------------------
! Arrays description
!...................
! rag            (1:3,npartint)      componenti x,y,z della distanza tra coppie di particelle
! nPartIntorno   (1:nag)             numero particelle intorno alla corrente;
! PartIntorno    (1:npartint)        indice particella che interagisce con quella corrente  
! PartKernel     (1,1:npartint)      Componente Gradiente ex function grad2w kernel
! PartKernel     (2,1:npartint)      grad2w / (r2+eta2)
! PartKernel     (3,1:npartint)      Componente Gradiente ex function grad2wacl kernel
! PartKernel     (4,1:npartint)      Componente Gradiente ex function w kernel
!
!----------------------------------------------------------------------------------------
!
!.. loop sul numero totale delle particelle presenti
!
  nPartIntorno = 0
!AA504 rm line
  ind_interfaces = 0

!AA406 start
!
!AA501 sub start
!AA601 sub
  if ((Domain%tipo == "bsph").and.(nag>0)) then
     nPartIntorno_fw = 0
     allocate (bounded(nag))
     bounded = 0
     allocate (dShep_old(nag))
   endif
!AA501 sub end
!
!AA406 end
!

!AA501b start
  if (n_bodies > 0) then 
     nPartIntorno_bp_f = 0
     nPartIntorno_bp_bp = 0
     aux2 = 0
  endif
!AA501b end

!AA601 sub
!AA504 sub omp directives
!$omp parallel do default(none) &
!$omp private(npi,nceli,irestocell,igridi,jgridi,kgridi,jgrid1,jgrid2,irang,jrang,krang,mm) &
!$omp private(npj,npartint,ncelj,ragtemp,rijtemp,rijtemp2,rij_su_h,rij_su_h_quad,denom) &
!$omp private(index_rij_su_h,gradmod,gradmodwacl,wu,contliq,fw,normal_int_abs,abs_vel,dis_fp_dbsph_inoutlet,dbsph_inoutlet_threshold) &
!$omp shared(nag,pg,Domain,Med,Icont,Npartord,NMAXPARTJ,rag,nPartIntorno,Partintorno,PartKernel) &
!$omp shared(ke_coef,kacl_coef,Doubleh,DoubleSquareh,squareh,nomsub,ncord,eta,eta2,nout,nscr,erosione) &
!$omp shared(ind_interfaces,DBSPH,pg_w,Icont_w,Npartord_w,rag_fw,nPartIntorno_fw,Partintorno_fw,kernel_fw,dShep_old,Granular_flows_options)
! 
  loop_nag: do npi = 1,nag
!
!AA406 start
!
!AA501
    if (Domain%tipo =="bsph") then
!
       pg(npi)%rhoSPH_old = pg(npi)%rhoSPH_new
       pg(npi)%rhoSPH_new = zero
       dShep_old(npi) = pg(npi)%dShep
       pg(npi)%dShep = zero 
       pg(npi)%sigma = zero 
       pg(npi)%FS = 0
!AA601 start       
       pg(npi)%DBSPH_inlet_ID = 0
       pg(npi)%DBSPH_outlet_ID = 0
!AA601 end

!AA601 start
      do npj=DBSPH%n_w+1,DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet  
         dis_fp_dbsph_inoutlet = dsqrt(dot_product(pg(npi)%coord-pg_w(npj)%coord,pg(npi)%coord-pg_w(npj)%coord)) 
         if (npj<=(DBSPH%n_w+DBSPH%n_inlet)) then
            dbsph_inoutlet_threshold = DBSPH%inlet_sections(npj-DBSPH%n_w,10)/dsqrt(2.d0) 
            else
               dbsph_inoutlet_threshold = DBSPH%outlet_sections(npj-DBSPH%n_w-DBSPH%n_inlet,7)/dsqrt(2.d0) 
         endif
         if (dis_fp_dbsph_inoutlet<=dbsph_inoutlet_threshold) then
            if (npj<=(DBSPH%n_w+DBSPH%n_inlet)) then
               pg(npi)%DBSPH_inlet_ID = npj
               else
                  pg(npi)%DBSPH_outlet_ID = npj    
            endif
         endif      
      enddo
!AA601 end

!AA501
    endif
!AA406 end

!AA504
    if (erosione) pg(npi)%normal_int(:) = 0.d0

    nceli = pg(npi)%cella
    if (nceli == 0) cycle
    irestocell = CellIndices (nceli,igridi,jgridi,kgridi)
!
    contliq = 0
    pg(npi)%indneighliqsol = 0
!AA504 start
    pg(npi)%ind_neigh_mix_bed = 0
    pg(npi)%ind_neigh_mob_for_granmob = 0
    pg(npi)%blt_flag = 0
!AA504 end    
    pg(npi)%rijtempmin = 99999.
!
!* indici inizio-fine scansione Y per casi 2d o 3d
    jgrid1 = jgridi - (ncord-2)
    jgrid2 = jgridi + (ncord-2)
!
    loop_jrang: do jrang = jgrid1,jgrid2    ! ---- a  loop sulle 9 celle 
      loop_irang: do irang = igridi-1,igridi+1    ! ---- b  loop sulle 9 celle  
        loop_krang: do krang = kgridi-1,kgridi+1    ! ---- c  loop sulle 9 celle  
!
          ncelj = CellNumber (irang,jrang,krang)
          if (ncelj == 0) cycle    ! cella fuori campo
!
!AA406rm
!          if (Icont(ncelj+1) <= Icont(ncelj)) cycle
!
          loop_mm: do mm = Icont(ncelj),Icont(ncelj+1)-1  ! +++ loop sulle particelle di una cella

             if (nPartIntorno(npi) >= NMAXPARTJ) then
                 write (nout,'(1x,a,i12,a,i12,a)')   ' The computational particle ',npi,' has reached ',NMAXPARTJ,' neighbouring particles.'
                 write (nout,'(1x,a,3f15.10)') '        Coordinate: ',pg(npi)%coord(:)
                 write (nscr,'(1x,a,i12,a,i12,a)')   ' The computational particle ',npi,' has reached ',NMAXPARTJ,' neighbouring particles.'
                 write (nscr,'(1x,a,3f15.10)') '        Coordinate: ',pg(npi)%coord(:)
                 cycle
             end if
            
!AA406
            if (Icont(ncelj+1) <= Icont(ncelj)) cycle
            npj = NPartOrd(mm)
!
!.. calcola le conponenti deltaX, deltaY e deltaZ della distanza tra le due particelle (corrente ed adiacente)
            ragtemp(1:3) = pg(npi)%coord(1:3) - pg(npj)%coord(1:3)
!
!AA406 rm
!.. se una sola delle componenti supera la sfera di influenza di una particella, scarta la coppia (non interagiscono)
!            if (abs(ragtemp(1)) > doubleh .or. abs(ragtemp(2)) > doubleh .or. abs(ragtemp(3)) > doubleh ) cycle
!
!.. calcola la distanza effettiva (considera direttamente il quadrato per aumentare l'accuratezza)
!
            rijtemp = ragtemp(1)*ragtemp(1) + ragtemp(2)*ragtemp(2) + ragtemp(3)*ragtemp(3)
!
!.. se il valore della distanza tra le particelle e' maggiore del raggio di influenza, scarta la coppia
            if (rijtemp > doublesquareh) cycle
!
!.. incremento il contatore totale delle particelle che interagiscono 
!.. incremento il contatore delle particelle che interagiscono con la particella corrente npi-esima
            nPartIntorno(npi) = nPartIntorno(npi) + 1
            npartint = (npi-1) * NMAXPARTJ + nPartIntorno(npi)
!
!AA406 sub
!.. controllo di non superare le dimensioni ipotizzate per i vettori e nel caso the simalution is killed
!

            if (nPartIntorno(npi) > NMAXPARTJ) then
              write (nout,'(1x,a,i12,a)')   ' ERROR! The particle ',npi,' has too many surround particles.'
              write (nout,'(1x,a,3f15.10)') '        Coordinate: ',pg(npi)%coord(:)
              write (nscr,'(1x,a,i12,a)')   ' ERROR! The particle ',npi,' has too many surround particles.'
              write (nscr,'(1x,a,3f15.10)') '        Coordinate: ',pg(npi)%coord(:)
              call diagnostic (arg1=10,arg2=1,arg3=nomsub)
            end if

!.. memorizzo l'indice della particella interagente con quella corrente, le distanze deltax,deltay e deltaz e la 
!.. distanza totale rij
!
            PartIntorno(npartint) = npj
            rijtemp2 = rijtemp
            rijtemp = Dsqrt(rijtemp)
            rij_su_h = rijtemp / Domain%h
            rij_su_h_quad = rijtemp2 / squareh
            index_rij_su_h = int(rij_su_h)
            denom = one / (rijtemp + eta)
            rag(1:3,npartint) = ragtemp(1:3)
!
!................................................................................................
!................................................................................................
!.. calcolo il gradiente e memorizza i relativi coseni direttori con segno in PartKernel
!................................................................................................
!.. partKernel(1) contiene la parte topologica della ex funzione grad2w (a meno di rag(1:3))
! s=r/h
! if (s<=1.)then
!    ssuh=s/h
!    q = -2.*ssuh + 1.5*s*ssuh
! else if(s>1. .and. s<2.) then
!    dms=2.-s
!    q = -dms*dms/(2*h)
! else
!    q = 0.
! end if
!! coef=0.682093/(h*h)
! grad2w = q *coef
!
!................................................................................................
!.. partKernel(2) contiene Partkernel(1)/(R2+eta2) per visc e diff 
!................................................................................................
!.. partKernel(3) contiene la parte topologica della ex funzione grad2wacl (a meno di rag(1:3))
! s = r/h
! if ( s <= 2. )then
!    dms = 2.0 - s
!    q   =-3.0 * dms * dms / h
!   !coef= 0.099472 / (h*h)    ! coeff 2D
!    grad2wacl = q * coef
! else
!    grad2wacl=0.
! end if
!................................................................................................
!.. partKernel(4) contiene la parte topologica della ex funzione w (a meno di rag(1:3))
! s = r/h
! if ( s<=1. )then
!    q = 0.666666667- s*s + (s*s*s)/2.
! else if(s>1. .and. s<2.) then
!    dms = 2.-s
!    q = dms*dms*dms/ 6.
! else
!    q = 0.
! end if
!! coef = 0.682093/(h*h)
! w = q *coef
!................................................................................................
            gradmod = zero
            gradmodwacl = zero
            wu = zero
            PartKernel(1:4,npartint) = zero
!
            if (index_rij_su_h >= 2) cycle
!
            gradmod = -two * rij_su_h + 1.5d0 * rij_su_h_quad
            gradmodwacl = -12.0d0 - 3.0d0 * rij_su_h_quad + 12.0d0 * rij_su_h 
            wu = 0.666666667d0 + rij_su_h_quad * (rij_su_h * half - one)
!
            if (index_rij_su_h > 0) then
              gradmod = -gradmod + rij_su_h_quad - two
              wu = (two - rij_su_h) * (two - rij_su_h) * (two - rij_su_h) * 0.166666666667d0
            end if 
!
            gradmod = gradmod * ke_coef
            gradmodwacl = gradmodwacl * kacl_coef
!
            PartKernel(1,npartint) = gradmod * denom 
            PartKernel(2,npartint) = PartKernel(1,npartint) / (rijtemp2 + eta2)
            PartKernel(3,npartint) = gradmodwacl * denom
            PartKernel(4,npartint) = wu * Domain%coefke
!AA406test 
!WendlandC4
!            PartKernel(1,npartint) = (3./(4.*PIGRECO*(2.*Domain%h**3))) * ((1-rij_su_h/2.)**5) &
!                                    * (-280. * (rij_su_h/2.)**2 - 56. * (rij_su_h/2.))
!            if (rij_su_h /= 0.) PartKernel(1,npartint) = PartKernel(1,npartint) * denom
!            PartKernel(4,npartint) = (3./(4.*PIGRECO*(Domain%h**2))) * &
!                                     ((1-rij_su_h/2.)**6) * (35. * ((rij_su_h/2.)**2) + 18.*(rij_su_h/2.) + 3.)
! Gallati
!               kernel_fw(2,npartint) = (5./(16.*PIGRECO*Domain%h**2)) * ((2-rij_su_h)**3) 
!
!AA406 start
             if (Domain%tipo == "bsph") then
                pg(npi)%sigma = pg(npi)%sigma + pg(npj)%mass * PartKernel(4,npartint) / pg(npj)%dens 
                pg(npi)%rhoSPH_new = pg(npi)%rhoSPH_new + pg(npj)%mass * PartKernel(4,npartint)
             endif
!AA406 end
!
!.. calcolo interfaccia liquida e solida per la particella e per la colonna
!

!AA504 sub start
! Searching for the nearest fluid/granular SPH particle 
            if (erosione) then
              if (Med(pg(npi)%imed)%tipo /= Med(pg(npj)%imed)%tipo) then
                 if ( (rijtemp<pg(npi)%rijtempmin(1)) .or. ( (rijtemp==pg(npi)%rijtempmin(1)).and.(npj>pg(npi)%indneighliqsol) ) ) then   
                    if ( (index(Med(pg(npi)%imed)%tipo,"liquid")>0) .and. (pg(npi)%coord(3)>pg(npj)%coord(3)) ) then
                       pg(npi)%indneighliqsol = npj
                       pg(npi)%rijtempmin(1) = rijtemp
                       else if ((index(Med(pg(npi)%imed)%tipo,"granular") > 0) .and. (pg(npi)%coord(3)<pg(npj)%coord(3))) then
                          pg(npi)%indneighliqsol = npj
                          pg(npi)%rijtempmin(1) = rijtemp
                    end if
                 end if
              end if
!AA504 sub end
!AA504 start
! Searching for the nearest mobile/fixed particle (if cycle apparently less efficient than the previous one because we have to compute the interface normal)              
              if (pg(npi)%state /= pg(npj)%state) then
!test                  
!                 if ((pg(npi)%state == "flu").and.(pg(npi)%coord(3)>pg(npj)%coord(3))) then
                    if (pg(npi)%state == "flu") then
                    if (pg(npi)%imed==Granular_flows_options%ID_granular) pg(npi)%normal_int(:) = pg(npi)%normal_int(:) &
                                                                              - (pg(npj)%coord(:)-pg(npi)%coord(:)) * PartKernel(4,npartint) 
                    if ( (rijtemp<pg(npi)%rijtempmin(2)) .or. ( (rijtemp==pg(npi)%rijtempmin(2)).and.(npj>pg(npi)%ind_neigh_mix_bed) ) ) then  
                       pg(npi)%ind_neigh_mix_bed = npj
                       pg(npi)%rijtempmin(2) = rijtemp
                    endif   
!test                    
!                    else if ((pg(npi)%state == "sol").and.(pg(npi)%coord(3)<pg(npj)%coord(3))) then
                    else if (pg(npi)%state == "sol") then
                       pg(npi)%normal_int(:) = pg(npi)%normal_int(:) + (pg(npj)%coord(:)-pg(npi)%coord(:)) * PartKernel(4,npartint)  
                    if ( (rijtemp<pg(npi)%rijtempmin(2)) .or. ( (rijtemp==pg(npi)%rijtempmin(2)).and.(npj>pg(npi)%ind_neigh_mix_bed) ) ) then 
                       pg(npi)%ind_neigh_mix_bed = npj
                       pg(npi)%rijtempmin(2) = rijtemp
                    end if
                 endif
              endif
! Searching for the nearest mobile granular/mobile particle
              if ((index(Med(pg(npi)%imed)%tipo,"granular")>0).and.(pg(npi)%state=="flu").and.(pg(npj)%state=="flu")) then
                 if ( (rijtemp<pg(npi)%rijtempmin(3)) .or. ( (rijtemp==pg(npi)%rijtempmin(3)).and.(npj>pg(npi)%ind_neigh_mob_for_granmob) ) ) then 
                    if (pg(npi)%coord(3)<pg(npj)%coord(3)) then
                        pg(npi)%ind_neigh_mob_for_granmob = npj
                        pg(npi)%rijtempmin(3) = rijtemp
                    endif
                 end if
              end if              
           end if  ! erosione
!AA504 end            

          end do loop_mm      ! +++ loop sulle part di una cella
!
!AA406 start
! Loop over the neighbouring wall particles in the cell
          if ((Domain%tipo == "bsph").and.(DBSPH%n_w > 0)) then
          loop_fw: do fw = Icont_w(ncelj),Icont_w(ncelj+1)-1
            npj = NPartOrd_w(fw)
! Relative positions and distances
            ragtemp(1:3) = pg(npi)%coord(1:3) - pg_w(npj)%coord(1:3)
            rijtemp = ragtemp(1)*ragtemp(1) + ragtemp(2)*ragtemp(2) + ragtemp(3)*ragtemp(3)
! Distance check
            if (rijtemp > doublesquareh) cycle
! counter incerements 
            nPartIntorno_fw(npi) = nPartIntorno_fw(npi) + 1
            npartint = (npi-1) * NMAXPARTJ + nPartIntorno_fw(npi)
! Savings
            PartIntorno_fw(npartint) = npj
            rijtemp2 = rijtemp
            rij_su_h = Dsqrt(rijtemp) / Domain%h
            rij_su_h_quad = rijtemp2 / squareh
            index_rij_su_h = int(rij_su_h)
            rag_fw(1:3,npartint) = ragtemp(1:3)
! Kernel computation
            if (index_rij_su_h >= 2) cycle
            wu = 0.666666667d0 + rij_su_h_quad * (rij_su_h * half - one)
            if (index_rij_su_h > 0) then
              wu = (two - rij_su_h) * (two - rij_su_h) * (two - rij_su_h) * 0.166666666667d0
            end if 
            kernel_fw(1,npartint) = wu * Domain%coefke
!
!AA406 start
            pg(npi)%dShep = pg(npi)%dShep + kernel_fw(1,npartint) * pg_w(npj)%weight  *  &
                            (  pg_w(npj)%normal(1) * (pg_w(npj)%vel(1)-pg(npi)%var(1)) +  &
                               pg_w(npj)%normal(2) * (pg_w(npj)%vel(2)-pg(npi)%var(2)) +  &
                               pg_w(npj)%normal(3) * (pg_w(npj)%vel(3)-pg(npi)%var(3))    ) 
!
!AA406test correction che dava il segmentation fault
            if ((rij_su_h*Domain%h) <= (1.3*Domain%dd + pg_w(npj)%weight/2.)) then
!
               pg_w(npj)%wet = 1
            endif
            denom = one / (Dsqrt(rijtemp) + eta)
            gradmod = -two * rij_su_h + 1.5d0 * rij_su_h_quad
            if (index_rij_su_h > 0) then
              gradmod = -gradmod + rij_su_h_quad - two
            end if 
            kernel_fw(2,npartint) = gradmod * ke_coef * denom
!
!AA501test
!            if (pg_w(npj)%wet == 1) then
!
            pg(npi)%sigma = pg(npi)%sigma +  pg_w(npj)%mass * kernel_fw(1,npartint) / pg_w(npj)%dens 
            pg(npi)%rhoSPH_new = pg(npi)%rhoSPH_new + pg_w(npj)%mass * kernel_fw(1,npartint)
!
!AA501test
!            endif
!
!
!AA406test Gallati
!            kernel_fw(3,npartint) = (5./(16.*PIGRECO*Domain%h**2)) * ((2-rij_su_h)**3) 
!            kernel_fw(4,npartint) = (-12.0d0 - 3.0d0 * rij_su_h_quad + 12.0d0 * rij_su_h) * kacl_coef * denom
!
!AA406test WendlandC4
!             kernel_fw(1,npartint) = (3./(4.*PIGRECO*(Domain%h**2))) * &
!                                      ((1-rij_su_h/2.)**6) * (35. * ((rij_su_h/2.)**2) + 18.*(rij_su_h/2.) + 3.)
!             kernel_fw(2,npartint) = (3./(4.*PIGRECO*(2.*Domain%h**3))) * ((1-rij_su_h/2.)**5) &
!                                    * (-280. * (rij_su_h/2.)**2 - 56. * (rij_su_h/2.))
!             if (rij_su_h /= 0.) kernel_fw(2,npartint) = kernel_fw(2,npartint) * denom
!
!AA406test Gallati
!            kernel_fw(2,npartint) = (5./(16.*PIGRECO*Domain%h**2)) * ((2-rij_su_h)**3) 
!AA406test
!            kernel_fw(2,npartint) = pg_w(npj)%dens 
!            kernel_fw(3,npartint) = pg_w(npj)%pres
!
         end do loop_fw
         endif
!AA406 end
!

        end do loop_krang   ! ---- c  loop sulle 9 celle    
      end do loop_irang   ! ---- b  loop sulle 9 celle    
    end do loop_jrang   ! ---- a  loop sulle 9 celle
!

    if (erosione) then
!.. trovo il pelo libero sulla colonna
      if (index(Med(pg(npi)%imed)%tipo,"liquid") > 0) then
!AA504 sub start 
!$omp critical (free_surface_detection)
        if (ind_interfaces(igridi,jgridi,1)==0) then
           ind_interfaces(igridi,jgridi,1) = npi   
           else
           if (pg(npi)%coord(3)>pg(ind_interfaces(igridi,jgridi,1))%coord(3)) then
              ind_interfaces(igridi,jgridi,1) = npi 
              elseif (pg(npi)%coord(3)==pg(ind_interfaces(igridi,jgridi,1))%coord(3)) then
              if (npi>ind_interfaces(igridi,jgridi,1)) ind_interfaces(igridi,jgridi,1) = npi
           endif
        end if
!$omp end critical (free_surface_detection)        
!AA504 sub end        
!
!!! 25gen2011
!!!      else if (index(Med(pg(npi)%imed)%tipo,"granular")) then
!!!! ............................... 10.41 .........................
!!!        pesoj = pesoj * contsol / nPartIntorno(npi) 
!!!        if (pesoj > zero .and. pesoj <= 0.10) pg(npi)%punta = .true.    !prova 10.42 - prova 10.41 con pesoj<=0.25
      
!AA504 start
         else
! In case of no free surface in the column and no erosion criterion, updating the miwture - fixed bed interface       
            if (Granular_flows_options%erosion_flag<2) then
!$omp critical (fixed_bed_detection)
               abs_vel = dsqrt(dot_product(pg(npi)%vel,pg(npi)%vel))
               if (abs_vel<=Granular_flows_options%velocity_fixed_bed) then
                  if (ind_interfaces(igridi,jgridi,4)==0) then
                     ind_interfaces(igridi,jgridi,4) = npi   
                     else
                        if (pg(npi)%coord(3)>pg(ind_interfaces(igridi,jgridi,4))%coord(3)) then
                           ind_interfaces(igridi,jgridi,4) = npi 
                           elseif (pg(npi)%coord(3)==pg(ind_interfaces(igridi,jgridi,4))%coord(3)) then
                              if (npi>ind_interfaces(igridi,jgridi,4)) ind_interfaces(igridi,jgridi,4) = npi
                        endif
                  end if
               endif
!$omp end critical (fixed_bed_detection) 
           endif
!AA504 end
      
      end if                                     
!
      
    end if
!
    
!AA504 start 
! Normalization of the interface normal between the granular mixture and the fixed bed (granular flows)
  if ((erosione).and.(pg(npi)%imed==Granular_flows_options%ID_granular)) then 
      normal_int_abs = dsqrt(dot_product(pg(npi)%normal_int,pg(npi)%normal_int)) 
      if (normal_int_abs>zero) then
         pg(npi)%normal_int(:) = pg(npi)%normal_int(:)/normal_int_abs
      endif    
  endif
!AA504 end  
  
!AA504 start
  if (erosione) then
!$omp critical (interface_definition)  
!Update the local position of the upper interface of the bed load transport region
  if ((index(Med(pg(npi)%imed)%tipo,"liquid")>0).and.(pg(npi)%indneighliqsol.ne.0)) then 
     if (ind_interfaces(igridi,jgridi,3)==0) then
        ind_interfaces(igridi,jgridi,3) = pg(npi)%indneighliqsol   
        elseif (pg(pg(npi)%indneighliqsol)%coord(3)>pg(ind_interfaces(igridi,jgridi,3))%coord(3)) then
           ind_interfaces(igridi,jgridi,3) = pg(npi)%indneighliqsol 
           elseif (pg(pg(npi)%indneighliqsol)%coord(3)==pg(ind_interfaces(igridi,jgridi,3))%coord(3)) then
              if (pg(npi)%indneighliqsol>ind_interfaces(igridi,jgridi,3)) ind_interfaces(igridi,jgridi,3) = pg(npi)%indneighliqsol
     endif
  endif        
!Update the local position of the upper interface of the bed load transport region (liquid side), only for Shields 2D and Mohr criteria   
  if ((Granular_flows_options%ID_erosion_criterion.ne.1).and.(index(Med(pg(npi)%imed)%tipo,"granular")>0).and.(pg(npi)%indneighliqsol.ne.0)) then   
     if (ind_interfaces(igridi,jgridi,2)==0) then
        ind_interfaces(igridi,jgridi,2) = pg(npi)%indneighliqsol   
        elseif (pg(pg(npi)%indneighliqsol)%coord(3)>pg(ind_interfaces(igridi,jgridi,2))%coord(3)) then
           ind_interfaces(igridi,jgridi,2) = pg(npi)%indneighliqsol 
           elseif (pg(pg(npi)%indneighliqsol)%coord(3)==pg(ind_interfaces(igridi,jgridi,2))%coord(3)) then
              if (pg(npi)%indneighliqsol>ind_interfaces(igridi,jgridi,2)) ind_interfaces(igridi,jgridi,2) = pg(npi)%indneighliqsol
     endif
  endif
!Update the local position of the fixed bed    
  if ((pg(npi)%state=="flu").and.(pg(npi)%ind_neigh_mix_bed.ne.0)) then  
     if (ind_interfaces(igridi,jgridi,4)==0) then
        ind_interfaces(igridi,jgridi,4) = pg(npi)%ind_neigh_mix_bed    
        elseif (pg(pg(npi)%ind_neigh_mix_bed)%coord(3)>pg(ind_interfaces(igridi,jgridi,4))%coord(3)) then
           ind_interfaces(igridi,jgridi,4) = pg(npi)%ind_neigh_mix_bed
           elseif (pg(pg(npi)%ind_neigh_mix_bed)%coord(3)==pg(ind_interfaces(igridi,jgridi,4))%coord(3)) then
              if (pg(npi)%ind_neigh_mix_bed>ind_interfaces(igridi,jgridi,4)) ind_interfaces(igridi,jgridi,4) = pg(npi)%ind_neigh_mix_bed
     endif
  endif 
!$omp end critical (interface_definition)   
  endif
!AA504 end

  end do loop_nag
!
!$omp end parallel do

!AA504 start
! Compute the interface flags (granular flows)
 if (erosione) then
 do npi=1,nag
    nceli = ParticleCellNumber(pg(npi)%coord)
    irestocell = CellIndices(nceli,igridi,jgridi,kgridi)
    if (pg(npi)%imed==Granular_flows_options%ID_granular) then
       if (ind_interfaces(igridi,jgridi,3)==0) then
          ind_interfaces(igridi,jgridi,3) = -npi   
          else
             if (ind_interfaces(igridi,jgridi,3)<0) then
                if (pg(npi)%coord(3)>pg(-ind_interfaces(igridi,jgridi,3))%coord(3)) then
                   ind_interfaces(igridi,jgridi,3) = -npi 
                      elseif (pg(npi)%coord(3)==pg(-ind_interfaces(igridi,jgridi,3))%coord(3)) then
                         if (npi>(-ind_interfaces(igridi,jgridi,3))) ind_interfaces(igridi,jgridi,3) = -npi
                endif
             end if
       endif
    endif
 end do
!$omp parallel do default(none) shared(Grid,pg,ind_interfaces,nout) private(i_grid,j_grid)
 do i_grid=1,Grid%ncd(1)
    do j_grid=1,Grid%ncd(2)
       if (ind_interfaces(i_grid,j_grid,3)<0) ind_interfaces(i_grid,j_grid,3) = -ind_interfaces(i_grid,j_grid,3)
       if ( (ind_interfaces(i_grid,j_grid,3).ne.0) .and. (ind_interfaces(i_grid,j_grid,1).ne.0) ) then
          if (pg(ind_interfaces(i_grid,j_grid,3))%coord(3)>pg(ind_interfaces(i_grid,j_grid,1))%coord(3)) then
             ind_interfaces(i_grid,j_grid,3) = ind_interfaces(i_grid,j_grid,1)
          endif
       endif
       if (ind_interfaces(i_grid,j_grid,4).ne.0) pg(ind_interfaces(i_grid,j_grid,4))%blt_flag = 3
       if (ind_interfaces(i_grid,j_grid,3).ne.0) pg(ind_interfaces(i_grid,j_grid,3))%blt_flag = 2
       if (ind_interfaces(i_grid,j_grid,1).ne.0) pg(ind_interfaces(i_grid,j_grid,1))%blt_flag = 1
    enddo
 enddo
!$omp end parallel do
 endif
!AA504 end 

!!! allocate (PartIntornoAppo(1:NMAXPARTJ*PARTICLEBUFFER), stat = ier)
!!! if (ier /= 0) then
!!!   write (nout,'(1x,a,i2)') "   Arrays PartIntornoAppo not allocated. Error code: ",ier
!!!   call diagnostic (arg1=4,arg3=nomsub)
!!! end if
!!! PartIntornoAppo = Partintorno
!
!!!!!ordinamento
!!!!  do npi = 1, nag
!!!!    do i = 1, nPartIntorno(npi) - 1
!!!!      npartint = (npi-1) * NMAXPARTJ + i
!!!!      irang = npartint
!!!!      do j = (i + 1) , nPartIntorno(npi)
!!!!        npartintj = (npi-1) * NMAXPARTJ + j
!!!!        if (Partintorno(npartintj) < Partintorno(irang)) then
!!!!          irang = npartintj
!!!!        end if
!!!!      end do
!!!!      AppoPartintorno = Partintorno(npartint)
!!!!      Partintorno(npartint) = Partintorno(irang)
!!!!      Partintorno(irang) = AppoPartintorno
!!!!      Apporag(1:3) = rag(1:3,npartint)
!!!!      rag(1:3,npartint) = rag(1:3,irang)
!!!!      rag(1:3,irang) = Apporag(1:3)
!!!!      AppoPartKernel(1:4) = PartKernel(1:4,npartint)
!!!!      PartKernel(1:4,npartint) = PartKernel(1:4,irang)
!!!!      PartKernel(1:4,irang) = AppoPartKernel(1:4)
!!!!    end do
!!!!  end do
!
! inverto ordine
!!!!  do npi = 1, nag
!!!!    jrang = nPartIntorno(npi)
!!!!    irang = 1
!!!!    do while (irang < jrang)
!!!!      npartint = (npi-1) * NMAXPARTJ + irang
!!!!      npartintj = (npi-1) * NMAXPARTJ + jrang
!!!!      AppoPartintorno = Partintorno(npartint)
!!!!      Partintorno(npartint) =  Partintorno(npartintj)
!!!!      Partintorno(npartintj) =  AppoPartintorno
!!!!      Apporag(1:3) = rag(1:3,npartint)
!!!!      rag(1:3,npartint) =  rag(1:3,npartintj)
!!!!      rag(1:3,npartintj) =  Apporag(1:3)
!!!!      AppoPartKernel(1:4) = PartKernel(1:4,npartint)
!!!!      PartKernel(1:4,npartint) =  PartKernel(1:4,npartintj)
!!!!      PartKernel(1:4,npartintj) =  AppoPartKernel(1:4)
!!!!      irang = irang + 1
!!!!      jrang = jrang - 1
!!!!    end do
!!!!  end do
!
!!! do npi = 1,nag
!!!   write (999,*) npi,tempo,'----------'
!!!   write (998,*) npi,tempo,'----------'
!!!   do irang = 1,NMAXPARTJ
!!!     npartint = (npi-1) * NMAXPARTJ + irang
!!!     write (999,*) irang,Partintorno(npartint)
!!!     write (998,*) irang,PartintornoAppo(npartint)
!!!   end do
!!! end do
!!! deallocate (PartIntornoAppo)
!!! stop
!
!
!!!if (tempo > 0.5) then
!!!  open (unit=98, file="File_Icont_v330.txt", form="formatted", status="unknown")
!!!  write(98,*) '  '
!!!  write(98,*) ' tempo = ',tempo,' -------------------- '
!!!  write(98,'(a,i8,a,i8)') ' Icont(grid%nmax+1) = ',Icont(grid%nmax+1)
!!!  do mm = 1,grid%nmax
!!!    write(98,'(a,i8,a,i8,a,i8)') ' mm = ',mm,'  Icont(mm) = ',Icont(mm),'  Icont(mm+1) = ',Icont(mm+1)
!!!    do mm1 = Icont(mm),Icont(mm+1)-1
!!!      write(98,'(a,i8,a,i8)') ' mm1 = ',mm1,'  NPartOrd(mm1) = ',NPartOrd(mm1)
!!!    end do
!!!  end do
!!!  close(98)
!!!!
!!!  open (unit=99, file="File_Intorno_v330.txt", form="formatted", status="unknown")
!!!!  npi = 2153  !2D
!!!  npi = 10023   !3D
!!!  nceli = pg(npi)%cella
!!!  irestocell = CellIndices (nceli,igridi,jgridi,kgridi)
!!!  write(99,*) '  '
!!!  write(99,*) ' tempo = ',tempo,' -------------------- '
!!!  write(99,'(a,i8,a,i8,a,a3)') '  particella         = ',npi,'  cellai   = ',nceli,'  status = ',pg(npi)%state
!!!  write(99,'(a,3f16.6)')  '  coordinate         = ',pg(npi)%coord
!!!  write(99,'(a,3f16.6)')  '  coordinate cellai  = ',(igridi-1)*grid%dcd(1)+grid%extr(1,1),(jgridi-1)*grid%dcd(2)+grid%extr(2,1),(kgridi-1)*grid%dcd(3)+grid%extr(3,1)
!!!  write(99,'(a,3f16.6)')  '  coordinate cellai  = ',igridi*grid%dcd(1)+grid%extr(1,1),jgridi*grid%dcd(2)+grid%extr(2,1),kgridi*grid%dcd(3)+grid%extr(3,1)
!!!  write(99,'(a,3f16.6)')  '  coordinate cellai  = ',(igridi+1)*grid%dcd(1)+grid%extr(1,1),(jgridi+1)*grid%dcd(2)+grid%extr(2,1),(kgridi+1)*grid%dcd(3)+grid%extr(3,1)
!!!  write(99,'(a,i8)')      '  nPartIntorno(i)    = ',nPartIntorno(i)
!!!  niniz = (i-1) * NMAXPARTJ + 1
!!!  nfine = (i-1) * NMAXPARTJ + nPartIntorno(i)
!!!  do ijkgrid = niniz,nfine
!!!    npj = PartIntorno(ijkgrid)
!!!    ncelj = pg(npj)%cella
!!!    write(99,*) '  '
!!!    write(99,'(a,i8,a,i8,a,a3)') '  PartIntorno per npi   = ',PartIntorno(ijkgrid),'  cellaj   = ',ncelj,'  status = ',pg(npj)%state
!!!    write(99,'(a,3f16.6)')  '  -> coordinatej      = ',pg(npj)%coord
!!!    write(99,'(a,3f16.8)')  '  rag(1:3,...)        = ',rag(1:3,ijkgrid)
!!!    write(99,'(a,4f16.6)')  '  PartKernel(1:4,...) = ',PartKernel(1:4,ijkgrid)
!!!!    jgrid1 = (j-1) * NMAXPARTJ + 1
!!!!    jgrid2 = (j-1) * NMAXPARTJ + nPartIntorno(npj)
!!!!    write(99,'(a,i8)')      '  -> nPartIntorno(npj)    = ',nPartIntorno(npj)
!!!!    do mm = jgrid1,jgrid2
!!!!      write(99,'(a,i8)')      '   --> PartIntorno per npj    = ',PartIntorno(mm)
!!!!      write(99,'(a,3f16.8)')  '   --> rag(1:3,...) = ',rag(1:3,mm)
!!!!      write(99,'(a,4f16.6)')  '   --> PartKernel(1:4,...) = ',PartKernel(1:4,mm)
!!!!    end do
!!!  end do
!!!  close(99)
!!! stop
!!!end if
!
!AA501
  if (Domain%tipo == "bsph") then
!AA406 start
  do npi=1,nag
!Gamma intialization for not inlet conditions
!AA601 sub
      if (it_corrente == -2) then
!
!AA406test
         if (nPartIntorno_fw(npi) == 0) then 
            pg(npi)%Gamma = one
         else
               pg(npi)%Gamma = pg(npi)%sigma
               pg(npi)%Gamma = min (pg(npi)%Gamma,one)    
         endif
!
!AA406test
!AA601 rm
      endif
!AA406!!!test
!      do contj=1,nPartIntorno_fw(npi)
!         npartint = (npi-1)* NMAXPARTJ + contj
!         npj = PartIntorno_fw(npartint)
!         if (pg_w(npj)%wet == 1) bounded(npi) = 1 
!      end do 
!
!AA601
      if (it_corrente>-2) then
!AA406test
!      min_sigma_Gamma = min((pg(npi)%sigma+0.005),pg(npi)%Gamma)
        min_sigma_Gamma = min((pg(npi)%sigma+0.05),pg(npi)%Gamma)
!
!AA406!!!test
      if (min_sigma_Gamma /= pg(npi)%Gamma) then
         pg(npi)%Gamma_last_active = zero
         pg(npi)%FS = 1
!AA406test
!         if ((nPartIntorno_fw(npi) > 0).and.(bounded(npi)==0)) then
!            pg(npi)%FS = 2
!            else
               pg(npi)%uni = pg(npi)%sigma
!         endif
         else
!AA406test
!            if ((nPartIntorno_fw(npi) > 0).and.(bounded(npi)==0)) then
!               pg(npi)%FS = 3
!               if (pg(npi)%Gamma_last_active == zero) then
!                  pg(npi)%Gamma_last_active = pg(npi)%Gamma - dShep_old(npi) * dt
!                  pg(npi)%uni = pg(npi)%Gamma_last_active 
!               endif   
!               else
!AA406test
!      if (it_corrente == it_start) then
!
                  pg(npi)%uni = pg(npi)%Gamma
!AA406test
!      endif
!
                  pg(npi)%Gamma_last_active = zero
!            endif     
      endif
!AA406test
!                  pg(npi)%uni = pg(npi)%Gamma
!
!AA406test
!      if (it_corrente == it_start) then
!         if (pg(npi)%FS == 1) then
!            pg(npi)%dens = pg(npi)%rhoSPH_new / pg(npi)%sigma
!            else
!               pg(npi)%dens = pg(npi)%rhoSPH_new / pg(npi)%Gamma
!         endif
!      endif
!AA406test
!       if (pg(npi)%rhoSPH_old == zero) pg(npi)%DensShep = pg(npi)%dens * pg(npi)%Gamma
        if (pg(npi)%rhoSPH_old == zero) then
           pg(npi)%DensShep = pg(npi)%rhoSPH_new * pg(npi)%Gamma
           if (pg(npi)%FS ==1) then
               pg(npi)%dens = pg(npi)%rhoSPH_new / pg(npi)%sigma
              else
!
!AA406test
!Benchmark2
!             if (it_corrente>1000) then
                pg(npi)%dens = pg(npi)%rhoSPH_new / pg(npi)%Gamma
!                else 
!                pg(npi)%dens = pg(npi)%rhoSPH_new / pg(npi)%sigma
!             endif
!
           endif
        endif
!AA601
      endif
  end do
!AA601 sub start
  if (allocated(bounded)) deallocate (bounded)
  if (allocated(dShep_old)) deallocate (dShep_old)
!AA601 sub end  
!
!AA501
  endif 
!
!
!AA406test start density from SPH approximation
!  do npi=1,nag 
!      if (pg(npi)%uni >0.1) then
!         pg(npi)%dens =  dens_app(npi) / pg(npi)%uni
!         else
!            pg(npi)%dens = 1000.
!      endif
!  end do
!  deallocate (dens_app)
!AA406 end

!AA501b start
! Parameters for the body dynamics
 if (n_bodies > 0) then
! Loop over the body particles
!AA501btest
!$omp parallel do default(none) &
!$omp private(npi,nceli,irestocell,igridi,jgridi,kgridi,jgrid1,jgrid2,irang,jrang,krang,bp_f) &
!$omp private(npj,npartint,ncelj,ragtemp,rijtemp,rijtemp2,rij_su_h,rij_su_h_quad,denom) &
!$omp private(index_rij_su_h,gradmod,bp,gradmodwacl,i,aux2) &
!$omp shared(NMAXPARTJ,ke_coef,kacl_coef,DoubleSquareh,squareh,ncord,Domain) &
!$omp shared (Icont_bp,NPartOrd_bp,bp_arr,nPartIntorno_bp_bp,PartIntorno_bp_bp,rag_bp_bp,n_body_part,Icont) &
!$omp shared(n_bodies,nPartIntorno_bp_f,PartIntorno_bp_f,rag_bp_f,KerDer_bp_f_cub_spl,KerDer_bp_f_Gal,NPartOrd,eta,pg)
    do npi=1,n_body_part

!AA501btest
!Computation of the ID of the surface body particles
       i = 0
       aux2 = 0
       do while (i<npi) 
          i = i+1 
          if (bp_arr(i)%area > 0.) aux2 = aux2+1
       enddo
!        aux2 = aux2+1 

       nceli = bp_arr(npi)%cell
       if (nceli == 0) cycle
       irestocell = CellIndices (nceli,igridi,jgridi,kgridi)
! indici inizio-fine scansione Y per casi 2d o 3d
       jgrid1 = jgridi - (ncord-2)
       jgrid2 = jgridi + (ncord-2)
! Loop over the neighbouring cells
       do jrang = jgrid1,jgrid2        ! ---- a  loop sulle 9 celle 
       do irang = igridi-1,igridi+1    ! ---- b  loop sulle 9 celle  
       do krang = kgridi-1,kgridi+1    ! ---- c  loop sulle 9 celle  
          ncelj = CellNumber (irang,jrang,krang)
          if (ncelj == 0) cycle    ! cella fuori campo

! Parameters for body particle - fluid particle  interactions (for body and fluid dynamics)
! Loop over the neighbouring body particles in the cell
          loop_bp_f: do bp_f = Icont(ncelj),Icont(ncelj+1)-1
             npj = NPartOrd(bp_f)
! Relative positions and distances
! sign inversion because body particle acts here as a computational particle
             ragtemp(1:3) = - pg(npj)%coord(1:3) + bp_arr(npi)%pos(1:3)  
             rijtemp = ragtemp(1)*ragtemp(1) + ragtemp(2)*ragtemp(2) + ragtemp(3)*ragtemp(3)
! Distance check
             if (rijtemp > doublesquareh) cycle
! neighbouring lists 
             nPartIntorno_bp_f(npi) = nPartIntorno_bp_f(npi) + 1
             npartint = (npi-1) * NMAXPARTJ + nPartIntorno_bp_f(npi)
             PartIntorno_bp_f(npartint) = npj
! Relative distance
             rijtemp2 = rijtemp
             rijtemp = Dsqrt(rijtemp)
             rij_su_h = rijtemp / Domain%h
             rij_su_h_quad = rijtemp2 / squareh
             index_rij_su_h = int(rij_su_h)
             denom = one / (rijtemp + eta)
             rag_bp_f(1:3,npartint) = ragtemp(1:3)  
             if (ncord == 2) rag_bp_f(2,npartint) = 0.
! Kernel gradients (cubic spline)
             gradmod = zero
             KerDer_bp_f_cub_spl(npartint) = zero
             if (index_rij_su_h >= 2) cycle
             gradmod = -two * rij_su_h + 1.5d0 * rij_su_h_quad
             if (index_rij_su_h > 0) then
                gradmod = -gradmod + rij_su_h_quad - two
             end if 
             gradmod = gradmod * ke_coef
             KerDer_bp_f_cub_spl(npartint) = gradmod * denom 
! Kernel gradients (Gallati's kernel derivative)
             gradmodwacl = -12.0d0 - 3.0d0 * rij_su_h_quad + 12.0d0 * rij_su_h 
             gradmodwacl = gradmodwacl * kacl_coef
             KerDer_bp_f_Gal(npartint) = gradmodwacl * denom 
          end do loop_bp_f
! End Loop over the neighbouring body particles in the cell
          
! Loop over the neighbouring body particles in the cell
          loop_bp: do bp = Icont_bp(ncelj),Icont_bp(ncelj+1)-1
             npj = NPartOrd_bp(bp)
! Only neighbours belonging to a surface of another body
             if ( (bp_arr(npi)%area>0.) .and. (bp_arr(npj)%area>0.) .and. (bp_arr(npi)%body /= bp_arr(npj)%body) ) then
! Relative positions and distances
                ragtemp(1:3) = bp_arr(npi)%pos(1:3) - bp_arr(npj)%pos(1:3)  
                rijtemp = ragtemp(1)*ragtemp(1) + ragtemp(2)*ragtemp(2) + ragtemp(3)*ragtemp(3)
! Distance check
                if (rijtemp > doublesquareh) cycle
! neighbouring lists 
                nPartIntorno_bp_bp(aux2) = nPartIntorno_bp_bp(aux2) + 1
                npartint = (aux2-1) * NMAXPARTJ + nPartIntorno_bp_bp(aux2)
                PartIntorno_bp_bp(npartint) = npj
! Relative distance
                rag_bp_bp(1:3,npartint) = ragtemp(1:3)  
                if (ncord == 2) rag_bp_bp(2,npartint) = 0.   
             endif
          end do loop_bp
          
       end do 
       end do 
       end do 
! End Loop over the neighbouring body particles in the cell       

    end do
!$omp end parallel do    
!AA501btest

 endif
!AA501b end

return
end subroutine CalcVarLength
!---split

