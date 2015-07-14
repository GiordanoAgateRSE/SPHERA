!cfile AggDens.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name    : AggDens
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
!
!************************************************************************************
! Module purpose : 
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D, time_integration
!
! Called routines: 
!
!************************************************************************************
!
subroutine AggDens
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
integer(4) :: npi,ii
double precision :: tirhoc, tirhow, vdens      !,tirhoc1,tirhow1,tirhoc2,tirhow2
double precision :: a1, b1, c1     !,gam,qut,qu,appo, 
!double precision :: dop, denom, del, a, b, vol_temp, controllo,controllo2    !,gam,qut,qu,
!
!.. Executable Statements ..
!appo,
!$omp parallel do default(none) &
!$omp private(npi,vdens,tirhoc,tirhow,a1,b1,c1,ii) &  
!$omp shared(nag,Pg,Med,Domain,dt,it_corrente,ncord,indarrayFlu,Array_Flu)
!
!!!  do npi = 1,nag
!!!!
!!!    if (pg(npi)%koddens /= 0) cycle
!!!!
!!!!!    if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!!!!!gio12mar2009
!!!    if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" .or. pg(npi)%state == "sol") cycle
!$$$$$
  do ii = 1,indarrayFlu
    npi = Array_Flu(ii)
    if (pg(npi)%koddens /= 0) cycle
!$$$$$
!
    vdens  = pg(npi)%dden / pg(npi)%dens
!
!§
    if(it_corrente <= 1) then
      pg(npi)%tiroc = pg(npi)%rhoc * pg(npi)%VolFra
      tirhoc = pg(npi)%tiroc + dt * (pg(npi)%tiroc * vdens + pg(npi)%rhoc * pg(npi)%diffu - Med(2)%settlingcoef)
    else    
      tirhoc = pg(npi)%tiroc + dt * (pg(npi)%tiroc * vdens + pg(npi)%rhoc * pg(npi)%diffu - Med(2)%settlingcoef)
      pg(npi)%tiroc = tirhoc     
    end if
!§
!    tirhoc = pg(npi)%rhoc * pg(npi)%VolFra           ! rhograintilde     
!    tirhoc = tirhoc + dt * (tirhoc * vdens + pg(npi)%rhoc * pg(npi)%diffu - med(2)%settlingcoef)      ! new rhograintilde
!
    if (tirhoc >= (pg(npi)%dens / pg(npi)%VolFra)) then
      tirhoc = pg(npi)%dens
      tirhow = zero
      pg(npi)%VolFra = VFmx
      pg(npi)%rhoc = pg(npi)%dens
      pg(npi)%mass = pg(npi)%dens * (Domain%dd**ncord)
    else if (tirhoc <= zero) then
      tirhoc = zero
      tirhow = pg(npi)%dens
      pg(npi)%rhow = pg(npi)%dens
      pg(npi)%VolFra = VFmn
      pg(npi)%mass = pg(npi)%dens * (Domain%dd**ncord)
    else
!
!§
!      if(pg(npi)%imed==1)then
!        tirhow = (Med(2)%den0 * VFmn + Med(1)%den0 * (1-VFmn)) - tirhoc
!      else if(pg(npi)%imed==2)then
!        tirhow = (Med(2)%den0 * VFmx + Med(1)%den0 * (1-VFmx)) - tirhoc
!      end if
!§
      tirhow = pg(npi)%dens - tirhoc
!
      a1 = med(2)%den0 * med(2)%celerita*med(2)%celerita - med(1)%den0 * med(1)%celerita*med(1)%celerita
      b1 = - (med(2)%celerita*med(2)%celerita) * (tirhoc + med(2)%den0) - (med(1)%celerita*med(1)%celerita) * (tirhow - med(1)%den0)
      c1 = (med(2)%celerita*med(2)%celerita) * tirhoc 
!      appo = b1 * b1 - 4.d0 * a1 * c1
!
!      vol_temp = pg(npi)%VolFra
!if(( b1 * b1 - 4. * a1 * c1) <= zero) then
!continue
!end if
      pg(npi)%VolFra = (- b1 - Dsqrt( b1 * b1 - 4.d0 * a1 * c1)) / (two * a1)   ! nuovo calcolo della frazione di volume
!      pg(npi)%VolFra = vol_temp + (pg(npi)%VolFra - vol_temp) / 10
!
!§
!      if (pg(npi)%VolFra >= one) then
      if (pg(npi)%VolFra >= VFmx) then
        pg(npi)%VolFra  = VFmx
!§
        pg(npi)%rhoc = pg(npi)%dens 
        pg(npi)%rhow = Med(1)%den0 
        pg(npi)%mass = pg(npi)%dens * (Domain%dd**ncord)
!§
!      else if(pg(npi)%VolFra <= zero)then
      else if(pg(npi)%VolFra <= VFmn)then
        pg(npi)%VolFra = VFmn
!§
        pg(npi)%rhoc = Med(2)%den0 
        pg(npi)%rhow = pg(npi)%dens
        pg(npi)%mass = pg(npi)%dens * (Domain%dd**ncord)
      else
        pg(npi)%rhoc = tirhoc / pg(npi)%VolFra            ! new dens grain
        pg(npi)%rhow = tirhow / (one-pg(npi)%VolFra)      ! new dens wat
        pg(npi)%mass = (pg(npi)%VolFra * pg(npi)%rhoc + (one-pg(npi)%VolFra)*pg(npi)%rhow) * (Domain%dd**ncord)
      end if
!
    end if
!
  end do       
!
!$omp end parallel do
!
return
end subroutine AggDens
!---split


!cfile calc_pelo.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : calc_pelo
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
! Module purpose : Module to write results for free surface
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
!
subroutine calc_pelo
!
!.. assign modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Scalars ..
integer(4)       :: i,j, ncel,npartcel,mm,jj  !!!!!,icord,n 
!!!!!double precision :: ro1, ro2, romezzi
!!!!!logical          :: trovato
!
!.. Local Arrays ..
integer(4),      dimension(1)  :: minpos1,minpos2
double precision,dimension(3)  :: ragtemp
double precision,dimension(3,nlines)  :: pelolib
integer(4),      dimension(:),allocatable :: PartCelnum
double precision,dimension(:),allocatable :: PartCel
 
!
!.. Executable Statements ..
!
!!!!! pelolib = zero
!!!!! romezzi = Med(Domain%ipllb_md)%den0 * half
!!!!! trovato = .false.
!!!!!
!!!!! do i = 1,nlines
!!!!!
!!!!!    POINTS_LOOP: do j = control_lines(i)%Icont(1)+1,control_lines(i)%Icont(2)
!!!!!!       do n = 1,ncord
!!!!!!          icord = icoordp(n,ncord-1)
!!!!!!          pelolib(n,i) = control_points(j-1)%coord(icord)
!!!!!!       end do
!!!!!!       ro1 = control_points(j-1)%dens * control_points(j-1)%uni
!!!!!!       ro2 = control_points(j  )%dens * control_points(j  )%uni
!!!!!       ro1 = control_points(j-1)%dens * half
!!!!!       ro2 = control_points(j  )%dens * half
!!!!!
!!!!!       if ( ro2 < romezzi ) then
!!!!!          trovato = .true.
!!!!!          do n = 1,ncord
!!!!!             icord = icoordp(n,ncord-1)
!!!!!             pelolib(n,i) = control_points(j-1)%coord(icord) + &
!!!!!                            ( romezzi - ro1 ) / ( ro2 - ro1 ) * &
!!!!!                            ( control_points(j)%coord(icord) - control_points(j-1)%coord(icord) )
!!!!!          end do
!!!!!          exit POINTS_LOOP
!!!!!       end if
!!!!!    end do POINTS_LOOP
!!!!! end do
!!!!!
!!!!! if (trovato) write (nplb,'(30g14.7)') tempo,pelolib

 allocate (PartCelnum(NMAXPARTJ), PartCel(NMAXPARTJ))
!.. loop su tutte le linee
 Pelolib = 0
 do i = 1,nlines
!.. loop sui punti della linea
   POINTS_LOOP: do j = control_lines(i)%Icont(1)+1,control_lines(i)%Icont(2)
     ncel = control_points(j)%cella
     if (ncel == 0) cycle    ! cella fuori campo
!.. se la cella a cui appartiene il punto non contiene particelle esco
     if (Icont(ncel+1) <= Icont(ncel)) cycle
!.. altrimenti trovo le due particelle piu' vicine al punto e faccio la media della posizione
     nPartCel = 0
     PartCelnum = 0
     PartCel = 99999
     do mm = Icont(ncel),Icont(ncel+1)-1  ! +++ loop sulle part di una cella
       jj = NPartOrd(mm)
!.. calcola le conponenti deltaX, deltaY e deltaZ della distanza tra il punto e la particella
       ragtemp(1:3) = control_points(j)%coord(1:3) - pg(jj)%coord(1:3)
!.. calcola la distanza effettiva (considera direttamente il quadrato per aumentare l'accuratezza)
       npartcel = npartcel + 1
       PartCelnum(npartcel) = jj
       PartCel(npartcel) = ragtemp(1)*ragtemp(1) + ragtemp(2)*ragtemp(2) + ragtemp(3)*ragtemp(3)
     end do  ! +++ loop sulle part di una cella
!.. cerco le particelle piu' vicine al punto
     minpos1 = minloc(PartCel,1)
     PartCel(minpos1) = 9999.
     minpos2 = minloc(PartCel,1)
 !.. posizione del pelo libero come media della posizione delle due particelle piu' vicine al punto
     pelolib(1:3,i) = (pg(PartCelnum(minpos1(1)))%coord(1:3) + pg(PartCelnum(minpos2(1)))%coord(1:3)) * half
   end do POINTS_LOOP
 end do
 write (nplb,'(30g14.7)') tempo,pelolib
 deallocate (PartCelnum,PartCel)

return
end subroutine calc_pelo
!---split

!cfile calcpre.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : calcpre
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
! Module purpose : Module calculation of pressure in every particles of the general 
!                  field
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D, time_integration
!
! Called routines: 
!
!************************************************************************************
!
  subroutine CalcPre
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
  integer(4)       :: npi
  double precision :: rhorif,c2,crhorif,wrhorif,wc2,cc2 !,maxpExpl   !,presc,presw
!
!.. Executable Statements ..
!
!  maxpExpl = max_negative_number
!
!.. modello bifluido
  if (diffusione) then
!
!$omp parallel do default(none) private(npi,crhorif,wrhorif,wc2,cc2) shared(nag,Pg,Med)
!
    do npi = 1,nag
!
      if (pg(npi)%koddens /= 0) cycle
!
      if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!
      crhorif = Med(2)%den0
      wrhorif = Med(1)%den0
      wc2     = Med(1)%celerita * Med(1)%celerita
      cc2     = Med(2)%celerita * Med(2)%celerita
!§
      if (pg(npi)%imed == 1)then
!        presc = cc2 * (pg(npi)%rhoc * pg(npi)%VolFra - crhorif * VFmn)
!        presw = wc2 * (pg(npi)%rhow  * (one-pg(npi)%VolFra) - wrhorif * (one-VFmn))
        pg(npi)%pres = wc2 * (pg(npi)%dens - (crhorif * VFmn + wrhorif * (1-VFmn))) !ALTERNATIVO
      else if (pg(npi)%imed == 2)then
!        presc = cc2 * (pg(npi)%rhoc * pg(npi)%VolFra - crhorif * VFmx)
!        presw = wc2 * (pg(npi)%rhow  * (one-pg(npi)%VolFra) - wrhorif * (one-VFmx))
        pg(npi)%pres = cc2 * (pg(npi)%dens - (crhorif * VFmx + wrhorif * (1-VFmx))) !ALTERNATIVO
      end if
!      pg(npi)%pres = presc + presw
!§
!      presc = cc2 * (pg(npi)%rhoc - crhorif)
!      presw = wc2 * (pg(npi)%rhow - wrhorif)
!      pg(npi)%pres = pg(npi)%VolFra * presc + (one - pg(npi)%VolFra) * presw
!
    end do
!
!$omp end parallel do
!
  else
!
    if (it_corrente > 0) then  !inutile ??
!
!AA501 sub
!$omp parallel do default(none) private(npi,rhorif,c2) shared(nag,pg,Domain,Med,esplosione) !,maxpExpl)
!
!.. loops on all the particles
!
      do npi = 1,nag
!
!.. skips the outgone particles
!.. skips the particles with velocity type different from "standard"
!
        if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!
        if (esplosione) then
!...................................... 2011 mar 08
!.. modify for Specific Internal Energy
          pg(npi)%pres = (Med(pg(npi)%imed)%gamma - one) * pg(npi)%IntEn * pg(npi)%dens
!          maxpExpl = max(maxpExpl,pg(npi)%pres)
!......................................
        else
!.. evaluates the new pressure condition for the particle
!ç
!!!!        if (it_corrente > 0) then  !.and. pg(npi)%state /= 'sol') then
          rhorif     = Med(pg(npi)%imed)%den0
          c2         = Med(pg(npi)%imed)%eps / rhorif
!
!AA406
           if ((pg(npi)%dens - rhorif) /= 0.) pg(npi)%pres = c2 * (pg(npi)%dens - rhorif)
!!!!        end if
!
!AA406test
!
!AA601rm
        end if
!
      end do
!
!......................................
!!!!      if (maxpExpl < pg(??)%pres) esplosione = .false.
!      if (maxpExpl < 0.99e+9) esplosione = .false.
!......................................
!
!$omp end parallel do
!
    end if
!
  end if
!
  return
  end subroutine CalcPre
!---split

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

!cfile CalcVarp.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : CalcVarp
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
! Module purpose : Module to calculate p on point (xpu,ypu)
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: GetVarPart
!
!************************************************************************************
!
subroutine CalcVarp 
!* calcola grandezze nei punti di monitor
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
integer(4) :: i,ii,jj,kk,pointcellnumber
double precision xp,yp,zp
type (TyCtlPoint) :: pglocal
!
!.. External routines ..
  integer(4), external   :: CellNumber
!
!.. Executable Statements ..
!
  if (Npointst < 1) return

  do i = 1, Npointst

    pglocal%coord(:) = control_points(i)%coord(:)
    pglocal%pres     = zero
    pglocal%dens     = zero
    pglocal%vel(:)   = zero
    pglocal%uni      = zero
    xp = pglocal%coord(1) - Grid%extr(1,1)
    yp = pglocal%coord(2) - Grid%extr(2,1)
    zp = pglocal%coord(3) - Grid%extr(3,1)
    ii = ceiling(xp / Grid%dcd(1))
    jj = ceiling(yp / Grid%dcd(2))
    kk = ceiling(zp / Grid%dcd(3)) 
    pointcellnumber = CellNumber(ii, jj, kk)
    pglocal%cella   = pointcellnumber

    call GetVarPart (pglocal)

    control_points(i)%pres   = pglocal%pres
    control_points(i)%dens   = pglocal%dens
    control_points(i)%vel(:) = pglocal%vel(:)
    control_points(i)%uni    = pglocal%uni
    control_points(i)%cella  = pglocal%cella

  end do

return
end subroutine CalcVarp
!---split


!cfile CellIndices.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
  integer(4) function CellIndices (nc, i, j, k)

!.Returns indices i, j, k, of the cell nc
!..in a 3D box frame with ni, nj, nk, cells on each side
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
  integer(4),intent(IN)  :: nc
  integer(4),intent(OUT) :: i,j,k
!
!.. Local Scalars ..
  integer(4) :: ncij, nucellsij,ni,nj   !,nk
!
!.. Executable Statements ..
!
  ni = Grid%ncd(1)
  nj = Grid%ncd(2)
!  nk = Grid%ncd(3)
!
  nucellsij = ni * nj
!
!.. grid index in the Z direction
!
  k = int((nc - 1) / nucellsij) + 1
  ncij = nc - nucellsij * (k - 1)
!
!.. grid index in the Y direction
!
  j = int((ncij - 1) / ni) + 1
!
!.. grid index in the X direction
!
  i = ncij - ni * (j - 1)
!
  CellIndices = ncij
!
  End function CellIndices
!---split

!cfile CellNumber.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
  integer(4) function CellNumber (i, j, k)

!Returns the number of the cell sited at integer coordinates i, j, k
!in a 3D box frame with ni, nj, nk cells on each side
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
  integer(4),intent(IN) :: i,j,k
!
!.. scalars
  integer(4) :: ni,nj,nk
!
!.. Executable Statements ..
!
  ni = Grid%ncd(1)
  nj = Grid%ncd(2)
  nk = Grid%ncd(3)
!
  if ( i<1 .or. i>ni .or. j<1 .or. j>nj .or. k<1 .or. k>nk) then
!
!.. the cell is outside the grid limits
!
    CellNumber = 0  !-1
!
  else
!
!.. return the cell number
!
    CellNumber = ((k - 1) * nj + (j - 1)) * ni + i
!
  end if
!
  End Function CellNumber
!---split

!cfile contrmach.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
!subroutine contrmach (k,celmax)
!
!use GLOBAL_MODULE
!use AdM_USER_TYPE
!use ALLOC_MODULE
!
!implicit none
!
!integer(4) :: k,
!
!double precision    :: celmax,amaxnmach,amaxnmach2,amachnumb2
!!double precision   :: vmod2,cel2,amachnumb
!
! amaxnmach  = 0.2d0
! cel2       = celmax*celmax
! amaxnmach2 = amaxnmach*amaxnmach
!
! controllo numero mach - scrive su nout
! vmod2 = pg(k)%vel(1)*pg(k)%vel(1) + pg(k)%vel(3)*pg(k)%vel(3)
! amachnumb2 = vmod2 / cel2
!
! if ( amachnumb2 > amaxnmach2 )then    ! -------------------------------
!    amachnumb = Dsqrt(amachnumb2)
!!   pg(k)%vel(:)=(amaxnmach*celmax/Dsqrt(vmod2))*pg(k)%vel(:) 
! end if        ! -------------------------------
!
!return
!end
!---split

!cfile CreaGrid.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : CreaGrid
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
! Module purpose : Module to create a grid for sorting on the whole domain
!
! Calling routine: Gest_Input
!
! Called routines: diagnostic
!
!************************************************************************************
!
  subroutine CreaGrid
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
!  double precision, parameter    :: epsi = 0.001d0
!
!.. local variables
  integer(4) :: ier
  double precision       :: epsi
  character(len=lencard) :: nomsub = "CreaGrid"
!
!.. local arrays
  double precision, dimension(3) :: dextr
!
!.. Executable Statements ..
!
  epsi = 0.01d0 * Domain%dd
!!  epsi = 0.1d0 * Domain%dd
!
!.. set the vertices of the virtual grid on the base of the domain dimensions and the particle smooth,
!.. considering a tolerance epsi on the coordinates for numerical reasons
!
  Grid%extr(:,1) = Domain%coord(:,1) - doubleh - epsi
!
  Grid%extr(:,2) = Domain%coord(:,2) + doubleh + epsi
!
!.. evaluates the dimensions of the grid in all the xyz directions
! 
  dextr(:)       = Grid%extr(:,2) - Grid%extr(:,1)
!
!.. evaluates the number of grid cells in all the directions, assuming that the cell is a cube having 
!.. the side length equal to the double of the smooth length
! 
  Grid%ncd(:)    = NINT(dextr(:) / doubleh)
!
!.. evaluates the real cell dimensions in all the directions, transforming the cubic cell into a real 
!.. hesaedric regular cell
!
  Grid%dcd(:)    = dextr(:) / Grid%ncd(:)
!
!.. in 2D calculation, the number of cells in the Y direction is forced to 1
!
  if ( ncord == 2 ) Grid%ncd(2) = 1
!
!.. evaluates the maximum number of virtual cells in the grid covering the domain parallelepiped
!
  Grid%nmax = Grid%ncd(1) * Grid%ncd(2) * Grid%ncd(3)
!
  write (nout,'(1x,a)') " "
  write (nout,'(1x,a,3i8)') " Number of grid in x, y, z directions : ",Grid%ncd(1),Grid%ncd(2),Grid%ncd(3)
  write (nout,'(1x,a,i10)') " Number of total grid : ",Grid%nmax
  write (nout,'(1x,a)') " "
!
!.. allocazione matrice 2d per calcolo pelolibero (caso erosione)
!
!AA504 sub the fifth element is removed all over the code
  allocate (ind_interfaces(Grid%ncd(1),Grid%ncd(2),4), stat = ier)
  if (ier /= 0) then
!AA504 sub      
    write (nout,'(1x,a,i2)') "    Array ind_interfaces not allocated. Error code: ",ier
    call diagnostic (arg1=4,arg3=nomsub)
  else
!AA504 sub      
    write (nout,'(1x,a)') "    Array ind_interfaces successfully allocated "
  end if
!
  return
  end subroutine CreaGrid
!---split

!cfile CreateSectionPoints.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine  CreateSectionPoints ( vp, values, opt, seccor )
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. implicit declarations
  implicit none
!
!.. dummy arguments
integer(4)   :: seccor
character(1) :: opt
double precision, dimension(3)   :: vp
double precision, dimension(3,2) :: values
!
!.. local scalars
integer(4)   :: i,j,k,n, npse
integer(4),      dimension(3)   :: Nmesh
double precision,dimension(3)   :: Cc,CcStart

character(80), external :: lcase
!
!.. executable statements
!
!  creazione mesh di lato dd
 Nmesh = 1
 npse = Control_Sections(seccor)%icont(1)
!
 CcStart(:) = vp(:) - Domain%dd
 do n = 1,SPACEDIM
   if ( n == 1 .AND. lcase(opt) == "x" ) cycle
   if ( n == 2 .AND. lcase(opt) == "y" ) cycle
   if ( n == 3 .AND. lcase(opt) == "z" ) cycle
   Nmesh(n)   = nint ( (values(n,2)-values(n,1)) / Domain%dd )
   CcStart(n) = values(n,1) - Domain%dd * half
 end do

 Cc(1) = CcStart(1)

 do i = 1, Nmesh(1)
   Cc(1) = Cc(1) + Domain%dd
   Cc(2) = CcStart(2)
   do j = 1, Nmesh(2)
     Cc(2) = Cc(2) + Domain%dd
     Cc(3) = CcStart(3)
     do k = 1, Nmesh(3)
       Cc(3) = Cc(3) + Domain%dd
       npse = npse + 1
       Control_Points(npse)%coord(:) = Cc(:)
     end do
   end do
 end do
 
return
end subroutine  CreateSectionPoints
!---split

!cfile defcolpartzero.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : defcolpartzero
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
! Module purpose : Module for particle color definition as from input file
!
! Calling routine: GenerateSourceParticles_2D, GenerateSourceParticles_3D
!                  SetParticles
!
! Called routines: 
!
!************************************************************************************
!
subroutine defcolpartzero (ir,partz,pg)
!* definisce colore particella nello stato zero
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
integer(4),       intent(IN)                         :: ir
type (TyZone),    intent(IN),   dimension(NPartZone) :: partz
type (TyParticle),intent(INOUT)                      :: pg
!
!.. Local Scalars ..
integer(4)       :: nbande, numbanda
double precision :: aldx
!
!.. Local Arrays ..
integer(4), dimension(5) :: iclnumb
!
!.. Executable Statements ..
!
 iclnumb(1)=1
 iclnumb(2)=2
 iclnumb(3)=4
 iclnumb(4)=5
 iclnumb(5)=6

 if ( partz(ir)%bend == "u")then        ! color uniform
    pg%icol = partz(ir)%icol

 else if ( partz(ir)%bend == "o")then   ! color optional (solo su opzione esterna)
    pg%icol = partz(ir)%icol

 else if(partz(ir)%bend == "b")then     ! bande verticali
    nbande = partz(ir)%icol
    aldx   =( partz(ir)%coordMM(1,2)-partz(ir)%coordMM(1,1))/nbande
    numbanda = int((pg%coord(1)-partz(ir)%coordMM(1,1))/aldx)+1
    numbanda = min(nbande,numbanda)
    numbanda = max(0,numbanda)
    pg%icol = iclnumb(numbanda)
 end if

return
end subroutine defcolpartzero
!---split

!cfile diagnostic.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : diagnostic
!
! Last updating : May 02, 2014
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Agate             08/05/2012     Limited and licensed version
! 04  Agate             27/03/2013     Limited and licensed version for single and multi users
! 05  Agate             02/05/14       Add more modules check in license
!
!************************************************************************************
! Module purpose : Diagnostic management
!
! Calling routine: all
!
! Called routines: memo_results
!                  print_results
!                  result_converter
!
!************************************************************************************
!
  subroutine diagnostic (ierr,ivalue,avalue)
!
!.. assign modules
  use FILES_ENTITIES
  use GLOBAL_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
  integer(4), intent(in)                      :: ierr
  integer(4), intent(in), optional            :: ivalue
  character(LEN=lencard), intent(in),optional :: avalue
!
!.. Local Scalars ..
  character(LEN=2)       :: error
  integer(4)             :: dato,la,it1,it2,it3
  double precision       :: dtvel
  character(len=lencard) :: stringa
  logical                :: print_out
!
!.. Executable Statements ..
!
  print_out = .false.
  if (present(avalue)) then
    stringa = adjustl(avalue)
    la = len_trim(stringa)
  else
    stringa = ' '
    la = 0
  end if
  if (present(ivalue)) then
    dato = ivalue
  else
    dato = -9999
  end if
!
  write (error,'(i2)') ierr
  write (nscr,'(a,a2,a)') ' >> Diagnostic ',error,' detected * Run ended unsuccessfully'
  if (ierr /= 1) then
    write (nout,'(1x,100(''=''))')
    write (nout,'(31x,5a)') ' *  ',acode,'  version  ',version,'  * '
    write (nout,*)
    write (nout,'(a,a2,a)') 'Diagnostic ',error,' detected * Run ended unsuccessfully'
    write (nout,*)
    flush(nout)
  end if
!
  select case (ierr)
    case ( 1)    ! Errors in arguments
      write (nscr,*) "-----------------------------------"
      write (nscr,*) "ERROR in command line argument !!!."
      write (nscr,*) "-----------------------------------"
      write (nscr,*) " Correct command sequence-> <executablename> <casename> [euristic]|[original]"
!
    case ( 2)    ! Errors in files
      write (nout,'(1x,a)') 'Checkfile function detected an error:'
      select case (dato)
        case (1); write (nout,'(1x,a)') '  1 - Output file cannot be opened.'
        case (2); write (nout,'(1x,a)') '  2 - Input file does not exist.'
        case (3); write (nout,'(1x,a)') '  3 - Input file cannot be opened.'
      end select
!
    case ( 3)    ! Errors in arrays
      write (nout,'(1x,2a)') 'Array deallocation failed.  Routine -> ',stringa(1:la)
      write (nout,'(1x,a)') 'Run is stopped.'
!
    case ( 4)    ! Errors in arrays
      write (nout,'(1x,2a)') 'Array allocation/deallocation failed.  Routine -> ',stringa(1:la)
      write (nout,'(1x,a)') 'Run is stopped.'
!
    case ( 5)    ! Errors in input
      write (nout,'(1x,a)') 'An error was detected analyzing the input data deck for the line:'
      write (nout,'(1x,a)') '->'//stringa(1:la)//'<-'
      write (nout,'(1x,a)') 'Specific code error is:'
      select case (dato)
        case (1); write (nout,'(1x,a)')   '    1 - Unknown Input Option'
        case (2); write (nout,'(1x,a)')   "    2 - Input file doesn't match program versions"
        case (3); write (nout,'(1x,a)')   '    3 - Unknown Domain Type. Routine -> '
        case (4); write (nout,'(1x,a)')   '    4 - Wrong type data or unsufficient number of data'
        case (5); write (nout,'(1x,a)')   '    5 - Type of Domain = rect is possible only with 2D problems '
        case (6); write (nout,'(1x,a)')   '    6 - Type of Erosion model is not correct (shields, mohr) '
        case (101); write (nout,'(1x,a)') '  101 - Unrecognised boundary type'
        case (102); write (nout,'(1x,a)') '  102 - PERIMETER boundary type: different first and last vertex'
        case (103); write (nout,'(1x,a)') '  103 - TAPIS boundary type: 2 vertices are requested'
        case (104); write (nout,'(1x,a)') '  104 - FACE definition: face already defined'
        case (201); write (nout,'(1x,a)') '  201 - Restart file cannot be opened'
        case (203); write (nout,'(1x,a)') '  203 - Restart step location failed'
        case (204); write (nout,'(1x,a)') '  204 - Restart '//stringa(1:la)//' fails at the record reported above'
        case (205); write (nout,'(1x,a)') '  205 - Restart option is unknown'
        case (301); write (nout,'(1x,a)') '  301 - Ascii input option is unknown'
        case (302); write (nout,'(1x,a)') '  302 - Ascii input version does not match code version'
        case (304); write (nout,'(1x,a)') '  304 - Ascii input '//stringa(1:la)//' fails at the record reported above'
        case (305); write (nout,'(1x,a)') '  305 - Ascii input '//stringa(1:la)//' fails at the record for level option'
        case (306); write (nout,'(1x,a)') '  306 - Ascii input '//stringa(1:la)//' fails at the record for erosion Shields Mohr module'
        case (307); write (nout,'(1x,a)') '  307 - Ascii input '//stringa(1:la)//' fails at the record for diffusion module'
        case (308); write (nout,'(1x,a)') '  308 - Ascii input '//stringa(1:la)//' fails at the record for esplosion module'
        case (309); write (nout,'(1x,a)') '  309 - Ascii input '//stringa(1:la)//' fails at the record for multi fluid module'
        case (310); write (nout,'(1x,a)') '  310 - Ascii input '//stringa(1:la)//' fails at the record for more fluid module'
        case (311); write (nout,'(1x,a)') '  311 - Ascii input '//stringa(1:la)//' fails at the record for Granular flux module'
        case (312); write (nout,'(1x,a)') '  312 - Ascii input '//stringa(1:la)//' fails at the record for erosion Shields Van Rijn Seminara module'
        case (320); write (nout,'(1x,a)') '  320 - Ascii input '//stringa(1:la)//' fails at the record for temporal scheme module'
        case (330); write (nout,'(1x,a)') '  330 - Ascii input '//stringa(1:la)//' fails at the record for body dynamics module'
        case (340); write (nout,'(1x,a)') '  340 - Ascii input '//stringa(1:la)//' fails at the record for DB-SPH module'
        case default; write (nout,'(1x,a)') '  ??? - Unknown error type'
      end select
!
    case ( 6)    ! Errors in Source
      write (nout,'(1x,a)') 'An error was detected analyzing the source of particles'
      write (nout,'(1x,a)') 'Specific code error is:'
      select case (dato)
        case (1); write (nout,'(1x,2a)') ' 1 - Number of particles nag > PARTICLEBUFFER (a).  Routine ->',stringa(1:la)
        case (2); write (nout,'(1x,2a)') ' 2 - Number of particles nag > PARTICLEBUFFER (b).  Routine ->',stringa(1:la)
      end select
!
    case ( 7); write (nout,'(1x,2a)') 'Number of cels in grid: NumCellmax < Grid%nmax.  Routine ->',stringa(1:la)
!
    case ( 8)    ! Errors in boundaries
      write (nout,'(1x,a)') 'An error was detected analyzing the boundaries'
      write (nout,'(1x,a)') 'Specific code error is:'
      select case (dato)
        case (1);  write (nout,'(1x,a,i10,2a)') ' 1 - New Max num particles*BoundaryCloseSides the new value is: MaxNcbs = ', &
                                                MaxNcbs,'.  Routine ->',stringa(1:la)
        case (2);  write (nout,'(1x,a,i10,2a)') ' 2 - New Max num particles*BoundaryCloseFaces the new value is: MaxNcbf = ', &
                                                MaxNcbf,'.  Routine ->',stringa(1:la)
        case (3);  write (nout,'(1x,2a)')       ' 3 - The kernel table is not implemented for 2D case.  Routine ->',stringa(1:la)
        case (4);  write (nout,'(1x,2a)')       ' 4 - Sides are not consecutive.  Routine ->',stringa(1:la)
        case (5);  write (nout,'(1x,3a)')       ' 5 - Number of cell calculated is different from the number stored for the ', &
                                               'current particle.  Routine ->',stringa(1:la)
        case (6);  write (nout,'(1x,2a)')       ' 6 - Number of close boundary faces ncbf > MAXCLOSEBOUNDFACES (a).  Routine ->', &
                                                stringa(1:la)
        case (7);  write (nout,'(1x,2a)')       ' 7 - Number of close boundary faces ncbf > MAXCLOSEBOUNDFACES (b).  Routine ->', &
                                                stringa(1:la)
        case (8);  write (nout,'(1x,2a)')       ' 8 - Number of close boundary sides Ncbs > MAXCLOSEBOUNDSIDES.  Routine ->', &
                                                stringa(1:la)
        case (9);  write (nout,'(1x,2a)')       ' 9 - Number of close boundary sides Ncbs > LIMCLOSEBOUNDSIDES.  Routine ->', &
                                                stringa(1:la)
        case (10); write (nout,'(1x,2a)')       '10 - Number of convex edges NumBEdges > MAXNUMCONVEXEDGES.  Routine ->', &
                                                stringa(1:la)
      end select
!
    case ( 9)    ! Errors in Loop
      write (nout,'(1x,a)') 'An error was detected during the loop'
      write (nout,'(1x,a)') 'Specific code error is:'
      select case (dato)
        case (1); write (nout,'(1x,a,i10,2a)') ' 1 - Number of fluid particles array Array_Flu greater than max ',PARTICLEBUFFER, &
                                               ' (a).  Routine ->',stringa(1:la)
        case (2); write (nout,'(1x,a,i10,2a)') ' 2 - Number of fluid particles array Array_Flu greater than max ',PARTICLEBUFFER, &
                                               ' (b).  Routine ->',stringa(1:la)
        case (3); write (nout,'(1x,2a)')       ' 3 - Increase parameter MAXCLOSEBOUNDFACES.  Routine ->',stringa(1:la)
      end select
      print_out = .true.
!
    case (10)    ! Errors in Tools
      write (nout,'(1x,a)') 'An error was detected during the tools execution'
      write (nout,'(1x,a)') 'Specific code error is:'
      select case (dato)
        case (1); write (nout,'(1x,a,i10,2a)') ' 1 - Number of surrounding particles of current particle is greater than max ', &
                                                NMAXPARTJ,' (a).  Routine ->',stringa(1:la)
        case (2); write (nout,'(1x,2a)')       ' 2 - Array BoundaryVertex bounds exceed.  Routine ->',stringa(1:la)
        case (3); write (nout,'(1x,2a)')       ' 3 - Array Vertice bounds exceed.  Routine ->',stringa(1:la)
        case (4); write (nout,'(1x,2a)')       ' 4 - Number of particles nag > PARTICLEBUFFER.  Routine ->',stringa(1:la)
        case (5); write (nout,'(1x,2a)')       ' 5 - Unknown Domain Type. Routine -> ',stringa(1:la)
        case (6); write (nout,'(1x,2a)')       ' 6 - Number of open sides NumOpenSides > MAXOPENSIDES. Routine -> ',stringa(1:la)
        case (7); write (nout,'(1x,2a)')       ' 7 - Number of open faces NumOpenFaces > MAXOPENFACES. Routine -> ',stringa(1:la)
        case (8); write (nout,'(1x,2a)')       ' 8 - Free level condition is not found. Routine -> ',stringa(1:la)
        case (88); write (nout,'(1x,2a)')      '88 - Error deltat law calculation for "law" particles. Routine -> ',stringa(1:la)
      end select
      print_out = .true.
!
    case (11)    ! Errors in Erosion Model
      write (nout,'(1x,a)') 'An error was detected during execution of erosion model'
      write (nout,'(1x,a)') 'Specific code error is:'
      select case (dato)
        case (1); write (nout,'(1x,2a)') ' 1 - The Ustar in Shields erosion model is NaN.  Routine ->',stringa(1:la)
        case (2); write (nout,'(1x,3a)') ' 2 - It is Impossible to find liquid interface and solid interface for the ', &
                                         'current particle.  Routine ->',stringa(1:la)
      end select
      print_out = .true.
!
    case (73)    ! Errors in License file
      write (nout,'(1x,a)')'License is not correct. Code error:'
      select case (dato)
        case (1) ; write (nout,'(1x,a)')  '1  - The license file has not been provided with the installation'
        case (2) ; write (nout,'(1x,3a)') '2  - Running machine name ',stringa(1:la),' does not match the license machine name'
        case (3) ; write (nout,'(1x,3a)') '3  - Running machine address ',stringa(1:la),' does not match the license machine address'
        case (4) ; write (nout,'(1x,3a)') '4  - Code release ',stringa(1:la),' does not match the license release'
        case (5) ; write (nout,'(1x,3a)') '5  - The user ',stringa(1:la),' does not match the license user'
        case (6) ; write (nout,'(1x,3a)') '6  - No modules are licensed: at least one module must be licensed'
        case (7) ; write (nout,'(1x,a)')  '7  - License has been expired',stringa(1:la),' days ago'
        case (8) ; write (nout,'(1x,a)')  '8  - License file is corrupted'
        case (9) ; write (nout,'(1x,a)')  '9  - The master node of the Linux cluster has not been setted in the proper environmental variable (MASTER_NAME)'
        case (10); write (nout,'(1x,a)')  '10 - Error in <whoami> call. License key cannot be verified'
        case (11); write (nout,'(1x,a)')  '11 - Error in <',stringa(1:la),'> call. License key cannot be verified'
        case (12); write (nout,'(1x,a)')  '12 - Error in the name of the executable file: the extension must be "exe" (Windows) or "x" (Linux/Unix)'
        case default 
          write (nout,*) 'Generic error reading the license. Please contact the assistance.'
       end select
!
    case default 
      write (nout,*) 'Unpredictable error'
      print_out = .true.
  end select
  if (ierr /= 1) then
    write (nout,'(1x,100(''=''))')
    flush(nout)
  end if
!
  if (print_out) then
    write (nout,'(1x,(a))') ' ---------------------------'
    write (nout,'(1x,(a))') ' Print results at last step.'
    write (nout,'(1x,(a))') ' ---------------------------'
    write (nscr,'(1x,(a))') ' ---------------------------'
    write (nscr,'(1x,(a))') ' Print results at last step.'
    write (nscr,'(1x,(a))') ' ---------------------------'
!
!=== SCRITTURA SU FILE RISULTATI ######################
!
    it1 = 0
    it2 = 0
    it3 = 0
    dtvel = 99999.0
!Salvataggio risultati e restart ULTIMO STEP (sempre?)
    if (nout > 0) then
      call print_results ( it1, it2, 'fine__' )
    end if
    if ( nres > 0 ) then
      call memo_results ( it1, it2, it3, dtvel, 'fine__' )
    end if
!
!=== SCRITTURA SU FILE FORMATO VTK ######################
!
    if ( vtkconv ) then
      call result_converter ('fine__')
    end if
!
  end if
!
  stop
!
  end subroutine diagnostic
!---split

!cfile diffumorris.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : diffumorris
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
! Module purpose : Module to 
! 
! Calling routine: inter_EqCont_2d, inter_EqCont_3d
!
! Called routines: 
!
!************************************************************************************
!
subroutine diffumorris (npi,npj,npartint,dervol,factdiff,rvw)
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
 integer(4)       :: npi,npj,npartint
 double precision :: dervol,factdiff,rvw
!
!.. Local Scalars ..
! double precision        :: gradd
! double precision        :: dvar(3),cuei,eta
 double precision        :: anuitilde,rhotilde,amassj,coei,coej
! double precision        :: coeffveli,coeffvelj,dist,distj,disti
! double precision        :: dvi(3),dvj(3)
!
!.. Executable Statements ..
!
! coeffvelj = (Dsqrt(pg(npj)%dv(1)*pg(npj)%dv(1)+pg(npj)%dv(3)*pg(npj)%dv(3)))
! coeffveli = (Dsqrt(pg(npi)%dv(1)*pg(npi)%dv(1)+pg(npi)%dv(3)*pg(npi)%dv(3)))
! coej=two*pg(npj)%coefdif * coeffvelj
! coei=two*pg(npi)%coefdif * coeffveli

 pg(npj)%coefdif = med(pg(npj)%imed)%codif
 pg(npi)%coefdif = med(pg(npi)%imed)%codif
 
 coej = pg(npj)%coefdif
 coei = pg(npi)%coefdif

!! Schema armonico
 amassj    = pg(npj)%mass
 anuitilde = 4.0d0 * (coei * coej) / (coei + coej + 0.00001d0)
 rhotilde  = pg(npj)%dens

 if (pg(npi)%visc == med(2)%numx .or. pg(npj)%visc == med(2)%numx) then
   anuitilde = zero
 end if
 
 if ( pg(npj)%vel_type /= "std" ) then          !non part fix o altro
    rhotilde  = pg(npi)%dens
    amassj    = pg(npi)%mass
    anuitilde = zero
 end if

!.. formula utilizzata nel nov 2005 da rapporto Bon/Gatti/Zuccala/DiMonaco/Gallati
! factdiff = amassj * anuitilde / rhotilde
!! eta      = 0.01d0 * hh
! gradd    = gradmod / (rij+eta)
! rvw      =-gradd * dervol

! disti = pg(npi)%coord(1) * pg(npi)%coord(1) + pg(npi)%coord(3) * pg(npi)%coord(3)
! distj = pg(npj)%coord(1) * pg(npj)%coord(1) + pg(npj)%coord(3) * pg(npj)%coord(3)
! dist  = distj - disti
! cuei  = -two * gradd * dist * amassj / pg(npj)%dens

!.. formula utilizzata nel nov 2007 come da letteratura
 factdiff = amassj * anuitilde / rhotilde
! eta     = 0.01d0 * hh * hh
 rvw     = -dervol * PartKernel(2,npartint) * &
            (rag(1,npartint)*rag(1,npartint) + rag(2,npartint)*rag(2,npartint) + rag(3,npartint)*rag(3,npartint)) 

return
end subroutine diffumorris
!---split

!cfile FindFrame.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
  subroutine FindFrame (Xmin, Xmax, Nt)
!
!Finds extremes of the rectangular frame which contains the boundary mib
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_User_Type
  use ALLOC_MODULE
!
!.. implicit declarations
  implicit none
!
!.. dummy arguments
  integer(4),intent(IN) :: Nt
  double precision,intent(INOUT),dimension(SPACEDIM,NumFacce) :: Xmin, Xmax
!
!.. local scalars
  integer(4) :: i, n, iv, nf, nod
!
!.. executable statements
!
  do iv = Tratto(Nt)%iniface, Tratto(Nt)%iniface + Tratto(Nt)%numvertices - 1
    nf = BFaceList(iv)
    do n = 1, 4
      nod = BoundaryFace(nf)%Node(n)%name
      if ( nod <= 0 ) cycle
      do i = 1, Ncord
        if ( Vertice(i,nod) < Xmin(i,Nt) ) Xmin(i,Nt) = Vertice(i,nod)
        if ( Vertice(i,nod) > Xmax(i,Nt) ) Xmax(i,Nt) = Vertice(i,nod)
      end do
    end do
  end do
!
  return
  end subroutine FindFrame
!---split

!cfile FindLine.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
  subroutine FindLine (Xmin, Xmax, Nt)
!
!..Finds extremes of the rectangular frame which contains the boundary mib
!
!.. assign modules ..
  use FILES_ENTITIES
  use GLOBAL_MODULE
  use AdM_User_Type
  use ALLOC_MODULE
  use DIAGNOSTIC_MODULE
!
!.. implicit declarations ..
  implicit none
!
!.. dummy arguments ..
  integer(4),      intent(IN)                                   :: Nt
  double precision,intent(INOUT), dimension(SPACEDIM,NumTratti) :: Xmin, Xmax
!
!.. local Scalars ..
  integer(4) :: i, iv, nf
  character(len=lencard)  :: nomsub = "FindLine"
!
!.. executable statements
!
!.. loops on the segments of the line
!
  do iv = Tratto(Nt)%inivertex, Tratto(Nt)%inivertex + Tratto(Nt)%numvertices - 1
!
!.. checks for the maximum boundary lines storage
!
    if (iv > NumBVertices) then
      call diagnostic (arg1=10,arg2=2,arg3=nomsub)
    end if
!
!.. evaluates the minimum and maximum coordinates with respect the current vertex nf in the line nt
!
    nf = BoundaryVertex(iv)
!
    do i = 1, SPACEDIM
      if ( Vertice(i,nf) < Xmin(i,Nt) ) Xmin(i,Nt) = Vertice(i,nf)
      if ( Vertice(i,nf) > Xmax(i,Nt) ) Xmax(i,Nt) = Vertice(i,nf)
    end do
  end do
!
  return
  end subroutine FindLine
!---split

!cfile GeneratePart.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : GeneratePart
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!AA504
! 03  Amicarelli        08Apr14        (v5.04) Modifications for IC (reservoir + dam) from PV file 
!
!************************************************************************************
! Module purpose : Module for geometry definition and initial condition definition
!                  as from input file
!
! Calling routine: Gest_Input
!
! Called routines: FindLine
!                  FindFrame
!                  random_seed
!                  SetParticles
!AA504
!                  point_inout_polygone
!                  SetParticleParameters
!
!************************************************************************************
!
!AA504sub
  subroutine GeneratePart(IC_loop)
!
!.. assign modules
  use FILES_ENTITIES
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
!AA504  
  use DIAGNOSTIC_MODULE
  
!
!.. Implicit Declarations ..
  implicit none
!

!AA504  
!Dummy variables
  integer(4),intent(IN) :: IC_loop

!Local Scalars
!AA504 sub  
 integer(4)       :: Nt,Nz,Mate,IsopraS,NumParticles,i,j,k,dimensioni,NumPartPrima,aux_factor,i_vertex,j_vertex,test_xy,test_z,n_levels,nag_aux,i_face,j_node,npi,ier,test_face
 integer(4)       :: test_dam,test_xy_2
!AA504
 double precision :: distance_hor,z_min,aux_scal,rnd
 double precision :: aux_vec(3)
!
!.. local Arrays ..
  integer(4),      dimension(SPACEDIM)         :: Npps
  double precision,dimension(SPACEDIM)         :: MinOfMin,XminReset
  double precision,dimension(:,:), allocatable :: Xmin,Xmax
!AA504 start
  double precision,dimension(:),allocatable    :: z_aux
  type(TyParticle),dimension(:),allocatable    :: pg_aux
  character(len=lencard)  :: nomsub = "GeneratePart"  
!AA504 end
  
!
!.. Executable Statements ..
!
!.. initializations
!
  call random_seed()

!AA504 sub start 
  test_z = 0
  if (IC_loop == 2) then
      do Nz=1,NPartZone
         if (Partz(Nz)%IC_source_type==2) test_z = 1
      end do 
  endif
  if (test_z==0) nag = 0 
!AA504 sub end
  
  NumParticles = 0
  
  if (ncord == 2) then
    dimensioni = NumTratti
  else
    dimensioni = NumFacce
  end if
  allocate (xmin(spacedim,dimensioni),xmax(spacedim,dimensioni))

  MinOfMin = max_positive_number
!
!.. loops on all the zones in order to set the initial minimum and maximum subzones coordinates
!
  first_cycle: do Nz = 1, NPartZone
!
!.. searches for the maximum and minimum partzone coordinates
!
    Partz(Nz)%coordMM(1:Spacedim,1) =  max_positive_number
    Partz(Nz)%coordMM(1:Spacedim,2) =  max_negative_number
!
!.. loops on the different regions belonging to a same zone
!
    do Nt = Partz(Nz)%Indix(1), Partz(Nz)%Indix(2)
!
      Xmin(1:SPACEDIM,nt) =  max_positive_number
      Xmax(1:SPACEDIM,nt) =  max_negative_number
!
!.. searches for the minimum and maximum coordinates of the current zone in 2D dimensions
!   
      if (ncord == 2) then
!
        call FindLine (Xmin, Xmax, Nt)
!
!.. searches for the minimum and maximum coordinates of the current subzone in 3D dimensions
!   
      else
!
        call FindFrame (Xmin, Xmax, Nt)
!
      end if
!
!.. evaluates the minimum and maximum coordinates of the zone checking for all the subzones
!
      do i = 1, Spacedim
!
        if ( Xmin(i,nt) < Partz(Nz)%coordMM(i,1) ) Partz(Nz)%coordMM(i,1) = Xmin(i,nt)
        if ( Xmax(i,nt) > Partz(Nz)%coordMM(i,2) ) Partz(Nz)%coordMM(i,2) = Xmax(i,nt)
!
!.. evaluates the minimum of subzones minimum coordinates
!
        MinOfMin(i) = min(MinOfMin(i),Xmin(i,nt))
!
      end do
!
    end do
!
  end do first_cycle
!
!.. loops on the zone in order to set the particle locations
!
  second_cycle: do Nz = 1, NPartZone
!
!.. skips the boundaries not having types "perimeter", "pool"
!
    if ( Partz(Nz)%tipo /= "peri" .AND. Partz(Nz)%tipo /= "pool" ) cycle second_cycle
!
!.. "perimeter", "pool" conditions are set
!
    Mate = Partz(Nz)%Medium
    Partz(Nz)%limit(1) = NumParticles + 1
!
    
!AA504 start
    if ((Partz(Nz)%tipo == "peri").and.(Partz(Nz)%IC_source_type == 2)) then
!During the first IC loop
       if (IC_loop==1) then
!Update number of points/vertices of the zone
          Partz(Nz)%npoints = Partz(Nz)%ID_last_vertex - Partz(Nz)%ID_first_vertex + 1
          aux_factor = anint(Partz(Nz)%dx_CartTopog/Domain%dd) 
!Allocate the auxiliary arrays 
          allocate(z_aux(Partz(Nz)%npoints))
!Loops over Cartesian topography points
          z_min = max_positive_number
          do i_vertex=Partz(Nz)%ID_first_vertex,Partz(Nz)%ID_last_vertex
             z_min = min(z_min,(Vertice(3,i_vertex)))
          end do 
          nag_aux = Partz(Nz)%nag_aux
!test          nag_aux = Partz(Nz)%npoints * aux_factor * (int((Partz(Nz)%H_res-z_min)/Domain%dd)+1)
          allocate(pg_aux(nag_aux))
          NumParticles = NumParticles + 1
!Loops over Cartesian topography points
          do i_vertex=Partz(Nz)%ID_first_vertex,Partz(Nz)%ID_last_vertex
! Check if the vertex is inside the plan_reservoir     
             call point_inout_polygone(Vertice(1:2,i_vertex),Partz(Nz)%plan_reservoir_points,Partz(Nz)%plan_reservoir_pos(1,1:2), &
                                       Partz(Nz)%plan_reservoir_pos(2,1:2),Partz(Nz)%plan_reservoir_pos(3,1:2),Partz(Nz)%plan_reservoir_pos(4,1:2),test_xy)
             if (test_xy==1) then
!Set minimum topography height among the closest 9 points
                z_aux(i_vertex) = max_positive_number 
                do j_vertex=1,Partz(Nz)%npoints
                   distance_hor = dsqrt( (Vertice(1,i_vertex)-Vertice(1,j_vertex))**2 + (Vertice(2,i_vertex)-Vertice(2,j_vertex))**2)
                   if (distance_hor<=(1.5*Partz(Nz)%dx_CartTopog)) z_aux(i_vertex)=min(z_aux(i_vertex),(Vertice(3,j_vertex)))
                end do
!Generating daughter points of the Cartesian topography, according to dx
                n_levels=int((Partz(Nz)%H_res-z_aux(i_vertex))/Domain%dd) !maybe ceiling instead of int would be better
                do i=1,aux_factor
                   do j=1,aux_factor
                      do k=1,n_levels
!Setting coordinates                       
!We avoid particles to be located on the very diagonals of topography (otherwise some particles could be non-removed even if below the topography; 
!                                                                      it is not clear why this rare eventuality would happen, but no problem remains anymore)
                         pg_aux(NumParticles)%coord(1) = Vertice(1,i_vertex) - aux_factor/2.*Domain%dd  + (i-1+0.501)*Domain%dd
                         pg_aux(NumParticles)%coord(2) = Vertice(2,i_vertex) - aux_factor/2.*Domain%dd  + (j-1+0.5)*Domain%dd
                         pg_aux(NumParticles)%coord(3) = (Partz(Nz)%H_res-Domain%dd/2.)-(k-1)*Domain%dd
!Test if the particle is below the reservoir
!Loop over boundaries 
                         test_z = 0
                         test_face = 0
!$omp parallel do default(none) shared(NumFacce,Tratto,BoundaryFace,Partz,Nz,pg_aux,test_face,test_z,NumParticles) private(i_face,aux_vec,aux_scal,test_xy)                         
                         do i_face=1,NumFacce
                            if (test_face==1) cycle
                            if (Tratto(BoundaryFace(i_face)%stretch)%zone==Partz(Nz)%Car_top_zone) then  
!Test if the point lies inside the plan projection of the face     
                               call point_inout_polygone(pg_aux(NumParticles)%coord(1:2),BoundaryFace(i_face)%nodes,BoundaryFace(i_face)%Node(1)%GX(1:2), &
                                    BoundaryFace(i_face)%Node(2)%GX(1:2),BoundaryFace(i_face)%Node(3)%GX(1:2),BoundaryFace(i_face)%Node(4)%GX(1:2),test_xy)
                               if (test_xy==1) then
!No need for a critical section to update test_face: only one face/process is interested (no conflict); 
!                                                    in particular cases there could be a conflict, but no face has a preference and test_face cannot come back to zero (no problem)
                                  test_face = 1
!Test if the point is invisible (test_z=1, below the topografy) or visible (test_z=0) to the face 
                                  aux_vec(:) = pg_aux(NumParticles)%coord(:) - BoundaryFace(i_face)%Node(1)%GX(:)
                                  aux_scal = dot_product(BoundaryFace(i_face)%T(:,3),aux_vec)
!No need for a critical section to update test_z: in particular cases there could be a conflict, but the "if" condition would always provide the same result
                                  if (aux_scal<=0.) test_z=1
                               endif
                            endif
                         end do 
!$omp end parallel do 

!Check the presence of a dam zone
                         test_dam = 0
                         if (Partz(Nz)%dam_zone_ID>0) then
!Test if the point lies inside the plan projection of the dam zone
                            call point_inout_polygone(pg_aux(NumParticles)%coord(1:2),Partz(Nz)%dam_zone_n_vertices,Partz(Nz)%dam_zone_vertices(1,1:2), &
                                                      Partz(Nz)%dam_zone_vertices(2,1:2),Partz(Nz)%dam_zone_vertices(3,1:2),Partz(Nz)%dam_zone_vertices(4,1:2),test_xy)
                            if (test_xy==1) then
                               test_face = 0
!Loop on the faces of the dam in the dam zone
!$omp parallel do default(none) shared(NumFacce,test_dam,Tratto,BoundaryFace,Partz,Nz,pg_aux,NumParticles,test_face) private(i_face,aux_vec,aux_scal,test_xy_2)
                               do i_face=1,NumFacce
                                  if (test_face==1) cycle
!Check if the face belongs to the dam_zone boundary and in particular to its top                                  
                                  if ((Tratto(BoundaryFace(i_face)%stretch)%zone==Partz(Nz)%dam_zone_ID).and.(BoundaryFace(i_face)%T(3,3)<0.)) then 
!Test if the particle horizontal coordinates ly inside the horizontal projection of the face
                                     call point_inout_polygone(pg_aux(NumParticles)%coord(1:2),BoundaryFace(i_face)%nodes,BoundaryFace(i_face)%Node(1)%GX(1:2), &
                                             BoundaryFace(i_face)%Node(2)%GX(1:2),BoundaryFace(i_face)%Node(3)%GX(1:2),BoundaryFace(i_face)%Node(4)%GX(1:2),test_xy_2)
                                     if (test_xy_2==1) then
!Note: even if a particle simultanously belongs to 2 faces test_dam will provide the same results: no matter the face order
!$omp critical
                                        test_dam = 1  
                                        test_face = 1
!Test if the particle position is invisible (test_dam=1) or visible (test_dam=0) to the dam. Visibility to the dam means invisibility to the current dam top face. 
                                        aux_vec(:) = pg_aux(NumParticles)%coord(:) - BoundaryFace(i_face)%Node(1)%GX(:)
                                        aux_scal = dot_product(BoundaryFace(i_face)%T(:,3),aux_vec)
                                        if (aux_scal<0.) then 
                                            test_dam=0
                                        endif
!$omp end critical 
                                     endif   
                                  endif
                               end do
!$omp end parallel do                            
                            endif   
                         endif
!Update fluid particle counter
!Check if the fluid particle is visible both to the bottom reservoir and the eventual dam
                         if ((test_z==0).and.(test_dam==0)) then
                            NumParticles = NumParticles + 1 
!Check the storage for the reached number of fluid particles
                            if (NumParticles>nag_aux) call diagnostic (arg1=10,arg2=4,arg3=nomsub)  
                         endif
                      end do
                   end do
                end do
             endif
          end do 
!Allocate fluid particle array       
          NumParticles = NumParticles - 1
          PARTICLEBUFFER = NumParticles * Domain%COEFNMAXPARTI
          nag_reservoir_CartTopog = NumParticles
          allocate(pg(PARTICLEBUFFER),stat=ier)
          if (ier/=0) then
             write (nout,'(1x,a,i2)') "    Array PG not allocated. Error code: ",ier
             call diagnostic (arg1=4,arg3=nomsub)
             else
                 write (nout,'(1x,a)') "    Array PG successfully allocated "
                 pg(:) = PgZero
          end if
!Loop over the auxiliary particle array
!$omp parallel do default(none) shared(pg,pg_aux,NumParticles) private(npi)
          do npi=1,NumParticles
!Copy fluid particle array from the corresponding auxiliary array
             pg(npi) = pg_aux(npi)
          end do
!$omp end parallel do  
!Deallocate auxiliary arrays
          deallocate(z_aux)
          deallocate(pg_aux)
!Second IC loop
          else
!Loops over fluid particles 
!$omp parallel do default(none) shared(nag_reservoir_CartTopog,Nz,Mate,Domain,pg) private(npi,rnd)  
             do npi=1,nag_reservoir_CartTopog
!Set particle parameters
                call SetParticleParameters(npi,Nz,Mate)   
!Modify random coordinates 
                if (Domain%RandomPos == 'r') then
                   call random_number(rnd)
                   pg(npi)%coord(1) = pg(npi)%coord(1) + (2. * rnd - 1.) * 0.1 * Domain%dd
                   call random_number(rnd)
                   pg(npi)%coord(2) = pg(npi)%coord(2) + (2. * rnd - 1.) * 0.1 * Domain%dd
                   call random_number(rnd)
                   pg(npi)%coord(3) = pg(npi)%coord(3) + (2. * rnd - 1.) * 0.1 * Domain%dd
                end if
             end do
!$omp end parallel do  
          NumParticles = nag_reservoir_CartTopog
       endif
    else
!AA504 end

!.. loops on the different regions belonging to a same zone
!
    do Nt = Partz(Nz)%Indix(1), Partz(Nz)%Indix(2)
!
!.. evaluates the number of particles for the subdomain Npps in all the directions
!
      if ((Xmin(1,nt) ==  max_positive_number) .or. (Xmax(1,nt) ==  max_negative_number)) then
!.. the zone is declared but is not used
        Npps = -1
      else 
        do i = 1, spacedim
          NumPartPrima = Nint((Xmin(i,nt) - MinOfMin(i))/Domain%dd)
          Xminreset(i) = MinOfMIn(i) + NumPartPrima * Domain%dd
          Npps(i)   = NInt((Xmax(i,nt) - XminReset(i)) / Domain%dd)  
        end do
      end if
!
!.. force almost one particle if the domain width in a direction is less than the dd particle dimension
!
      if (ncord == 2) where ( Npps == 0 ) Npps = 1
!
!.. considers the "pool" condition: set Isopra = 1 if it crosses the "virtual" side
!
      if ( Partz(Nz)%tipo == "pool" ) then
        Xmax(Partz(Nz)%ipool,Nt) = Partz(Nz)%pool
        IsopraS = 1   
      else
        IsopraS = 0
      end if     
!
!.. set the particles in the zone
!
      Call SetParticles ( Nt, Nz, Mate, XminReset, Npps, NumParticles, IsopraS )
!
    end do
    
!
!.. set the upper pointer for the particles in the nz-th zone
!
    Partz(Nz)%limit(2) = NumParticles
!

!AA504 
    endif

  end do second_cycle
!
  
!AA504 sub start 
  test_z = 0
  if (IC_loop == 2) then
      do Nz=1,NPartZone
         if (Partz(Nz)%IC_source_type==2) test_z = 1
      end do 
  endif
  if (test_z==0)   nag = NumParticles
!AA504 sub end

  nagpg = nag
!
  deallocate (xmin,xmax)
!
  return
  end subroutine GeneratePart
!---split

!cfile gest_dealloc.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : gest_dealloc
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Amicarelli/Agate  30Nov11        BSPH parameters
!
!AA501b
! 04  Amicarelli-Agate  23Jul12        Body dynamics
!AA601
! 05  Amicarelli        26Jan15        New memory managemnt for DBSPH 
!
!************************************************************************************
! Module purpose : Module for array deallocation
!
! Calling routine: Sphera, Gest_Input
!
! Called routines: diagnostic
!
!************************************************************************************
!
subroutine Gest_Dealloc ( nomsub )
!
!.. assign modules
use FILES_ENTITIES
use GLOBAL_MODULE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
!character(*),intent(IN)    :: nomsub
character(LEN=lencard), intent(IN) :: nomsub
!
!.. local scalars

!AA501b modified
integer(4)   :: ier,i

!
logical      :: check
!
!.. execution statements
!
 write (nout,'(1x,a)') ">> Storage deallocation in routine "//trim(nomsub)
 check = .true.
!
 if ( allocated(Vertice) ) then
    deallocate ( Vertice, stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: VERTICE not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: VERTICE successfully deallocated "
    end if
 end if

 if ( allocated(Tratto) ) then
    deallocate ( Tratto, stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: TRATTO not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: TRATTO successfully deallocated "
    end if
 end if

 if ( allocated(BoundaryFace) ) then
    deallocate ( BoundaryFace, stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: BoundaryFace not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: BoundaryFace successfully deallocated "
    end if
 end if

 if ( allocated(BFaceList) ) then 
    deallocate ( BFaceList, stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: BFACELIST not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: BFACELIST successfully deallocated "
    end if
 end if 

 if ( allocated(BoundaryVertex) ) then
    deallocate ( BoundaryVertex, stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: BOUNDARYVERTEX not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: BOUNDARYVERTEX successfully deallocated "
    end if
 end if

 if ( allocated(BoundarySide) ) then
    deallocate ( BoundarySide, stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: BOUNDARYSIDE not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: BOUNDARYSIDE successfully deallocated "
    end if
 end if

 if ( allocated(Partz) ) then
    deallocate ( Partz, stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: PARTZ not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: PARTZ successfully deallocated "
    end if
 end if

 if ( allocated(Med) ) then
    deallocate ( Med, stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: MED not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: MED successfully deallocated "
    end if
 end if

 if ( allocated(OpCount) ) then
    deallocate ( OpCount, stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: OpCount not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: OpCount successfully deallocated "
    end if
 end if

 if ( allocated(SpCount) ) then
    deallocate ( SpCount, stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: SpCount not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: SpCount successfully deallocated "
    end if
 end if

 if ( allocated(EpCount) ) then
    deallocate ( EpCount, stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: EpCount not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: EpCount successfully deallocated "
    end if
 end if

 if ( allocated(EpOrdGrid) ) then
    deallocate ( EpOrdGrid, stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: EpOrdGrid not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: EpOrdGrid successfully deallocated "
    end if
 end if

 if ( allocated(Control_Sections) ) then
    deallocate ( Control_Sections, stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: CONTROL_SECTION not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: CONTROL_SECTION successfully deallocated "
    end if
 end if

 if ( allocated(Control_Points) ) then
    deallocate ( Control_Points, stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: CONTROL_POINTS not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: CONTROL_POINTS successfully deallocated "
    end if
 end if

 if ( allocated(Control_Lines) ) then
    deallocate ( Control_Lines, stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: CONTROL_LINES not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: CONTROL_LINES successfully deallocated "
    end if
 end if
!
 if ( allocated(Pg) ) then
    deallocate ( Pg , stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: PG not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: PG successfully deallocated "
    end if
 end if
!
 if ( allocated(ts0_pg) ) then
    deallocate ( ts0_pg , stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: ts0_pg not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: ts0_pg successfully deallocated "
    end if
 end if
 !
 if ( allocated(Section_Points) ) then
    deallocate ( Section_Points, stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: SECTION_POINTS not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: SECTION_POINTS successfully deallocated "
    end if
 end if
!
 if ( allocated(NPartOrd) ) then
    deallocate ( NPartOrd , stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: NPARTORD not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: NPARTORD successfully deallocated "
    end if
 end if
!
 if ( allocated(Icont) ) then
    deallocate ( Icont , stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: ICONT not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: ICONT successfully deallocated "
    end if
 end if
!
 if ( allocated(nPartIntorno) ) then
    deallocate ( nPartIntorno , stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: nPartIntorno not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: nPartIntorno successfully deallocated "
    end if
 end if

 if ( allocated(PartIntorno) ) then
    deallocate ( PartIntorno , stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: PartIntorno not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: PartIntorno successfully deallocated "
    end if
 end if

 if ( allocated(PartKernel) ) then
    deallocate ( PartKernel , stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: PartKernel not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: PartKernel successfully deallocated "
    end if
 end if

 if ( allocated(rag) ) then
    deallocate ( rag , stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: rag not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: rag successfully deallocated "
    end if
 end if

 if ( allocated(BoundaryDataTab) ) then
    deallocate ( BoundaryDataTab , stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: BoundaryDataTab not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: BoundaryDataTab successfully deallocated "
    end if
 end if

 if ( allocated(BoundaryDataPointer) ) then
    deallocate ( BoundaryDataPointer , stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: BoundaryDataPointer not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: BoundaryDataPointer successfully deallocated "
    end if
 end if

 if ( allocated(Array_Flu) ) then
    deallocate ( Array_Flu , stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: Array_Flu not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: Array_Flu successfully deallocated "
    end if
 end if

! if ( allocated(Array_Sol) ) then
!    deallocate ( Array_Sol , stat = ier )
!    if (ier /= 0) then
!      write (nout,'(1x,a,i2)') "   Array: Array_Sol not deallocated with error code: ",ier
!      check = .false.
!    else
!      write (nout,'(1x,a)') "   Array: Array_Sol successfully deallocated "
!    end if
! end if
!
!.. final control
!
!AA406 start
 if (Domain%tipo == "bsph") then
 if ( allocated(NPartOrd_w) ) then
    deallocate ( NPartOrd_w , stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: NPARTORD_w not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: NPARTORD_w successfully deallocated "
    end if
 end if
!
 if ( allocated(Icont_w) ) then
    deallocate ( Icont_w , stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: ICONT_w not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: ICONT_w successfully deallocated "
    end if
 end if
!
 if ( allocated(nPartIntorno_fw) ) then
    deallocate ( nPartIntorno_fw , stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: nPartIntorno_fw not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: nPartIntorno_fw successfully deallocated "
    end if
 end if

 if ( allocated(PartIntorno_fw) ) then
    deallocate ( PartIntorno_fw , stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: PartIntorno_fw not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: PartIntorno_fw successfully deallocated "
    end if
 end if

 if ( allocated(kernel_fw) ) then
    deallocate ( kernel_fw , stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: kernel_fw not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: kernel_fw successfully deallocated "
    end if
 end if

 if ( allocated(rag_fw) ) then
    deallocate ( rag_fw , stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: rag_fw not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: rag_fw successfully deallocated "
    end if
 end if
 endif
!AA406 end
!

!AA501b start
 do i=1,n_bodies
 if (allocated(body_arr(i)%body_kinematics)) then
    deallocate (body_arr(i)%body_kinematics,stat =ier)
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: body_arr(i)%body_kinematics not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: body_arr(i)%body_kinematics successfully deallocated "
    end if
 end if
 end do
 if (allocated(body_arr)) then
    deallocate (body_arr,stat =ier)
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: body_arr not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: body_arr successfully deallocated "
    end if
 end if
 if (allocated(bp_arr)) then
    deallocate (bp_arr,stat =ier)
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: bp_arr not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: bp_arr successfully deallocated "
    end if
 end if
 if (allocated(Icont_bp)) then
    deallocate (Icont_bp,stat =ier)
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: Icont_bp not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: Icont_bp successfully deallocated "
    end if
 end if
 if (allocated(NPartOrd_bp)) then
    deallocate (NPartOrd_bp,stat =ier)
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: NPartOrd_bp not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: NPartOrd_bp successfully deallocated "
    end if
 end if
 if (allocated(nPartIntorno_bp_f)) then
    deallocate (nPartIntorno_bp_f,stat =ier)
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: nPartIntorno_bp_f not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: nPartIntorno_bp_f successfully deallocated "
    end if
 end if
 if (allocated(PartIntorno_bp_f)) then
    deallocate (PartIntorno_bp_f,stat =ier)
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: PartIntorno_bp_f not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: PartIntorno_bp_f successfully deallocated "
    end if
 end if
 if (allocated(KerDer_bp_f_cub_spl)) then
    deallocate (KerDer_bp_f_cub_spl,stat =ier)
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: KerDer_bp_f_cub_spl not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: KerDer_bp_f_cub_spl successfully deallocated "
    end if
 end if
 if (allocated(KerDer_bp_f_Gal)) then
    deallocate (KerDer_bp_f_Gal,stat =ier)
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: KerDer_bp_f_Gal not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: KerDer_bp_f_Gal successfully deallocated "
    end if
 end if 
 if (allocated(rag_bp_f)) then
    deallocate (rag_bp_f,stat =ier)
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: rag_bp_f not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: rag_bp_f successfully deallocated "
    end if
 end if
 if (allocated(surf_body_part)) then
    deallocate (surf_body_part,stat =ier)
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: surf_body_part not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: surf_body_part successfully deallocated "
    end if
 end if 
 if (allocated(nPartIntorno_bp_bp)) then
    deallocate (nPartIntorno_bp_bp,stat =ier)
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: nPartIntorno_bp_bp not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: nPartIntorno_bp_bp successfully deallocated "
    end if
 end if
 if (allocated(PartIntorno_bp_bp)) then
    deallocate (PartIntorno_bp_bp,stat =ier)
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: PartIntorno_bp_bp not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: PartIntorno_bp_bp successfully deallocated "
    end if
 end if
 if (allocated(rag_bp_bp)) then
    deallocate (rag_bp_bp,stat =ier)
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: rag_bp_bp not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: rag_bp_bp successfully deallocated "
    end if
 end if
 if (allocated(impact_vel)) then
    deallocate (impact_vel,stat =ier)
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: impact_vel not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: impact_vel successfully deallocated "
    end if
 end if
!AA501b end

!AA504 start
 if (allocated(BoundaryConvexEdge)) then
    deallocate (BoundaryConvexEdge,stat =ier)
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: BoundaryConvexEdge not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: BoundaryConvexEdge successfully deallocated "
    end if
 end if
 if (allocated(GCBFVector)) then
    deallocate (GCBFVector,stat =ier)
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "   Array: GCBFVector not deallocated with error code: ",ier
      check = .false.
    else
      write (nout,'(1x,a)') "   Array: GCBFVector successfully deallocated "
    end if
 end if 
 if(allocated(Q_sections%section)) then
   do i=1,Q_sections%n_sect
      if (allocated(Q_sections%section(i)%flow_rate)) deallocate(Q_sections%section(i)%flow_rate) 
   end do    
   deallocate(Q_sections%section)
 endif 
 if(allocated(Granular_flows_options%lines)) then
   deallocate(Granular_flows_options%lines)
 endif
!AA504 end 
!AA601 start
 if (allocated(DBSPH%kinematics)) then
    deallocate(DBSPH%kinematics,STAT=ier)
    if (ier/=0) then
       write(nout,*) 'Deallocation of aux_array in GestDealloc failed; the program terminates here.'
       call diagnostic (arg1=5,arg2=340)       
       else
          write (nout,'(1x,a)') "Deallocation of aux_array in GestDealloc successfully completed."
    endif
 endif
!AA601 start
 if (allocated(DBSPH%inlet_sections)) then
    deallocate(DBSPH%inlet_sections,STAT=ier)
    if (ier/=0) then
       write(nout,*) 'Deallocation of DBSPH%inlet_sections in GestDealloc failed; the program terminates here.'
       call diagnostic (arg1=5,arg2=340)       
       else
          write (nout,'(1x,a)') "Deallocation of DBSPH%inlet_sections in GestDealloc successfully completed."
    endif
 endif
 if (allocated(DBSPH%outlet_sections)) then
    deallocate(DBSPH%outlet_sections,STAT=ier)
    if (ier/=0) then
       write(nout,*) 'Deallocation of DBSPH%outlet_sections in GestDealloc failed; the program terminates here.'
       call diagnostic (arg1=5,arg2=340)       
       else
          write (nout,'(1x,a)') "Deallocation of DBSPH%outlet_sections in GestDealloc successfully completed."
    endif
 endif
!AA601 end

  if (check) return
  call diagnostic(arg1=3,arg3=nomsub)
!
return
end subroutine Gest_dealloc
!---split

!cfile gest_input.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name      : gest_input
!
! Last updating : Nov 14, 2012
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Amicarelli/Agate  30Nov11        BSPH approach
! 04  Amicarelli-Agate  13nov12        (AA501b) Body dynamics 
!AA504
! 05  Amicarelli        08Apr14        (v5.04) Several modifications (mainly for granular flows and 3D erosion criterion)
!AA601
! 06  Amicarelli        26Jan15        DBSPH-input (AA601). New DBSPH input management.
!
!************************************************************************************
! Module purpose : Module for input check and management
!
! Calling routine: sphera
!
! Called routines: diagnostic
!                  CompleteBoundaries3D
!                  CreaGrid
!                  DefineBoundaryFaceGeometry3D
!                  DefineBoundarySideGeometry2D
!                  FindBoundaryConvexEdges3D
!                  GeneratePart
!                  Gest_Dealloc
!                  Init_Arrays
!                  ModifyFaces
!                  OrdGrid1
!                  ReadInput
!                  ReadRestartFile
!                  SubCalcPreIdro
!                  Input_Body_Dynamics (AA501b)
!AA601 start
!                  Import_ply_surface_meshes 
!                  DBSPH_IC_surface_elements
!                  DBSPH_find_close_faces
!                  semi_particle_volumes
!AA601 end
!
!************************************************************************************
!
subroutine Gest_Input
!
!.. assign modules
use FILES_ENTITIES
use GLOBAL_MODULE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
integer(4), parameter :: ner0 = 0
!
!AA406!!!test
!real :: dx_w,xmin,xmax,ymin,ymax,zmin,zmax,xmin_obs,xmax_obs,ymin_obs,ymax_obs,zmin_obs,zmax_obs
!integer(4) :: nx,ny,nz,nx_obs,ny_obs,nz_obs   
!
!.. Local Scalars ..
!
!AA406sub
integer(4) :: npi, ier,i,n,isi,nfc,nt !,j
!
!AA504sub
integer(4) :: nrecords,IC_loop
integer(4) :: InputErr
character(len=lencard) :: nomsub = "GEST_INPUT"
character(len=lencard) :: ainp,msg_err
logical :: OnlyTriangle != .FALSE. ! .TRUE. !Opzione trasformazione facce 4 in 2 da 3 lati
!AA504
integer(4) :: machine_Julian_day,machine_hour,machine_minute,machine_second 
!
!.. Local Arrays ..
integer(4), dimension(20)  :: NumberEntities 
!
! External functions and subrotuines
character(10), external :: ltrim
character(80), external :: lcase
!
!AA406 test
!integer(4) :: ParticleCellNumber
!
!.. executable statements
!
  write (nout,'(1x,a)') ">> Input data management starts... "
  write (nscr,'(1x,a)') ">> Input data management starts... "
!
!.. deallocation of arrays
!
  call Gest_Dealloc ( nomsub )
!
!.. initializations
!
  ncord = 0                         ! spatial coordinates
  NumberEntities  = 0
  Domain%istart   = 0               ! initialization for first execution
  Domain%start    = zero            ! initialization for first execution
  Domain%file     = " "             ! initialization for first execution 
  Domain%NormFix  = .FALSE.         !Calcolo SI/NO normali particelle fisse con moto
  Domain%Slip     = .FALSE.
  OnlyTriangle    = .FALSE.
  diffusione      = .FALSE.
  esplosione      = .FALSE.
  erosione        = .FALSE.
  Restart         = .FALSE.
!
!.. time loop parameters are initialized
!
  tempo    = zero
  dt       = zero
  it_start = 0
!
  
!AA504 start
 if (exetype == "linux") then
    call system("date +%j%H%M%S > date_0.txt")
    open (unit_time_elapsed,file='date_0.txt',status="unknown",form="formatted")
    read (unit_time_elapsed,'(i3,i2,i2,i2)') machine_Julian_day,machine_hour,machine_minute,machine_second
    close (unit_time_elapsed)
    Domain%t0 = machine_Julian_day*24*60*60+machine_hour*60*60+machine_minute*60+machine_second
 endif
!AA504 end  
  
!.. allocations of temporary arrays
!
  write (nout,'(1x,a)') ">> Temporary storage allocation in routine "//trim(nomsub)
  allocate ( Vertice(SPACEDIM,1), BoundaryFace(1), Tratto(1), BoundaryVertex(1), stat = ier )
  if (ier /= 0) then
    write (nout,'(1x,a,i2)') "    Arrays VERTICE, BoundaryFace, TRATTO, BOUNDARYVERTEX not allocated. Error code: ",ier
    call diagnostic (arg1=4,arg3=nomsub)
  else
    write (nout,'(1x,a)') "    Arrays VERTICE, BoundaryFace, TRATTO, BOUNDARYVERTEX successfully allocated "
  end if
!
  allocate ( Partz(1), Med(1), Control_Sections(0:1), Control_Points(1), Control_Lines(1), stat = ier )
  if (ier /= 0) then
    write (nout,'(1x,a,i2)') &
    "    Arrays PARTZ, MED, Control_Sections, Control_Points, Control_Lines not allocated. Error code: ",ier
    call diagnostic (arg1=4,arg3=nomsub)
  else
    write (nout,'(1x,a)') "    Arrays PARTZ, MED, Control_Sections, Control_Points, Control_Lines successfully allocated "
  end if
!
!.. input data reading
!.. read for the first time the input file in order to have the parametrs for the arrays dimension
!
  NumVertici   = 1
  NumFacce     = 1
  NumTratti    = 1
  NumBVertices = 1
  NPartZone    = 1
  NMedium      = 1
  NSections    = 1
  NPoints      = 1
  NLines       = 1
!
  call ReadInput( NumberEntities,OnlyTriangle,InputErr,ainp )
!.. an error was detected in the input data deck. Execution fails
!
  msg_err = trim("dimensioning")
  if ( InputErr /= 0 ) then
    InputErr = InputErr + 300
    call diagnostic (arg1=5,arg2=InputErr,arg3=msg_err)
  end if
!
  write (nout,'(1x,a)') " > Data are read from an ASCII input file in the routine ReadInput"
  write (nscr,'(1x,a)') " > Data are read from an ASCII input file in the routine ReadInput"
!
!.. deallocations of temporary arrays
!
  write (nout,'(1x,a)') ">> Deallocation of temporary arrays "
  deallocate ( Vertice, BoundaryFace, Tratto, BoundaryVertex, stat = ier )
  if (ier /= 0) then
    write (nout,'(1x,a,i2)') "    Arrays VERTICE, BoundaryFace, TRATTO, BOUNDARYVERTEX not deallocated. Error code: ",ier
    call diagnostic (arg1=4,arg3=nomsub)
  else
    write (nout,'(1x,a)') "    Arrays VERTICE, BoundaryFace, TRATTO, BOUNDARYVERTEX successfully deallocated "
  end if
!
  deallocate ( Partz , Med , Control_Sections, Control_Points , Control_Lines , stat = ier )
  if (ier /= 0) then
    write (nout,'(1x,a,i2)') &
    "    Arrays PARTZ, MED, Control_Sections, Control_Points, Control_Lines not deallocated. Error code: ",ier
    call diagnostic (arg1=4,arg3=nomsub)
  else
    write (nout,'(1x,a)') "    Arrays PARTZ, MED, Control_Sections, Control_Points, Control_Lines successfully deallocated "
  end if
!
!.. --------------------------
!
!.. --------------------------
!.. a restart procedure has been invoked: restart positioning (step number / step time)
!.. --------------------------
!
  if ( Domain%istart > 0 .or. Domain%start > zero ) then
!
    Restart = .True.
!
!.. open the restart file from which restart data will be restored
!
    open ( unit=nsav, file=trim(Domain%file), form="unformatted", status="old", iostat=ier)
    if ( ier /= 0 ) then
      ainp = Domain%file
      call diagnostic (arg1=5,arg2=201,arg3=trim(ainp))
    else
      write (nout,'(1x,a)') " > Data are read from the restart file"//trim(Domain%file)//" in the routine ReadRestartFile"
      write (nscr,'(1x,a)') " > Data are read from the restart file in the routine ReadRestartFile"
    end if
!
!.. restore the data from the restart file
!
    call ReadRestartFile ( trim("heading"), ier, nrecords)

    msg_err = trim("heading")
    if ( ier /= 0 ) then
      ier = ier + 200
      call diagnostic (arg1=5,arg2=ier,arg3=msg_err)
    end if
!
    NPoints        = NumberEntities(4)
    NPointsl       = NumberEntities(6)
    NPointst       = NumberEntities(4) + NumberEntities(6) + NumberEntities(13)
    NLines         = NumberEntities(5)
    NSections      = NumberEntities(12)
    NPointse       = NumberEntities(13)
!
! .. two fluids case and diffusion are considered
!
!    do i = 1, NMedium
!      if (Med(i)%codif /= zero) diffusione = .TRUE.
!      if (Med(i)%Gamma /= zero) esplosione = .TRUE.
!      if ((index(Med(i)%tipo,"granular") > 0)) then
!        erosione = .TRUE.
!        modelloerosione = Med(i)%modelloerosione
!      end if
!    end do
!
!.. -----------------------------------
!.. no restart. Standard initialization
!.. -----------------------------------
!
  else
!
    ncord          = NumberEntities(1)
    nmedium        = NumberEntities(2)
    NPartZone      = NumberEntities(3)
    NPoints        = NumberEntities(4)
    NPointsl       = NumberEntities(6)
    NPointst       = NumberEntities(4) + NumberEntities(6) + NumberEntities(13)
    NLines         = NumberEntities(5)
    NumVertici     = NumberEntities(7)
    NumTratti      = NumberEntities(8)
    NumBVertices   = NumberEntities(9)
    NumBSides      = NumberEntities(10)
    NumFacce       = NumberEntities(11)
    if ( OnlyTriangle ) NumFacce = NumFacce + NumberEntities(18)
    NSections      = NumberEntities(12)
    NPointse       = NumberEntities(13)
    if ( NumberEntities(19) == 1 )  Domain%Slip     = .TRUE.
    if ( NumberEntities(20) == 1 )  Domain%NormFix = .TRUE.
!
  end if
!
!.. --------------------------
!
!.. allocate the arrays for the calculation depending on the quantities read or restarted
!
  write (nout,'(1x,a)') ">> Final storage allocation in routine "//trim(nomsub)
  allocate ( Vertice       (SPACEDIM,max(1,NumVertici)) , &
             BoundaryFace  (max(1,NumFacce)), &
             Tratto        (max(1,NumTratti)), &
             BoundaryVertex(max(1,NumBVertices)), & 
             BoundarySide  (max(1,NumBSides)), &
             stat = ier )
  if (ier /= 0) then
    write (nout,'(1x,a,i2)') &
    "    Arrays VERTICE, BoundaryFace, TRATTO, BOUNDARYVERTEX, BOUNDARYSIDE not allocated. Error code: ",ier
    call diagnostic (arg1=4,arg3=nomsub)
  else
    write (nout,'(1x,a)') "    Arrays VERTICE, BoundaryFace, TRATTO, BOUNDARYVERTEX, BOUNDARYSIDE successfully allocated "
  end if
!
  allocate ( Partz(NPartZone), &
             Med(NMedium), OpCount(NMedium), SpCount(NMedium), EpCount(NMedium), EpOrdGrid(NMedium), &
             Control_Sections(0:NSections+1), &
             Control_Points(max(1,NPointst)), &
             Control_Lines(max(1,NLines)), &
             Section_Points(1), &
             stat = ier )
  if (ier /= 0) then
    write (nout,'(1x,a,i2)') &
    "    Arrays PARTZ, MED, Control_Sections, Control_Points, Control_Lines, Section_Points not allocated. Error code: ",ier
    call diagnostic (arg1=4,arg3=nomsub)
  else
    write (nout,'(1x,a)') &
    "    Arrays PARTZ, MED, Control_Sections, Control_Points, Control_Lines, Section_Points successfully allocated "
  end if
!
  rewind (ninp)
!
!NumberEntities = 0
!
!.. array initializations
!
  call Init_Arrays
!
!.. ----------------------------
!.. no restart. Data acquisition
!.. ----------------------------
!
  if ( .not. Restart ) then
!
    call ReadInput( NumberEntities,OnlyTriangle,InputErr,ainp )
!
    msg_err = trim("readinput")
    if ( InputErr /= 0 ) then
      InputErr = InputErr + 300
      call diagnostic (arg1=5,arg2=InputErr,arg3=msg_err)
    end if
!
! .. two fluids case and diffusion are considered
!
    do i = 1, NMedium
      if (Med(i)%codif /= zero) diffusione = .TRUE.
      if (Med(i)%Gamma /= zero) esplosione = .TRUE.
      if ((index(Med(i)%tipo,"granular") > 0)) then
        erosione = .TRUE.
        modelloerosione = Med(i)%modelloerosione
      end if
    end do
!
    close (ninp)
!
!.. set the domain parameters
!
    nag    = 0
!
!.. 2D layout
!
    if (ncord == 2)then
!    
!AA406 sub
      Domain%coefke    = 0.682093d0 / squareh  ! 10 / (7 * pigreco) *(3./2.)
!
      Domain%coefkacl = 0.099472d0 / squareh   ! 5 / (16 * pigreco)
! calcolo volume particella
      Domain%PVolume = Domain%dd * Domain%dd
!
!.. 3D layout
!
    else if (ncord == 3)then
      Domain%coefke    = 0.477464d0 / cubich   ! 1 / pigreco
      Domain%coefkacl = 0.074603d0 / cubich    ! 15 / (64 * pigreco)
! calcolo volume particella
      Domain%PVolume = Domain%dd * Domain%dd * Domain%dd
    end if
!
    Control_Sections(NSections+1)%XYZRange(1:3,1) = Domain%coord(1:3,1)
    Control_Sections(NSections+1)%XYZRange(1:3,2) = Domain%coord(1:3,2)
!
!.. an irregular domain is considered (standard or semianalytical)
!
!AA406 sub
    if ( (Domain%tipo == "semi") .or. (Domain%tipo == "bsph") ) then
!
!.. 2D parameter setting
!
      if ( ncord == 2 ) then
!
!.. define the boundary structure
!
        call DefineBoundarySideGeometry2D
!
!.. 3D parameter setting
!
      else if ( ncord == 3 ) then
!
!.. modifies four-sided structures into three-sided structures
!
        if ( OnlyTriangle ) call ModifyFaces ( NumberEntities )
!
        allocate ( BFaceList(NumFacce), stat = ier )
        if (ier /= 0) then
          write (nout,'(1x,a,i2)') "    Array BFACELIST not allocated. Error code: ",ier
          call diagnostic (arg1=4,arg3=nomsub)
        else
          write (nout,'(1x,a)') "    Array BFACELIST successfully allocated "
        end if
!
!.. define the boundary structure
!
        call CompleteBoundaries3D
!
        call DefineBoundaryFaceGeometry3D
!
        
!AA504 start
        allocate (BoundaryConvexEdge(1:Domain%MAXNUMCONVEXEDGES), stat = ier)
        if (ier /= 0) then
           write (nout,'(1x,a,i2)') "   Array BoundaryConvexEdge not allocated. Error code: ",ier
           call diagnostic (arg1=4,arg3=nomsub)
           else
              write (nout,'(1x,a)') "   Array BoundaryConvexEdge successfully allocated "
        end if
!AA504 end

        call FindBoundaryConvexEdges3D
!
      end if
!
!.. creation of the field particles in order to allocate the storage
!
      nagpg = 0

!AA504 start
      if (Granular_flows_options%ID_erosion_criterion > 0) then
          Med(Granular_flows_options%ID_granular)%den0_s = Med(Granular_flows_options%ID_granular)%den0
          Med(Granular_flows_options%ID_granular)%den0 = (1.d0-Med(Granular_flows_options%ID_granular)%gran_vol_frac_max) * &
          Med(Granular_flows_options%ID_main_fluid)%den0 + Med(Granular_flows_options%ID_granular)%gran_vol_frac_max * Med(Granular_flows_options%ID_granular)%den0_s 
          Med(Granular_flows_options%ID_granular)%eps = Med(Granular_flows_options%ID_granular)%eps * &
                                                            Med(Granular_flows_options%ID_granular)%den0/Med(Granular_flows_options%ID_granular)%den0_s
      endif
!AA504 end      
      
!AA504sub
      IC_loop = 1
      call GeneratePart(IC_loop)

      
!AA504
      if (.not.(allocated(pg))) then      

!.. evaluates the number of field particles. Total number of particles is allocated depending on the nag value
!
      if (nag < 100) then
!.. initial domain empty and with source

!AA504 sub
        PARTICLEBUFFER = INIPARTICLEBUFFER * Domain%COEFNMAXPARTI

      else
          
!AA504sub          
        PARTICLEBUFFER = nag * Domain%COEFNMAXPARTI

      end if
!
     
!AA406 sub
      if ( ((Domain%tipo == "semi") .or. (Domain%tipo == "bsph")) ) then
!
        allocate ( pg(PARTICLEBUFFER), stat = ier )
      else
        call diagnostic (arg1=10,arg2=5,arg3=nomsub)
      end if   
      if (ier /= 0) then
        write (nout,'(1x,a,i2)') "    Array PG not allocated. Error code: ",ier
        call diagnostic (arg1=4,arg3=nomsub)
      else
        write (nout,'(1x,a)') "    Array PG successfully allocated "
        Pg(:) = PgZero
      end if
      
!AA504 sub
      endif
      
!
!AA402 start
      if ( Domain%RKscheme > 1 ) then 
        if ( Domain%tipo == "semi" ) then 
          allocate ( ts0_pg(PARTICLEBUFFER), stat = ier )
        else
          call diagnostic (arg1=10,arg2=5,arg3=nomsub)
        end if   
        if (ier /= 0) then
          write (nout,'(1x,a,i2)') "    Array ts0_pg not allocated. Error code: ",ier
          call diagnostic (arg1=4,arg3=nomsub)
        else
          write (nout,'(1x,a)') "    Array ts0_pg successfully allocated "
          ts0_pg(:) = ts_pgZero
        end if
      end if
!AA402 end
!
!.. virtual spatial grid is generated on all the domain
!
      call CreaGrid
!
!.. initial areas are discretized and the particles are created and initialized
!

!AA504sub
      IC_loop = 2
      call GeneratePart(IC_loop)

!AA406test
!
!AA501 sub
!AA601 rm
    else
      call diagnostic (arg1=10,arg2=5,arg3=nomsub)
    end if
!
!.. --------------------------
!.. a restart option is active
!.. --------------------------
!
  else
!
    if (nag < 100) then
!.. initial domain empty and with source
      PARTICLEBUFFER = INIPARTICLEBUFFER * Domain%COEFNMAXPARTI
    else
      PARTICLEBUFFER = nag * Domain%COEFNMAXPARTI
    end if
!
    if ( Domain%tipo == "semi" ) then   
      allocate ( pg(PARTICLEBUFFER), stat = ier )  
    else
      call diagnostic (arg1=10,arg2=5,arg3=nomsub)
    end if   
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "    Array PG not allocated. Error code: ",ier
      call diagnostic (arg1=4,arg3=nomsub)
    else
      write (nout,'(1x,a)') "    Array PG successfully allocated "
      Pg(:) = PgZero
    end if
!
!AA402 start
    if (Domain%RKscheme > 1) then
      if ( Domain%tipo == "semi" ) then   
        allocate ( ts0_pg(PARTICLEBUFFER), stat = ier )  
      else
        call diagnostic (arg1=10,arg2=5,arg3=nomsub)
      end if   
      if (ier /= 0) then
        write (nout,'(1x,a,i2)') "    Array ts0_pg not allocated. Error code: ",ier
        call diagnostic (arg1=4,arg3=nomsub)
      else
        write (nout,'(1x,a)') "    Array ts0_pg successfully allocated "
        ts0_pg(:) = ts_pgZero
      end if
    end if
!AA402 end
!
    allocate ( BFaceList(NumFacce), stat = ier )
    if (ier /= 0) then
      write (nout,'(1x,a,i2)') "    Array BFACELIST not allocated. Error code: ",ier
      call diagnostic (arg1=4,arg3=nomsub)
    else
      write (nout,'(1x,a)') "    Array BFACELIST successfully allocated "
    end if
!
    call ReadRestartFile ( trim("reading"), ier, nrecords)
    msg_err = trim("reading")
    if ( ier /= 0 ) then
      ier = ier + 200
      call diagnostic (arg1=5,arg2=ier,arg3=msg_err)
    end if
!
    close (nsav)
!
    call ReadInput( NumberEntities,OnlyTriangle,InputErr,ainp )
!
    msg_err = trim("restart reading?")
    if ( InputErr /= 0 ) then
      InputErr = InputErr + 300
      call diagnostic (arg1=5,arg2=InputErr,arg3=msg_err)
    end if
!
! .. two fluids case and diffusion are considered
!
!$$$$$$$$$$$$$$$
!esplosione = .FALSE.
!do i=1,nag
!  pg(i)%pres = zero
!  pg(i)%IntEn = zero
!  pg(i)%EnVar = zero
!  pg(i)%dEdT = zero
!end do
!$$$$$$$$$$$$$$$
    do i = 1, NMedium
      if (Med(i)%codif /= zero) diffusione = .TRUE.
!      if (Med(i)%Gamma /= zero) esplosione = .TRUE.
      if ((index(Med(i)%tipo,"granular") > 0)) then
        erosione = .TRUE.
        modelloerosione = Med(i)%modelloerosione
      end if
    end do
!
!.. save current time for result_converter
!
    val_time = tempo
!
    close (ninp)
!
  end if
!
!.. --------------------------
!
!
!.. final prints
!
!  write (nout,*) 
!  write (nout,*) "Number of particles          NAG= ",nag
!
! scrittura su file di ouput delle particelle di campo
  if ( Domain%ioutopt < 0 ) then
    write (nout,*) 
    write (nout,*) "======== PARTICLES COORDINATES =========="
    do n = 1, NPartZone
      write (nout,*) 
      write (nout,"(a,i5,3x,a)") "ZONE",n,Partz(n)%label
      do npi = Partz(n)%limit(1), Partz(n)%limit(2)
        write (nout,"(i10,4f14.5)") npi, pg(npi)%coord, pg(npi)%tstop  !, pg(npi)%vel
      end do
    end do
  end if
!

!AA501b start
! Management of body dynamics input
 if (n_bodies > 0) then
    call Input_Body_Dynamics
 endif
!AA501b end

!allocazione della memoria per i vettori di ordinamento delle particelle
!
!AA406 sub
  if ( (Domain%tipo == "semi") .or. (Domain%tipo == "bsph") ) then
!
    allocate ( NPartOrd(PARTICLEBUFFER),Icont(grid%nmax+1), stat = ier ) 
  else
    call diagnostic (arg1=10,arg2=5,arg3=nomsub)
  end if    
!
  if (ier /= 0) then
    write (nout,'(1x,a,i2)') "    Array NPARTORD,ICONT not allocated. Error code: ",ier
    call diagnostic (arg1=4,arg3=nomsub)
  else
    write (nout,'(1x,a)') "    Array NPARTORD,ICONT successfully allocated "
    NPartOrd(:) = 0
    Icont(:) = 0
  end if
!AA601 rm
!AA501b start
  if (n_bodies > 0.) then
     allocate (NPartOrd_bp(n_body_part),Icont_bp(grid%nmax+1),stat = ier) 
     if (ier /= 0) then
       write (nout,'(1x,a,i2)') "    Arrays NPARTORD_bp and ICONT_bp not allocated. Error code: ",ier
       call diagnostic (arg1=4,arg3=nomsub)
       else
          write (nout,'(1x,a)') "    Arrays NPARTORD_bp and ICONT_bp successfully allocated "
          NPartOrd_bp(:) = 0
          Icont_bp(:) = 0
     end if
  endif
!AA501b end
!AA601 start
  if (Domain%tipo == "bsph") then
     call Import_ply_surface_meshes
     call DBSPH_IC_surface_elements
     if (.not.allocated(NPartOrd_w)) then
        allocate(NPartOrd_w(DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet),Icont_w(grid%nmax+1),stat=ier) 
        if (ier/=0) then
           write(nout,*) 'Error! Allocation of NPartOrd_w or Icont_w Gest_Input failed.'           
           call diagnostic (arg1=5,arg2=340)
           else
              write (nout,'(1x,a)') "Arrays NPARTORD_w and ICONT_w successfully allocated."
              NPartOrd_w(:) = 0
              Icont_w(:) = 0
        end if
     endif
  endif
!AA601 end
  call OrdGrid1 ( nout )
!AA601 start
  if (Domain%tipo == "bsph") then
     call DBSPH_find_close_faces 
     call semi_particle_volumes
  endif
!AA601 end
!.. initialize the pressure and density fields
!
  if (.not. Restart) call SubCalcPreIdro
!
!.. evaluates and stores the opened boundary sides for 2D calculation
!
  if ( ncord == 2 ) then     
!
!.. searches for the opened boundary sides
!
    NumOpenSides = 0
!
    do isi = 1, NumBSides
      if (BoundarySide(isi)%tipo == "leve" .OR. BoundarySide(isi)%tipo == "velo" .OR. BoundarySide(isi)%tipo == "flow" .OR. &
          BoundarySide(isi)%tipo == "crit" .OR. BoundarySide(isi)%tipo == "open") then  
        NumOpenSides = NumOpenSides + 1
        if (NumOpenSides > MAXOPENSIDES) call diagnostic (arg1=10,arg2=6,arg3=nomsub)
        OpenSide(NumOpenSides) = isi
      end if
    end do
!
!.. evaluates and stores the opened boundary sides for 3D calculation
!
  else
!
!.. searches for the boundary faces having opened conditions
!
    NumOpenFaces = 0
!
    do nfc = 1, NumFacce
      nt = BoundaryFace(nfc)%stretch
      if (Tratto(nt)%tipo == "leve" .OR. Tratto(nt)%tipo == "velo" .OR. Tratto(nt)%tipo == "flow" .OR. &
          Tratto(nt)%tipo == "crit" .OR. Tratto(nt)%tipo == "open") then
        NumOpenFaces = NumOpenFaces + 1
        if (NumOpenFaces > MAXOPENFACES ) call diagnostic (arg1=10,arg2=7,arg3=nomsub)
        OpenFace(NumOpenFaces) = nfc
      end if
    end do
!
!.. searches for the maximum zone index ("maxzone") among them identifying volumes(?) of particles in the domain
!.. the sources (water material) must have currently the maximum zone index because the particle array is partitioned 
!.. by zones, i.e. the particles belonging to the first zone is loaded at first, after that those of the second one and so on.
!.. Therefore, the source particles can be added only in the last zone (last array section) without restructuring all the 
!.. zone pointers in the array (Part(zone)%limit(1) is the first particle of the zone-th, Part(zone)%limit(2) is the last
!.. one)
!
!    call SearchforParticleZone_3D(maxzone)      
!
  end if
!
  OpCount = 0
  EpCount = 0
  EpOrdGrid = 0

  return
  end subroutine Gest_Input
!---split

!cfile gest_trans.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : gest_trans
!
! Last updating : May 08, 2012
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Amicarelli/Agate  30Nov11        BSPH approach
! 04  Amicarelli/Agate  13nov12        (AA501) bBody dynamics 
! 05  Amicarelli/Agate  18apr13        add maximum number of neighbours read in input
!AA504
! 06  Amicarelli        08Apr14        (v5.04) Modifications for maximum depth and specific flow rate
!
!************************************************************************************
! Module purpose : Module for loop management
!
! Calling routine: sphera
!
! Called routines: diagnostic
!                  ComputeBoundaryIntegralTab
!                  ComputeKernelTable
!                  EvaluateBER_TimeStep
!                  GridCellBoundaryFacesIntersections3D
!                  loop_irre_2D
!                  loop_irre_3D
!                  s_ctime
!                  start_and_stop
!
!************************************************************************************
!
subroutine Gest_Trans 
!
!.. assign modules
!AA503
use AdM_USER_TYPE
use FILES_ENTITIES
use GLOBAL_MODULE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
integer(4),parameter :: ner0 = 0
!
!.. Local Scalars ..
integer(4)             :: npi,NumCellmax,i,ier,k,kk,k1,k2,nlinee,nvalori
character(len=lencard) :: nomsub = "GEST_TRANS"
character(len=lencard) :: filename,stringa,prefix,filevtk
!AA504 sub
character(len=200)     :: cargo
!
!.. Local Arrays ..
double precision,dimension(:),allocatable :: verticecolore
!
!.. External functions and subroutines ..
!
character(80), external :: lcase
logical,       external :: check_files2
integer(4),    external :: stepdata
!
!.. Executable Statements ..
!
 write ( nout,"(/)" )
 write (nout,*) "Initial number of particles      NAG= ",nag
 write ( nout,"(/)" )

!! scrittura su file di ouput delle particelle di campo
 if ( Domain%ioutopt < 0 ) then
   write (nout,*) 
   write (nout,*) "======== PARTICLES COORDINATES =========="
   do npi = 1,nag
     write (nout,"(i10,4f14.5)") npi, pg(npi)%coord, pg(npi)%tstop  !, pg(npi)%vel
   end do

 end if
!
!.. initialization arrays for vtkconverter
 nblocchi = 0
 blocchi = 0
 block = -1
 Time_Block = zero

!.. allocation arrays to store close boundaries and integrals for the loop
!!!! MaxNcbs = int(COEFNMAXPARTJ * MAXCLOSEBOUNDSIDES * Nag)
!!!! MaxNcbf = int(COEFNMAXPARTJ * MAXCLOSEBOUNDFACES * Nag)
!! MaxNcbs = int(COEFNMAXPARTJ * MAXCLOSEBOUNDSIDES * PARTICLEBUFFER)
!! MaxNcbf = int(COEFNMAXPARTJ * MAXCLOSEBOUNDFACES * PARTICLEBUFFER)
 MaxNcbs = int(MAXCLOSEBOUNDSIDES * PARTICLEBUFFER)
!AA504sub 
 MaxNcbf = int(Domain%MAXCLOSEBOUNDFACES * PARTICLEBUFFER)

!AA406
 if ((Domain%tipo == "semi") .or. (Domain%tipo == "bsph"))  then
!
 allocate (BoundaryDataPointer(1:3,1:PARTICLEBUFFER), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays BoundaryDataPointer not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays BoundaryDataPointer successfully allocated "
 end if
!
 if (Ncord == 2) then
   write(nout,'(a,i15)') "     Max num of close boundary sides: MaxNcbs = ",MaxNcbs
   allocate (BoundaryDataTab(1:MaxNcbs), stat = ier)
   if (ier /= 0) then
     write (nout,'(1x,a,i2)') "   Arrays BoundaryDataTab not allocated. Error code: ",ier
     call diagnostic (arg1=4,arg3=nomsub)
   else
     write (nout,'(1x,a)') "   Arrays BoundaryDataTab successfully allocated "
   end if
 else
   write(nout,'(a,i15)') "     Max num of close boundary faces: MaxNcbf = ",MaxNcbf
   allocate (BoundaryDataTab(1:MaxNcbf), stat = ier)
   if (ier /= 0) then
     write (nout,'(1x,a,i2)') "   Arrays BoundaryDataTab not allocated. Error code: ",ier
     call diagnostic (arg1=4,arg3=nomsub)
   else
     write (nout,'(1x,a)') "   Arrays BoundaryDataTab successfully allocated "
   end if
 end if
!
!AA406
 endif
!
!.. allocation arrays to store values for the loop
 NMAXPARTJ = Domain%COEFNMAXPARTJ * (Domain%h * four / Domain%dd) ** Ncord
 write(nout,'(a,i15)') "     Max num particles surrounding the current particle: NMAXPARTJ = ",NMAXPARTJ
!
 allocate (Array_Flu(1:PARTICLEBUFFER), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays Array_Flu not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays Array_Flu successfully allocated "
 end if
!
! allocate (Array_Sol(1:PARTICLEBUFFER), stat = ier)
! if (ier /= 0) then
!   write (nout,'(1x,a,i2)') "   Arrays Array_Sol not allocated. Error code: ",ier
!   call diagnostic (arg1=4,arg3=nomsub)
! else
!   write (nout,'(1x,a)') "   Arrays Array_Sol successfully allocated "
! end if
!
 allocate (nPartIntorno(1:PARTICLEBUFFER), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays NPARTINTORNO not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays NPARTINTORNO successfully allocated "
 end if
!
 allocate (PartIntorno(1:NMAXPARTJ*PARTICLEBUFFER), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays PARTINTORNO not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays PARTINTORNO successfully allocated "
 end if
 allocate (PartKernel(1:4,1:NMAXPARTJ*PARTICLEBUFFER), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays PARTKERNEL not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays PARTKERNEL successfully allocated "
 end if
 allocate (rag(1:3,1:NMAXPARTJ*PARTICLEBUFFER), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays RAG not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays RAG successfully allocated "
 end if
!
!AA406 start
 if (Domain%tipo == "bsph") then
 allocate (nPartIntorno_fw(1:PARTICLEBUFFER), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays NPARTINTORNO_fw not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays NPARTINTORNO_fw successfully allocated "
 end if
!
 allocate (PartIntorno_fw(1:NMAXPARTJ*PARTICLEBUFFER), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays PARTINTORNO_fw not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays PARTINTORNO_fw successfully allocated "
 end if
!
!AA406test
!AA501 sub
 if (Domain%tipo == "bsph") allocate (kernel_fw(2,1:NMAXPARTJ*PARTICLEBUFFER), stat = ier)
!
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays kernel_fw not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays kernel_fw successfully allocated "
 end if
 allocate (rag_fw(1:3,1:NMAXPARTJ*PARTICLEBUFFER), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays RAG_fw not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays RAG_fw successfully allocated "
 end if
 endif
!AA406 end
!

!AA501b start
 if (n_bodies > 0) then
 allocate (nPartIntorno_bp_f(n_body_part), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays nPartIntorno_bp_f not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays nPartIntorno_bp_f successfully allocated "
 end if
 allocate (PartIntorno_bp_f(NMAXPARTJ*n_body_part), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays PartIntorno_bp_f not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays PartIntorno_bp_f successfully allocated "
 end if
 allocate (KerDer_bp_f_cub_spl(NMAXPARTJ*n_body_part), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays KerDer_bp_f_cub_spl not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays KerDer_bp_f_cub_spl successfully allocated "
 end if
 allocate (KerDer_bp_f_Gal(NMAXPARTJ*n_body_part), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays KerDer_bp_f_Gal not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays KerDer_bp_f_Gal successfully allocated "
 end if 
 allocate (rag_bp_f(3,NMAXPARTJ*n_body_part), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays rag_bp_f not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays rag_bp_f successfully allocated "
 end if
 allocate (nPartIntorno_bp_bp(n_surf_body_part), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays nPartIntorno_bp_bp not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays nPartIntorno_bp_bp successfully allocated "
 end if
 allocate (PartIntorno_bp_bp(n_surf_body_part*NMAXPARTJ), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays PartIntorno_bp_bp not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays PartIntorno_bp_bp successfully allocated "
 end if
 allocate (rag_bp_bp(3,n_surf_body_part*NMAXPARTJ), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Arrays rag_bp_bp not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Arrays rag_bp_bp successfully allocated "
 end if
 if (ncord==2) allocate (impact_vel(n_surf_body_part,(n_bodies+NumBSides)), stat = ier)
 if (ncord==3) allocate (impact_vel(n_surf_body_part,(n_bodies+NumFacce)), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   Array impact_vel not allocated. Error code: ",ier
   call diagnostic (arg1=4,arg3=nomsub)
 else
   write (nout,'(1x,a)') "   Array impact_vel successfully allocated "
 end if
 impact_vel = 0.
 endif
!AA501b end

 write (nout,'(1x,a)') "..."
 write (nout,'(a,i15)') " Max number of particles  : PARTICLEBUFFER = ",PARTICLEBUFFER
 write (nout,*) " Size # of elements in array pg                  : ",size(pg)
 write (nout,*) " Size # of elements in array BoundaryDataTab     : ",size(BoundaryDataTab)
 write (nout,*) " Size # of elements in array BoundaryDataPointer : ",size(BoundaryDataPointer)
 write (nout,*) " Size # of elements in array Array_Flu           : ",size(Array_Flu)
! write (nout,*) " Size # of elements in array Array_Sol           : ",size(Array_Sol)
 write (nout,*) " Size # of elements in array nPartIntorno        : ",size(nPartIntorno)
 write (nout,*) " Size # of elements in array PartIntorno         : ",size(PartIntorno)
 write (nout,*) " Size # of elements in array PartKernel          : ",size(PartKernel)
 write (nout,*) " Size # of elements in array rag                 : ",size(rag)
!
!AA406 start
 if ((Domain%tipo == "bsph").and.(DBSPH%n_w > 0)) then
 write (nout,*) " Size # of elements in array pg_w                : ",size(pg_w)
 write (nout,*) " Size # of elements in array nPartIntorno_fw     : ",size(nPartIntorno_fw)
 write (nout,*) " Size # of elements in array PartIntorno_fw      : ",size(PartIntorno_fw)
 write (nout,*) " Size # of elements in array kernel_fw           : ",size(kernel_fw)
 write (nout,*) " Size # of elements in array rag_fw              : ",size(rag_fw)
 endif
!AA406 end

!AA501b start
 if (n_bodies > 0) then
 write (nout,*) " Size # of elements in array nPartIntorno_bp_f   : ",size(nPartIntorno_bp_f)
 write (nout,*) " Size # of elements in array PartIntorno_bp_f    : ",size(PartIntorno_bp_f)
 write (nout,*) " Size # of elements in array KerDer_bp_f_cub_spl : ",size(KerDer_bp_f_cub_spl)
 write (nout,*) " Size # of elements in array KerDer_bp_f_Gal     : ",size(KerDer_bp_f_Gal) 
 write (nout,*) " Size # of elements in array rag_bp_f            : ",size(rag_bp_f)
 write (nout,*) " Size # of elements in array surf_body_part      : ",size(surf_body_part) 
 write (nout,*) " Size # of elements in array nPartIntorno_bp_bp  : ",size(nPartIntorno_bp_bp)
 write (nout,*) " Size # of elements in array PartIntorno_bp_bp   : ",size(PartIntorno_bp_bp)
 write (nout,*) " Size # of elements in array rag_bp_bp           : ",size(rag_bp_bp)
 endif
!AA501b end

!
 write (nout,'(1x,a)') "..."
 write (nout,*) " Size in bytes of array pg                       : ",sizeof(pg)
 write (nout,*) " Size in bytes of array BoundaryDataTab          : ",sizeof(BoundaryDataTab)
 write (nout,*) " Size in bytes of array BoundaryDataPointer      : ",sizeof(BoundaryDataPointer)
 write (nout,*) " Size in bytes of array Array_Flu                : ",sizeof(Array_Flu)
! write (nout,*) " Size in bytes of array Array_Sol                : ",sizeof(Array_Sol)
 write (nout,*) " Size in bytes of array nPartIntorno             : ",sizeof(nPartIntorno)
 write (nout,*) " Size in bytes of array PartIntorno              : ",sizeof(PartIntorno)
 write (nout,*) " Size in bytes of array PartKernel               : ",sizeof(PartKernel)
 write (nout,*) " Size in bytes of array rag                      : ",sizeof(rag)
!
!AA406 start
 if ((Domain%tipo == "bsph").and.(DBSPH%n_w > 0)) then
 write (nout,*) " Size in bytes of array pg_w                     : ",sizeof(pg_w)
 write (nout,*) " Size in bytes of array nPartIntorno_fw          : ",sizeof(nPartIntorno_fw)
 write (nout,*) " Size in bytes of array PartIntorno_fw           : ",sizeof(PartIntorno_fw)
 write (nout,*) " Size in bytes of array kernel_fw                : ",sizeof(kernel_fw)
 write (nout,*) " Size in bytes of array rag_fw                   : ",sizeof(rag_fw)
 endif
!AA406 end

!AA501b start
 if (n_bodies > 0.) then
 write (nout,*) " Size in bytes of array nPartIntorno_bp_f        : ",sizeof(nPartIntorno_bp_f)
 write (nout,*) " Size in bytes of array PartIntorno_bp_f         : ",sizeof(PartIntorno_bp_f)
 write (nout,*) " Size in bytes of array KerDer_bp_f_cub_spl      : ",sizeof(KerDer_bp_f_cub_spl)
 write (nout,*) " Size in bytes of array KerDer_bp_f_Gal          : ",sizeof(KerDer_bp_f_Gal)
 write (nout,*) " Size in bytes of array rag_bp_f                 : ",sizeof(rag_bp_f)
 write (nout,*) " Size in bytes of array surf_body_part           : ",sizeof(surf_body_part) 
 write (nout,*) " Size in bytes of array nPartIntorno_bp_bp       : ",sizeof(nPartIntorno_bp_bp)
 write (nout,*) " Size in bytes of array PartIntorno_bp_bp        : ",sizeof(PartIntorno_bp_bp)
 write (nout,*) " Size in bytes of array rag_bp_bp                : ",sizeof(rag_bp_bp)
 endif
!AA501b end

!
 write (nout,'(1x,a)') "..."
!
 write (nout,'(1x,a)') " "
 write (nout,'(1x,a)') "   end allocation step. "
 write (nout,'(1x,a)') " "
!
!scrittura su file di output del numero cella e delle  particelle che gravano su di essa
 if ( Domain%ioutopt < 0 ) then
   write (nout,*) 
   write (nout,*) "Number of cells           NCELLS= ",grid%nmax
   write (nout,*) 
   write (nout,*) "======== CELLS AND RELATED PARTICLES =========="
   do i = 1,grid%nmax   
     if ( Icont(i+1) > Icont(i) ) then  !20051230
       write (nout,"(3(a,i5),a)") " cell", i," from",Icont(i)," to",Icont(i+1),"  particles:"
       write (nout,"(5i8)") NPartOrd(Icont(i):Icont(i+1)-1)  !20051230
     end if
   end do
 end if
!
 write (nout,*) 
 write (nout,*) 
 call s_ctime( nout )
 write (nout,*) 
 write (nscr,*) "Transient loop begins..."
 write (nout,*) "Transient loop begins..."
 write (nout,*)
!
!inizializza file post-processing
 if ( Domain%imemo_fr > 0 .OR. Domain%memo_fr > zero ) then
   open ( nres, file=nomefile(2) &
        , status="unknown"      &
        , access="sequential"   &
        , form  ="unformatted"  )
 else
   nres = -nres
 end if
!
!inizializza file restart in output
! if ( Domain%irest_fr > 0 .OR. Domain%rest_fr > zero ) then
!   open ( nsav, file=nomefile(3) &
!        , status="unknown"      &
!        , access="sequential"   &
!        , form  ="unformatted"  )
! else
!   nsav = -nsav
! end if
!
!inizializza file pelo libero in output
 if ( Domain%ipllb_fr > 0 .OR. Domain%pllb_fr > zero ) then
   open ( nplb, file=nomefile(4) &
        , status="unknown"      &
        , access="sequential"   &
        , form  ="formatted"  )
   write (nplb,"(a)") "tempo         free_surface_quota"
 else
   nplb = -nplb
 end if
!
!inizializza file fronte in output
 if ( Domain%imemo_fr > 0 .OR. Domain%memo_fr > zero ) then
   open ( nfro, file=nomefile(5) &
        , status="unknown"      &
        , access="sequential"   &
        , form  ="formatted"  )
   write (nfro,"(a)") "tempo         x fronte      (y fronte)    z fronte"
 else
   nfro = -nfro
 end if
!
! write (nres) menouno,const_m_9999,nag,ncord
! write (nres) domain, grid, &
!              NPartZone, NMedium, NPointst, NLines, NSections, &
!              NumVertici, NumFacce, NumTratti, NumBVertices, NumBSides, &
!              Vertice(1:SPACEDIM,1:NumVertici), &
!              BoundaryFace(1:NumFacce), &
!              Tratto (1:NumTratti),  &
!              BoundaryVertex(1:NumBVertices), & 
!              BoundarySide  (1:NumBSides), &
!              Partz(1:NPartZone), &
!              Med(1:nmedium), &
!              Control_Sections(0:NSections+1), &
!              Control_Points(1:NPointst), &
!              Control_Lines(1:NLines)

!
!.. create file for vtk to design the boundaries of the domain
!
 if (vtkconv) then
!
   prefix = nomecaso
   filevtk = "VTKConverter_"//prefix(1:len_trim(prefix))//"_domain.vtk"
   open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',status='unknown')
   write(unitvtk,'(a)') '# vtk DataFile Version 2.0'
   write(unitvtk,'(a,a)') 'Domain limits for the case:',prefix(1:len_trim(prefix))
   write(unitvtk,'(a)') 'ASCII'
   write(unitvtk,'(a)') 'DATASET POLYDATA'
   write(unitvtk,'(a,i8,a)') 'POINTS ',numvertici,' float'

   do i = 1,numvertici,4
     k1 = i
     k2 = k1 + 3
     if (k2 > numvertici) k2 = numvertici
     write (stringa,'(4(3(e12.5,1x)))') (vertice(1,k),vertice(2,k),vertice(3,k),k=k1,k2)
     stringa = adjustl(trim(stringa))
     write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
   end do
!
!.. 2D case
!
   if (ncord == 2) then
!
     nlinee = 0
     nvalori = 0
     do i = 1,numbsides
       if (boundaryside(i)%TIPO == 'peri') cycle
       nlinee = nlinee + 1
       nvalori = nvalori + 1
!!??       nvalori = nvalori + (boundaryside(i)%VERTEX(2) - boundaryside(i)%VERTEX(1) + 1)
       nvalori = nvalori + 2
     end do
     write(unitvtk,'(a,2i8)') 'LINES ', nlinee, nvalori
!
     allocate (verticecolore(numvertici))
     verticecolore = zero
     do i = 1,numbsides
       if (boundaryside(i)%TIPO == 'peri') cycle
       stringa = ' '
       k1 = boundaryside(i)%VERTEX(1) - 1
       k2 = boundaryside(i)%VERTEX(2) - 1
       write (stringa,'(a,2i10)') '2',k1,k2
       stringa = adjustl(trim(stringa))
       write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
       if (boundaryside(i)%TIPO == 'velo' .or. boundaryside(i)%TIPO == 'flow' .or. boundaryside(i)%TIPO == 'open') then
         verticecolore(k1+1) = 2.0
         verticecolore(k2+1) = 2.0
       else if (boundaryside(i)%TIPO == 'sour') then
         verticecolore(k1+1) = 1.0
         verticecolore(k2+1) = 1.0
       end if
     end do
!
     write(unitvtk,'(a,i8)') 'POINT_DATA ', numvertici
     write(unitvtk,'(a,a,a)') 'SCALARS ', prefix(1:len_trim(prefix)),' float 1'
     write(unitvtk,'(a)') 'LOOKUP_TABLE mytable'
!
     do i = 1,numvertici,8
       k1 = i
       k2 = k1 + 7
       if (k2 > numvertici) k2 = numvertici
       write (stringa,'(8(f10.3,1x))') (verticecolore(k),k=k1,k2)
       stringa = adjustl(trim(stringa))
       write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
     end do
     deallocate (verticecolore)
!
     write(unitvtk,'(a)') 'LOOKUP_TABLE mytable 3'
     write(unitvtk,'(a)') '0.0 0.0 0.0 1.0'
     write(unitvtk,'(a)') '0.0 0.0 1.0 1.0'
     write(unitvtk,'(a)') '1.0 0.0 0.0 1.0'
!
!.. 3D case
!
   else
!
     nlinee = 0
     nvalori = 0
     do i = 1,numfacce
       if (tratto(BoundaryFace(i)%stretch)%tipo == 'peri') cycle
       nlinee = nlinee + 1
       nvalori = nvalori + 1 + BoundaryFace(i)%nodes + 1
     end do
     write(unitvtk,'(a,2i8)') 'LINES ', nlinee, nvalori
!
     allocate (verticecolore(numvertici))
     verticecolore = 0.0
     do i = 1,numfacce
       if (tratto(BoundaryFace(i)%stretch)%tipo == 'peri') cycle
       stringa = ' '
       do k = 1,BoundaryFace(i)%nodes,8
         k1 = k 
         k2 = k1 + 7
         if (k2 > BoundaryFace(i)%nodes) k2 = BoundaryFace(i)%nodes
         write (stringa,'(10i10)') BoundaryFace(i)%nodes+1,(BoundaryFace(i)%node(kk)%name-1,kk=k1,k2),BoundaryFace(i)%node(1)%name-1
         stringa = adjustl(trim(stringa))
         write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
       end do
       if (tratto(BoundaryFace(i)%stretch)%tipo == 'velo' .or. tratto(BoundaryFace(i)%stretch)%tipo == 'flow' .or. &
           tratto(BoundaryFace(i)%stretch)%tipo == 'open') then
         do kk = k1,k2
           verticecolore(BoundaryFace(i)%node(kk)%name) = 2.0
           verticecolore(BoundaryFace(i)%node(kk)%name) = 2.0
         end do
       else if (tratto(BoundaryFace(i)%stretch)%tipo == 'sour') then
         do kk = k1,k2
           verticecolore(BoundaryFace(i)%node(kk)%name) = 1.0
           verticecolore(BoundaryFace(i)%node(kk)%name) = 1.0
         end do
       end if
     end do
!
     write(unitvtk,'(a,i8)') 'POINT_DATA ', numvertici
     write(unitvtk,'(a,a,a)') 'SCALARS ', prefix(1:len_trim(prefix)),' float 1'
     write(unitvtk,'(a)') 'LOOKUP_TABLE mytable'
!
     do i = 1,numvertici,8
       k1 = i
       k2 = k1 + 7
       if (k2 > numvertici) k2 = numvertici
       write (stringa,'(8(f10.3,1x))') (verticecolore(k),k=k1,k2)
       stringa = adjustl(trim(stringa))
       write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
     end do
     deallocate (verticecolore)
!
     write(unitvtk,'(a)') 'LOOKUP_TABLE mytable 3'
     write(unitvtk,'(a)') '0.0 0.0 0.0 1.0'
     write(unitvtk,'(a)') '0.0 0.0 1.0 1.0'
     write(unitvtk,'(a)') '1.0 0.0 0.0 1.0'
!
   end if
!
   flush(unitvtk)
   close (unitvtk)
 end if
!
!..
!
!AA406 sub
   if ( (Domain%tipo == "semi") .or. (Domain%tipo == "bsph") ) then
!
!.. calcolo DTmin reazione elastica boundary
   call EvaluateBER_TimeStep
!
   if (ncord==2) then
!       
     call start_and_stop(2,5)
     call loop_irre_2D 
     call start_and_stop(3,5)
!
   else if (ncord==3) then
!
     NumCellmax = Grid%nmax
     allocate ( GCBFPointers(NumCellmax,2), stat = ier )
     if (ier /= 0) then
       write (nout,'(1x,a,i2)') "   Array GCBFPointers not allocated. Error code: ",ier
       call diagnostic (arg1=4,arg3=nomsub)
     else
       write (nout,'(1x,a)') "   Array GCBFPointers successfully allocated "
     end if
!
!AA406
     if (Domain%tipo == "semi") then
!AA504 rm spart         
!.. To select the grid cells intercepting a boundary faces (Boundaries.f90)
     call GridCellBoundaryFacesIntersections3D (NumCellmax)
!
!.. To compute the local coordinates, solid angle and solid normal, relative to a boundary element (semi-analytic approach; Boundaries.f90)
     call ComputeBoundaryIntegralTab
!
!.. Computation of the boundary contributions for the continuity equation (Boundaries.f90)
     call ComputeKernelTable
!
!AA406 
     endif
     
!AA504 start
!Allocation and initialization of the Z_fluid_max array
!loop over the zones
     do i=1,NPartZone
        if (Partz(i)%IC_source_type == 2) then
            allocate(Z_fluid_max(Grid%ncd(1)*Grid%ncd(2)))
            Z_fluid_max = -999.
            allocate(q_max(Grid%ncd(1)*Grid%ncd(2)))
            q_max = 0.d0
            exit
        endif
     end do     
!AA504 end
     call start_and_stop(2,5)
!.. main loop
     call loop_irre_3D
     call start_and_stop(3,5)

!AA504 start 
!Writing the h_max array and deallocation of the Z_fluid_max array
     if (allocated(Z_fluid_max)) then
        call write_h_max        
        if (allocated(Z_fluid_max)) deallocate(Z_fluid_max)
     endif
!AA504 end
     
   end if 
!
 else
   call diagnostic (arg1=10,arg2=5,arg3=nomsub)
 end if
!
!!!!! deallocate ( NPartOrd,Icont, stat = ier )
!!!!! if (ier /= 0) then
!!!!!   write (nout,'(1x,a,i2)') "   Arrays NPartOrd,Icont  not deallocated. Error code: ",ier
!!!!!   call diagnostic (arg1=4,arg3=nomsub)
!!!!! else
!!!!!   write (nout,'(1x,a)') "   Arrays NPartOrd,Icont successfully deallocated "
!!!!! end if
!
!.. create file for vtk
!
  if (vtkconv) then
    prefix = nomecaso
    filevtk = "VTKConverter_"//prefix(1:len_trim(prefix))//".pvd"
    open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',status='unknown')
!
    write(unitvtk,'(a)') '<?xml version="1.0"?>'
    write(unitvtk,'(a)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
    write(unitvtk,'(a)') '  <Collection>'
    do i = 1,nblocchi
      stringa = ' '
      write(cargo,'(i10)') blocchi(i)
      cargo = adjustl(trim(cargo))
      filename = "VTKConverter_"//prefix(1:len_trim(prefix))//"_block_"//cargo(1:len_trim(cargo))//".vtu"
      stringa = '  <DataSet timestep="'
!AA504 sub     
      write(cargo,'(f15.6)')  Time_Block(i)
      cargo = adjustl(trim(cargo))
      stringa = stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" group="" part="'
      write(cargo,'(i10)') i
      cargo = adjustl(trim(cargo))
      stringa = stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" file="'//filename(1:len_trim(filename))//'"/>'
      write(unitvtk,'(a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '  </Collection>'
    write(unitvtk,'(a)') '</VTKFile>'
!
    flush(unitvtk)
    close (unitvtk)
  end if
!
!AA406 start
!.. create file for vtk (wall elements)
!
 if ((Domain%tipo == "bsph").and.(DBSPH%n_w > 0)) then
  if (vtkconv) then
    filevtk = "VTKConverter_wall_"//prefix(1:len_trim(prefix))//".pvd"
    open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',status='unknown')
    write(unitvtk,'(a)') '<?xml version="1.0"?>'
    write(unitvtk,'(a)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
    write(unitvtk,'(a)') '  <Collection>'
    do i = 1,nblocchi
      stringa = ' '
      write(cargo,'(i10)') blocchi(i)
      cargo = adjustl(trim(cargo))
      filename = "VTKConverter_"//prefix(1:len_trim(prefix))//"_block_wall_"//cargo(1:len_trim(cargo))//".vtu"
      stringa = '  <DataSet timestep="'
!AA504 sub       
      write(cargo,'(f15.6)')  Time_Block(i)
      cargo = adjustl(trim(cargo))
      stringa = stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" group="" part="'
      write(cargo,'(i10)') i
      cargo = adjustl(trim(cargo))
      stringa = stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" file="'//filename(1:len_trim(filename))//'"/>'
      write(unitvtk,'(a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '  </Collection>'
    write(unitvtk,'(a)') '</VTKFile>'
    flush(unitvtk)
    close (unitvtk)
  end if
  endif
!AA406 end

!AA501b start
 if (n_bodies > 0) then
 
! creation of the .pvd file for body particles
  if (vtkconv) then
    filevtk = "VTKConverter_body-part_"//prefix(1:len_trim(prefix))//".pvd"
    open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',status='unknown')
    write(unitvtk,'(a)') '<?xml version="1.0"?>'
    write(unitvtk,'(a)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
    write(unitvtk,'(a)') '  <Collection>'
    do i = 1,nblocchi
      stringa = ' '
      write(cargo,'(i10)') blocchi(i)
      cargo = adjustl(trim(cargo))
      filename = "VTKConverter_"//prefix(1:len_trim(prefix))//"_block_body-part_"//cargo(1:len_trim(cargo))//".vtu"
      stringa = '  <DataSet timestep="'
!AA504 sub 
      write(cargo,'(f15.6)')  Time_Block(i)
      cargo = adjustl(trim(cargo))
      stringa = stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" group="" part="'
      write(cargo,'(i10)') i
      cargo = adjustl(trim(cargo))
      stringa = stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" file="'//filename(1:len_trim(filename))//'"/>'
      write(unitvtk,'(a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '  </Collection>'
    write(unitvtk,'(a)') '</VTKFile>'
    flush(unitvtk)
    close (unitvtk)
  end if

! creation of the .pvd file for bodies
    if (vtkconv) then
    filevtk = "VTKConverter_body_"//prefix(1:len_trim(prefix))//".pvd"
    open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',status='unknown')
    write(unitvtk,'(a)') '<?xml version="1.0"?>'
    write(unitvtk,'(a)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
    write(unitvtk,'(a)') '  <Collection>'
    do i = 1,nblocchi
      stringa = ' '
      write(cargo,'(i10)') blocchi(i)
      cargo = adjustl(trim(cargo))
      filename = "VTKConverter_"//prefix(1:len_trim(prefix))//"_block_body_"//cargo(1:len_trim(cargo))//".vtu"
      stringa = '  <DataSet timestep="'
!AA504 sub       
      write(cargo,'(f15.6)')  Time_Block(i)
      cargo = adjustl(trim(cargo))
      stringa = stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" group="" part="'
      write(cargo,'(i10)') i
      cargo = adjustl(trim(cargo))
      stringa = stringa(1:len_trim(stringa))//cargo(1:len_trim(cargo))//'" file="'//filename(1:len_trim(filename))//'"/>'
      write(unitvtk,'(a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '  </Collection>'
    write(unitvtk,'(a)') '</VTKFile>'
    flush(unitvtk)
    close (unitvtk)
  end if
  
  endif
!AA501b end

!
return
end subroutine Gest_Trans
!---split

!cfile GetToken.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
character(80) function GetToken(ainp,itok,ioerr)
implicit none
!
character(*) :: ainp
integer(4)   :: itok, ioerr
integer(4)   :: n
!
integer(4)   :: number_token
integer(4), dimension(2,20) :: index_token
logical      :: blank 
!
 number_token = 0
 index_token  = 0
 blank        = .TRUE.

 do n = 1, len_trim(ainp)
 
   !if ( ainp(n:n) == "!" ) exit  !Commento dopo i dati
    if ( blank .AND. (ainp(n:n) /= " ") ) then
       number_token = number_token + 1
       index_token(1,number_token) = n
       index_token(2,number_token) = n
       blank = .FALSE.
    else if ( .NOT.blank .AND. (ainp(n:n) /= " ") ) then
       index_token(2,number_token) = n
    else if ( ainp(n:n) == " " ) then 
       blank = .TRUE.
    end if

 end do

 if ( itok <= number_token ) then
    ioerr = 0
    GetToken = ainp( index_token(1,itok):index_token(2,itok) )
 else
    ioerr = itok
    GetToken = ""
 end if
!
return
end function GetToken
!---split

!cfile GetVarPart.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : GetVarPart
!
! Last updating : May 08, 2012
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!
!************************************************************************************
! Module purpose : Module to act on pl(0) where are pre-fixed x,y,z and define the
!                  particle variables
!
! Calling routine: CalcVarp
!
! Called routines: 
!
!************************************************************************************
!
subroutine GetVarPart (pglocal)
!* agisce sulla pl(0) di cui sono prefissati x,y,h
!* e definisce le variabili della particella
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
!
!AA406 sub
integer(4)        :: nceli,ncel,npar,igridi,kgridi,jgridi,irang,krang,jrang,fw
!
integer(4)        :: mm,npj,irestocell
double precision  :: rijlocal,uni, pesoj,plocal,ro   
!
!.. Local Arrays ..
double precision, dimension (3) :: raglocal,vel
type (TyCtlPoint) :: pglocal

integer(4),      external :: CellIndices, CellNumber
double precision,external :: w
!
!.. Executable Statements ..
!
 nceli = pglocal%cella
 if (nceli == 0) return
 irestocell = CellIndices(nceli,igridi,jgridi,kgridi)
!
 uni    = zero
 plocal = zero
 ro     = zero
 vel(:) = zero
!
 npar = 0
 do jrang = jgridi-1,jgridi+1    ! ---- a  loop sulle 9 celle 
   do irang = igridi-1,igridi+1    ! ---- b  loop sulle 9 celle  
     do krang = kgridi-1,kgridi+1    ! ---- c  loop sulle 9 celle  
!
       ncel = CellNumber (irang,jrang,krang)
       if ( ncel == 0 ) cycle    ! cella fuori campo
!
       if (Icont(ncel+1) <= Icont(ncel)) cycle
       do mm = Icont(ncel),Icont(ncel+1)-1     ! loop sulle part di una cella  !20051230
         npj   = NPartOrd(mm)
!
         if ( pg(npj)%vel_type == "fix" ) cycle
!
         raglocal(:) = abs ( pglocal%coord(:)-pg(npj)%coord(:) )
         if ( raglocal(1) >= doubleh ) cycle
         if ( raglocal(2) >= doubleh ) cycle
         if ( raglocal(3) >= doubleh ) cycle
         rijlocal = raglocal(1)*raglocal(1) + raglocal(2)*raglocal(2) + raglocal(3)*raglocal(3)
         if ( rijlocal > doublesquareh ) cycle
         npar = npar+1
!
!AA406test
         rijlocal = dsqrt(rijlocal)       
!
         pesoj = pg(npj)%mass * w(rijlocal,Domain%h,Domain%coefke) / pg(npj)%dens
!
!.. calcolo den
         uni   = uni + pesoj
!
!.. calcolo vars
         plocal  = plocal  + pg(npj)%pres * pesoj        ! p locale
         ro = ro + pg(npj)%dens * pesoj                  ! ro VA FATTO SOLO SULLO STESSO MEZZO
         vel(:) = vel(:) + pg(npj)%vel(:) * pesoj        ! vx
!
       end do
!    
!AA406 start
       if ((Domain%tipo == "bsph").and.(DBSPH%n_w > 0)) then
! Loop over the neighbouring wall particles in the cell
          do fw = Icont_w(ncel),Icont_w(ncel+1)-1
             npj = NPartOrd_w(fw)
! Relative positions and distances
             raglocal(1:3) = pglocal%coord(:) - pg_w(npj)%coord(1:3)
             rijlocal = raglocal(1)*raglocal(1) + raglocal(2)*raglocal(2) + raglocal(3)*raglocal(3)
! Distance check
             if (rijlocal > doublesquareh) cycle
             rijlocal = dsqrt(rijlocal)
             pesoj = pg_w(npj)%mass * w(rijlocal,Domain%h,Domain%coefke) / pg_w(npj)%dens
             uni   = uni + pesoj
             plocal  = plocal  + pg_w(npj)%pres * pesoj
             ro = ro + pg_w(npj)%dens * pesoj                  
             vel(:) = vel(:) + pg_w(npj)%vel(:) * pesoj        
         end do 
      endif
!AA406 end
       
!
     end do
   end do
 end do
!
 pglocal%uni = uni
!!!! if ( npar > 0 ) then
! prova correzione per escludere punti che stanno al di fuori del campo
!
!AA501 sub
 if ( npar > 0 .and. uni >= 0.1) then
!
    pglocal%pres   = plocal/uni
    pglocal%dens   = ro/uni
    pglocal%vel(:) = vel(:)/uni
 end if
!
return
end subroutine GetVarPart
!---split

!cfile inidt2.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : inidt2
!
! Last updating : November 23, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Agate             2011           varie
! 04  Amicarelli        23/11/2011     multiple inlet
!
!************************************************************************************
! Module purpose : Module to act on pl(0) where are pre-fixed x,y,z and define the
!                  particle variables
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
  subroutine inidt2 

!Versione AdM del 19/11/05 - Da usarsi INSIEME a RUNDT2

!Calcola il passo dt iniziale secondo la procedura euristica proposta da AdM
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
!double precision,parameter :: k1 = 2.33333333333333d0  ! 7/3
!double precision,parameter :: Skmaxq = 0.55d0 * 0.55d0
!double precision,parameter :: Ckmax = 1.33333333333333d0  ! 4/3
!
!.. Local Scalars ..
!AA401
!  integer(4)       :: j
!  double precision :: celmedq, viscmed, viscequi, dtmed, dtmin, dtdiff1, diffmax, Denom, TermB, TermC  !, hh, viscterm
  integer(4)       :: j,npi,mate,ii
  double precision :: dtmin,dt_CFL,dt_dif,dt_vis,diffmax,U
!
!.. Executable Statements ..
!
!.. initializations
! 
  dtmin   = 1.0d+30                                                        
  diffmax = zero
!
!.. evaluates the time step minimum value depending on the different celerities in the different materials
!
!  do j = 1,nmedium
!!
!    celmedq  = Med(j)%celerita * Med(j)%celerita
!    viscmed  = max ( Med(j)%visc, Med(j)%numx )
!!    viscequi = 0.75d0 * Med(j)%alfaMon * Med(j)%celerita * Domain%h + 2.3333333d0 * viscmed
!!    viscterm = 3.287311d0 * viscequi / celmedq
!!    dtmed    = Dsqrt( viscterm * viscterm + 3.698224d0 * Domain%h * Domain%h / celmedq ) - viscterm
!    viscequi = Med(j)%alfaMon * Med(j)%celerita * Domain%h + k1 * viscmed
!    Denom    = Skmaxq * celmedq
!    termB    = Ckmax * viscequi / (Denom + Denom)
!    termC    = squareh / Denom
!    dtmed    = Dsqrt(termB * termB + termC) - termB
!!
!    diffmax  = max ( diffmax,Med(j)%codif)
!    dtdiff1  = half * squareh / (diffmax+0.000000001d0)
!    dtmin    = min ( dtmed, dtmin , dtdiff1)
!!
!  end do
!!
!!.. evaluates the time step current value applying the time coefficient to the found value
!!
!!  dt = Domain%cote * dtmin
!  dt = Domain%CFL * dtmin
!
!AA401 loop to compute the time step, according to 3 conditions:
!1) the CFL condition: dt_CFL=min(2h/(c+U))
!2) viscous stability condition dt_vis=min(rho*h**2/(0.5*mu)) 
!3) interface diffusion condition dt_diff=(h**2/2*teta)
  do ii = 1,indarrayFlu
    npi = Array_Flu(ii)
    mate = pg(npi)%imed
    U = sqrt(pg(npi)%vel(1)**2+pg(npi)%vel(2)**2+pg(npi)%vel(3)**2)
    dt_CFL = 2.*Domain%h/(Med(mate)%celerita+U)
!    dt_vis = pg(npi)%dens*Domain%h**2/(0.5*pg(npi)%mu)
    dt_vis = pg(npi)%dens*squareh/(half*pg(npi)%mu)
    dtmin  = min(dtmin,dt_CFL,dt_vis)
  end do
  do j = 1,nmedium
    diffmax  = max(diffmax,Med(j)%codif)
  end do
  dt_dif  = half * squareh / (diffmax+0.000000001d0)
  dtmin   = min (dtmin,dt_dif)
!AA401 end
!.. evaluates the time step current value applying the time coefficient to the found value
!
!AA405 
!initial dt for a jet
  if (indarrayFlu == 0) dtmin = 2.*Domain%h/(Med(1)%celerita)
!
!AA401
! CFL is used as a constant for every condition (the CFL, the viscous and the diffusive one)
  dt = Domain%CFL * dtmin
!
!.. initializes the medium time value for step to the initial value
!                             
  dt_average = dt
!
  return
  end subroutine inidt2
!---split

!cfile Init_arrays.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : Init_arrays
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
! Module purpose : Initialize the storage arrays for calculation
!
! Calling routine: Gest_input
!
! Called routines: none
!
!************************************************************************************

subroutine Init_Arrays
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
integer(4) :: i, j, n, i1
!
!.. Executable Statements ..
!
!Initialize Arrays

 do i = 1,size(Med)
   Med(i)%tipo         = "Empty   "
   Med(i)%modelloerosione = "        "
   Med(i)%index        = 0
   Med(i)%NIterSol     = 0
   Med(i)%den0         = zero
   Med(i)%eps          = zero
   Med(i)%celerita     = zero
   Med(i)%alfaMon      = zero
   Med(i)%betaMon      = zero
   Med(i)%visc         = zero
   Med(i)%coes         = zero
   Med(i)%numx         = zero
   Med(i)%mumx         = zero
   Med(i)%taucri       = zero
   Med(i)%cuin         = zero
   Med(i)%phi          = zero
   Med(i)%cons         = zero
   Med(i)%Cs           = zero
   Med(i)%RoughCoef    = zero
   Med(i)%D50          = zero
   Med(i)%SettlingCoef = zero
   Med(i)%Codif        = zero
   Med(i)%Gamma        = zero
   Med(i)%InitialIntEn = zero
 end do


 do i = 1,size(Partz)
    Partz(i)%label    = "        "
    Partz(i)%tipo     = "    "
    Partz(i)%shape    = " "
    Partz(i)%bend     = " "
    Partz(i)%pressure = "  "
    Partz(i)%move     = "   "
    Partz(i)%slip     = " "

    Partz(i)%ipool    = 0
    Partz(i)%npoints  = 0
    Partz(i)%icol     = 0
    Partz(i)%Medium   = 0
    Partz(i)%npointv  = 0
    Partz(i)%Indix(1) = 0
    Partz(i)%Indix(2) = 0
    Partz(i)%limit(1) = 1
    Partz(i)%limit(2) = 0
    Partz(i)%pool           = zero
    do j = 1,3
      Partz(i)%coordMM(j,1) = zero
      Partz(i)%coordMM(j,2) = zero
      Partz(i)%vel(j)       = zero
    end do
    do j = 0,3
      do n = 1,MAXPOINTSVLAW
        Partz(i)%vlaw(j,n) = zero
      end do
    end do
    Partz(i)%trampa = zero
    Partz(i)%valp = zero
 end do


 do i = 1,size(Control_Points)
   Control_Points(i)%cella    = 0
   Control_Points(i)%coord(1) = zero
   Control_Points(i)%coord(2) = zero
   Control_Points(i)%coord(3) = zero
   Control_Points(i)%vel  (1) = zero
   Control_Points(i)%vel  (2) = zero
   Control_Points(i)%vel  (3) = zero
   Control_Points(i)%pres     = zero
   Control_Points(i)%dens     = zero
   Control_Points(i)%uni      = zero
   Control_Points(i)%dist     = zero
 end do


 do i = 1,size(Section_Points)
   Section_Points(i)%cella    = 0
   Section_Points(i)%coord(1) = zero
   Section_Points(i)%coord(2) = zero
   Section_Points(i)%coord(3) = zero
   Section_Points(i)%vel  (1) = zero
   Section_Points(i)%vel  (2) = zero
   Section_Points(i)%vel  (3) = zero
   Section_Points(i)%pres     = zero
   Section_Points(i)%dens     = zero
   Section_Points(i)%uni      = zero
   Section_Points(i)%dist     = zero
 end do


 Control_Sections(0            )%Label = "Domain  "
 Control_Sections(1:NSections  )%Label = "        "
 Control_Sections(  NSections+1)%Label = "One more"
 do i1 = 1,size(Control_Sections)
   i = i1 - 1
   Control_Sections(i)%Tipo          = "  "
   Control_Sections(i)%Icont(1)      = 0
   Control_Sections(i)%Icont(2)      = 0
   Control_Sections(i)%ColorCode     = 1
   Control_Sections(i)%Constant(1)   = zero
   Control_Sections(i)%Constant(2)   = zero
   Control_Sections(i)%Constant(3)   = zero
   Control_Sections(i)%XYZRange(1,1) = zero
   Control_Sections(i)%XYZRange(2,1) = zero
   Control_Sections(i)%XYZRange(3,1) = zero
   Control_Sections(i)%XYZRange(1,2) = zero
   Control_Sections(i)%XYZRange(2,2) = zero
   Control_Sections(i)%XYZRange(3,2) = zero
   do j = 1,SPACEDIM
     do n = 1,SPACEDIM
       Control_Sections(i)%TGLsection(j,n) = zero
     end do
   end do
   Control_Sections(i)%TGLsection(1,2) = one
   Control_Sections(i)%TGLsection(2,1) =-one
   Control_Sections(i)%TGLsection(3,3) = one
 end do


 do i = 1,size(Control_Lines)
   Control_Lines (i)%label    = "Empty   "
   Control_Lines (i)%icont(1) = 0
   Control_Lines (i)%icont(2) = 0
 end do


 do j = 1,SPACEDIM
   Vertice(j,1:NumVertici) = zero
 end do

 do i = 1,size(BoundaryFace)
    do j = 1,MAXFACENODES
       BoundaryFace(i)%Node(j)%name = 0
       do n = 1,SPACEDIM
          BoundaryFace(i)%Node(j)%GX(n) = zero
          BoundaryFace(i)%Node(j)%LX(n) = zero
       end do
    end do
    BoundaryFace(i)%nodes   = 0
    BoundaryFace(i)%stretch  = 0
    BoundaryFace(i)%CloseParticles = 0
    BoundaryFace(i)%CloseParticles_maxQuota = const_m_9999
    BoundaryFace(i)%area    = zero
    do j = 1,SPACEDIM
      do n = 1,SPACEDIM
        BoundaryFace(i)%T(j,n)    = zero
        BoundaryFace(i)%RPsi(j,n) = zero
        BoundaryFace(i)%RFi(j,n)  = zero
      end do
      BoundaryFace(i)%velocity(j) = zero
    end do
 end do


 if (allocated(BFaceList)) then
   do i = 1,size(BFaceList)
     BFaceList(i) = 0
   end do
 end if


 do i = 1,size(BoundaryVertex)
   BoundaryVertex(i) = 0
 end do


 do i = 1,size(Tratto)
   Tratto(i)%tipo         = "    "
   Tratto(i)%ColorCode    = 0
   Tratto(i)%numvertices  = 0
   Tratto(i)%inivertex    = 0
   Tratto(i)%iniside      = 0
   Tratto(i)%iniface      = 0
   Tratto(i)%medium       = 0
   Tratto(i)%zone         = 0
   Tratto(i)%NormVelocity = zero
   Tratto(i)%trampa       = zero
   Tratto(i)%ShearCoeff   = zero
   do j = 1,SPACEDIM
     Tratto(i)%velocity(j) = zero
     Tratto(i)%PsiCoeff(j) = zero
     Tratto(i)%FiCoeff(j)  = zero
   end do
 end do


 do i = 1,size(BoundarySide)
   BoundarySide(i)%tipo           = "    "
   BoundarySide(i)%stretch        = 0
   BoundarySide(i)%previous_side  = 0
   BoundarySide(i)%vertex(1)      = 0
   BoundarySide(i)%vertex(2)      = 0
   BoundarySide(i)%CloseParticles = 0
   BoundarySide(i)%length         = zero
   BoundarySide(i)%CloseParticles_maxQuota = const_m_9999
   do n = 1,SPACEDIM
     BoundarySide(i)%T(n,1:SPACEDIM)  = zero
     BoundarySide(i)%R(n,1:SPACEDIM)  = zero
     BoundarySide(i)%RN(n,1:SPACEDIM) = zero
   end do
   BoundarySide(i)%angle = zero
   do j = 1,SPACEDIM
     BoundarySide(i)%velocity(J) = zero
   end do
 end do

! caso di restart non azzero domian e grid
if (Restart) return

 do j = 1,3
   Grid%ncd(j)    = 0
   Grid%dcd(j)    = zero
   Grid%extr(j,1) = zero
   Grid%extr(j,2) = zero
 end do
 Grid%nmax = 0


 Domain%tipo      = "semi"
 Domain%file      = "                                                                                "
 Domain%Psurf     = " "
 Domain%RandomPos = " "
 Domain%iplot_fr  = 0
 Domain%imemo_fr  = 0
 Domain%irest_fr  = 0
 Domain%icpoi_fr  = 0
 Domain%ipllb_fr  = 0
 Domain%ipllb_md  = 0
 Domain%istart    = 0
 Domain%ioutopt   = 0
 Domain%itmax     = 0
 do j = 1,3
   Domain%coord(j,1) = zero
   Domain%coord(j,2) = zero
   Domain%grav(j)    = zero
 end do
 Domain%tmax     = zero
 Domain%dd       = zero
 Domain%trunc    = zero
 Domain%coefke   = zero
 Domain%coefkacl = zero
! Domain%cote     = zero
 Domain%CFL      = zero
 Domain%prif     = zero
 Domain%plot_fr  = zero
 Domain%memo_fr  = zero
 Domain%rest_fr  = zero
 Domain%cpoi_fr  = zero
 Domain%pllb_fr  = zero
! Domain%TetaX    = zero
 Domain%TetaP    = zero
 Domain%TetaV    = zero
 Domain%pre      = zero
 Domain%h        = zero
 Domain%start    = zero
 Domain%NormFix  = .false.
 Domain%Slip     = .false.

return
end subroutine Init_Arrays
!---split

!cfile IsParticleInternal2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : IsParticleInternal2D
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
! Module purpose : Module check if a particle is internal to the domain 2D
!
! Calling routine: SetParticles
!
! Called routines: diagnostic
!
!************************************************************************************
!
  logical function IsParticleInternal2D (Nt, PX)
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
!.. Formal Arguments ..
  integer(4),intent(IN)                           :: Nt
  double precision,intent(IN),dimension(SPACEDIM) :: PX
!
!.. Local Scalars ..
  integer(4)       :: inizio,fine,iv, ni,n,n2
  double precision :: xa,za, xba,zba,  xi, zs
  character(len=lencard)  :: nomsub = "IsParticleInternal2D"
!
!.. Executable Statements ..
!
!.. initializations
!
  IsParticleInternal2D = .FALSE.
  ni = 0
  inizio = Tratto(Nt)%inivertex
  fine   = Tratto(Nt)%inivertex + Tratto(Nt)%numvertices - 2
!
!.. loops on the line vertices
!
  do iv = inizio, fine
!
!. checks for the storage limits
!
    if (iv > NumBVertices) then
      call diagnostic (arg1=10,arg2=2,arg3=nomsub)
    end if
!
!.. set the pointer to the current vertex
!
    n = BoundaryVertex(iv)
!
    if (n > NumVertici) then
      call diagnostic (arg1=10,arg2=3,arg3=nomsub)
    end if
!
!.. set the coordinates of the current vertex
!
    xa = Vertice(1,n)
    za = Vertice(3,n)
!
!. set the next vertex and its coordinates
!
    n2 = BoundaryVertex(iv+1)    
    xba = Vertice(1,n2)
    zba = Vertice(3,n2)
!
!.. the current segment is not vertical
!
    if ( abs(xa-xba) >= xyz_tolerance ) then
!
!.. the current segment is not horizontal
!
      if ( abs(za-zba) >= xyz_tolerance ) then
!
!.. evaluates the x coordinate of the particle projection on the segment along X
!
        xi = xa + (xba - xa) * (px(3) - za)/(zba - za)
!
      else
!
!.. the segment is horizontal: the X value of the mean point is assumed
!
          xi = half * (xa+xba)
!
      end if
    else
!
!.. the segment is vertical: the X value of the vertices is assumed
!
       xi = xa
!
    end if
!
!.. order the vertices in order to have the first one with the lower Z value
!
    if ( za > zba ) then
      zs  = za
      za  = zba
      zba = zs
    end if
!
!.. the Z value of the particle is inside the segment Z values
!
    if ( (PX(3)- za) > xyz_tolerance .AND. (PX(3) - zba) < xyz_tolerance ) then
!
      if ( xi > PX(1) ) then
        ni = ni + 1
      end if
    end if
!
  end do
!
  if ( MOD(ni,2) == 1 ) then
    IsParticleInternal2D = .TRUE.
  end if
!
  return
  end function IsParticleInternal2D
!---split

!cfile IsParticleInternal3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : IsParticleInternal3D
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
! Module purpose : Module check if a particle is internal to the domain 3D
!
! Calling routine: SetParticles
!
! Called routines: LocalNormalCoordinates
!                  IsPointInternal
!
!************************************************************************************
!
Logical Function IsParticleInternal3D ( mib, PX, IsopraS )
!
!Checks if point Px() is internal to the perimeter mib;
!in the affirmative returns 'true'; otherwise returns 'false'
!The perimeter can be both convex and concave !!
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
integer(4),parameter       :: intxy = 3
double precision,parameter :: eps = 0.001d0
!
!.. Formal Arguments ..
integer(4),      intent(IN)   :: mib
double precision,intent(IN),dimension(SPACEDIM) :: PX
integer(4),      intent(IN)   :: IsopraS
!
!.. Local Scalars ..
integer(4)       :: kf, nf, i, j, sd, nnodes, norig
integer(4)       :: Nints, IntSotto, IntSopra, fkod
double precision :: tpar
double precision,dimension(SPACEDIM) :: P1, Pint, LPint
double precision,dimension(3)        :: csi

!Dynamic Array
double precision,dimension(Tratto(mib)%numvertices) :: XYInts
!
!.. External Routines ..
logical, external :: IsPointInternal
!
!.. Executable Statements ..
!
 Nints    = 0
 IntSotto = 0
 IntSopra = 0
 IntSopra = IsopraS
 IsParticleInternal3D = .FALSE.

 do kf = Tratto(mib)%iniface, Tratto(mib)%iniface + Tratto(mib)%numvertices - 1

    nf     = BFaceList(kf)
    nnodes = 4
    if ( BoundaryFace(nf)%Node(4)%name <= 0 ) nnodes = 3
    norig  = nnodes                                      !nodo origine del sistema locale

    do sd = 1, SPACEDIM
       P1(sd) = Vertice(sd,BoundaryFace(nf)%Node(norig)%name)
    end do
    tpar = zero
    do sd = 1, SPACEDIM
       tpar = tpar + BoundaryFace(nf)%T(sd, 3) * (P1(sd) - PX(sd))
    end do

    if ( Abs(BoundaryFace(nf)%T(3, 3)) > eps ) then

        tpar = tpar / BoundaryFace(nf)%T(3, 3)

        do sd = 1, SPACEDIM                          !Pint()= coordinate globali del punto di intersezione
           Pint(sd) = PX(sd)
        end do
        Pint(3) = Pint(3) + tpar

        LPint = zero
        do sd = 1, PLANEDIM                          !LPint()= coordinate locali del punto di intersezione
           LPint(sd) = zero
           do j = 1, SPACEDIM
              LPint(sd) = LPint(sd) + BoundaryFace(nf)%T(j, sd) * (Pint(j) - P1(j))
           end do
        end do

        call LocalNormalCoordinates ( LPint, csi, nf )

        fkod = nnodes - 2
        if ( IsPointInternal ( fkod, csi ) ) then    !Il punto di intersezione è interno alla faccia;
                                                     !cioè la faccia interseca la verticale per Px.
                                                     !Memorizza coordinata z di intersezione
            Nints = Nints + 1
            XYInts(Nints) = Pint(3)
        end if
    end if

 end do

 if ( Nints > 0 ) then
    do i = 1, Nints
       if ( XYInts(i) <= PX(intxy) ) then
          IntSotto = IntSotto + 1
       Else
          IntSopra = IntSopra + 1
       end if
    end do
    if ( Mod(IntSotto,2) == 1 .AND. Mod(IntSopra,2) == 1 ) then
        IsParticleInternal3D = .TRUE.
    end if
 end if

return
End Function IsParticleInternal3D
!---split

!cfile IsPointInternal.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
Logical Function IsPointInternal ( fk, csi )

!Checks wheather a point with local normal coordinates csi() is
!internal to the face whose code is fk (=1 triangle, =2 parallelogram)

!.. assign modules
use GLOBAL_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Scalars ..
integer(4) :: i
integer(4) :: fk
!
!.. Local Arrays ..
double precision, dimension(1:SPACEDIM) :: csi
!
!.. Executable statements ..
!
 IsPointInternal = .FALSE.
!
 if ( fk == 1 ) then            !triangle
   do i = 1, 3
     if ( csi(i) < zero ) return
     if ( csi(i) > one ) return
   end do
   IsPointInternal = .TRUE.
 else if ( fk == 2 ) then        !parallelogram
   do i = 1, 2
     if ( csi(i) < zero ) return
     if ( csi(i) > one ) return
   end do
   IsPointInternal = .TRUE.
 end if
!
return
End Function IsPointInternal
!---split

!cfile language_writime2.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
module english_writime2
 character(24) :: cpulbl = "CPU     TIME (SECONDS):"
 character(24) :: elalbl = "ELAPSED TIME (SECONDS):"
 character(12) :: totlbl =  "(TOTAL)"
 character(12) :: przlbl =  "(PARTIAL)"
end module

module italiano_writime2
 character(24) :: cpulbl = "TEMPO CPU     (SECONDI):"
 character(24) :: elalbl = "TEMPO ELAPSED (SECONDI):"
 character(12) :: totlbl =  "(TOTALE)"
 character(12) :: przlbl =  "(PARZIALE)"
end module

module language_writime2
 use english_writime2
end module
!---split

!cfile lcase.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
character(80) function lcase(ainp)
implicit none

character(*) :: ainp
integer(4)   :: n, ia

 do n = 1,80
    lcase(n:n) = " "
 end do

 do n = 1, len_trim(ainp)
    ia = iachar(ainp(n:n))
!    write(6,*) n,ia,ainp(n:n); flush(6)
    if ( ia >= 65 .and. ia <= 90 ) then
        lcase(n:n) = char(ia+32)
    else
        lcase(n:n) = ainp(n:n)
    end if
 end do

return
end function lcase
!---split

!cfile LocalNormalCoordinates.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : LocalNormalCoordinates
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
! Module purpose : Module that Given the local coordinates PX(1 to 2) of a point P
!                  laing on the plane of the boundary face nf, the procedure
!                  computes in csi(1 to 3) the normal coordinates of the point Q
!                  corresponding to P in the inverse linear tranformation
!
! Calling routine: AddElasticBoundaryReaction_3D, ComputeBoundaryVolumeIntegrals_P0,
!                  FindCloseBoundaryFaces3D, CancelOutgoneParticles_3D, IsParticleInternal3D
! 
! Called routines: 
!
!************************************************************************************
!
subroutine LocalNormalCoordinates ( PX, csi, nf )
!Given the local coordinates PX(1 to 2) of a point P laing on the plane
!of the boundary face nf, the procedure
!computes in csi(1 to 3) the normal coordinates of the point Q
!corresponding to P in the inverse linear tranformation
!
!.. assign modules
use GLOBAL_MODULE
use AdM_User_Type
use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
double precision,intent(IN),   dimension(1:SPACEDIM) :: PX
double precision,intent(INOUT),dimension(1:SPACEDIM) :: csi
integer(4),      intent(IN)                          :: nf
!
!.. Local Scalars ..
integer(4)       :: i, j, k, nodes, fkod
double precision :: AA, BB, CC, DueArea, UsuDueArea, xj, yj, xk, yk
!
!.. Local Arrays ..
integer(4), dimension(3)   :: iseg = (/ 2,3,1 /)
integer(4), dimension(2,3) :: mainod ! = (/ 1,1, 2,3, 3,4 /) ! modifica per compatibilita xlf90
!
!.. Executable Statements ..
!
! modifica per compatibilita xlf90
  mainod(1, 1) = 1
  mainod(1, 2) = 2
  mainod(1, 3) = 3
  mainod(2, 1) = 1
  mainod(2, 2) = 3
  mainod(2, 3) = 4
! fine modifica
!
  nodes = 4
  if ( BoundaryFace(nf)%Node(4)%name <= 0 ) nodes = 3
  fkod = nodes - 2                    !=1 (triangolo), =2 (parallelogramma)
  DueArea = (3 - fkod) * BoundaryFace(nf)%Area
  UsuDueArea = one / DueArea
  do i = 1, 2
    j = iseg(i)
    k = iseg(j)
    xj = BoundaryFace(nf)%Node(mainod(fkod, j))%LX(1)
    yj = BoundaryFace(nf)%Node(mainod(fkod, j))%LX(2)
    xk = BoundaryFace(nf)%Node(mainod(fkod, k))%LX(1)
    yk = BoundaryFace(nf)%Node(mainod(fkod, k))%LX(2)
!
    AA = xj * yk - xk * yj
    BB = yj - yk
    CC = xk - xj
!
    csi(i) = (AA + BB * PX(1) + CC * PX(2)) * UsuDueArea
  end do
!
  csi(3) = one - (csi(1) + csi(2))
!
return
end subroutine LocalNormalCoordinates
!---split

!cfile ltrim.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
character(10) function ltrim(txt)
implicit none
character(*) :: txt

integer(4)   :: i,l,n

 l = len_trim(txt)
 do n = 1,l
    if ( txt(n:n) /= " " ) then
       txt(1:l-n+1) = txt(n:l)
       do i = l-n+2, l
          txt(i:i) = " "
       end do
       exit
    end if
 end do
 ltrim = trim(txt)

return
end function ltrim
!---split

!cfile MatrixProduct.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : MatrixProduct
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
! Module purpose : Module that Returns in CC() the product between matrices
!                  AA() and BB()
!
! Calling routine: BoundaryMassForceMatrix3D, BoundaryPressureGradientMatrix3D,
!AA501b
!                  vector_rotation
!
! Called routines: 
!
!************************************************************************************
!
subroutine MatrixProduct ( AA, BB, CC, nr, nrc, nc )
!Returns in CC() the product between matrices AA() and BB()
!nr =   number of rows of AA() and CC()
!nc =   number of columns of BB() and CC()
!nrc =  number of columns of AA() = number of rows of BB()
!
!.. assign modules
use GLOBAL_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
integer(4),      intent(IN) :: nr
integer(4),      intent(IN) :: nrc
integer(4),      intent(IN) :: nc
double precision,intent(IN),   dimension(nr,nrc) :: AA
double precision,intent(IN),   dimension(nrc,nc) :: BB
double precision,intent(INOUT),dimension(nr, nc) :: CC
!
!.. Local Scalars ..
integer(4)       :: i,j,k
double precision :: sum
!
!.. Executable Statements ..
!
 do i = 1, nr
    do j = 1, nc
       sum = zero
       do k = 1, nrc
          sum = sum + AA(i, k) * BB(k, j)
       end do
       CC(i, j) = sum
    end do
 end do
!
return
end subroutine MatrixProduct
!---split

!cfile MatrixTransposition.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : MatrixTransposition
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
! Module purpose : Module that Returns in AAT(n, m) the transponse of AA(m, n)
!
! Calling routine: BoundaryMassForceMatrix3D, BoundaryPressureGradientMatrix3D
!
! Called routines: 
!
!************************************************************************************
!
subroutine MatrixTransposition (AA, AAT, m, n )
!Returns in AAT(n, m) the transponse of AA(m, n)
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
integer(4),      intent(IN)   :: m
integer(4),      intent(IN)   :: n
double precision,intent(IN),   dimension(m,n) :: AA
double precision,intent(INOUT),dimension(n,m) :: AAT
!
!.. Local Scalars ..
integer(4) :: i,j
!
!.. Executable Statements ..
!
 do i = 1, n
    do j = 1, m
        AAT(i, j) = AA(j, i)
    end do
 end do

return
end subroutine MatrixTransposition
!---split

!cfile Memo_Ctl.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : Memo_Ctl
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
! Module purpose : Module to write results in the control lines and control points
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
!
subroutine Memo_Ctl
!
!.. assign modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Scalars ..
!
integer(4)     :: i,j
character(255) :: nomefilectl
!
!.. Executable Statements ..
!
!.. stampa control points
!
if (Npoints > 0) then
  write(nomefilectl,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_',it_corrente,".cpt"
  open ( ncpt, file=nomefilectl, status="unknown", form="formatted" )
  write (ncpt,*) "Control points "
  write (ncpt,'(1x,2(a,10x),3(a,8x),3(a,5x),a,7x,a)') &
  " Time","Iter","X Coord","Y Coord","Z Coord","X Velocity","Y Velocity","Z Velocity"," Pressure","Density "
  flush(ncpt)
  do i = 1,Npoints
    if (control_points(i)%cella == 0) then
      write (ncpt,'(a,i10,a,3g14.7)') "control point ",i," is outside. Coord=",Control_Points(i)%coord(:)
    else
      write (ncpt,'(g14.7,i14,8(1x,g14.7))') tempo, it_corrente &
                             ,Control_Points(i)%coord(:) &
                             ,Control_Points(i)%vel(:)   &
                             ,Control_Points(i)%pres     &
                             ,Control_Points(i)%dens
    end if
  end do
  close (ncpt)
end if
!
!.. stampa control lines
!
if (Nlines > 0) then
  write(nomefilectl,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_',it_corrente,".cln"
  open ( ncpt, file=nomefilectl, status="unknown", form="formatted" )
  write (ncpt,*) "Control lines "
  write (ncpt,'(1x,2(a,10x),3(a,8x),3(a,5x),a,7x,a)') &
  " Time","Iter","X Coord","Y Coord","Z Coord","X Velocity","Y Velocity","Z Velocity"," Pressure","Density "
  flush(ncpt)
  do i = 1,Nlines
    write (ncpt,*) "line #", i,"    Label ",Control_Lines(i)%label
    do j = Control_Lines(i)%icont(1),Control_Lines(i)%icont(2)
      if (control_points(j)%cella == 0) then
        write (ncpt,'(a,i10,a,g14.7)') "control point ",j," is outside. Coord=",Control_Points(j)%coord(:)
      else
        write (ncpt,'(g14.7,i14,8(1x,g14.7))') tempo, it_corrente &
                               ,Control_Points(j)%coord(:) &
                               ,Control_Points(j)%vel(:)   &
                               ,Control_Points(j)%pres     &
                               ,Control_Points(j)%dens
      end if
    end do
  end do
  close (ncpt)
end if
!AA601 rm part
!.. stampa control sections
!
if (Nsections > 0) then
  write(nomefilectl,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_',it_corrente,".csc"
  open ( ncpt, file=nomefilectl, status="unknown", form="formatted" )
  write (ncpt,*) "Control sections "
  write (ncpt,'(1x,2(a,10x),3(a,8x),3(a,5x),a,7x,a)') &
  " Time","Iter","X Coord","Y Coord","Z Coord","X Velocity","Y Velocity","Z Velocity"," Pressure","Density "
  flush(ncpt)
  do i = 1,Nsections
    write (ncpt,*) "section #", i,"    Label ",Control_sections(i)%label,"    Type ",Control_sections(i)%Tipo
    do j = Control_sections(i)%icont(1),Control_sections(i)%icont(2)
      if (control_points(j)%cella == 0) then
        write (ncpt,'(a,i10,a,g14.7)') "control point ",j," is outside. Coord=",Control_Points(j)%coord(:)
      else
        write (ncpt,'(g14.7,i14,8(1x,g14.7))') tempo, it_corrente &
                               ,Control_Points(j)%coord(:) &
                               ,Control_Points(j)%vel(:)   &
                               ,Control_Points(j)%pres     &
                               ,Control_Points(j)%dens
      end if
    end do
  end do
  close (ncpt)
end if

return

end subroutine Memo_Ctl
!---split

!cfile Memo_Results.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : Memo_Results
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
! Module purpose : Module to write results for post-processor
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D, diagnostic
!
! Called routines: 
!
!************************************************************************************
!
subroutine Memo_Results ( it, it_memo, it_rest, dtvel, str)
!
!.. assign modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
integer(4),      intent(IN)    :: it
integer(4),      intent(INOUT) :: it_memo
integer(4),      intent(INOUT) :: it_rest
double precision,intent(IN)    :: dtvel
character(6),    intent(IN)    :: str
!
!.. Local Scalars ..
integer(4) :: nrecords, restartcode
!
!
!.. Executable Statements ..
!
restartcode = 0
!
!=== SCRITTURA SU FILE RISULTATI e di RESTART ######################
!
! caso step a scrittura forzata (iniziale o finale)
! if ( it == it_memo ) then
 if ( index(str,'inizio')/=0 ) then
!
! memorizzazione informazioni fisse
    nrecords = 5
    if ( NumVertici   > 0 ) nrecords = nrecords + 1
    if ( NumFacce     > 0 ) nrecords = nrecords + 1
    if ( NumFacce     > 0 ) nrecords = nrecords + 1
    if ( NumTratti    > 0 ) nrecords = nrecords + 1
    if ( NPartZone    > 0 ) nrecords = nrecords + 1
    if ( NumBVertices > 0 ) nrecords = nrecords + 1
    if ( NumBSides    > 0 ) nrecords = nrecords + 1
!    if ( NPointst     > 0 ) nrecords = nrecords + 1
!    if ( NLines       > 0 ) nrecords = nrecords + 1
!    if ( NSections    > 0 ) nrecords = nrecords + 1
    write (nres) version,nrecords
    write (nres) Ncord, Nag, NMedium, NPartZone, &
                 NumVertici, NumFacce, NumTratti, NumBVertices, NumBSides, &
                 NPointst,NPoints,NPointsl,NPointse, NLines, NSections, GCBFVecDim, &
                 doubleh
    write (nres) domain
    write (nres) grid
    write (nres) Med(1:NMedium)
    if ( NumVertici   > 0 ) write (nres) Vertice(1:SPACEDIM,1:NumVertici)
    if ( NumFacce     > 0 ) write (nres) BoundaryFace(1:NumFacce)
    if ( NumFacce     > 0 ) write (nres) BFaceList(1:NumFacce)
    if ( NumTratti    > 0 ) write (nres) Tratto(1:NumTratti)
    if ( NPartZone    > 0 ) write (nres) Partz(1:NPartZone)
    if ( NumBVertices > 0 ) write (nres) BoundaryVertex(1:NumBVertices)
    if ( NumBSides    > 0 ) write (nres) BoundarySide(1:NumBSides)
!    if ( NPointst     > 0 ) write (nres) control_points(1:NPointst)
!    if ( NLines       > 0 ) write (nres) control_lines(1:NLines)
!    if ( NSections    > 0 ) write (nres) Control_Sections(0:NSections+1)
    flush(nres)
!
    write (nout,'(a,i10,a,f15.5)') " ----------------------------------------------------------------------------"
    write (nout,'(a,i10,a,f15.5)') " Results and restart heading saved   step: ",it,"   time: ",tempo
    write (nout,'(a,i10,a,f15.5)') " ----------------------------------------------------------------------------"
    write (nscr,'(a,i10,a,f15.5)') " ----------------------------------------------------------------------------"
    write (nscr,'(a,i10,a,f15.5)') " Results and restart heading saved   step: ",it,"   time: ",tempo
    write (nscr,'(a,i10,a,f15.5)') " ----------------------------------------------------------------------------"
!
  end if
!  else 
! memorizzazione informazioni variabili nel tempo
! caso memorizzazione nel loop
! 
! caso delta step con restart 
    if (Domain%irest_fr > 0) then
     if ( mod(it,Domain%irest_fr) == 0 ) then
       it_rest = it
     end if
! caso delta time con restart 
    else if (Domain%rest_fr > zero) then
     if ( it > 1 .and. mod(tempo,Domain%rest_fr) <= dtvel ) then
       it_rest = it
     end if
    end if
!
! caso delta step con solo risultati 
    if (Domain%imemo_fr > 0) then
     if ( mod(it,Domain%imemo_fr) == 0 ) then
        it_memo = it
      end if
! caso delta time con solo risultati 
     else if (Domain%memo_fr > zero) then
      if ( it > 1 .and. mod(tempo,Domain%memo_fr) <= dtvel ) then
        it_memo = it
      end if
    end if
!
!    if (it_rest == it) then
    if (it_rest == it .or. index(str,'inizio') /= 0 .or. index(str,'fine') /= 0) then
! memorizzo per restart
      restartcode = 1 ! se restartcode=1 memorizzo tutto il vettore pg
      write (nres) it,tempo,dt,nag,ncord,restartcode
      write (nres) pg(1:nag)
      flush(nres)
!
      if ( index(str,'inizio') == 0 ) then
        write (nout,'(a,i10,a,f15.5)') " --------------------------------------------------------------------"
        write (nout,'(a,i10,a,f15.5)') " Results and restart saved   step: ",it,"   time: ",tempo
        write (nout,'(a,i10,a,f15.5)') " --------------------------------------------------------------------"
        write (nscr,'(a,i10,a,f15.5)') " --------------------------------------------------------------------"
        write (nscr,'(a,i10,a,f15.5)') " Results and restart saved   step: ",it,"   time: ",tempo
        write (nscr,'(a,i10,a,f15.5)') " --------------------------------------------------------------------"
      end if
!
    else if (it_memo == it) then
! memorizzo per risultati
      restartcode = 0 ! se restartcode=0 memorizzo il vettore pg solo per visualizzazione 
      write (nres) it,tempo,dt,nag,ncord,restartcode
      write (nres) pg(1:nag)%coord(1),pg(1:nag)%coord(2),pg(1:nag)%coord(3), &
                   pg(1:nag)%vel(1),  pg(1:nag)%vel(2),  pg(1:nag)%vel(3), &
                   pg(1:nag)%pres,   &
                   pg(1:nag)%dens,   &
                   pg(1:nag)%mass,   &
                   pg(1:nag)%visc,   &
                   pg(1:nag)%IntEn,  &
                   pg(1:nag)%VolFra, &
                   pg(1:nag)%imed,   &
                   pg(1:nag)%icol
      flush(nres)
!
      if ( index(str,'inizio') == 0 ) then
        write (nout,'(a,i10,a,f15.5)') " --------------------------------------------------------"
        write (nout,'(a,i10,a,f15.5)') " Results saved   step: ",it,"   time: ",tempo
        write (nout,'(a,i10,a,f15.5)') " --------------------------------------------------------"
        write (nscr,'(a,i10,a,f15.5)') " --------------------------------------------------------"
        write (nscr,'(a,i10,a,f15.5)') " Results saved   step: ",it,"   time: ",tempo
        write (nscr,'(a,i10,a,f15.5)') " --------------------------------------------------------"
      end if
    end if
!
! end if
!
  return
  end subroutine Memo_Results
!---split

!cfile ModifyFaces.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ModifyFaces
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
! Module purpose : Module 
!
! Calling routine: Gest_Input
!
! Called routines: 
!
!************************************************************************************
!
subroutine ModifyFaces ( NumberEntities )
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
integer(4),intent(IN),dimension(20) :: NumberEntities
!
!.. Local Scalars ..
integer(4)       :: n,i,new
double precision :: d13, d24
!
!.. Executable Statements ..
!
!genera i nuovi triangoli, divisione lungo diagonale minima
 new = NumberEntities(11)
 do n = 1, NumberEntities(11)
    if ( BoundaryFace(n)%Node(4)%name == 0 ) cycle
    new = new + 1
    d13 = zero
    d24 = zero
    do i =1, SPACEDIM
       d13 = d13 + (Vertice(i,BoundaryFace(n)%Node(1)%name)-Vertice(i,BoundaryFace(n)%Node(3)%name)) * &
                   (Vertice(i,BoundaryFace(n)%Node(1)%name)-Vertice(i,BoundaryFace(n)%Node(3)%name))
       d24 = d24 + (Vertice(i,BoundaryFace(n)%Node(2)%name)-Vertice(i,BoundaryFace(n)%Node(4)%name)) * &
                   (Vertice(i,BoundaryFace(n)%Node(2)%name)-Vertice(i,BoundaryFace(n)%Node(4)%name))
    end do
    if ( d13 < d24 ) then
       BoundaryFace(new) = BoundaryFace(n)
       BoundaryFace(n)%Node(4)%name   =-BoundaryFace(n)%Node(4)%name
       BoundaryFace(new)%Node(2)%name = BoundaryFace(new)%Node(3)%name
       BoundaryFace(new)%Node(3)%name = BoundaryFace(new)%Node(4)%name
       BoundaryFace(new)%Node(4)%name =-99999999
    else
       BoundaryFace(new) = BoundaryFace(n)
       i                 = BoundaryFace(n)%Node(1)%name
       BoundaryFace(n)%Node(1)%name   = BoundaryFace(n)%Node(2)%name
       BoundaryFace(n)%Node(2)%name   = BoundaryFace(n)%Node(3)%name
       BoundaryFace(n)%Node(3)%name   = BoundaryFace(n)%Node(4)%name
       BoundaryFace(n)%Node(4)%name   =-i
       BoundaryFace(new)%Node(3)%name = BoundaryFace(new)%Node(4)%name
       BoundaryFace(new)%Node(4)%name =-99999999
    end if
 end do

return
end subroutine ModifyFaces
!---split

!cfile NormFix.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine NormFix

use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE

implicit none

integer(4)                    :: npi
double precision              :: unity
double precision,dimension(3) :: appo
!
!.. Executable statements
!
 if ( Ncord == 2 ) then
!
!$omp parallel do default(none) private(npi,appo,unity) shared(nag,Pg)
!
    do npi = 1,nag

       if (pg(npi)%cella == 0 .or. pg(npi)%vel_type == "std") cycle
!
       call InterFix(npi,appo,unity)
!
      !Componenti della normale generica 
       pg(npi)%mno = Dsqrt( (appo(1)*appo(1)) + (appo(3)*appo(3)) )
       pg(npi)%zer(1) = appo(1)/(pg(npi)%mno+0.0001d0)
       pg(npi)%zer(2) = zero
       pg(npi)%zer(3) = appo(2)/(pg(npi)%mno+0.0001d0)
!
       pg(npi)%ang    = (ATAN2 (pg(npi)%zer(1),pg(npi)%zer(3)))
!
       if ( abs(pg(npi)%zer(1)) < 0.5 .AND. abs(pg(npi)%zer(3)) < 0.5 ) then
          pg(npi)%zer(1) = zero
          pg(npi)%zer(2) = zero
          pg(npi)%zer(3) = zero
          pg(npi)%mno    = zero
          pg(npi)%ang    = zero
       end if
!
    end do
!
!$omp end parallel do
!
 else if ( Ncord == 3 ) then
!
!$omp parallel do default(none) private(npi,appo,unity) shared(nag,Pg)
!
    do npi = 1,nag

       if (pg(npi)%cella == 0 .or. pg(npi)%vel_type == "std") cycle

       call InterFix(npi,appo,unity)

       pg(npi)%mno = Dsqrt( (appo(1)*appo(1)) +(appo(2)*appo(2)) + (appo(3)*appo(3)) )
       pg(npi)%zer(:) = appo(:)/(pg(npi)%mno+0.0001d0)

!      pg(npi)%ang  =(ATAN2 (pg(npi)%zer(1),pg(npi)%zer(3)))
!      if ( abs(pg(npi)%zer(1)) < 0.5 .AND. abs(pg(npi)%zer(3)) < 0.5 ) then
   !3D    Decidere che verifica fare in 3D
!         pg(npi)%zer(1) = zero
!         pg(npi)%zer(2) = zero
!         pg(npi)%zer(3) = zero
!         pg(npi)%mno    = zero
!         pg(npi)%ang    = zero
!      end if

    end do
!
!$omp end parallel do
!
 end if
!
return
end subroutine NormFix
!---split

!cfile NumberSectionPoints.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
integer(4) function NumberSectionPoints ( values, opt )
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. implicit declarations
  implicit none
!
!.. dummy arguments
double precision,dimension(3,2) :: values
character(1)                    :: opt
!
!.. local Scalars ..
integer(4) :: n
!
!.. local Arrays ..
integer(4),dimension(3) :: Nmesh
!
!.. executable statements
!
!creazione mesh di lato dd
 Nmesh = 1
 do n = 1,SPACEDIM  !Ncord
    if ( n == 1 .AND. opt == "x" ) cycle
    if ( n == 2 .AND. opt == "y" ) cycle
    if ( n == 3 .AND. opt == "z" ) cycle
    Nmesh(n)   = nint ( (values(n,2)-values(n,1)) / Domain%dd )
 end do

 NumberSectionPoints = Product(Nmesh)

return
end function NumberSectionPoints
!---split

!cfile OrdGrid1.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : OrdGrid1
!
! Last updating : April 18, 2013
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Amicarelli/Agate  30Nov11        BSPH: wall element ordering
! 04  Amicarelli/Agate  13nov12        (AA501b) Body dynamics
! 05  Amicarelli/Agate  18apr13        add check on body_part_reorder
!
!************************************************************************************
! Module purpose : Module for sorting on the particles grid
!
! Calling routine: Gest_Input, Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: ncella
!
!************************************************************************************
!
subroutine OrdGrid1 ( nout )
!* ordinamento su griglia delle particelle
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
integer(4),intent(IN) :: nout
!
!.. Local Scalars ..
integer(4) :: npi,ncel,i
!
!.. Local Arrays ..
integer(4),dimension(Grid%nmax) :: numpartincelgiaaposto
!
!.. External Routines ..
integer(4),external :: ParticleCellNumber
!
!.. Executable Statements ..
!
!* annullamento ordine preesistente
 numpartincelgiaaposto = 0
 Icont                 = 0
!
!* 1 passata
!* si trova la cella in cui cade la particella e si contano le
!* particelle che ci sono in ogni cella
!
!!!!!!  do npi = 1,nag
!!!!!!    if ( pg(npi)%cella == 0 ) cycle
!!!!!!    ncel = ParticleCellNumber(npi)
!!!!!!    if ( ncel <= 0 .or. ncel > Grid%nmax) then
!!!!!!!
!!!!!!!.. the particle has been detected out of the grid and it is removed from the particle array
!!!!!!!
!!!!!!      write(nout,'(a,i7,a,i7,3x,3e15.8)') "ORDGRID1 particle #",npi, "   cell:",ncel,pg(npi)%coord(:)
!!!!!!!ag      EpOrdGrid(pg(npi)%imed) = EpOrdGrid(pg(npi)%imed) + 1
!!!!!!      pg(npi)%cella = 0
!!!!!!!ag      pg(npi) = PgZero
!!!!!!    else
!!!!!!!
!!!!!!!.. count the particles in each cell and the total number of particles (location nmax+1)
!!!!!!!
!!!!!!      Icont(ncel) = Icont(ncel) + 1
!!!!!!      pg(npi)%cella = ncel
!!!!!!      Icont(Grid%nmax+1) = Icont(Grid%nmax+1) + 1
!!!!!!    end if
!!!!!!  end do
!
!ag
!.. inserire eliminazione della particella e fare loop con while
!.. togliere EpOrdGrid(pg(npi)%imed) = EpOrdGrid(pg(npi)%imed) + 1    e   pg(npi) = PgZero   sopra
  npi = 0
  do while (nag > 0 .AND. npi < nag)
    npi = npi + 1
    ncel = ParticleCellNumber(pg(npi)%coord)
    if (pg(npi)%cella <= 0) then
      pg(npi) = pg(nag)
      pg(nag) = PgZero
!
!AA402
     if (Domain%RKscheme > 1) ts0_pg(nag) = ts_pgZero
!
      nag = nag - 1
      npi = npi - 1
!
    else if (ncel <= 0 .or. ncel > Grid%nmax) then
!
!.. the particle has been detected out of the grid and it is removed from the particle array
      write(nout,'(a,i7,a,i7,3x,3e15.8,a,i7)') &
            "ORDGRID1 particle #",npi, "   cell:",ncel,pg(npi)%coord(:),"  Total number of particles is = ",nag-1
      EpOrdGrid(pg(npi)%imed) = EpOrdGrid(pg(npi)%imed) + 1
      pg(npi) = pg(nag)
      pg(nag) = PgZero
!
!AA402
      if (Domain%RKscheme > 1) ts0_pg(nag) = ts_pgZero
!
      nag = nag - 1
      npi = npi - 1
!
    else
!
!.. count the particles in each cell and the total number of particles (location nmax+1)
      Icont(ncel) = Icont(ncel) + 1
      pg(npi)%cella = ncel
      Icont(Grid%nmax+1) = Icont(Grid%nmax+1) + 1
    end if
  end do
!
!ag
!
!* 2 passata
!* si definisce il puntatore d'inizio cella
  do i = Grid%nmax, 1, -1
    Icont(i) = Icont(i+1) - Icont(i)
  end do
!
!* 2 passata bis
!* si incrementa il contatore di 1
  do i = 1, Grid%nmax+1
    Icont(i) = Icont(i) + 1
  end do
!
!* 3 passata finale
!* si mettono le particelle al loro posto nel vettore NPartOrd
!
 do npi = 1,nag
    ncel = pg(npi)%cella
    if ( ncel == 0 ) cycle
    NPartOrd ( Icont(ncel) + numpartincelgiaaposto(ncel) ) = npi
    numpartincelgiaaposto(ncel) = numpartincelgiaaposto(ncel) + 1
 end do
 !
!AA406 start
!AA503sub 
  if ((DBSPH%n_w > 0) .and. ((it_corrente == it_start) .or. (Domain%body_part_reorder == 1)) ) then
!
    Icont_w = 0
    numpartincelgiaaposto = 0
!
!* 1 passata
!* si trova la cella in cui cade la particella e si contano le
!* particelle che ci sono in ogni cella
!
    npi = 0
!AA601 sub    
    do while (npi < (DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet))
       npi = npi + 1
       pg_w(npi)%cella = ParticleCellNumber(pg_w(npi)%coord)
!.. count the particles in each cell and the total number of particles (location nmax+1)
       Icont_w(pg_w(npi)%cella) = Icont_w(pg_w(npi)%cella) + 1
       Icont_w(Grid%nmax+1) = Icont_w(Grid%nmax+1) + 1
    end do
!
!* 2 passata
!* si definisce il puntatore d'inizio cella
    do i = Grid%nmax, 1, -1
       Icont_w(i) = Icont_w(i+1) - Icont_w(i)
    end do
!
!* 2 passata bis
!* si incrementa il contatore di 1
    do i = 1, Grid%nmax+1
       Icont_w(i) = Icont_w(i) + 1
    end do
!
!* 3 passata finale
!* si mettono le particelle al loro posto nel vettore NPartOrd
!AA601
    do npi = 1,(DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet)
       NPartOrd_w(Icont_w(pg_w(npi)%cella)+numpartincelgiaaposto(pg_w(npi)%cella)) = npi
       numpartincelgiaaposto(pg_w(npi)%cella) = numpartincelgiaaposto(pg_w(npi)%cella) + 1
!AA501
!       pg_w(npi)%wet = 1
!
    end do
! 
 endif
!AA406 end
!

!AA501b start
 if (n_bodies > 0) then
    Icont_bp = 0
    numpartincelgiaaposto = 0

!* 1 passata
!* si trova la cella in cui cade la particella e si contano le
!* particelle che ci sono in ogni cella
!
    npi = 0
    do while (npi < n_body_part)
       npi = npi + 1
       bp_arr(npi)%cell = ParticleCellNumber(bp_arr(npi)%pos)
!.. count the particles in each cell and the total number of particles (location nmax+1)
       Icont_bp(bp_arr(npi)%cell) = Icont_bp(bp_arr(npi)%cell) + 1
       Icont_bp(Grid%nmax+1) = Icont_bp(Grid%nmax+1) + 1
    end do

!* 2 passata
!* si definisce il puntatore d'inizio cella
    do i = Grid%nmax, 1, -1
       Icont_bp(i) = Icont_bp(i+1) - Icont_bp(i)
    end do
!
!* 2 passata bis
!* si incrementa il contatore di 1
    do i = 1, Grid%nmax+1
       Icont_bp(i) = Icont_bp(i) + 1
    end do

!* 3 passata finale
!* si mettono le particelle al loro posto nel vettore NPartOrd
    do npi = 1,n_body_part
       NPartOrd_bp(Icont_bp(bp_arr(npi)%cell)+numpartincelgiaaposto(bp_arr(npi)%cell)) = npi
       numpartincelgiaaposto(bp_arr(npi)%cell) = numpartincelgiaaposto(bp_arr(npi)%cell) + 1
    end do
! 
 endif
!AA501b end

return
end subroutine OrdGrid1
!---split

!cfile ParticleCellNumber.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ParticleCellNumber
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
! Module purpose : Module Returns the number of the grid cell where particle np is 
!                  at the present instant 
!                  if particle is outside of the reference grids returns -1
!
! Calling routine: 
!
! Called routines: CellNumber
!
!************************************************************************************
!
!AA406 sub
  Integer(4) function ParticleCellNumber (coord)

!.. Returns the number of the grid cell where particle np is at the present instant
!.. if particle is outside of the reference grids returns -1
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
!
!AA406 sub
  double precision,intent(in)  :: coord(3)
!
!.. Local Scalars ..
  integer(4)       :: i, j, k
  double precision :: xp, yp, zp
!
!.. External Routines ..
  integer(4), external :: CellNumber
!
!.. Executable Statements ..
!
!.. evaluates the progressive cell pointer in all the directions with respect the lower grid coordinates
!
!AA406 sub
  xp = coord(1) - Grid%extr(1,1)
  yp = coord(2) - Grid%extr(2,1)
  zp = coord(3) - Grid%extr(3,1)
!
  i = ceiling(xp / Grid%dcd(1))
  k = ceiling(zp / Grid%dcd(3)) 
!
  if (ncord == 3) then
    j = ceiling(yp / Grid%dcd(2))
  else
    j = 1
  end if
!
!.. return the cell number (-1 if the particle is out of the grid)
!
  ParticleCellNumber = CellNumber(i, j, k)
!
  return
  End Function ParticleCellNumber
!---split

!cfile Print_Results.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : Print_Results
!
! Last updating : November 23, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Agate             2011           varie
! 05  Amicarelli        23/11/2011     multiple inlet
! 06  Amicarelli/Agate  30Nov11        BSPH: wall element parameters
!AA501b comment
! 07  Amicarelli-Agate  13nov12        Body dynamics
!AA504
! 08  Amicarelli        08apr14        (v5.04) Modifications for 3D erosion criterion and elapsed time
!
!************************************************************************************
! Module purpose : Module to write results on output ascii file
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D, diagnostic
!
! Called routines: 
!
!************************************************************************************
!
  subroutine Print_Results (it, it_print, str)
!
!.. assign modules
  use FILES_ENTITIES
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
!
!.. implicit declarations
  implicit none
!
!.. dummy arguments
  integer(4),  intent(IN)    :: it
  integer(4),  intent(INOUT) :: it_print
  character(6),intent(IN)    :: str
!
!.. local Scalars ..
  integer(4)        :: npi,i,codice,dummy,OpCountot,SpCountot
  integer(4)        :: minlocvelo,maxlocvelo
  integer(4)        :: minlocvelx,maxlocvelx
  integer(4)        :: minlocvely,maxlocvely
  integer(4)        :: minlocvelz,maxlocvelz
  integer(4)        :: minlocpres,maxlocpres
  integer(4)        :: minlocdens,maxlocdens
  integer(4)        :: minlocvisc,maxlocvisc
  integer(4)        :: minloccodi,maxloccodi
  integer(4)        :: minlocInEn,maxlocInEn
  double precision  :: minvelo,maxvelo
  double precision  :: minvelx,maxvelx
  double precision  :: minvely,maxvely
  double precision  :: minvelz,maxvelz
  double precision  :: minpres,maxpres
  double precision  :: mindens,maxdens
  double precision  :: minvisc,maxvisc
  double precision  :: mincodi,maxcodi
  double precision  :: minInEn,maxInEn
  double precision  :: modvel
  character(len=42) :: fmt100="(a,i10,a,e18.9,a,e18.9,a,i  ,a,i  ,a,i  ))"
  character(len=47) :: fmt101="(a,2(1x,f10.4,1x,a,1x,i6,1x,a,3(1x,f8.2,1x,a)))"
  character(len=47) :: fmt102="(a,2(1x,f10.1,1x,a,1x,i6,1x,a,3(1x,f8.2,1x,a)))"
  character(len=47) :: fmt103="(a,2(1x,f10.6,1x,a,1x,i6,1x,a,3(1x,f8.2,1x,a)))"
  character(len=47) :: fmt104="(a,2(1x,g10.4,1x,a,1x,i6,1x,a,3(1x,f8.2,1x,a)))"
!AA504  
  character(len=47) :: fmt105="(a,2(1x,e10.4,1x,a,1x,i6,1x,a,3(1x,f8.2,1x,a)))"
  character(len=12) :: stringa
  character(len=2)  :: coppia
!
!AA406 start
  integer(4)        :: minlocvelo_w,maxlocvelo_w,minlocpres_w,maxlocpres_w
  double precision  :: minvelo_w,maxvelo_w,minpres_w,maxpres_w
!AA406 end

!AA501b start
  integer(4)        :: minlocvelo_bp,maxlocvelo_bp,minlocpres_bp,maxlocpres_bp
  double precision  :: minvelo_bp,maxvelo_bp,minpres_bp,maxpres_bp
  integer(4)        :: minlocvelo_body,maxlocvelo_body,minlocomega_body,maxlocomega_body,nbi
  double precision  :: minvelo_body,maxvelo_body,minomega_body,maxomega_body,modomega
!AA501b end

!AA504 start
  integer(4)        :: minloctau_tauc,maxloctau_tauc,minlock_BetaGamma,maxlock_BetaGamma,minlocu_star,maxlocu_star
  integer(4)        :: machine_Julian_day,machine_hour,machine_minute,machine_second
  double precision  :: mintau_tauc,maxtau_tauc,mink_BetaGamma,maxk_BetaGamma,minu_star,maxu_star,time_elapsed_tot_est
!AA504 end

!
!.. local Arrays ..
  integer(4),  dimension(1) :: pos
!
!.. executable statements
!
!.. detects the output code for prints
!
  codice = abs(Domain%ioutopt)
  if ( index(str,'inizio')/=0 .or. index(str,'fine')/=0) then
    if ( codice == 0 ) return
  else
    if ( codice == 0 .or. mod(it,codice) /= 0 ) return
  end if
!
!.. evaluates the proper formats
!
  dummy = count(pg(1:nag)%cella >0)
  write(stringa,'(i12)') dummy
  i = len_trim(adjustl(stringa))
  write (coppia,'(i2)') i
  fmt100(27:28) = adjustl(coppia)
!
  OpCountot = 0
  SpCountot = 0
  do i = 1,NMedium
    OpCountot = OpCountot + OpCount(i)
    SpCountot = SpCountot + SpCount(i)
  end do
  write(stringa,'(i12)') OpCountot
  i = len_trim(adjustl(stringa))
  write (coppia,'(i2)') i
  fmt100(33:34) = adjustl(coppia)
!
  write(stringa,'(i12)') SpCountot
  i = len_trim(adjustl(stringa))
  write (coppia,'(i2)') i
  fmt100(39:40) = adjustl(coppia)
!
!AA406 
  if (nag>0) then
!
!.. searches for the minimum and maximum values of all the distributions and 
!.. the corresponding particle identifiers
!
!.. module of the velocity and their components
!
  minvelo = max_positive_number
  maxvelo = max_negative_number
!
!AA406 start
 if ((Domain%tipo == "bsph").and.(DBSPH%n_w > 0)) then
     minvelo_w = max_positive_number
     maxvelo_w = max_negative_number
  endif
!AA406 end
!

!AA501b start
 if (n_bodies > 0) then
     minvelo_bp = max_positive_number
     maxvelo_bp = max_negative_number
     minvelo_body = max_positive_number
     maxvelo_body = max_negative_number
     minomega_body = max_positive_number
     maxomega_body = max_negative_number
  endif
!AA501b end

!AA504 start
 if ((Granular_flows_options%erosion_flag.ne.1).and.(Granular_flows_options%ID_erosion_criterion==1)) then
    mintau_tauc = max_positive_number
    maxtau_tauc = max_negative_number
    minu_star = max_positive_number
    maxu_star = max_negative_number
    if (Granular_flows_options%ID_erosion_criterion == 1) then    
       mink_BetaGamma = max_positive_number
       maxk_BetaGamma = max_negative_number
    endif
 endif
!AA504 end

  do npi = 1,nag
!   
    if (pg(npi)%cella == 0) cycle
    modvel = pg(npi)%vel(1)*pg(npi)%vel(1)+pg(npi)%vel(2)*pg(npi)%vel(2)+pg(npi)%vel(3)*pg(npi)%vel(3)
    if (modvel < minvelo) then
      minvelo = modvel
      minlocvelo = npi
    end if
    if (modvel > maxvelo) then
      maxvelo = modvel
      maxlocvelo = npi
    end if
 !
  end do
  minvelo = Dsqrt(minvelo)
  maxvelo = Dsqrt(maxvelo)
!
  minvelx = minval(pg(1:nag)%vel(1),mask=pg(1:nag)%cella /= 0)
  maxvelx = maxval(pg(1:nag)%vel(1),mask=pg(1:nag)%cella /= 0)
  pos = minloc(pg(1:nag)%vel(1),mask=pg(1:nag)%cella /= 0)
  minlocvelx = pos(1)
  pos = maxloc(pg(1:nag)%vel(1),mask=pg(1:nag)%cella /= 0)
  maxlocvelx = pos(1)
!
  minvely = minval(pg(1:nag)%vel(2),mask=pg(1:nag)%cella /= 0)
  maxvely = maxval(pg(1:nag)%vel(2),mask=pg(1:nag)%cella /= 0)
  pos = minloc(pg(1:nag)%vel(2),mask=pg(1:nag)%cella /= 0)
  minlocvely = pos(1)
  pos = maxloc(pg(1:nag)%vel(2),mask=pg(1:nag)%cella /= 0)
  maxlocvely = pos(1)
!
  minvelz = minval(pg(1:nag)%vel(3),mask=pg(1:nag)%cella /= 0)
  maxvelz = maxval(pg(1:nag)%vel(3),mask=pg(1:nag)%cella /= 0)
  pos = minloc(pg(1:nag)%vel(3),mask=pg(1:nag)%cella /= 0)
  minlocvelz = pos(1)
  pos = maxloc(pg(1:nag)%vel(3),mask=pg(1:nag)%cella /= 0)
  maxlocvelz = pos(1)
!
!.. pressure
!
  minpres = minval(pg(1:nag)%pres,mask=pg(1:nag)%cella /= 0)
  maxpres = maxval(pg(1:nag)%pres,mask=pg(1:nag)%cella /= 0)
  pos = minloc(pg(1:nag)%pres,mask=pg(1:nag)%cella /= 0)
  minlocpres = pos(1)
  pos = maxloc(pg(1:nag)%pres,mask=pg(1:nag)%cella /= 0)
  maxlocpres = pos(1)
!
!.. density
!
  mindens = minval(pg(1:nag)%dens,mask=pg(1:nag)%cella /= 0)
  maxdens = maxval(pg(1:nag)%dens,mask=pg(1:nag)%cella /= 0)
  pos = minloc(pg(1:nag)%dens,mask=pg(1:nag)%cella /= 0)
  minlocdens = pos(1)
  pos = maxloc(pg(1:nag)%dens,mask=pg(1:nag)%cella /= 0)
  maxlocdens = pos(1)
!
!.. viscosity
!
  minvisc = minval(pg(1:nag)%visc,mask=pg(1:nag)%cella /= 0)
  maxvisc = maxval(pg(1:nag)%visc,mask=pg(1:nag)%cella /= 0)
  pos = minloc(pg(1:nag)%visc,mask=pg(1:nag)%cella /= 0)
  minlocvisc = pos(1)
  pos = maxloc(pg(1:nag)%visc,mask=pg(1:nag)%cella /= 0)
  maxlocvisc = pos(1)
!
!.. diffusion coefficient
!
  if (diffusione) then
    mincodi = minval(pg(1:nag)%coefdif,mask=pg(1:nag)%cella /= 0)
    maxcodi = maxval(pg(1:nag)%coefdif,mask=pg(1:nag)%cella /= 0)
    pos = minloc(pg(1:nag)%coefdif,mask=pg(1:nag)%cella /= 0)
    minloccodi = pos(1)
    pos = maxloc(pg(1:nag)%coefdif,mask=pg(1:nag)%cella /= 0)
    maxloccodi = pos(1)
  end if
!
!.. esplosion coefficient
!
  if (esplosione) then
    minInEn = minval(pg(1:nag)%IntEn,mask=pg(1:nag)%cella /= 0)
    maxInEn = maxval(pg(1:nag)%IntEn,mask=pg(1:nag)%cella /= 0)
    pos = minloc(pg(1:nag)%IntEn,mask=pg(1:nag)%cella /= 0)
    minlocInEn = pos(1)
    pos = maxloc(pg(1:nag)%IntEn,mask=pg(1:nag)%cella /= 0)
    maxlocInEn = pos(1)
  end if
!
!
!AA406 start
!Wall parameter limits
 if ((Domain%tipo == "bsph").and.(DBSPH%n_w > 0)) then
! Wall pressure 
     minpres_w = minval(pg_w(1:DBSPH%n_w)%pres,mask=pg_w(1:DBSPH%n_w)%cella /= 0)
     maxpres_w = maxval(pg_w(1:DBSPH%n_w)%pres,mask=pg_w(1:DBSPH%n_w)%cella /= 0)
     pos = minloc(pg_w(1:DBSPH%n_w)%pres,mask=pg_w(1:DBSPH%n_w)%cella /= 0)
     minlocpres_w = pos(1)
     pos = maxloc(pg_w(1:DBSPH%n_w)%pres,mask=pg_w(1:DBSPH%n_w)%cella /= 0)
     maxlocpres_w = pos(1)
! Wall velocity
     do npi = 1,DBSPH%n_w
        if (pg_w(npi)%cella == 0) cycle
        modvel = pg_w(npi)%vel(1)*pg_w(npi)%vel(1)+pg_w(npi)%vel(2)*pg_w(npi)%vel(2)+pg_w(npi)%vel(3)*pg_w(npi)%vel(3)
        if (modvel < minvelo_w) then
           minvelo_w = modvel
           minlocvelo_w = npi
        end if
        if (modvel > maxvelo_w) then
           maxvelo_w = modvel
           maxlocvelo_w = npi
        end if
     end do
     minvelo_w = Dsqrt(minvelo_w)
     maxvelo_w = Dsqrt(maxvelo_w)
  endif
!AA406 end

!AA501b start
!limits for body particles and bodies
 if (n_bodies > 0) then
! pressure (body particles)
     minpres_bp = minval(bp_arr(1:n_body_part)%pres,mask=bp_arr(1:n_body_part)%cell /= 0)
     maxpres_bp = maxval(bp_arr(1:n_body_part)%pres,mask=bp_arr(1:n_body_part)%cell /= 0)
     pos = minloc(bp_arr(1:n_body_part)%pres,mask=bp_arr(1:n_body_part)%cell /= 0)
     minlocpres_bp = pos(1)
     pos = maxloc(bp_arr(1:n_body_part)%pres,mask=bp_arr(1:n_body_part)%cell /= 0)
     maxlocpres_bp = pos(1)
! velocity (body particles)
     do npi = 1,n_body_part
        if (bp_arr(npi)%cell == 0) cycle
        modvel = bp_arr(npi)%vel(1)*bp_arr(npi)%vel(1)+bp_arr(npi)%vel(2)*bp_arr(npi)%vel(2)+bp_arr(npi)%vel(3)*bp_arr(npi)%vel(3)
        if (modvel < minvelo_bp) then
           minvelo_bp = modvel
           minlocvelo_bp = npi
        end if
        if (modvel > maxvelo_bp) then
           maxvelo_bp = modvel
           maxlocvelo_bp = npi
        end if
     end do
     minvelo_bp = Dsqrt(minvelo_bp)
     maxvelo_bp = Dsqrt(maxvelo_bp)
! velocity (bodies)
     do nbi = 1,n_bodies
        modvel = body_arr(nbi)%u_CM(1)*body_arr(nbi)%u_CM(1)+body_arr(nbi)%u_CM(2)*body_arr(nbi)%u_CM(2)+ &
                 body_arr(nbi)%u_CM(3)*body_arr(nbi)%u_CM(3)
        if (modvel < minvelo_body) then
           minvelo_body = modvel
           minlocvelo_body = nbi
        end if
        if (modvel > maxvelo_body) then
           maxvelo_body = modvel
           maxlocvelo_body = nbi
        end if
     end do
     minvelo_body = Dsqrt(minvelo_body)
     maxvelo_body = Dsqrt(maxvelo_body)
! angular velocity (bodies)
     do nbi = 1,n_bodies
        modomega = dsqrt(dot_product(body_arr(nbi)%omega,body_arr(nbi)%omega))
        if (modomega < minomega_body) then
           minomega_body = modomega
           minlocomega_body = nbi
        end if
        if (modomega > maxomega_body) then
           maxomega_body = modomega
           maxlocomega_body = nbi
        end if
     end do
     minomega_body = Dsqrt(minomega_body)
     maxomega_body = Dsqrt(maxomega_body)
  endif
!AA501b end

!AA504 start
!Limits for supplementary bed load transport parameters 
 if ((Granular_flows_options%erosion_flag.ne.1).and.(Granular_flows_options%ID_erosion_criterion==1)) then
! tau_tauc 
     mintau_tauc = minval(pg(1:nag)%tau_tauc,mask=pg(1:nag)%cella /= 0)
     maxtau_tauc = maxval(pg(1:nag)%tau_tauc,mask=pg(1:nag)%cella /= 0)
     pos = minloc(pg(1:nag)%tau_tauc,mask=pg(1:nag)%cella /= 0)
     minloctau_tauc = pos(1)
     pos = maxloc(pg(1:nag)%tau_tauc,mask=pg(1:nag)%cella /= 0)
     maxloctau_tauc = pos(1)
! u_star  
     minu_star = minval(pg(1:nag)%u_star,mask=pg(1:nag)%cella /= 0)
     maxu_star = maxval(pg(1:nag)%u_star,mask=pg(1:nag)%cella /= 0)
     pos = minloc(pg(1:nag)%u_star,mask=pg(1:nag)%cella /= 0)
     minlocu_star = pos(1)
     pos = maxloc(pg(1:nag)%u_star,mask=pg(1:nag)%cella /= 0)
     maxlocu_star = pos(1) 
! k_BetaGamma 
     if (Granular_flows_options%ID_erosion_criterion == 1) then   
        mink_BetaGamma = minval(pg(1:nag)%k_BetaGamma,mask=pg(1:nag)%cella /= 0)
        maxk_BetaGamma = maxval(pg(1:nag)%k_BetaGamma,mask=pg(1:nag)%cella /= 0)
        pos = minloc(pg(1:nag)%k_BetaGamma,mask=pg(1:nag)%cella /= 0)
        minlock_BetaGamma = pos(1)
        pos = maxloc(pg(1:nag)%k_BetaGamma,mask=pg(1:nag)%cella /= 0)
        maxlock_BetaGamma = pos(1)      
     endif   
 endif
!AA504 end

!AA504 start
 if (exetype == "linux") then
    if (Domain%tmax>0.d0) then
       call system("date +%j%H%M%S > date_now.txt")
       open (unit_time_elapsed,file='date_now.txt',status="unknown",form="formatted")
       read (unit_time_elapsed,'(i3,i2,i2,i2)') machine_Julian_day,machine_hour,machine_minute,machine_second
       close (unit_time_elapsed)
       call system("rm -f date_now.txt")
    endif   
 endif 
!AA504 end
 
!
!.. final prints
!
!AA406 rm (moved some lines before)
!AA405 
!  if (nag>0) then
!
  write (nout,'(128("."))') 
  write (nout,fmt100) &
  " Print at:     | step: ",it," | time: ",tempo," | Dt: ",dt," | Particles: inside ",dummy," gone out ",OpCountot, &
  " gone in ",SpCountot
!AA504 start  
  if (exetype == "linux") then
    time_elapsed_tot_est = ((Domain%t_pre_iter-Domain%t0) + &
                            (machine_Julian_day*24*60*60+machine_hour*60*60+machine_minute*60+machine_second - Domain%t_pre_iter) * 1.d0) / (3600.0d0)  
    if (time_elapsed_tot_est<0.d0) time_elapsed_tot_est = time_elapsed_tot_est + 366.d0*24.d0*60.d0*60.d0   
    write (nout,'(a,g12.5,a,g12.5,a)') "Elapsed time: ",time_elapsed_tot_est," hours = ",time_elapsed_tot_est/24.d0," days."
    time_elapsed_tot_est = ((Domain%t_pre_iter-Domain%t0) + &
                            (machine_Julian_day*24*60*60+machine_hour*60*60+machine_minute*60+machine_second - Domain%t_pre_iter) &
                             * (Domain%tmax/tempo)) / (3600.0d0)  
    if (time_elapsed_tot_est<0.d0) time_elapsed_tot_est = time_elapsed_tot_est + 366.d0*24.d0*60.d0*60.d0  
    write (nout,'(a,g12.5,a,g12.5,a)') "Elapsed time (at the end of the simulation, real time estimation): ",time_elapsed_tot_est," hours = ", &
                                        time_elapsed_tot_est/24.d0," days."
  end if
!AA504 end  
  write (nout,fmt101) &
  " ............. |  Min. val. |Min.loc.| X coord. | Y coord. | Z coord. ||  Max. val. |Max.loc.| X coord. | Y coord. | Z coord. |"
  write (nout,fmt101) &
  "  Tot velocity |",minvelo,"|",minlocvelo,"|",pg(minlocvelo)%coord(1),"|",pg(minlocvelo)%coord(2),"|",pg(minlocvelo)%coord(3), &
  "||",maxvelo,"|",maxlocvelo,"|",pg(maxlocvelo)%coord(1),"|",pg(maxlocvelo)%coord(2),"|",pg(maxlocvelo)%coord(3),"|"
  write (nout,fmt101) &
  "  velocity x   |",minvelx,"|",minlocvelx,"|",pg(minlocvelx)%coord(1),"|",pg(minlocvelx)%coord(2),"|",pg(minlocvelx)%coord(3), &
  "||",maxvelx,"|",maxlocvelx,"|",pg(maxlocvelx)%coord(1),"|",pg(maxlocvelx)%coord(2),"|",pg(maxlocvelx)%coord(3),"|"
  write (nout,fmt101) &
  "  velocity y   |",minvely,"|",minlocvely,"|",pg(minlocvely)%coord(1),"|",pg(minlocvely)%coord(2),"|",pg(minlocvely)%coord(3), &
  "||",maxvely,"|",maxlocvely,"|",pg(maxlocvely)%coord(1),"|",pg(maxlocvely)%coord(2),"|",pg(maxlocvely)%coord(3),"|"
  write (nout,fmt101) &
  "  velocity z   |",minvelz,"|",minlocvelz,"|",pg(minlocvelz)%coord(1),"|",pg(minlocvelz)%coord(2),"|",pg(minlocvelz)%coord(3), &
  "||",maxvelz,"|",maxlocvelz,"|",pg(maxlocvelz)%coord(1),"|",pg(maxlocvelz)%coord(2),"|",pg(maxlocvelz)%coord(3),"|"
  if (esplosione) then
    write (nout,fmt104) &
  "  pressure     |",minpres,"|",minlocpres,"|",pg(minlocpres)%coord(1),"|",pg(minlocpres)%coord(2),"|",pg(minlocpres)%coord(3), &
  "||",maxpres,"|",maxlocpres,"|",pg(maxlocpres)%coord(1),"|",pg(maxlocpres)%coord(2),"|",pg(maxlocpres)%coord(3),"|"
    write (nout,fmt104) &
  "  Int.Energy   |",minInEn,"|",minlocInEn,"|",pg(minlocInEn)%coord(1),"|",pg(minlocInEn)%coord(2),"|", &
  pg(minlocInEn)%coord(3),"||",maxInEn,"|",maxlocInEn,"|",pg(maxlocInEn)%coord(1),"|",pg(maxlocInEn)%coord(2),"|", &
  pg(maxlocInEn)%coord(3),"|"
    write (nout,fmt104) &
  "  density      |",mindens,"|",minlocdens,"|",pg(minlocdens)%coord(1),"|",pg(minlocdens)%coord(2),"|",pg(minlocdens)%coord(3), &
  "||",maxdens,"|",maxlocdens,"|",pg(maxlocdens)%coord(1),"|",pg(maxlocdens)%coord(2),"|",pg(maxlocdens)%coord(3),"|"
  else
    write (nout,fmt102) &
  "  pressure     |",minpres,"|",minlocpres,"|",pg(minlocpres)%coord(1),"|",pg(minlocpres)%coord(2),"|",pg(minlocpres)%coord(3), &
  "||",maxpres,"|",maxlocpres,"|",pg(maxlocpres)%coord(1),"|",pg(maxlocpres)%coord(2),"|",pg(maxlocpres)%coord(3),"|"
    write (nout,fmt102) &
  "  density      |",mindens,"|",minlocdens,"|",pg(minlocdens)%coord(1),"|",pg(minlocdens)%coord(2),"|",pg(minlocdens)%coord(3), &
  "||",maxdens,"|",maxlocdens,"|",pg(maxlocdens)%coord(1),"|",pg(maxlocdens)%coord(2),"|",pg(maxlocdens)%coord(3),"|"
  end if
!AA504 sub  
  write (nout,fmt105) &
  "  viscosity    |",minvisc,"|",minlocvisc,"|",pg(minlocvisc)%coord(1),"|",pg(minlocvisc)%coord(2),"|",pg(minlocvisc)%coord(3), &
  "||",maxvisc,"|",maxlocvisc,"|",pg(maxlocvisc)%coord(1),"|",pg(maxlocvisc)%coord(2),"|",pg(maxlocvisc)%coord(3),"|"
  if (diffusione) then
    write (nout,fmt103) &
  "  coef.diff.   |",mincodi,"|",minloccodi,"|",pg(minloccodi)%coord(1),"|",pg(minloccodi)%coord(2),"|", &
  pg(minloccodi)%coord(3),"||",maxcodi,"|",maxloccodi,"|",pg(maxloccodi)%coord(1),"|",pg(maxloccodi)%coord(2),"|", &
  pg(maxloccodi)%coord(3),"|"
  end if
  write (nscr,'(128("."))') 
  write (nscr,fmt100) &
  " Print at:     | step: ",it," | time: ",tempo," | Dt: ",dt," | Particles: inside ",dummy," gone out ",OpCountot, &
  " gone in ",SpCountot
  write (nscr,fmt101) &
  " ............. |  Min. val. |Min.loc.| X coord. | Y coord. | Z coord. ||  Max. val. |Max.loc.| X coord. | Y coord. | Z coord. |"
  write (nscr,fmt101) &
  "  Tot velocity |",minvelo,"|",minlocvelo,"|",pg(minlocvelo)%coord(1),"|",pg(minlocvelo)%coord(2),"|",pg(minlocvelo)%coord(3), &
  "||",maxvelo,"|",maxlocvelo,"|",pg(maxlocvelo)%coord(1),"|",pg(maxlocvelo)%coord(2),"|",pg(maxlocvelo)%coord(3),"|"
  write (nscr,fmt101) &
  "  velocity x   |",minvelx,"|",minlocvelx,"|",pg(minlocvelx)%coord(1),"|",pg(minlocvelx)%coord(2),"|",pg(minlocvelx)%coord(3), &
  "||",maxvelx,"|",maxlocvelx,"|",pg(maxlocvelx)%coord(1),"|",pg(maxlocvelx)%coord(2),"|",pg(maxlocvelx)%coord(3),"|"
  write (nscr,fmt101) &
  "  velocity y   |",minvely,"|",minlocvely,"|",pg(minlocvely)%coord(1),"|",pg(minlocvely)%coord(2),"|",pg(minlocvely)%coord(3), &
  "||",maxvely,"|",maxlocvely,"|",pg(maxlocvely)%coord(1),"|",pg(maxlocvely)%coord(2),"|",pg(maxlocvely)%coord(3),"|"
  write (nscr,fmt101) &
  "  velocity z   |",minvelz,"|",minlocvelz,"|",pg(minlocvelz)%coord(1),"|",pg(minlocvelz)%coord(2),"|",pg(minlocvelz)%coord(3), &
  "||",maxvelz,"|",maxlocvelz,"|",pg(maxlocvelz)%coord(1),"|",pg(maxlocvelz)%coord(2),"|",pg(maxlocvelz)%coord(3),"|"
  if (esplosione) then
    write (nscr,fmt104) &
  "  pressure     |",minpres,"|",minlocpres,"|",pg(minlocpres)%coord(1),"|",pg(minlocpres)%coord(2),"|",pg(minlocpres)%coord(3), &
  "||",maxpres,"|",maxlocpres,"|",pg(maxlocpres)%coord(1),"|",pg(maxlocpres)%coord(2),"|",pg(maxlocpres)%coord(3),"|"
    write (nscr,fmt104) &
  "  Int.Energy   |",minInEn,"|",minlocInEn,"|",pg(minlocInEn)%coord(1),"|",pg(minlocInEn)%coord(2),"|",&
  pg(minlocInEn)%coord(3),"||",maxInEn,"|",maxlocInEn,"|",pg(maxlocInEn)%coord(1),"|",pg(maxlocInEn)%coord(2),"|", &
  pg(maxlocInEn)%coord(3),"|"
    write (nscr,fmt104) &
  "  density      |",mindens,"|",minlocdens,"|",pg(minlocdens)%coord(1),"|",pg(minlocdens)%coord(2),"|",pg(minlocdens)%coord(3), &
  "||",maxdens,"|",maxlocdens,"|",pg(maxlocdens)%coord(1),"|",pg(maxlocdens)%coord(2),"|",pg(maxlocdens)%coord(3),"|"
  else
    write (nscr,fmt102) &
  "  pressure     |",minpres,"|",minlocpres,"|",pg(minlocpres)%coord(1),"|",pg(minlocpres)%coord(2),"|",pg(minlocpres)%coord(3), &
  "||",maxpres,"|",maxlocpres,"|",pg(maxlocpres)%coord(1),"|",pg(maxlocpres)%coord(2),"|",pg(maxlocpres)%coord(3),"|"
    write (nscr,fmt101) &
  "  density      |",mindens,"|",minlocdens,"|",pg(minlocdens)%coord(1),"|",pg(minlocdens)%coord(2),"|",pg(minlocdens)%coord(3), &
  "||",maxdens,"|",maxlocdens,"|",pg(maxlocdens)%coord(1),"|",pg(maxlocdens)%coord(2),"|",pg(maxlocdens)%coord(3),"|"
  end if
!AA504sub
  write (nscr,fmt105) &
  "  viscosity    |",minvisc,"|",minlocvisc,"|",pg(minlocvisc)%coord(1),"|",pg(minlocvisc)%coord(2),"|",pg(minlocvisc)%coord(3), &
  "||",maxvisc,"|",maxlocvisc,"|",pg(maxlocvisc)%coord(1),"|",pg(maxlocvisc)%coord(2),"|",pg(maxlocvisc)%coord(3),"|"
  if (diffusione) then
    write (nscr,fmt103) &
  "  coef.diff.   |",mincodi,"|",minloccodi,"|",pg(minloccodi)%coord(1),"|",pg(minloccodi)%coord(2),"|",&
  pg(minloccodi)%coord(3),"||",maxcodi,"|",maxloccodi,"|",pg(maxloccodi)%coord(1),"|",pg(maxloccodi)%coord(2),"|", &
  pg(maxloccodi)%coord(3),"|"
  end if
!
!AA406 start
!AA601 sub
 if ((Domain%tipo == "bsph").and.(DBSPH%n_w > 0)) then
     write (nout,fmt101) &
  "  Wall velocity|",minvelo_w,"|",minlocvelo_w,"|",pg_w(minlocvelo_w)%coord(1),"|",pg_w(minlocvelo_w)%coord(2),"|",pg_w(minlocvelo_w)%coord(3), &
  "||",maxvelo_w,"|",maxlocvelo_w,"|",pg_w(maxlocvelo_w)%coord(1),"|",pg_w(maxlocvelo_w)%coord(2),"|",pg_w(maxlocvelo_w)%coord(3),"|"
     write (nout,fmt102) &
  "  Wall pressure|",minpres_w,"|",minlocpres_w,"|",pg_w(minlocpres_w)%coord(1),"|",pg_w(minlocpres_w)%coord(2),"|",pg_w(minlocpres_w)%coord(3), &
  "||",maxpres_w,"|",maxlocpres_w,"|",pg_w(maxlocpres_w)%coord(1),"|",pg_w(maxlocpres_w)%coord(2),"|",pg_w(maxlocpres_w)%coord(3),"|"
  endif
!AA406 end
!

!AA501b start
  if (n_bodies > 0) then
     write (nout,fmt101) &
  "Body part.vel. |",minvelo_bp,"|",minlocvelo_bp,"|",bp_arr(minlocvelo_bp)%pos(1),"|",bp_arr(minlocvelo_bp)%pos(2),"|", &
  bp_arr(minlocvelo_bp)%pos(3),"||",maxvelo_bp,"|",maxlocvelo_bp,"|",bp_arr(maxlocvelo_bp)%pos(1),"|", &
  bp_arr(maxlocvelo_bp)%pos(2),"|",bp_arr(maxlocvelo_bp)%pos(3),"|"
     write (nout,fmt102) &
  "Body part.pres.|",minpres_bp,"|",minlocpres_bp,"|",bp_arr(minlocpres_bp)%pos(1),"|",bp_arr(minlocpres_bp)%pos(2),"|", &
  bp_arr(minlocpres_bp)%pos(3),"||",maxpres_bp,"|",maxlocpres_bp,"|",bp_arr(maxlocpres_bp)%pos(1),"|", &
  bp_arr(maxlocpres_bp)%pos(2),"|",bp_arr(maxlocpres_bp)%pos(3),"|"
      write (nout,fmt101) &
  "Body velocity  |",minvelo_body,"|",minlocvelo_body,"|",body_arr(minlocvelo_body)%x_CM(1),"|", &
  body_arr(minlocvelo_body)%x_CM(2),"|",body_arr(minlocvelo_body)%x_CM(3),"||",maxvelo_body,"|",maxlocvelo_body,"|", &
  body_arr(maxlocvelo_body)%x_CM(1),"|",body_arr(maxlocvelo_body)%x_CM(2),"|",body_arr(maxlocvelo_body)%x_CM(3),"|"
      write (nout,fmt101) &
  "Body omega     |",minomega_body,"|",minlocomega_body,"|",body_arr(minlocomega_body)%x_CM(1),"|", &
  body_arr(minlocomega_body)%x_CM(2),"|",body_arr(minlocomega_body)%x_CM(3),"||",maxomega_body,"|",maxlocomega_body,"|", &
  body_arr(maxlocomega_body)%x_CM(1),"|",body_arr(maxlocomega_body)%x_CM(2),"|",body_arr(maxlocomega_body)%x_CM(3),"|"
  endif
!AA501b end

!AA504 start
  if ((Granular_flows_options%erosion_flag.ne.1).and.(Granular_flows_options%ID_erosion_criterion==1)) then
     write (nout,fmt101) &
  "tau_tauc       |",mintau_tauc,"|",minloctau_tauc,"|",pg(minloctau_tauc)%coord(1),"|",pg(minloctau_tauc)%coord(2),"|", &
  pg(minloctau_tauc)%coord(3),"||",maxtau_tauc,"|",maxloctau_tauc,"|",pg(maxloctau_tauc)%coord(1),"|", &
  pg(maxloctau_tauc)%coord(2),"|",pg(maxloctau_tauc)%coord(3),"|"
     write (nout,fmt101) &
  "u_star         |",minu_star,"|",minlocu_star,"|",pg(minlocu_star)%coord(1),"|",pg(minlocu_star)%coord(2),"|", &
  pg(minlocu_star)%coord(3),"||",maxu_star,"|",maxlocu_star,"|",pg(maxlocu_star)%coord(1),"|", &
  pg(maxlocu_star)%coord(2),"|",pg(maxlocu_star)%coord(3),"|"  
  if (Granular_flows_options%ID_erosion_criterion == 1) then       
     write (nout,fmt101) &
  "k_BetaGamma    |",mink_BetaGamma,"|",minlock_BetaGamma,"|",pg(minlock_BetaGamma)%coord(1),"|",pg(minlock_BetaGamma)%coord(2),"|", &
  pg(minlock_BetaGamma)%coord(3),"||",maxk_BetaGamma,"|",maxlock_BetaGamma,"|",pg(maxlock_BetaGamma)%coord(1),"|", &
  pg(maxlock_BetaGamma)%coord(2),"|",pg(maxlock_BetaGamma)%coord(3),"|"
  endif   
  endif
!AA504 end

!AA405 start
  else
     write (nout,'(128("."))') 
     write (nout,'(a)') "No particles inside the domain at the moment"
     write (nout,'(128("."))') 
  endif
!AA405 end
!
  it_print = it
!
  return
  end subroutine Print_Results
!---split

!cfile result_converter.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : result_converter
!
! Last updating : November 23, 2011
!
! Improvement traceback:
!
! 00  R. Guandalini  31/07/2007   Module creation
! 01  G. Agate       24/10/2007   Inserted in SPHERA
! 02  Amicarelli     23/11/2011   multiple inlet
! 03  Amicarelli/Agate 30Nov11    BSPH: wall element parameters
!AA501b comment
! 04  Amicarelli-Agate 13nov12    Body dynamics
!AA504
! 05  Amicarelli       08Apr14    (v5.04) Modifications for granular flows and 3D erosion criterion
!AA601
! 06  Amicarelli       26Jan15    DBSPH-input (AA601). New DBSPH PV output. 
!
!************************************************************************************
! Module purpose : Converts the results into VTK distributions 
!
! Calling modules: Loop_Irre_2D, Loop_Irre_3D, diagnostic
!
! Called modules : none
!      
!************************************************************************************

  subroutine result_converter (str)
!
!.. assign modules
  use FILES_ENTITIES
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Local Parameters
!  double precision, parameter    :: tol = 1.0d-04             ! tollerance
!
!.. Formal Arguments ..
  character(6),intent(IN) :: str
!
!.. Local Scalars ..
  character(len=256) :: stringa,header_string
  character(len=120) :: filevtk,prefix
  character(len=10)  :: cargo
  integer(4)         :: npi,i,j,k,k1,k2
  integer(4)         :: numcells,numpoints
  double precision   :: curtime
!
!.. Local Arrays ..
!  double precision, dimension(16)       :: appo
  integer(4), dimension(24)             :: iappo
  integer(4), dimension(:), allocatable :: finger

!.. executable statements
!
!.. check for time sampling is active
!
!!!    curtime = tempo
    curtime = tempo - MOD(tempo,abs(freq_time))  ! arrotondamento tempo corrente per animazioni
    if (curtime < val_time .and. index(str,'fine')==0) return
!    if (str == 'loop__') then
!      if (abs(val_time - curtime) > tol .and. abs(freq_time) > zero)  return
!      if (abs(val_time - curtime) > tol .and. freq_time < zero ) return
!    end if
!AA405
    if (nag > 0) then
!
    block = block + 1
    nblocchi = nblocchi + 1
    if (nblocchi > maxnumblock) then
      write (nscr,'(a)') ' ATTENZIONE !! nblocchi > maxnumblock in routine result_converter per file VTK.'
      write (nscr,'(a)') '               aumentare maxnumblock oppure diminuire frequenza memorizzazione per file VTK.'
      write (nout,'(a)') ' ATTENZIONE !! nblocchi > maxnumblock in routine result_converter per file VTK.'
      write (nout,'(a)') '               aumentare maxnumblock oppure diminuire frequenza memorizzazione per file VTK.'
      nblocchi = maxnumblock
    end if
    blocchi(nblocchi) = block
    Time_Block(nblocchi) = curtime
    prefix = nomecaso
!
!.. open the VTK formatted file VTKConverter_<casename>_<block>.vtk for the results storing
!
!!    write (cargo,'(f10.4)') curtime
!!    cargo = adjustl(cargo)
!!    i = scan(cargo,'.')
!!    cargo(i:(len_trim(cargo)-1)) = cargo(i+1:(len_trim(cargo)))
!!    filevtk = "VTKConverter_"//prefix(1:len_trim(prefix))//"_time_"//cargo(1:len_trim(cargo))//".vtk"
    write (cargo,'(i6)') block
    cargo = adjustl(cargo)
    filevtk = "VTKConverter_"//prefix(1:len_trim(prefix))//"_block_"//cargo(1:len_trim(cargo))//".vtu"
    write (nscr,'(a)') "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
    write (nout,'(a)') "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
    open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',status='unknown')
    rewind (unit=unitvtk)
!
!.. initialize the VTK formatted file
!
!.. cells are: 1 = vtk_poly_vertex including all the particle triplets in xyz cartesian system
!
    numcells = 1
!
!.. write the heading records
!
    write(unitvtk,'(a)') '<?xml version="1.0"?>'
    write(unitvtk,'(a)') &
    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
    write(unitvtk,'(a)') '<UnstructuredGrid>'
!
    write(cargo,'(i6)') int(curtime)
    cargo = adjustl(cargo)
    header_string = "case "//prefix(1:len_trim(prefix))//" * time "//cargo(1:len_trim(cargo))//" (s)"
!
!.. evaluates and write the point coordinates 
!
    numpoints = count(pg(1:nag)%cella > 0)
    allocate (finger(numpoints))
    k = 0
    do npi = 1,nag
      if (pg(npi)%cella == 0) cycle
      k = k + 1
      finger(k) = npi
    end do
!
    write(cargo,'(i6)') numpoints
    cargo = adjustl(trim(cargo))
    stringa = '  <Piece NumberOfPoints="'//cargo(1:len_trim(cargo))
    write(cargo,'(i6)') numcells
    cargo = adjustl(trim(cargo))
    stringa = stringa(1:len_trim(stringa))//'"        NumberOfCells="'//cargo(1:len_trim(cargo))//'"      >'
    write(unitvtk,'(a)') stringa(1:len_trim(stringa))
!
!.. write the coordinates of all the active particles
!
    write(unitvtk,'(a)') '    <Points>'
    write(unitvtk,'(a)') '      <DataArray type="Float32" NumberOfComponents="3" format="ascii" >'
    do i = 1,numpoints,6
      k1 = i
      k2 = k1 + 5
      if (k2 > numpoints) k2 = numpoints
      write (stringa,'(6(3(e12.5,1x)))') (pg(finger(k))%coord(1),pg(finger(k))%coord(2),pg(finger(k))%coord(3),k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </Points>'
!
!.. writes the topology 
!
!.. write the first cell: vtk_poly_vertex of type 2
!
    write(unitvtk,'(a)') '    <Cells>'
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="connectivity" format="ascii" >'
    do i = 0,numpoints-1,30
      k1 = i
      k2 = k1 + 29
      if (k2 > numpoints) k2 = numpoints
      write (stringa,'(30i8)') (k,k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the cell offsets
!
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="offsets" format="ascii" >'
    write(stringa,'(i8)') numpoints
    write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the cell types
!
    stringa = ' '
    write(unitvtk,'(a)') '      <DataArray type="UInt8" Name="types" format="ascii" >'
    stringa(1:6) = '     2'
    write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    write(unitvtk,'(a)') '      </DataArray>'
!
    write(unitvtk,'(a)') '    </Cells>'
!
!.. write the results on the VTK format file for post processing as point data
!
    write(unitvtk,'(a)') '<PointData>'
!
!AA501 rm start
!.. write the total velocity distribution
!
!      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Total velocity (m/s)" format="ascii" >'
!        do i = 1, numpoints,16
!          k1 = i
!          k2 = k1 + 15
!          if (k2 > numpoints) k2 = numpoints
!          j = 0
!          do k = k1,k2
!            j = j + 1
!            appo(j) = Dsqrt(pg(finger(k))%vel(1)*pg(finger(k))%vel(1) + pg(finger(k))%vel(2)*pg(finger(k))%vel(2) + &
!                      pg(finger(k))%vel(3)*pg(finger(k))%vel(3))
!          end do
!          write(unitvtk,'(8x,16(1x,e12.5))') (appo(k),k = 1,j)
!        end do
!      write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm end
!.. write the total velocity as vector distribution
!
      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Velocity vectors"  NumberOfComponents="3"  format="ascii" >'
        do i = 1, numpoints,6
          k1 = i
          k2 = k1 + 5
          if (k2 > numpoints) k2 = numpoints
          write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg(finger(k))%vel(1),pg(finger(k))%vel(2),pg(finger(k))%vel(3),k = k1,k2)
        end do
      write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm start
!.. write the X velocity component distribution
!
!      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="X velocity (m/s)" format="ascii" >'
!        do i = 1, numpoints,16
!          k1 = i
!          k2 = k1 + 15
!          if (k2 > numpoints) k2 = numpoints
!          write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%vel(1),k = k1,k2)
!        end do
!      write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the Y velocity component distribution
!
!      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Y velocity (m/s)" format="ascii" >'
!        do i = 1, numpoints,16
!          k1 = i
!          k2 = k1 + 15
!          if (k2 > numpoints) k2 = numpoints
!          write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%vel(2),k = k1,k2)
!        end do
!      write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the Z velocity component distribution
!
!      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Z velocity (m/s)" format="ascii" >'
!        do i = 1, numpoints,16
!          k1 = i
!          k2 = k1 + 15
!          if (k2 > numpoints) k2 = numpoints
!          write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%vel(3),k = k1,k2)
!        end do
!      write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm end
!
!.. write the pressure distribution
!
      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Pressure (Pa)" format="ascii" >'
        do i = 1, numpoints,16
          k1 = i
          k2 = k1 + 15
          if (k2 > numpoints) k2 = numpoints
          write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%pres,k = k1,k2)
        end do
      write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the density distribution
!
      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Density (kg/mc)" format="ascii" >'
        do i = 1, numpoints,16
          k1 = i
          k2 = k1 + 15
          if (k2 > numpoints) k2 = numpoints
          write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%dens,k = k1,k2)
        end do
      write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the mass distribution
!
      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Mass (kg)" format="ascii" >'
        do i = 1, numpoints,16
          k1 = i
          k2 = k1 + 15
          if (k2 > numpoints) k2 = numpoints
          write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%mass,k = k1,k2)
        end do
      write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the viscosity distribution
!
      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Viscosity (mq/s)" format="ascii" >'
        do i = 1, numpoints,16
          k1 = i
          k2 = k1 + 15
          if (k2 > numpoints) k2 = numpoints
          write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%visc,k = k1,k2)
        end do
      write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the Volume Fraction distribution
!
      if (diffusione) then
        write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Volume Fraction (mc/mc)" format="ascii" >'
        do i = 1, numpoints,16
          k1 = i
          k2 = k1 + 15
          if (k2 > numpoints) k2 = numpoints
          write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%VolFra,k = k1,k2)
        end do
        write(unitvtk,'(a)') '      </DataArray>'
      end if
!
!.. write the Internal Energy distribution
!
      if (esplosione) then
        write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Internal Energy (J/kg)" format="ascii" >'
        do i = 1, numpoints,16
          k1 = i
          k2 = k1 + 15
          if (k2 > numpoints) k2 = numpoints
          write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%IntEn,k = k1,k2)
        end do
        write(unitvtk,'(a)') '      </DataArray>'
      end if
!
!.. write the medium distribution
!
      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Medium" format="ascii" >'
        do i = 1, numpoints,24
          k1 = i
          k2 = k1 + 23
          if (k2 > numpoints) k2 = numpoints
          write(unitvtk,'(8x,24(1x,i6))') (pg(finger(k))%imed,k = k1,k2)
        end do
      write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the particle state distribution
!
      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Status" format="ascii" >'
        do i = 1, numpoints,24
          k1 = i
          k2 = k1 + 23
          if (k2 > numpoints) k2 = numpoints
          iappo = 0
          j = 0
          do k = k1,k2
            j = j + 1
            if ( (pg(finger(k))%state == 'flu') .and. (index(Med(pg(finger(k))%imed)%tipo,"liquid") > 0) ) then
              iappo(k-k1+1) = 1
            else if ( (pg(finger(k))%state == 'flu') .and. (index(Med(pg(finger(k))%imed)%tipo,"granular") > 0) ) then
              iappo(k-k1+1) = 2
            else if ( (pg(finger(k))%state == 'sol') .and. (index(Med(pg(finger(k))%imed)%tipo,"granular") > 0) ) then
              iappo(k-k1+1) = 3
            else if ( (pg(finger(k))%state == 'flu') .and. (index(Med(pg(finger(k))%imed)%tipo,"gas") > 0) ) then
              iappo(k-k1+1) = 4
            end if
          end do
          write(unitvtk,'(8x,24(1x,i6))') (iappo(k),k = 1,j)
        end do
      write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the particle position in the pg array
!
      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Finger" format="ascii" >'
        do i = 1, numpoints,24
          k1 = i
          k2 = k1 + 23
          if (k2 > numpoints) k2 = numpoints
          write(unitvtk,'(8x,24(1x,i6))') (finger(k),k = k1,k2)
        end do
      write(unitvtk,'(a)') '      </DataArray>'
!
!AA406 start
 if (Domain%tipo == "bsph") then
!Writing Shepard's coefficient
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Shepard coefficient" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%uni,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!Writing the discrete Shepard's coefficient
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="discrete Shepard" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%sigma,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!Writing the integral Shepard's coefficient
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="integral Shepard" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%Gamma,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!Writing the free surface flag
    write(unitvtk,'(a)') '      <DataArray type="UInt8" Name="free surface" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,i8))') (pg(finger(k))%FS,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm start
!
!Writing the absolute value of the density gradient
!    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Density gradient (kg/(m^4))" format="ascii" >'
!    do i = 1, numpoints,16
!       k1 = i
!       k2 = k1 + 15
!       if (k2 > numpoints) k2 = numpoints
!       j = 0
!       do k = k1,k2
!          j = j + 1
!          appo(j) = Dsqrt(pg(finger(k))%drho(1)*pg(finger(k))%drho(1) + pg(finger(k))%drho(2)*pg(finger(k))%drho(2) + &
!                          pg(finger(k))%drho(3)*pg(finger(k))%drho(3))
!       end do
!       write(unitvtk,'(8x,16(1x,e12.5))') (appo(k),k = 1,j)
!    end do
!    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm end
!
!Writing the vector of the density gradient
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="drho vectors"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg(finger(k))%drho(1),pg(finger(k))%drho(2),pg(finger(k))%drho(3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm start
!
!Writing the absolute value of the gradient of the x-component of velocity (d(u_x))
!    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="d(u_x) (1/(s))" format="ascii" >'
!    do i = 1, numpoints,16
!       k1 = i
!       k2 = k1 + 15
!       if (k2 > numpoints) k2 = numpoints
!       j = 0
!       do k = k1,k2
!          j = j + 1
!          appo(j) = Dsqrt(pg(finger(k))%dvel(1,1)*pg(finger(k))%dvel(1,1) + pg(finger(k))%dvel(1,2)*pg(finger(k))%dvel(1,2) + &
!                          pg(finger(k))%dvel(1,3)*pg(finger(k))%dvel(1,3))
!       end do
!       write(unitvtk,'(8x,16(1x,e12.5))') (appo(k),k = 1,j)
!    end do
!    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm end
!
!Writing the vector of the gradient of the x-component of velocity
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="d(u_x) vectors"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg(finger(k))%dvel(1,1),pg(finger(k))%dvel(1,2),pg(finger(k))%dvel(1,3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    if (ncord==3) then
!
!AA501 rm start
!
!Writing the absolute value of the gradient of the y-component of velocity (d(u_y))
!    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="d(u_y) (1/(s))" format="ascii" >'
!    do i = 1, numpoints,16
!       k1 = i
!       k2 = k1 + 15
!       if (k2 > numpoints) k2 = numpoints
!       j = 0
!       do k = k1,k2
!          j = j + 1
!          appo(j) = Dsqrt(pg(finger(k))%dvel(2,1)*pg(finger(k))%dvel(2,1) + pg(finger(k))%dvel(2,2)*pg(finger(k))%dvel(2,2) + &
!                          pg(finger(k))%dvel(2,3)*pg(finger(k))%dvel(2,3))
!       end do
!       write(unitvtk,'(8x,16(1x,e12.5))') (appo(k),k = 1,j)
!    end do
!    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm end
!
!Writing the vector of the of the gradient of the y-component of velocity
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="d(u_y) vectors"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg(finger(k))%dvel(2,1),pg(finger(k))%dvel(2,2),pg(finger(k))%dvel(2,3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    endif
!
!AA501 rm start
!
!Writing the absolute value of the gradient of the z-component of velocity (d(u_z))
!    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="d(u_z) (1/(s))" format="ascii" >'
!    do i = 1, numpoints,16
!       k1 = i
!       k2 = k1 + 15
!       if (k2 > numpoints) k2 = numpoints
!       j = 0
!       do k = k1,k2
!          j = j + 1
!          appo(j) = Dsqrt(pg(finger(k))%dvel(3,1)*pg(finger(k))%dvel(3,1) + pg(finger(k))%dvel(3,2)*pg(finger(k))%dvel(3,2) + &
!                          pg(finger(k))%dvel(3,3)*pg(finger(k))%dvel(3,3))
!       end do
!       write(unitvtk,'(8x,16(1x,e12.5))') (appo(k),k = 1,j)
!    end do
!    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm end
!
!Writing the vector of the of the gradient of the z-component of velocity
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="d(u_z) vectors"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg(finger(k))%dvel(3,1),pg(finger(k))%dvel(3,2),pg(finger(k))%dvel(3,3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    endif ! Domain%tipo == "bsph"
!AA406 end

!AA504 start
    if (Granular_flows_options%ID_erosion_criterion>=1) then
!Writing sigma_prime
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="sigma_prime" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%sigma_prime,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!Writing sec_inv
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="sqrt_I2_eij" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%secinv,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!Writing blt_flag
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="blt_flag" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,i8))') (pg(finger(k))%blt_flag,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'  
!Writing z-coordinate for top view 2D field of free surface
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="z(m)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%coord(3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>' 
!Writing Bingham number
    if (Granular_flows_options%viscosity_blt_formula==4) then
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Bn" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%Bn,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    endif
    if (Granular_flows_options%erosion_flag.ne.1) then    
!Writing Beta
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Beta(radians)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%Beta_slope,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>' 
!Writing the vector normal_int
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="normal_int vector"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg(finger(k))%normal_int(1),pg(finger(k))%normal_int(2),pg(finger(k))%normal_int(3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!Writing tau_tauc
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="tau_tauc" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%tau_tauc,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!Writing u_star
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="u_star(m/s)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%u_star,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'  
    if (Granular_flows_options%ID_erosion_criterion == 1) then  
    if (ncord==3) then    
!Writing C_D
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="C_D" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%C_D,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>' 
!Writing C_L
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="C_L" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%C_L,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!Writing Gamma
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Gamma(radians)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%Gamma_slope,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>' 
    endif
!Writing k_BetaGamma
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="k_BetaGamma" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%k_BetaGamma,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'  
    endif
    endif
    endif 
!AA504 end

    write(unitvtk,'(a)') '    </PointData>'
!.. write the distribution data on cells
!
    write(unitvtk,'(a)') '    <CellData>'
    write(unitvtk,'(a)') '    </CellData>'
!
!.. close the result VTK file
!
    write(unitvtk,'(a)') '  </Piece>'
    write(unitvtk,'(a)') '</UnstructuredGrid>'
    write(unitvtk,'(a)') '</VTKFile>'
!
!.. Flush of file contents is forced
!
    flush(unitvtk)
    close (unitvtk)
    deallocate (finger)
!
!AA406 rm !!! (moved later on)
!AA405
!    endif  !nag>0
!
!AA406 start
    if ((DBSPH%n_w > 0) .and. (Domain%tipo == "bsph")) then
!
! Open the VTK unstructured grid formatted file VTKConverter_<casename>_wall_<block>.vtk for the results storing
    write (cargo,'(i6)') block
    cargo = adjustl(cargo)
    filevtk = "VTKConverter_"//prefix(1:len_trim(prefix))//"_block_wall_"//cargo(1:len_trim(cargo))//".vtu"
    write (nscr,'(a)') "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
    write (nout,'(a)') "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
    open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',status='unknown')
    rewind (unit=unitvtk)
! Initialization of the VTK formatted file (1 cell formally groupes all the points)
! Writing the heading records
    write(unitvtk,'(a)') '<?xml version="1.0"?>'
    write(unitvtk,'(a)') &
    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
    write(unitvtk,'(a)') '<UnstructuredGrid>'
    write(cargo,'(i6)') int(curtime)
    cargo = adjustl(cargo)
    header_string = "case "//prefix(1:len_trim(prefix))//" * time "//cargo(1:len_trim(cargo))//" (s)"
! Evaluating and writing the point coordinates 
    numpoints = count(pg_w(1:DBSPH%n_w)%cella > 0)
    allocate (finger(numpoints))
    k = 0
    do npi = 1,DBSPH%n_w
      if (pg_w(npi)%cella == 0) cycle
      k = k + 1
      finger(k) = npi
    end do
    write(cargo,'(i6)') numpoints
    cargo = adjustl(trim(cargo))
    stringa = '  <Piece NumberOfPoints="'//cargo(1:len_trim(cargo))
    write(cargo,'(i6)') numcells
    cargo = adjustl(trim(cargo))
    stringa = stringa(1:len_trim(stringa))//'"        NumberOfCells="'//cargo(1:len_trim(cargo))//'"      >'
    write(unitvtk,'(a)') stringa(1:len_trim(stringa))
! Writing the coordinates of all the particles
    write(unitvtk,'(a)') '    <Points>'
    write(unitvtk,'(a)') '      <DataArray type="Float32" NumberOfComponents="3" format="ascii" >'
    do i = 1,numpoints,6
      k1 = i
      k2 = k1 + 5
      if (k2 > numpoints) k2 = numpoints
      write (stringa,'(6(3(e12.5,1x)))') (pg_w(finger(k))%coord(1),pg_w(finger(k))%coord(2),pg_w(finger(k))%coord(3),k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </Points>'
! Writing the topology 
! Writing the first cell: vtk_poly_vertex of type 2
    write(unitvtk,'(a)') '    <Cells>'
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="connectivity" format="ascii" >'
    do i = 0,numpoints-1,30
      k1 = i
      k2 = k1 + 29
      if (k2 > numpoints) k2 = numpoints
      write (stringa,'(30i8)') (k,k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing the cell offset
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="offsets" format="ascii" >'
    write(stringa,'(i8)') numpoints
    write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    write(unitvtk,'(a)') '      </DataArray>'
! Writing the cell type
    stringa = ' '
    write(unitvtk,'(a)') '      <DataArray type="UInt8" Name="types" format="ascii" >'
    stringa(1:6) = '     2'
    write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </Cells>'
! Writing the results on the VTK format file for post processing as point data
    write(unitvtk,'(a)') '<PointData>'
!
!AA501 rm start
!
! Writing the absolute value of velocity
!    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Velocity (m/s)" format="ascii" >'
!    do i = 1, numpoints,16
!       k1 = i
!       k2 = k1 + 15
!       if (k2 > numpoints) k2 = numpoints
!       j = 0
!       do k = k1,k2
!          j = j + 1
!          appo(j) = Dsqrt(pg_w(finger(k))%vel(1)*pg_w(finger(k))%vel(1) + pg_w(finger(k))%vel(2)*pg_w(finger(k))%vel(2) + &
!                    pg_w(finger(k))%vel(3)*pg_w(finger(k))%vel(3))
!       end do
!       write(unitvtk,'(8x,16(1x,e12.5))') (appo(k),k = 1,j)
!    end do
!    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm end
!
! Writing the velocity vector 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Velocity vectors"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg_w(finger(k))%vel(1),pg_w(finger(k))%vel(2),pg_w(finger(k))%vel(3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm start
!
! Writing the absolute value of the normal (it should be the unity)
!    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="normal norm" format="ascii" >'
!    do i = 1, numpoints,16
!       k1 = i
!       k2 = k1 + 15
!       if (k2 > numpoints) k2 = numpoints
!       j = 0
!       do k = k1,k2
!          j = j + 1
!          appo(j) = Dsqrt(pg_w(finger(k))%normal(1)*pg_w(finger(k))%normal(1) + pg_w(finger(k))%normal(2)*pg_w(finger(k))%normal(2) + &
!                    pg_w(finger(k))%normal(3)*pg_w(finger(k))%normal(3))
!       end do
!       write(unitvtk,'(8x,16(1x,e12.5))') (appo(k),k = 1,j)
!    end do
!    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm end
!
! Writing the normal vectors n
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Normal vectors"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg_w(finger(k))%normal(1),pg_w(finger(k))%normal(2),pg_w(finger(k))%normal(3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm start
!
! Writing n_x
!    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="n_x " format="ascii" >'
!    do i = 1, numpoints,16
!       k1 = i
!       k2 = k1 + 15
!       if (k2 > numpoints) k2 = numpoints
!       write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%normal(1),k = k1,k2)
!    end do
!    write(unitvtk,'(a)') '      </DataArray>'
! Writing n_y
!    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="n_y " format="ascii" >'
!    do i = 1, numpoints,16
!       k1 = i
!       k2 = k1 + 15
!       if (k2 > numpoints) k2 = numpoints
!       write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%normal(2),k = k1,k2)
!    end do
!    write(unitvtk,'(a)') '      </DataArray>'
! Writing n_z
!    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="n_z " format="ascii" >'
!    do i = 1, numpoints,16
!       k1 = i
!      k2 = k1 + 15
!      if (k2 > numpoints) k2 = numpoints
!      write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%normal(3),k = k1,k2)
!   end do
!    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm end
!
! Writing pressure
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Pressure (Pa)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%pres,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing weight
    if (ncord == 2) then
       write(unitvtk,'(a)') '      <DataArray type="Float32" Name="weight (m)" format="ascii" >'
       else
          write(unitvtk,'(a)') '      <DataArray type="Float32" Name="weight (m^2)" format="ascii" >'
    endif
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%weight,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing semi-particle mass 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="semi-particle mass (kg)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%mass,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!AA601 start
! Writing semi-particle k_d 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="semi-particle k_d" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%k_d,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="semi-particle volume (m^3)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%volume,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'    
!AA601 end
! Write the particle position in the pg_w array
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Finger" format="ascii" >'
    do i = 1, numpoints,24
       k1 = i
       k2 = k1 + 23
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,24(1x,i6))') (finger(k),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </PointData>'
! Writing the distribution data on cells
    write(unitvtk,'(a)') '    <CellData>'
    write(unitvtk,'(a)') '    </CellData>'
! Closing the VTK file
    write(unitvtk,'(a)') '  </Piece>'
    write(unitvtk,'(a)') '</UnstructuredGrid>'
    write(unitvtk,'(a)') '</VTKFile>'
! Flushing the unit explicitly (FORTRAN 2003)
    flush(unitvtk)
    close (unitvtk)
    deallocate (finger)
!   
    endif  ! DBSPH%n_w>0 and Domain%tipo == "bsph"
!

!AA501b start
!Body dynamics option
    if (n_bodies > 0) then

! Body particles
! Open the VTK unstructured grid formatted file VTKConverter_<casename>_body-part_<block>.vtk for the results storing
    write (cargo,'(i6)') block
    cargo = adjustl(cargo)
    filevtk = "VTKConverter_"//prefix(1:len_trim(prefix))//"_block_body-part_"//cargo(1:len_trim(cargo))//".vtu"
    write (nscr,'(a)') "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
    write (nout,'(a)') "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
    open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',status='unknown')
    rewind (unit=unitvtk)
! Initialization of the VTK formatted file (1 cell formally groupes all the points)
! Writing the heading records
    write(unitvtk,'(a)') '<?xml version="1.0"?>'
    write(unitvtk,'(a)') &
    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
    write(unitvtk,'(a)') '<UnstructuredGrid>'
    write(cargo,'(i6)') int(curtime)
    cargo = adjustl(cargo)
    header_string = "case "//prefix(1:len_trim(prefix))//" * time "//cargo(1:len_trim(cargo))//" (s)"
! Evaluating and writing the point coordinates 
    numpoints = count(bp_arr(1:n_body_part)%cell > 0)
    allocate (finger(numpoints))
    k = 0
    do npi = 1,n_body_part
      if (bp_arr(npi)%cell == 0) cycle
      k = k + 1
      finger(k) = npi
    end do
    write(cargo,'(i6)') numpoints
    cargo = adjustl(trim(cargo))
    stringa = '  <Piece NumberOfPoints="'//cargo(1:len_trim(cargo))
    write(cargo,'(i6)') numcells
    cargo = adjustl(trim(cargo))
    stringa = stringa(1:len_trim(stringa))//'"        NumberOfCells="'//cargo(1:len_trim(cargo))//'"      >'
    write(unitvtk,'(a)') stringa(1:len_trim(stringa))
! Writing the coordinates of all the particles
    write(unitvtk,'(a)') '    <Points>'
    write(unitvtk,'(a)') '      <DataArray type="Float32" NumberOfComponents="3" format="ascii" >'
    do i = 1,numpoints,6
      k1 = i
      k2 = k1 + 5
      if (k2 > numpoints) k2 = numpoints
      write (stringa,'(6(3(e12.5,1x)))') (bp_arr(finger(k))%pos(1),bp_arr(finger(k))%pos(2),bp_arr(finger(k))%pos(3),k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </Points>'
! Writing the topology 
! Writing the first cell: vtk_poly_vertex of type 2
    write(unitvtk,'(a)') '    <Cells>'
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="connectivity" format="ascii" >'
    do i = 0,numpoints-1,30
      k1 = i
      k2 = k1 + 29
      if (k2 > numpoints) k2 = numpoints
      write (stringa,'(30i8)') (k,k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing the cell offset
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="offsets" format="ascii" >'
    write(stringa,'(i8)') numpoints
    write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    write(unitvtk,'(a)') '      </DataArray>'
! Writing the cell type
    stringa = ' '
    write(unitvtk,'(a)') '      <DataArray type="UInt8" Name="types" format="ascii" >'
    stringa(1:6) = '     2'
    write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </Cells>'
! Writing the results on the VTK format file for post processing as point data
    write(unitvtk,'(a)') '<PointData>'
! Writing the velocity vector 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Velocity vectors"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (bp_arr(finger(k))%vel(1),bp_arr(finger(k))%vel(2),bp_arr(finger(k))%vel(3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing the normal vectors n
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Normal vectors"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (bp_arr(finger(k))%normal(1),bp_arr(finger(k))%normal(2), &
                                             bp_arr(finger(k))%normal(3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing pressure
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Pressure (Pa)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (bp_arr(finger(k))%pres,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing mass 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="mass (kg)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (bp_arr(finger(k))%mass,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing area 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="area(m^2)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (bp_arr(finger(k))%area,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Write the particle position in the bp_arr array
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="Finger" format="ascii" >'
    do i = 1, numpoints,24
       k1 = i
       k2 = k1 + 23
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,24(1x,i6))') (finger(k),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Write the body identifier 
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="Body" format="ascii" >'
    do i = 1, numpoints,24
       k1 = i
       k2 = k1 + 23
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,24(1x,i6))') (bp_arr(finger(k))%body,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </PointData>'
! Writing the distribution data on cells
    write(unitvtk,'(a)') '    <CellData>'
    write(unitvtk,'(a)') '    </CellData>'
! Closing the VTK file
    write(unitvtk,'(a)') '  </Piece>'
    write(unitvtk,'(a)') '</UnstructuredGrid>'
    write(unitvtk,'(a)') '</VTKFile>'
! Flushing the unit explicitly (FORTRAN 2003)
    flush(unitvtk)
    close (unitvtk)
    deallocate (finger)

! Bodies
! Open the VTK unstructured grid formatted file VTKConverter_<casename>_body_<block>.vtk for the results storing
    write (cargo,'(i6)') block
    cargo = adjustl(cargo)
    filevtk = "VTKConverter_"//prefix(1:len_trim(prefix))//"_block_body_"//cargo(1:len_trim(cargo))//".vtu"
    write (nscr,'(a)') "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
    write (nout,'(a)') "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
    open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',status='unknown')
    rewind (unit=unitvtk)
! Initialization of the VTK formatted file (1 cell formally groupes all the points)
! Writing the heading records
    write(unitvtk,'(a)') '<?xml version="1.0"?>'
    write(unitvtk,'(a)') &
    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
    write(unitvtk,'(a)') '<UnstructuredGrid>'
    write(cargo,'(i6)') int(curtime)
    cargo = adjustl(cargo)
    header_string = "case "//prefix(1:len_trim(prefix))//" * time "//cargo(1:len_trim(cargo))//" (s)"
! Evaluating and writing the point coordinates 
    numpoints = n_bodies
    allocate (finger(numpoints))
    k = 0
    do npi = 1,n_bodies
      k = k + 1
      finger(k) = npi
    end do
    write(cargo,'(i6)') numpoints
    cargo = adjustl(trim(cargo))
    stringa = '  <Piece NumberOfPoints="'//cargo(1:len_trim(cargo))
    write(cargo,'(i6)') numcells
    cargo = adjustl(trim(cargo))
    stringa = stringa(1:len_trim(stringa))//'"        NumberOfCells="'//cargo(1:len_trim(cargo))//'"      >'
    write(unitvtk,'(a)') stringa(1:len_trim(stringa))
! Writing the coordinates of all the particles
    write(unitvtk,'(a)') '    <Points>'
    write(unitvtk,'(a)') '      <DataArray type="Float32" NumberOfComponents="3" format="ascii" >'
    do i = 1,numpoints,6
      k1 = i
      k2 = k1 + 5
      if (k2 > numpoints) k2 = numpoints
      write (stringa,'(6(3(e12.5,1x)))') (body_arr(finger(k))%x_CM(1),body_arr(finger(k))%x_CM(2), &
                                          body_arr(finger(k))%x_CM(3),k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </Points>'
! Writing the topology 
! Writing the first cell: vtk_poly_vertex of type 2
    write(unitvtk,'(a)') '    <Cells>'
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="connectivity" format="ascii" >'
    do i = 0,numpoints-1,30
      k1 = i
      k2 = k1 + 29
      if (k2 > numpoints) k2 = numpoints
      write (stringa,'(30i8)') (k,k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing the cell offset
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="offsets" format="ascii" >'
    write(stringa,'(i8)') numpoints
    write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    write(unitvtk,'(a)') '      </DataArray>'
! Writing the cell type
    stringa = ' '
    write(unitvtk,'(a)') '      <DataArray type="UInt8" Name="types" format="ascii" >'
    stringa(1:6) = '     2'
    write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </Cells>'
! Writing the results on the VTK format file for post processing as point data
    write(unitvtk,'(a)') '<PointData>'
! Writing the velocity vector 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Velocity vectors"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (body_arr(finger(k))%u_CM(1),body_arr(finger(k))%u_CM(2), &
                                             body_arr(finger(k))%u_CM(3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing the angular velocity 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Angular velocity"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (body_arr(finger(k))%omega(1),body_arr(finger(k))%omega(2), &
                                             body_arr(finger(k))%omega(3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing the moment of inertia
! first row 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Ic(x,:)"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (body_arr(finger(k))%Ic(1,1),body_arr(finger(k))%Ic(1,2), &
                                             body_arr(finger(k))%Ic(1,3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! second row 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Ic(y,:)"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (body_arr(finger(k))%Ic(2,1),body_arr(finger(k))%Ic(2,2), &
                                             body_arr(finger(k))%Ic(2,3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! third row 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Ic(z,:)"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (body_arr(finger(k))%Ic(3,1),body_arr(finger(k))%Ic(3,2), &
                                             body_arr(finger(k))%Ic(3,3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing mass 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="mass (kg)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (body_arr(finger(k))%mass,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Write the particle position in the body_arr array
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="Finger" format="ascii" >'
    do i = 1, numpoints,24
       k1 = i
       k2 = k1 + 23
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,24(1x,i6))') (finger(k),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </PointData>'
! Writing the distribution data on cells
    write(unitvtk,'(a)') '    <CellData>'
    write(unitvtk,'(a)') '    </CellData>'
! Closing the VTK file
    write(unitvtk,'(a)') '  </Piece>'
    write(unitvtk,'(a)') '</UnstructuredGrid>'
    write(unitvtk,'(a)') '</VTKFile>'
! Flushing the unit explicitly (FORTRAN 2003)
    flush(unitvtk)
    close (unitvtk)
    deallocate (finger)
 
    endif  
!AA501b end 
    
    endif  !nag>0
!AA406 end    

!.. increase the time sampling if active
!
     if (freq_time > zero) then
       val_time = val_time + abs(freq_time)
!
!AA404 
! truncation for compatibility with "curtime"
       val_time = val_time - MOD(val_time,abs(freq_time))
       if (val_time==curtime) val_time = val_time+abs(freq_time)
!
     else if (curtime >= val_time .and. val_time > zero) then
       val_time = 1.e20
     end if
!
!!  end do listing_loop
!
  return
  end subroutine result_converter
!---split

!cfile rundt2.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : rundt2
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Agate/Amicarelli  2011           Introduction CFL and substitution of COTE
!AA504
! 04  Amicarelli        08Apr14        (v5.04) Modifications for Monaghan viscosity and management of low-velocity SPH granular particles
!
!************************************************************************************
! Module purpose : Module to calculate the delta t 
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
!
!AA401 sub start
  subroutine rundt2
! Computation of the time step  according to 3 conditions (CFL,viscosity,interface diffusion)
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
  integer(4)       :: j,npi,mate,ii
  double precision :: dtmin,dt_CFL,dt_dif,dt_vis,diffmax,U,celiq
!AA504  
  double precision :: visc_Mon_j,dt_Mon_j,dt_Mon
  
!
!.. Executable Statements ..
!
!.. initializations
! 
  dtmin   = 1.0d+30                                                        
  diffmax = zero
!AA504
  dt_Mon = max_positive_number 
  
!
!AA401
!.. loop to compute the time step, according to 3 conditions:
!.. 1) the CFL condition: dt_CFL<min(2h/(c+U))
!.. 2) viscous stability condition dt_vis<min(rho*h^2/(0.5*mu)) 
!.. 3) interface diffusion condition (h^2/2*teta)
  do ii = 1,indarrayFlu
     npi = Array_Flu(ii)
!AA504 start
     if (Granular_flows_options%ID_erosion_criterion==1) then
         if (pg(npi)%state=="sol") cycle
         if ( (pg(npi)%coord(1)<Granular_flows_options%x_min_dt) .or. (pg(npi)%coord(1)>Granular_flows_options%x_max_dt) .or. &
              (pg(npi)%coord(2)<Granular_flows_options%y_min_dt) .or. (pg(npi)%coord(2)>Granular_flows_options%y_max_dt) .or. &
              (pg(npi)%coord(3)<Granular_flows_options%z_min_dt) .or. (pg(npi)%coord(3)>Granular_flows_options%z_max_dt) ) then
             cycle
         endif    
     endif    
!AA504 end
     mate = pg(npi)%imed
     U = sqrt(pg(npi)%vel(1)**2+pg(npi)%vel(2)**2+pg(npi)%vel(3)**2)
     dt_CFL = 2.*Domain%h/(Med(mate)%celerita+U)
!     dt_vis = pg(npi)%dens*Domain%h**2/(0.5*pg(npi)%mu)
     dt_vis = pg(npi)%dens*squareh/(half*pg(npi)%mu)
     dtmin  = min(dtmin,dt_CFL,dt_vis)
     celiq  = Med(mate)%eps / pg(npi)%dens
     if ( celiq >= zero ) pg(npi)%Csound = Dsqrt(celiq)
   end do
   do j = 1,nmedium
       diffmax  = max(diffmax,Med(j)%codif)
!AA601
       if (dt_alfa_Mon == .true.) then
!AA504 start
          visc_Mon_j = Med(j)%alfaMon * Med(j)%celerita * Domain%h / Med(j)%den0
          dt_Mon_j = squareh/(half*visc_Mon_j)
          dt_Mon = min(dt_Mon,dt_Mon_j)
!AA504 end
!AA601
       endif
   end do
   dt_dif = half * squareh / (diffmax+0.000000001d0)
   dtmin  = min(dtmin,dt_dif)
!AA601 sub
   if (dt_alfa_Mon == .true.) dtmin = min(dtmin,dt_Mon)
  
!.. evaluates the time step current value applying the time coefficient to the found value
!AA401 sub end
!
!AA401 sub 
! CFL is used as a constant for every condition (the CFL, the viscous and the diffusive one)
  dt = (one - pesodt) * Domain%CFL * dtmin + pesodt * dt_average
!$$$$$$$$$$$$$$ assegnazione provvisoria dt costante per test $$$$$$$$$$$$$$$$$$$$$$$$$
!
!  dt = 0.005
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
  dt_average = (dt_average * (it_corrente - 1) + dt) / it_corrente
!
  return
!
!!!!!
!!!!  subroutine rundt2
!!!!!
!!!!!.. Computation of the time step  according to 3 conditions (CFL,viscosity,interface diffusion)
!!!!!Calcola il passo dt secondo la procedura euristica proposta da AdM
!!!!!
!!!!!.. assign modules
!!!!  use GLOBAL_MODULE
!!!!  use AdM_USER_TYPE
!!!!  use ALLOC_MODULE
!!!!!
!!!!!.. Implicit Declarations ..
!!!!  implicit none
!!!!!
!!!!!.. Local Parameters ..
!!!!double precision,parameter :: k1 = 2.33333333333333d0  ! 7/3
!!!!double precision,parameter :: Skmaxq = 0.55d0 * 0.55d0
!!!!double precision,parameter :: Ckmax = 1.33333333333333d0  ! 4/3
!!!!!
!!!!!.. Local Scalars ..
!!!!  integer(4)       :: npi, mate, ii
!!!!  double precision :: celiq, celi, viscequi, Denom, termB, termC, dti, dtmin
!!!!  double precision :: dtdiff, difflim
!!!!!
!!!!!.. Executable Statements ..
!!!!!
!!!!!.. initializations
!!!!!
!!!!  dtmin = 1.0d+30
!!!!  difflim  = -1000.d0   !.. modello diffusione
!!!!!
!!!!!.. loops on all the particles
!!!!!
!!!!!!!  do npi = 1,nag
!!!!!!!!
!!!!!!!!!    if ( pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" ) cycle
!!!!!!!    if ( pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" .or. pg(npi)%state == "sol") cycle
!!!!!$$$$$
!!!!  do ii = 1,indarrayFlu
!!!!    npi = Array_Flu(ii)
!!!!!$$$$$
!!!!!
!!!!    mate     = pg(npi)%imed
!!!!    celiq    = Med(mate)%eps / pg(npi)%dens
!!!!!
!!!!    if ( celiq < zero ) cycle    
!!!!!
!!!!    celi     = Dsqrt(celiq)
!!!!!
!!!!!........................................... 2011 mar 15
!!!!!    if (esplosione) then
!!!!!      celiq = Med(pg(npi)%imed)%gamma * (Med(pg(npi)%imed)%gamma - one) * pg(npi)%IntEn
!!!!!      pg(npi)%Csound = Dsqrt(Med(pg(npi)%imed)%gamma * (Med(pg(npi)%imed)%gamma - one) * pg(npi)%IntEn)
!!!!!      celi  = Dsqrt(celiq)
!!!!!    end if
!!!!    pg(npi)%Csound = celi
!!!!!!........................................... 2011 mar 15
!!!!!
!!!!!!!!    termB    = Domain%h / celi                             ! CFL condition
!!!!!!!!    termC    = Domain%h * Domain%h / (two * pg(npi)%visc)    ! viscous condition
!!!!!!!!    difflim  = max(difflim,pg(npi)%coefdif)
!!!!!!!!    dtmin    = min(dtmin,termB,termC)
!!!!!
!!!!    viscequi = Med(mate)%alfaMon * celi * Domain%h + k1 * pg(npi)%visc
!!!!    Denom    = Skmaxq * celiq
!!!!    termB    = Ckmax * viscequi / (Denom + Denom)
!!!!    termC    = squareh / Denom
!!!!    dti      = Dsqrt(termB * termB + termC) - termB
!!!!    difflim  = max(difflim,pg(npi)%coefdif)
!!!!    dtmin    = min(dti,dtmin)
!!!!!
!!!!  end do
!!!!!
!!!!!.. is the diffusion model
!!!!  if (diffusione) then
!!!!!    dtdiff   = 0.125d0 * squareh * Domain%cote /(difflim+0.00000001)
!!!!    dtdiff   = 0.125d0 * squareh * Domain%CFL /(difflim+0.00000001)
!!!!    dtmin    = min(dtmin,dtdiff)
!!!!  end if
!!!!!
!!!!!.. riduzione sperimentale del dtmin per geometrie 3D e per condizioni al contorno.
!!!!!!!!!  dtmin = dtmin * half
!!!!!
!!!!    
!!!!!.. controllo DTmin reazione elastica boundary
!!!!  dtmin = min (dtmin, DTminBER)
!!!!
!!!!!.. set the new time step value
!!!!!
!!!!!  dt = (one - pesodt) * Domain%cote * dtmin + pesodt * dt_average
!!!!  dt = (one - pesodt) * Domain%CFL * dtmin + pesodt * dt_average
!!!!!
!!!!!$$$$$$$$$$$$$$ assegnazione provvisoria dt costante per test $$$$$$$$$$$$$$$$$$$$$$$$$
!!!!!
!!!!!  dt = 0.005
!!!!!
!!!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!!!!
!!!!  dt_average = (dt_average * (it_corrente - 1) + dt) / it_corrente
!!!!!
!!!!  return
!!!!
  end subroutine rundt2
!---split

!cfile SearchforParticleZone_3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : SearchforParticleZone_3D
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
! Module purpose : Module to search the highest pointer of the fluid
!
! Calling routine: PreSourceParticles_3D, Gest_Input
!
! Called routines: LocalNormalCoordinates
!
!************************************************************************************
!
  subroutine SearchforParticleZone_3D (partizone)   

!Restituisce in 'partizone' l'indice di zona più alto cui corrisponde un volume fluido (di particelle)
!Nel caso non esista una tale volume (zona) fluido assegna partzone = sourzone (zona dellasorgente)
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
  integer(4),intent(INOUT)   :: partizone
!
!.. Local Scalars ..
  integer(4)  :: iz, sourzone, mate
  character(4) :: tipo
!
!.. Executable Statements ..
!
  partizone = 0
  sourzone = 0
!
!.. 
  do iz = NPartZone, 1, -1
    tipo = Partz(iz)%tipo
    if (tipo /= "sour") then
      mate = Partz(iz)%Medium
      if (mate > 0) then
        partizone = iz
        exit
      end if
    else
      sourzone = iz
    end if
  end do
  if (partizone == 0) partizone = sourzone

  return
  end subroutine SearchforParticleZone_3D
!---split

!AA504: all the subroutine is modified to treat only coordinates and to call SetParticleParameters
!cfile SetParticles.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : SetParticles
!
! Last updating : April 08, 2014
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!AA504
! 03  Amicarelli        08Apr14        All the subroutine is modified to treat only coordinates 
!                                      and then to call SetParticleParameters (note: all the comments of v5.03 have been removed)
!AA601
! 04  Amicarelli        26Jan15        DBSPH-input (AA601). DBSPH IC for particle distribution and inlet treatment. 
!
!************************************************************************************
! Module purpose : Module for creation and setting the particles uniformly into the
!                  initial perimeter
!
! Calling routine: GeneratePart
!
! Called routines: IsParticleInternal2D
!                  IsParticleInternal3D
!                  stoptime
!                  vellaw
!AA601
!                  wavy_inlet
!
!
!************************************************************************************

  subroutine SetParticles (Nt, Nz, mate, Xmin, npps, NumParticles, IsopraS)

!Creates and sets particles uniformly into the initial perimeter mib

!Using Modules
  use GLOBAL_MODULE
  use FILES_ENTITIES
  use AdM_USER_TYPE
  use ALLOC_MODULE
  use DIAGNOSTIC_MODULE

!Implicit Declarations
  implicit none

!.. Formal Arguments ..
  integer(4),      intent(IN)                      :: Nt, Nz, mate
  double precision,intent(IN), dimension(SPACEDIM) :: Xmin
  integer(4),      intent(IN), dimension(SPACEDIM) :: npps
  integer(4),      intent(INOUT)                   :: NumParticles, IsopraS

!Local variables
!AA504 sub
  integer(4)       :: i,j,k,iaux,test,Nz_aux,nag_aux
  double precision :: aux1,aux2,aux3,rnd,tstop
  logical          :: particellainterna
  character(len=lencard)  :: nomsub = "SetParticles"
  double precision, dimension(SPACEDIM) :: PX

!External routines
  logical, external    :: IsParticleInternal3D
  logical, external    :: IsParticleInternal2D

!Executable Statements

  if ( nagpg > 0 ) then
!calcolo time stop per particelle tipo 'law'
    call stoptime ( partz(Nz), tstop )
!calcolo velocita' per particelle tipo 'law'
    call vellaw   ( partz(Nz)%vlaw,Partz(Nz)%vel,Partz(Nz)%npointv)
  end if

!AA601 sub start  
  if (Domain%tipo=="bsph") then
     if (ncord==3) then
        aux1 = + 0.25d0*Domain%dd
        else
           aux1 = - 0.25d0*Domain%dd
     endif
     aux2 = - 0.25d0*Domain%dd 
     aux3 = - 0.25d0*Domain%dd 
     iaux = 0
     else
        iaux=0
        aux1 = - Domain%dd * half
        aux2 = - Domain%dd * half
        aux3 = - Domain%dd * half
endif
!AA601 sub end
     
  PX(1) = Xmin(1) +aux1

!In case the zone is declared but is not used
  if (npps(1) < 0) return

!Loops on the X direction

 do i = 1, (npps(1)-iaux)
    PX(1) = PX(1) + Domain%dd
    PX(2) = Xmin(2) + aux2

!Loops on the Y direction
    if (ncord==2) iaux=0
    do j = 1, (npps(2)-iaux)   
       PX(2) = PX(2) + Domain%dd
       PX(3) = Xmin(3) + aux3

!Loops on the Z direction
       do k = 1, (npps(3)-iaux)
          PX(3) = PX(3) + Domain%dd

!Checks if the particle falls inside the zone
          if (ncord == 2) then
             particellainterna = IsParticleInternal2D (Nt, PX)
             else 
                particellainterna = IsParticleInternal3D (Nt, PX, IsopraS)
          end if

!In case the particle is inside the domain
          if ( particellainterna ) then

!the zone counter is increased
             NumParticles = NumParticles + 1

!the total particle number is increased
             if ( nagpg == 0 ) cycle
             
!AA504 sub start 
             test = 0
             do Nz_aux=1,NPartZone
                if (Partz(Nz_aux)%IC_source_type==2) test = 1
             end do 
             if (test==0) then
                nag = nag + 1 
!Check the storage for the reached number of particles
                if (nag > PARTICLEBUFFER) call diagnostic (arg1=10,arg2=4,arg3=nomsub)
                nag_aux = nag 
                else
                   nag_aux = NumParticles 
             endif    
!Modify the coordinates, if random
             if (Domain%RandomPos == 'r') then
               call random_number(rnd)
               pg(nag_aux)%coord(1) = PX(1) + (two * rnd - one) * 0.1d0 * Domain%dd
               call random_number(rnd)
               pg(nag_aux)%coord(2) = PX(2) + (two * rnd - one) * 0.1d0 * Domain%dd
               call random_number(rnd)
               pg(nag_aux)%coord(3) = PX(3) + (two * rnd - one) * 0.1d0 * Domain%dd
               else
                 pg(nag_aux)%coord = PX
             end if
             pg(nag_aux)%CoordOld = pg(nag_aux)%coord
!AA504 sub end             

!Setting Particle Parameters
!AA504 sub start 
             if (test==0) then
                call SetParticleParameters(nag,Nz,mate)  
                else
                   call SetParticleParameters(NumParticles,Nz,mate)  
             endif    
!AA504 sub end             
!AA601 sub start
! Nz is the zone ID, PartZ the vector of zones
             if ((Domain%tipo=="bsph").and.(Partz(pg(nag_aux)%izona)%tipo=="sour")) call wavy_inlet(Partz(pg(nag_aux)%izona))
!AA601 sub end 
          end if !particellainterna
          
        end do
     end do
  end do

  return
  end subroutine SetParticles
!---split

!AA504 all the subroutine is adapted from SetParticles Spherav503 (note: all the v503 comments are removed)
!cfile SetParticleParameters.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : SetParticleParameters
!
! Versions:
! 01   Amicarelli 08Apr14   (v5.04) New subroutine adapted from "SetParticles" (without setting coordinates)
!
!************************************************************************************
! Module purpose : Setting initial particle parameters
!
! Calling routine: SetParticles
!
! Called routines: ParticleCellNumber,stoptime,vellaw,defcolpartzero
!
!************************************************************************************

  subroutine SetParticleParameters (npi,Nz,Mate)

!Using modules
  use GLOBAL_MODULE
  use FILES_ENTITIES
  use AdM_USER_TYPE
  use ALLOC_MODULE
  use DIAGNOSTIC_MODULE

!Implicit Declaration 
  implicit none

!Formal Arguments 
  integer(4),intent(IN) :: npi,Nz,Mate

!Local variables
  double precision :: tstop
  
!External routines
  integer(4), external :: ParticleCellNumber

!Executable Statements

!AA504rm zeroing pg is already performed when pg is allocated

  if (Domain%RKscheme>1) ts0_pg(npi) = ts_pgZero

!Set the zone identifier
  pg(npi)%izona = Nz      

!Particle mass.. 
  if (ncord == 2) then
     pg(npi)%mass = Domain%PVolume * Med(Mate)%den0 
     pg(npi)%coord(2) = zero              
     pg(npi)%CoordOld(2) = zero              
     else 
        pg(npi)%mass = Domain%PVolume * Med(Mate)%den0
  end if

!Current velocity and initial velocity
  pg(npi)%vel    = partz(Nz)%vel
  pg(npi)%vstart = partz(Nz)%vel

!Calcolo time stop per particelle tipo 'law'
   call stoptime (partz(Nz),tstop)
!Calcolo velocita' per particelle tipo 'law'
   call vellaw   (partz(Nz)%vlaw,partz(Nz)%vel,partz(Nz)%npointv)
!Stopping time for blocks in movement
   pg(npi)%tstop = tstop                    

!Material ID
  pg(npi)%imed = Mate          

!Viscosities
  pg(npi)%visc = Med(Mate)%visc
  pg(npi)%mu   = Med(Mate)%visc * Med(Mate)%den0

!DB-SPH parameters
  if (Domain%tipo == "bsph") then
     pg(npi)%Gamma = one 
     pg(npi)%rhoSPH_new = zero
     pg(npi)%uni = zero
     pg(npi)%sigma = zero
     pg(npi)%dShep = zero 
     pg(npi)%FS = 0 
     pg(npi)%Gamma_last_active = zero
!AA601 start       
     pg(npi)%DBSPH_inlet_ID = 0
     pg(npi)%DBSPH_outlet_ID = 0
!AA601 end     
  endif

!AA504 start  
 ! Mixture density for granular SPH particles (bed load transport) 
  if (Granular_flows_options%ID_erosion_criterion==1) then
     if (Med(pg(npi)%imed)%tipo=="granular") then
         pg(npi)%dens = Med(Granular_flows_options%ID_granular)%den0_s 
! IC viscosity from pure fluid, as I cannot calculate the right one at this stage         
         pg(npi)%mu = Med(Granular_flows_options%ID_main_fluid)%visc*Med(Granular_flows_options%ID_main_fluid)%den0
         pg(npi)%visc = Med(Granular_flows_options%ID_main_fluid)%visc
     endif
     call initialization_fixed_granular_particle(npi)
     pg(npi)%sigma_prime = 0.0d0
  endif
!AA504 end  
  
!Particle status, depending on the velocity components (fluid or solid)
  if (index(Med(Mate)%tipo,"liquid") > 0 .or. index(Med(Mate)%tipo,"smagorin") > 0) then
     pg(npi)%state = "flu" 
     else if (index(Med(Mate)%tipo,"granular") > 0 .or. index(Med(Mate)%tipo,"general") > 0) then
        pg(npi)%state = "sol"
        else if (index(Med(Mate)%tipo,"gas") > 0) then
           pg(npi)%state = "flu" 
  end if    
  if (index(Med(Mate)%tipo,"granular") > 0 .or. index(Med(Mate)%tipo,"general") > 0) then
     if (pg(npi)%vel(1)/=zero .or. pg(npi)%vel(2)/=zero .or. pg(npi)%vel(3)/=zero) pg(npi)%state = "flu"     
  end if

!Motion index
  pg(npi)%vel_type = partz(Nz)%move        
  if ( partz(Nz)%move /= "std" ) pg(npi)%visc = zero

!Boundary slip condition
  pg(npi)%slip = partz(Nz)%slip  

!Grid cell 
  pg(npi)%cella = ParticleCellNumber(pg(npi)%coord)

!Particle color definition as from input file
  call defcolpartzero (Nz,partz,pg(npi))

!Modulo diffusione
  if (diffusione) then
     if (pg(npi)%imed == 1) then
         pg(npi)%VolFra    = VFmn
     end if
     if (pg(npi)%imed == 2) then          
         pg(npi)%VolFra    = VFmx
     end if
     else
        pg(npi)%VolFra    = one
  end if

!Modulo esplosione(gas)
  if (esplosione) then
     pg(npi)%IntEn  = Med(pg(npi)%imed)%InitialIntEn
  end if
  
  return
  end subroutine SetParticleParameters
!---split

!cfile stoptime.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : stoptime
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
! Module purpose : Module for definition of the stop time
!
! Calling routine: SetParticles
!
! Called routines: diagnostic
!
!************************************************************************************
!
  subroutine stoptime ( partzlocal, tstop )
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
  type(TyZone),    intent(INOUT)  :: partzlocal
  double precision,intent(INOUT)  :: tstop
!
!.. Local Scalars ..
  integer(4)       :: k,n, icord
  double precision :: tstopc, acc, deltat, spo, dspo, rad
  logical     :: out
  character(len=lencard)  :: nomsub = "stoptime"
!
!.. local arrays
  double precision, dimension(3)   :: dxyz
  double precision, dimension(3,2) :: vlimits, tlimits
!
!.. Executable Statements ..
!
  tstop = max_positive_number
!
!.. checks if there are fixed particles
!
  if ( partzlocal%move == "fix" ) then
!
    do n = 1, ncord
!
       icord = icoordp(n,ncord-1)
       if ( partzlocal%vel(icord) > zero ) then
          tstopc = (Domain%coord(icord,2)-partzlocal%coordMM(icord,2)-one*Domain%dd)/partzlocal%vel(icord)
       else if ( partzlocal%vel(icord) < zero ) then
          tstopc = (Domain%coord(icord,1)-partzlocal%coordMM(icord,1)+one*Domain%dd)/partzlocal%vel(icord)
       else
          tstopc = Domain%tmax
       end if
       tstop = min(tstop,tstopc)
!
     end do
!
!.. checks if there are assigned motion laws for particles
!
   else if ( partzlocal%move == "law" ) then
!
!.. evaluates the paths
!
     vlimits = zero
     do n = 1, ncord
       icord = icoordp(n,ncord-1)
!!       vlimits(icord,1) = Domain%coord(icord,1) + Domain%dd
!!       vlimits(icord,2) = Domain%coord(icord,2) - Domain%dd
       vlimits(icord,1) = Domain%coord(icord,1)
       vlimits(icord,2) = Domain%coord(icord,2)
     end do
!
     tlimits = zero
     dxyz    = zero
     out     = .FALSE.
!
     LAW_ZONE_LOOP: do k = 2, partzlocal%npointv
!
       COORDS_LOOP: do n = 1, ncord
!
         icord = icoordp(n,ncord-1)
!
!.. evaluates the acceleration
!
         deltat = partzlocal%vlaw(0,k) - partzlocal%vlaw(0,k-1)
         acc    = ( partzlocal%vlaw(icord,k) - partzlocal%vlaw(icord,k-1) ) / deltat
!
!.. upgrade the path
!
         dspo = partzlocal%vlaw(icord,k-1) * deltat + acc * deltat * deltat * half
         spo  = dxyz(icord) + dspo
!
!.. checks if the minimum limit has been overridden and how much time has been required
!
         if ( (partzlocal%coordMM(icord,1)+spo) < vlimits(icord,1) ) then
!
           out  = .TRUE.
           dspo = vlimits(icord,1) - (partzlocal%coordMM(icord,1)+dxyz(icord))
           if ( acc == zero ) then
             if (partzlocal%vlaw(icord,k-1) == zero) then
               deltat = max_positive_number
             else
               deltat = dspo / partzlocal%vlaw(icord,k-1)
             end if
           else
             rad  = partzlocal%vlaw(icord,k-1)*partzlocal%vlaw(icord,k-1) - 4.0*0.5*acc*dspo
             if ( rad >= zero ) then
               rad = Dsqrt(rad)  
               deltat = ( partzlocal%vlaw(icord,k-1) + rad ) / ( 2*0.5*acc )
             else
               call diagnostic (arg1=10,arg2=88,arg3=nomsub)       
             end if
           end if
         end if
!
!.. add the interval time to the total time
!
         tlimits(icord,1) = tlimits(icord,1) + deltat  
!
!.. check if the maximum limit has been overriden and how much time has been required 
!
         if ( (partzlocal%coordMM(icord,2)+spo) > vlimits(icord,2) ) then 
           out  = .TRUE.
           dspo = vlimits(icord,2) - (partzlocal%coordMM(icord,2)+dxyz(icord))
           if ( acc == zero ) then
             if (partzlocal%vlaw(icord,k-1) == zero) then
               deltat = max_positive_number
             else
               deltat = dspo / partzlocal%vlaw(icord,k-1)
             end if
           else
             rad  = partzlocal%vlaw(icord,k-1)*partzlocal%vlaw(icord,k-1) - 4.0*0.5*acc*dspo
             if ( rad >= zero ) then
               rad = Dsqrt(rad)  
               deltat = ( -partzlocal%vlaw(icord,k-1) + rad ) / ( 2*0.5*acc )
             else
               call diagnostic (arg1=10,arg2=88,arg3=nomsub)       
             end if
           end if
         end if
!
!.. add the interval time to the total time
!
         tlimits(icord,2) = tlimits(icord,2) + deltat  ! somma il tempo dell'intervallo
!
!.. save the evaluated displacement along the path
!
         dxyz(icord) = dxyz(icord) + dspo
!
       end do COORDS_LOOP
!
!.. ends if it is gone out of the domain
!
       if ( out ) exit LAW_ZONE_LOOP

     end do LAW_ZONE_LOOP
!
!.. evaluates the minimum time
!
     do n = 1, ncord
       icord = icoordp(n,ncord-1)
       tstop = min(tstop,tlimits(icord,1),tlimits(icord,2))
     end do
     partzlocal%move = "fix"
!
   else
!
     tstop = Domain%tmax
!
   end if
!
  return
  end subroutine stoptime
!---split

!cfile subCalcPreIdro.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : subCalcPreIdro
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
! Module purpose : Module to calculate hydrostatic pressure (submerged)
!
! Calling routine: Gest_Input
!
! Called routines: 
!
!************************************************************************************
!
  subroutine SubCalcPreIdro
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
  use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Local Scalars ..
  integer(4)       :: idpel, Nz, pl_id
  integer(4)       :: npi,iappo,i,j,k,numcell,ncelcorr,m,nnlocal,igridi,jgridi,kgridi,nnsave
  double precision :: ZQuotaMediumCorr,ZQuotaColonna,ZQuotaSecondMedium,affond1,affond2
  double precision :: pl_quote,gravmod,coshor,senhor
  logical          :: foundcell
  character(len=lencard)  :: nomsub = "SubCalcPreIdro"
!
! External functions and subroutines
  integer(4), external :: CellNumber,ParticleCellNumber,CellIndices
!
!.. Executable Statements ..
!
!.. evaluates the gravity module and the versors of the Z reference direction with respect the physical vertical direction
!
  gravmod = Dsqrt(Domain%grav(1)*Domain%grav(1) + Domain%grav(2)*Domain%grav(2) + Domain%grav(3)*Domain%grav(3) )
  if ( gravmod > zero) then
    coshor= -Domain%grav(3)/gravmod
    senhor=  Domain%grav(1)/gravmod
  else
    coshor= zero
    senhor= zero
  end if
!
   pl_quote = max_negative_number
   pl_id = 0
!
!.. searches for free level conditions, if present ("pl" type)
!
!
!!! !_!$omp parallel do default(none) private(npi,nz) shared(nag,Pg,Partz,pl_id,pl_quote)
!
   do npi = 1,nag
!
     Nz = pg(npi)%izona
     if (partz(Nz)%pressure /= "pl" ) cycle 
     if (pl_id == 0) then
       pl_id = int(partz(Nz)%valp)
     else if (pl_id /= int(partz(Nz)%valp)) then
       call diagnostic (arg1=10,arg2=8,arg3=nomsub)
     end if
     pl_quote = max(pl_quote,pg(npi)%coord(3))
!
   end do 
!
!!! !_!$omp end parallel do
!
!.. loops on all the particles
!
!
!!!!   !!_!$omp parallel do default(private) shared(nag,Pg,Med,Domain,Grid,Partz,Icont,npartord,pl_id,pl_quote,diffusione)
!
  particle_loop:  do npi = 1,nag
!
!.. the pressure for the current particle is assigned ("pa" type), so the input value is forced
!
    Nz = pg(npi)%izona
    if (partz(Nz)%pressure == "pa" ) then  
!     
      pg(npi)%pres = partz(Nz)%valp
      pg(npi)%dens = med(pg(npi)%imed)%den0
      cycle particle_loop
!
    end if   
!
!.. detect the cell number for the current particle
!
    ncelcorr = ParticleCellNumber(pg(npi)%coord)
    iappo = CellIndices(ncelcorr,igridi,jgridi,kgridi)
    i = igridi
    j = jgridi
!
!.. check if there is a cell including a medium interface in the column of cells above it, 
!.. and if this is it evaluates the medium interface level
!
    idpel = 0
    foundcell = .FALSE.
    ZQuotaMediumCorr   = pg(npi)%coord(3)
    ZQuotaColonna      = pg(npi)%coord(3)
    ZQuotaSecondMedium = pg(npi)%coord(3)
!
!.. loops on the above cell column
!
    do k = kgridi,Grid%ncd(3)
!
      numcell = CellNumber(i,j,k)
      if (numcell == 0) cycle
!
!.. loops on the number of particles inside the cell considered
!
      do m = Icont(numcell),Icont(numcell+1)-1
!
        nnlocal = npartord(m)
!
!.. a particle of different medium is found, so the interface cell identifier is set
!
!!        if (med(pg(npi)%imed)%index /= med(pg(nnlocal)%imed)%index ) then
!        if (med(pg(npi)%imed)%index /= med(pg(nnlocal)%imed)%index .and. &
!            pg(nnlocal)%IntEn < 10000.0) then  ! modificare per riconoscere il mezzo gas esplosione
        if (med(pg(npi)%imed)%index /= med(pg(nnlocal)%imed)%index .and. &
            index(Med(pg(nnlocal)%imed)%tipo,"gas") == 0) then
!
!.. the interface cell is set
!
          if (.not. foundcell) then
            idpel = numcell
            foundcell = .TRUE.
            nnsave = nnlocal
          end if
!
!..the minimum level in the interface cell for the medium different from the current one is set
!
          if (numcell == idpel) then
            ZQuotaSecondMedium = min(ZQuotaSecondMedium,pg(nnlocal)%coord(3))
          end if
!
!.. increase in any case the reference level of the column 
!
          ZQuotaColonna = max(ZQuotaColonna,pg(nnlocal)%coord(3))
      
        else
!
!.. evaluates the reference levels: ZQuotaMediumCorr is the maximum level of the same medium of the current particle,
!.. while ZQuotaColonna is the maximum level of the other medium located above (maximum of the column of cells)
!
          ZQuotaMediumCorr = max(ZQuotaMediumCorr,pg(nnlocal)%coord(3))
          ZQuotaColonna    = max(ZQuotaColonna,pg(nnlocal)%coord(3))
!
        end if
      end do
    end do
!
!.. checks if the current particles is inside the intermediate cell but it is of the upper medium type
!
    if (abs(ZQuotaMediumCorr - ZQuotaColonna) < xyz_tolerance) foundcell = .false.
!
!.. set the reference pressure quote depending on the condition type
!
    if (partz(Nz)%pressure == "qp" ) ZQuotaColonna = partz(Nz)%valp
    if (partz(Nz)%pressure == "pl" ) ZQuotaColonna = pl_quote
!
!.. an upper medium interface cell has been found 
!
    if (foundcell) then
!
!.. the average quote of the medium interface in the cell is calculated
!
      ZQuotaMediumCorr  = (ZquotaMediumCorr + ZQuotaSecondMedium) * half
!
!.. the pressure and density are evaluated for the current particle accounting for the medium interface
!
      affond1 = (ZQuotaColonna - ZQuotaMediumCorr) * coshor + pg(nnsave)%coord(1) * senhor
      affond2 = (ZQuotaMediumCorr - pg(npi)%coord(3)) * coshor + pg(npi)%coord(1) * senhor 
!AA601 start
      if (Domain%tipo == "bsph") then   
         pg(npi)%pres = 0.d0 * (affond1 * med(pg(npi)%imed)%den0 * gravmod) * (1.d0-pg(npi)%coord(1)/0.5925)
         else
!AA601 end          
            pg(npi)%pres = (affond1 * Med(pg(nnsave)%imed)%den0 + affond2 * med(pg(npi)%imed)%den0) * gravmod
!AA401 dam break ad-hoc correction test
!       pg(npi)%pres = 0.1 * ((affond1 * Med(pg(nnsave)%imed)%den0 + affond2 * med(pg(npi)%imed)%den0) * gravmod)
!AA601
      endif
!.. the column includes only one of the media
!
    else
!
      affond1 = (ZQuotaColonna-pg(npi)%coord(3)) * coshor + pg(npi)%coord(1) * senhor 
      pg(npi)%pres = affond1 * med(pg(npi)%imed)%den0 * gravmod
!
!AA406sub start
      if (Domain%tipo == "bsph") then
!AA601 sub          
         pg(npi)%pres = 0.d0 * (affond1 * med(pg(npi)%imed)%den0 * gravmod) * (1.d0-pg(npi)%coord(1)/0.5925)
      else
!AA401 dam break ad-hoc correction test
!AA501 sub
             pg(npi)%pres = (affond1 * med(pg(npi)%imed)%den0 * gravmod)
!            pg(npi)%pres = 0.1 * (affond1 * med(pg(npi)%imed)%den0 * gravmod)
!
      end if
!AA406sub end
!
  endif
!
!.. evaluates the density value
!
      pg(npi)%dens = (one + pg(npi)%pres/med(pg(npi)%imed)%eps) * med(pg(npi)%imed)%den0
!
!AA406test
!      pg(npi)%mass = pg(npi)%mass*pg(npi)%dens/med(pg(npi)%imed)%den0
!
      pg(npi)%dden = zero
!
!.. Diffusion model
!
    if (diffusione) then
      if ( pg(npi)%VolFra == VFmx) then
!§
        pg(npi)%pres = ((ZQuotaColonna - ZQuotaMediumCorr) * (med(2)%den0 * VFmn + med(1)%den0 * (one - VFmn)) &
                         + (ZQuotaMediumCorr - pg(npi)%coord(3)) * (med(2)%den0 * VFmx + med(1)%den0 * (one - VFmx))) * gravmod 
!        pg(npi)%pres = ((ZQuotaColonna-0.1d0) * (med(1)%den0 * (1- VFmn) + med(2)%den0 * VFmn) &
!                          + (0.1d0-pg(npi)%coord(3)) * (med(1)%den0 * (1- VFmx) + med(2)%den0 * VFmx)) * gravmod 
        pg(npi)%rhoc = (one + pg(npi)%pres/med(pg(npi)%imed)%eps) * med(pg(npi)%imed)%den0 
        pg(npi)%rhow = (one + pg(npi)%pres/med(1)%eps) * med(1)%den0
        pg(npi)%dens = VFmx * pg(npi)%rhoc + (one - VFmx) * pg(npi)%rhow
!        pg(npi)%rhoc = pg(npi)%dens
!        pg(npi)%rhow = med(1)%den0
      else if (pg(npi)%VolFra == VFmn) then
        pg(npi)%pres = ((ZQuotaColonna-pg(npi)%coord(3)) * (med(2)%den0 * VFmn + med(1)%den0 * (one - VFmn))) * gravmod 
!        pg(npi)%pres = ((ZQuotaColonna-pg(npi)%coord(3)) * (med(1)%den0 * (one - VFmn) + med(2)%den0 * VFmn)) * gravmod !&
        pg(npi)%rhoc = (one + pg(npi)%pres/med(2)%eps) * med(2)%den0 
        pg(npi)%rhow = (one + pg(npi)%pres/med(pg(npi)%imed)%eps) * med(pg(npi)%imed)%den0
        pg(npi)%dens = VFmn * pg(npi)%rhoc + (one - VFmn) * pg(npi)%rhow
!        pg(npi)%rhoc = med(2)%den0
!        pg(npi)%rhow = pg(npi)%dens
!§
      end if
    end if
!
  end do particle_loop
!
!!!!   !!_!$omp end parallel do
!
  return
  end subroutine subCalcPreIdro
!---split

!cfile s_ctime.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine s_ctime ( nout )

implicit none

!include "use_dflib.f90"
!use DFPORT
!PANE_WINNT      integer time,time_

integer(4)    :: nout
character(24) :: str
integer(4),    external :: time
character(24), external :: ctime

!PANE_WINNT      time = time_()
!PANE_WINNT      call ctime_(str,time)
 str = ctime(time())
 write (nout,'(a,a)') 's_ctime routine --> ',str

return
end subroutine s_ctime
!---split

!cfile s_secon2.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine s_secon2(tempo) 
!----------------------------------------------------------------------
!
!     subroutine: S_SECON2 
!
!     DESCRIZIONE:
!       Fornisce il tempo di CPU user+system ed il tempo di ELAPSED.
!
!     AUTORE: F. Riccio, E. Bon
!
!     NOTE:   Si utilizza la funzione di sistema  ETIME    (Alliant)
!             Si utilizza la funzione di sistema  ETIME_   (RISC6000)
!             Si utilizza la funzione di sistema  TIMEF    (RISC6000)
!
!     PARAMETRI DI INPUT:
!
!     PARAMETRI DI OUTPUT:
!       TEMPO(1) : tempo di CPU  all'atto della chiamata (user + system)
!       TEMPO(2) : tempo elapsed all'atto della chiamata.
!
!     VARIABILI UTILIZZATE:
!
!     ROUTINES CHIAMATE:
!
!
!     DATA: 17-11-93
!
!     SUCCESSIVE MODIFICHE:
!
!     DATA    AUTORE      OGGETTO
!     20-12-94 E.Bon      Adattamento a RISC 6000
!
!
!----------------------------------------------------------------------
implicit none

double precision,dimension(2) :: tempo
double precision,dimension(2) :: t2
integer(4)                    :: time
double precision              :: etime

 tempo(1) = etime(t2) 
 tempo(2) = dfloat(time())

return 
end subroutine s_secon2
!---split

!cfile start_and_stop.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : start_and_stop
!
! Last updating : May 26, 2014
!
! Improvement traceback:
!
! 00  Agate/Guandalini  18/12/07       INITIAL
! 01  Agate             07/05/12       Limited and licensed version
! 02  Agate/Guandalini  27/03/13       modify the key system
!
!************************************************************************************
! Module purpose : management of the time statistics and machine checks
!
! Calling modules: sphera, gest_trans
!
! Called modules : getarg       :acquire the command line parameters
!                  system_clock :check the time independently on the machine system
!                  iargc        :detect the number of command line parameter
!                  hostname     :check the operating machine
!                  getcwd       :check the current directory
!                  getenv       : check the contents of the SWEET_DIR env variable
!      
!************************************************************************************
!
  subroutine start_and_stop (iset,itot)
!
!.. directive to compiler on Windows in order to get the machine identifier
!.. it must be removed for Linux or Unix compilers
!
!....... MS$ ATTRIBUTES ALIAS:'_HOSTNAM@8' :: hostnm  ! compilazione windows fortran 6.x
!
use time_usertype
use FILES_ENTITIES
!
!.. Implicits ..
  implicit none
!
!.. Formal Arguments ..
  integer(4), intent(in) :: iset,itot
!
!.. Local variables
  integer(4)          :: is_cwd,i,hours,days,minutes,seconds,ishost
  character(LEN=1024) :: pwd_name,sphera_dir
  character(LEN=256)  :: exe_name  !text, 
  character(LEN=8)    :: dat
  character(LEN=5)    :: zone
  character(LEN=10)   :: ct
  double precision    :: appo
!
!.. Local Arrays ..
  integer(4),      dimension(8)  :: dat_array
  character(LEN=3),dimension(12) :: mesi
!
!.. External Functions ..
  integer :: HOSTNM
  integer :: GETCWD
  double precision,external :: omp_get_wtime
!
!.. Data Assignments ..
  data mesi/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/
!
!.. Executable Statements ...
!
  routine='start_and_stop      '
!
!AA406
 ishost = 0 
! 
!.. start detecting the executable module and the project prefix. if not
!.. assigned, the excution aborts ..
!
  select case (iset)
    case (0)
!
!..detect the data and time when execution starts ..
!
!
!.. initialize the CPU time statistics ..
!
!.. tot_times(1) = total initialization time ???
!.. tot_times(2) = total Gest_Input time
!.. tot_times(3) = total Gest_Trans time
!.. tot_times(4) = total Loop_ghost time
!.. tot_times(5) = total Loop_Irre_2D/3D time
!.. tot_times(6) = total Storage array & Boundary integrals time
!.. tot_times(7) = total Motion Equation time
!.. tot_times(8) = total Velocity smoothing time
!.. tot_times(9) = total Outgone/Source/Ordgrid
!.. tot_times(10) = total Search for neighbourhood particles
!.. tot_times(11) = total Boundary integrals
!.. tot_times(12) = total Continuity equation time
!.. tot_times(13) = total State equation time
!.. tot_times(14) = total Pressure smoothing time
!.. tot_times(15) = total Apparent viscosity time
!.. tot_times(16) = total Diffusion model time
!.. tot_times(17) = RK time integration
!AA406
!.. tot_times(18) = wall parameter updates
!AA501b comment
!.. tot_times(19) = rigid body transport 
!..
!.. tot_times(.) = total ... time
!..
!.. tot_times(numb_subr) = total time
!
!
      tot_times = zero
      tot_call = 0
      call cpu_time(tot_times(numb_subr,1))
!$ tot_times(numb_subr,1) = omp_get_wtime()

      call DATE_AND_TIME(dat,ct,zone,dat_array)
      date_exec = mesi(dat_array(2))//" "//dat(7:8)//", "//dat(1:4)// &
                " at "//ct(1:2)//":"//ct(3:4)//":"//ct(5:10)//" "//zone//" GMT"
      case_data(1:12) = date_exec(1:12)
      case_hour(1:8) = date_exec(17:24)
      read (ct(1:2),'(i2)') itime_struct%ihr
      read (ct(3:4),'(i2)') itime_struct%imin
      read (ct(5:6),'(i2)') itime_struct%isec
!
!.. read the name of the machine from environment ..
!
      ishost = hostnm(host_name)
!
!.. read the pathname of the working directory ..
!
      is_cwd = getcwd(pwd_name)
!
!.. read the SPHERA_DIR environmental variable 
!
!      is_env = getenv ("SPHERA_DIR",sphera_dir)
      call getenv ("SPHERA_DIR",sphera_dir)
!      sphera_dir = "./"
      is_cwd = len_trim(sphera_dir)
!
!.. read the name of executable file
      call getarg (0,exe_name)
      exe_name = adjustl(exe_name)
!      exe_name = text(is_cwd+1:len_trim(text))
!
! .. Open the output listing file (file ASCII) and assign the other file names
!
!      nomefile_tempi =  TRIM(nomecaso)//"_tempi.sta"
!
!      open (unit=nout,file=nomefile_tempi,form='formatted',recl=124,status='REPLACE',err=997,iostat=ierr)
!
! 997 stop 'errore open file: stoppato programma'
!
! ..Write the calculation heading on the output listing ..
!
      write (nout,*)
      write (nout,'(a,a)') " >>> SPHERA execution started the ",TRIM(date_exec)
      write (nout,*)
      write (nout,'(a,a)') " >>> SPHERA execution module :",TRIM(exe_name)
      write (nout,*)
      write (nout,'(a,a)') " >>> SPHERA execution machine :",TRIM(host_name)
!
      write (nout,*)
      write (nout,'(a,a)') " >>> Project working directory :",TRIM(pwd_name)
      write (nout,*)
!      write (nout,'(a,a)') " >>> Project identifier :",TRIM(prefix)
!      write (nout,*)
      write (nout,'(a,a)') " >>> Case identifier :",TRIM(nomecaso)
!
!.. check for the license correctness
!

!AA601 test sub start (with no license)
!call KeyDecoderCheck (sphera_dir,exe_name,dat)
Erosion_Module_Shields_Seminara = .true.
Erosion_Module_Shields_Mohr = .true.
Diffusion_Module = .true.
Explosion_Module = .true.
TemporalScheme_Module = .true.
BodyDynamics_Module = .true.
DBSPH_Module  = .true.
MultiFluid_Module = .true.
MoreFluids_Module = .true.
Granular_flux = .true.
!AA601 test sub end (with no license)

!
!.. detect the current time and check the progressive time ..
!
    case (1)
!
!.. end the execution evaluating the CPU and elapsed time statistics
!
      call cpu_time(appo)
!$ appo = omp_get_wtime()
      tot_times(numb_subr,2) = appo - tot_times(numb_subr,1)


      call DATE_AND_TIME(dat,ct,zone,dat_array)
      date_exec = mesi(dat_array(2))//" "//dat(7:8)//", "//dat(1:4)// &
                " at "//ct(1:2)//":"//ct(3:4)//":"//ct(5:10)//" "//zone//" GMT"
!
      write (nout,'(a)')   "----------------------------------------------------------------------------------------"
      write (nout,'(a,a)') " SPHERA execution ended the ",TRIM(date_exec)
      write (nout,'(a)')   "----------------------------------------------------------------------------------------"
      write (nout,'(a)')   " "
      write (nout,'(a)')   " Summary of the execution times (s) and internal iteration statistics:"
      write (nout,'(a)')   "----------------------------------------------------------------------------------------"
      write (nout,'(a)')   "   Task                                                  Total          %      number of"
      write (nout,'(a)')   "                                                         Elapsed                 calls" 
      write (nout,'(a)')   "----------------------------------------------------------------------------------------"
      do i = 1,numb_subr
        if (tot_call(i) == 0) cycle
        write (nout,"(4x,i3,2x,a,(f15.5,3x,f6.2,1x),i14)") &
        i,tot_routines(i),tot_times(i,2),tot_times(i,2)*100./tot_times(numb_subr,2),tot_call(i)
      end do
      write (nout,'(a)')   "----------------------------------------------------------------------------------------"
      days = int(tot_times(numb_subr,2)/86400)
      hours = int((tot_times(numb_subr,2)-days*86400)/3600)
      minutes = int((tot_times(numb_subr,2)-days*86400-hours*3600)/60)
      seconds = ceiling(tot_times(numb_subr,2)-days*86400-hours*3600-minutes*60)
      write (nout,"(4x,a,f10.2,a,4(i2,a3))")  &
        "Total elapsed time    : ",tot_times(numb_subr,2)," s equal to ",days," d ",hours," h ",  &
        minutes," m ",seconds," s " 
      write (nout,'(a)')   "----------------------------------------------------------------------------------------"
!
! .. close the files
!
!      close (unit=nout)
!
! ..Detect the time for the initial call to subroutines ..
!
    case (2)
      tot_call(itot) = tot_call(itot) + 1
      tot_call(numb_subr) = tot_call(numb_subr) + 1

      call cpu_time(appo)
!$ appo = omp_get_wtime()
      tot_times(itot,1) = appo
!
!.. evaluates the incremental time for the different subroutines
!
    case (3)
      call cpu_time(appo)
!$ appo = omp_get_wtime()
      tot_times(itot,2) = tot_times(itot,2) + (appo - tot_times(itot,1))
!
!.. ends the execution of the binary restart file converter
!
  end select
!
!.. I/O Formats ..
!
  return
end subroutine start_and_stop
!---split

!cfile Vector_Product.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : Vector_Product
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
! Module purpose : Module to Compute and return in ww(1 to SPACEDIM) the components
!                  of the vector product of vectors uu(1 to SPACEDIM) and
!                  vv(1 to SPACEDIM)
!
!AA501b modified
!AA601 sub
! Calling routine: DefineLocalSystemVersors,RHS_body_dynamics,area_triangle
!
! Called routines: 
!
!************************************************************************************
!
subroutine Vector_Product ( uu, VV, ww, SPACEDIM )
!Computes and returns in ww(1 to SPACEDIM) the components of the vector product
!of vectors uu(1 to SPACEDIM) and vv(1 to SPACEDIM)
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
integer(4),      intent(IN)                        :: SPACEDIM
double precision,intent(IN),   dimension(SPACEDIM) :: uu
double precision,intent(IN),   dimension(SPACEDIM) :: VV
double precision,intent(INOUT),dimension(SPACEDIM) :: ww
!
!.. Local Scalars ..
integer(4) :: i, j, k
!
!.. Local Arrays ..
integer(4), dimension(3) :: iseg = (/ 2,3,1 /)
!
!.. Executable Statements ..
!
 do i = 1, SPACEDIM
   j = iseg(i)
   k = iseg(j)
   ww(i) = uu(j) * VV(k) - uu(k) * VV(j)
 end do
!
return
end subroutine Vector_Product
!---split

!cfile VelLaw.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : VelLaw
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
! Module purpose : Module for velocity law calculation
!
! Calling routine: Loop_Irre_2D, Loop_Irre3D, SetParticles
!
! Called routines: 
!
!************************************************************************************
!
subroutine VelLaw (vlaw,vel,np)
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
double precision,intent(IN), dimension(0:3,MAXPOINTSVLAW) :: vlaw
double precision,intent(OUT),dimension(3)                 :: vel
integer(4),      intent(IN)                               :: np
!
!.. Local Scalars ..
double precision :: fra
integer(4)       :: n
!
!.. Executable Statements ..
!
 if ( np <= 1 ) return
 do n = 2,np
   if ( tempo > vlaw(0,n) ) cycle
   fra = ( tempo - vlaw(0,n-1) ) / ( vlaw(0,n) - vlaw(0,n-1) )   
   vel(1:3) = vlaw(1:3,n-1) +  ( vlaw(1:3,n) - vlaw(1:3,n-1) ) * fra
   return
 end do
 vel(1:3) = vlaw(1:3,np)

return
end subroutine VelLaw
!---split

!cfile viscapp.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : viscapp
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
! Module purpose : Module to define for each particle the viscosity based on 
!                  reology model
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
!
subroutine viscapp
!* definisce per ogni particella la viscosita' in base 
!* al modello reologico di appartenenza
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
double precision,parameter :: alfa = 1.0d0
!
!.. Local Scalars ..
integer(4)       :: npi
double precision :: mu, mumax, secinv,cuin, smalen, smalenq, visc1, visc2 
!
!.. External routines ..
character(80),external :: lcase
!
!.. Executable Statements ..
!
if (.not. diffusione) then
!
!$omp parallel do default(none) &
!$omp private(npi,visc1,smalen,smalenq,visc2,secinv,cuin,mu,mumax) &
!$omp shared(nag,pg,Med,Domain,it_corrente)
!
  do npi = 1,nag
    if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle
!
    select case ( Med(pg(npi)%imed)%tipo )
!
    case ( "liquid  " )
       pg(npi)%visc = Med(pg(npi)%imed)%visc
!   
    case ( "gas     " )
       pg(npi)%visc = Med(pg(npi)%imed)%visc
!   
    case ( "smagorin" )  
       visc1   = Med(pg(npi)%imed)%visc 
       smalen  = Med(pg(npi)%imed)%Cs * Domain%h               
       smalenq = smalen * smalen  
       visc2   = smalenq * two * pg(npi)%secinv
       pg(npi)%visc = visc1 + visc2 
!
    case ( "general " )
       secinv = two * pg(npi)%secinv
       cuin   = Med(pg(npi)%imed)%cuin
       mu     = Med(pg(npi)%imed)%taucri/(secinv+0.0001d0) + Med(pg(npi)%imed)%cons*((secinv+0.0001d0)**(cuin-one))
       mumax  = Med(pg(npi)%imed)%mumx
       pg(npi)%visc = min(mumax,mu) / pg(npi)%dens
!
!!!!!VERIFICARE
       if (pg(npi)%state == "sol") then
!.. controllo prime NIterSol iterazioni mezzo granulare non si muove
         if ((pg(npi)%visc /= mumax/pg(npi)%dens .and. it_corrente > Med(pg(npi)%imed)%NIterSol) .or. pg(npi)%kodvel == 2) then
           pg(npi)%state = "flu"
         end if
       end if
!!!!!VERIFICARE
!
    case ( "granular" )
!       secinv = two * pg(npi)%secinv        
!!       cuin   = Med(pg(npi)%imed)%cuin
!       pre    = (max(zero,pg(npi)%pres))
!       coes   = Med(pg(npi)%imed)%coes
!       coeff  = sin (Med(pg(npi)%imed)%phi)
!       mu     = (coes * cos(Med(pg(npi)%imed)%phi) + pre * coeff) / (secinv+0.00001d0)    !11-4 coesione
!       mumax  = Med(pg(npi)%imed)%mumx
!       pg(npi)%visc = min(mumax,mu) / pg(npi)%dens
       pg(npi)%visc = pg(npi)%mu / pg(npi)%dens
!
!.. verifies the "state" of a granular particle
!
!ç       if (pg(npi)%state == "sol") then
!.. controllo prime NIterSol iterazioni mezzo granulare non si muove
!AGGR         if (pg(npi)%visc /= mumax/pg(npi)%dens .and. it_corrente > Med(pg(npi)%imed)%NIterSol) then
!AGGR         if ((pg(npi)%visc /= mumax/pg(npi)%dens .and. it_corrente > Med(pg(npi)%imed)%NIterSol) .or. pg(npi)%kodvel == 2) then
!ç         if ((pg(npi)%visc /= mumax/pg(npi)%dens .and. it_corrente > Med(pg(npi)%imed)%NIterSol) .or. pg(npi)%kodvel == 2) then
!!!!         if (it_corrente > Med(pg(npi)%imed)%NIterSol) then
!ç           pg(npi)%state = "flu"
!ç         end if
!ç       end if
!
    case default
!
    end select
!
  end do
!
!$omp end parallel do
!
!
else if (diffusione) then
!!! modello bifluido
!
!$omp parallel do default(none) private(npi) shared(nag,pg,Med)
!
  do npi = 1,nag
    if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std") cycle

    if (index(Med(1)%tipo,"liquid") > 0 .and. index(Med(2)%tipo,"liquid") > 0) then

      pg(npi)%visc = pg(npi)%VolFra*Med(2)%visc + (one-pg(npi)%VolFra)*Med(1)%visc

    else if (index(Med(1)%tipo,"liquid") > 0 .and. index(Med(2)%tipo,"granular") > 0) then

!     if (pg(npi)%VolFra > 0.1d0) then
!§
      if (pg(npi)%VolFra >= VFmn .and. pg(npi)%VolFra <= VFmx) then
!       secinv = two * pg(npi)%secinv        
!       cuin   = Med(2)%cuin
!       pre    = (max(zero,pg(npi)%pres))
!       coes   = Med(2)%coes
!       coeff  = sin (Med(2)%phi)
!       mu     = pg(npi)%VolFra * (coes * cos(Med(2)%phi) + pre*coeff) / (secinv+0.00001d0) + &
!                (one-pg(npi)%VolFra) * Med(1)%visc
!       mumax  = Med(2)%mumx
!       pg(npi)%visc = min(mumax,mu) / pg(npi)%dens
        pg(npi)%visc = pg(npi)%VolFra*Med(2)%visc + (one-pg(npi)%VolFra)*Med(1)%visc
     
      else if (pg(npi)%VolFra < VFmn) then

        pg(npi)%visc = Med(1)%visc

      else if (pg(npi)%VolFra < VFmx) then

        pg(npi)%visc = Med(2)%visc
       
      end if
!§
!
!§
!!!!!VERIFICARE
!     if (pg(npi)%state == "sol") then
!       if (pg(npi)%visc /= (mumax/pg(npi)%dens) .or. pg(npi)%VolFra <= VFmx) then
!.. controllo prime NIterSol iterazioni mezzo granulare non si muove
!         if ( it_corrente > Med(pg(npi)%imed)%NIterSol) then
!           pg(npi)%state = "flu"
!         end if
!       end if
!     end if
!!!!!VERIFICARE
!§
    end if
!
  end do
!
!$omp end parallel do
!
end if
!
return
end subroutine viscapp
!---split

!cfile viscomon.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : viscomon
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!AA504
! 03  Amicarelli        08Apr14        (v5.04) Modifications for Monaghan's viscosisty (active also for separating elements) and 
!                                      neglection of the molecular viscosity depending on the velocity divergence (momentum equation)
!
!************************************************************************************
! Module purpose : Module to 
! 
! Calling routine: inter_EqMoto
!
! Called routines: 
!
!************************************************************************************
!
subroutine viscomon (npi,npj,npartint,dervel,rvwalfa,rvwbeta)
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
integer(4),      intent(IN)  :: npi
integer(4),      intent(IN)  :: npj
integer(4),      intent(IN)  :: npartint
double precision,intent(IN)  :: dervel(3)
double precision,intent(OUT) :: rvwalfa(3), rvwbeta(3)
!
!.. Local Scalars ..
double precision :: celtilde,rhotilde,amassj,vrij,TermMon
!
!.. Executable Statements ..
!
 vrij = -dervel(1) * rag(1,npartint) - dervel(2) * rag(2,npartint) - dervel(3) * rag(3,npartint)
!
!AA504 removed part: Monaghan's term is always active even if particle detach themselves 
   if ( pg(npj)%vel_type /= "std" ) then          !non part fix o altro
     amassj   = pg(npi)%mass
     rhotilde = two * pg(npi)%dens
     celtilde = Med(pg(npi)%imed)%celerita + Med(pg(npi)%imed)%celerita
   else
     amassj   = pg(npj)%mass
     rhotilde = pg(npi)%dens + pg(npj)%dens
!........................................... 2011 mar 15
     if (esplosione) then  !introdotto per stabilita' flushing
       celtilde = pg(npi)%Csound + pg(npj)%Csound
     else
       celtilde = Med(pg(npi)%imed)%celerita + Med(pg(npj)%imed)%celerita
     end if
!........................................... 2011 mar 15
   end if
!
!.. formula utilizzata nel nov 2007 come da letteratura
   TermMon = Med(pg(npi)%imed)%alfaMon * celtilde * Domain%h / rhotilde

!AA504 Molecular viscosity term depending on velocity divergence (compressible flows) can be neglected (and may causes several problems when activated)
   rvwalfa(1:3) = amassj * TermMon * &
                  vrij * rag(1:3,npartint) * PartKernel(2,npartint)

!.. rvwbeta e' da verificare
!   rvwbeta(1:3)= - amassj * two * Med(pg(npi)%imed)%betaMon * vrij*vrij * squareh / rhotilde & 
!                 * (PartKernel(2,npartint)/PartKernel(1,npartint)*(PartKernel(2,npartint)/PartKernel(1,npartint)))
   rvwbeta(1:3) = zero

!
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!if (it_corrente == 2 .and. npi > 16900 .and. npi < 17200) then
!  write (99,*) ' it = ',it_corrente,' n. particella = ',npi,npj
!  write (99,'(a,4d23.15)') ' viscomon   vrij,amassj,rhotilde,celtilde ',vrij,amassj,rhotilde,celtilde
!  write (99,'(a,3d23.15)') ' viscomon   rvwalfa ',rvwalfa
!end if
!$$$$$$$$$$$$$$$$!prova_arrotondamento
!
return
end subroutine viscomon
!---split

!cfile viscomorris.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : viscomorris
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
! Module purpose : Module to 
! 
! Calling routine: inter_EqMoto
!
! Called routines: 
!
!************************************************************************************
!
subroutine viscomorris (npi,npj,npartint,dervel,rvw)
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
integer(4),      intent(IN)  :: npi
integer(4),      intent(IN)  :: npj
integer(4),      intent(IN)  :: npartint
double precision,intent(IN)  :: dervel(3)
double precision,intent(OUT) :: rvw(3)
!
!.. Local Scalars ..
double precision :: amassj,rhotilde,anuitilde,factivis
!
!.. Executable Statements ..
!
 if ( pg(npj)%vel_type /= "std" ) then          !non part fix o altro
   amassj    = pg(npi)%mass
   rhotilde  = pg(npi)%dens
   anuitilde = two * pg(npi)%visc         !CONTROLLARE
 else
!   anuitilde = two * (pg(npi)%visc + pg(npj)%visc)
!   rhotilde  = (pg(npj)%dens + pg(npi)%dens)
   amassj = pg(npj)%mass
   rhotilde  = (pg(npi)%visc * pg(npi)%dens + pg(npj)%visc * pg(npj)%dens + 0.001d0)
   anuitilde = 4.0d0 * (pg(npi)%visc * pg(npj)%visc)      ! ATTENZIONE!! viscosita' cinematica
 end if
!
!.. formula utilizzata nel nov 2005 da rapporto Bon/Gatti/Zuccala/DiMonaco/Gallati
! eta     = 0.01d0 * hh
! gradd   = gradmod / (rij+eta)
! rvw(:)  =-gradd * dervel(:)
!.. formula utilizzata nel nov 2007 come da letteratura
! eta     = 0.01d0 * hh * hh
 factivis = amassj * anuitilde / rhotilde
 
 rvw(1:3) = factivis * &
            (-dervel(1:3) * PartKernel(2,npartint) * (rag(1,npartint) * rag(1,npartint) +  &
             rag(2,npartint) * rag(2,npartint) + rag(3,npartint) * rag(3,npartint)) )
!
return
end subroutine viscomorris
!---split

!cfile writime2.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine writime2 (ti,tf,nout)
!
use language_writime2
!
implicit none
!
double precision, dimension(2) :: ti,tf
integer(4)                     :: nout
!
! local arrays and scalars
logical, save          :: first=.true.
double precision       :: telat, tcput, tcpup, telap
double precision, save :: tcpus, telas
!
if ( first ) then
   tcpus = ti(1)
   telas = ti(2) !- ti(1)  ! toglie tcpus per aggiustare tempo elapsed (granularita' in secondi)
   tcpup = ti(1)
   telap = ti(2)
end if
first = .false.
!
tcput = tf(1)
telat = tf(2) - telas
!
write(nout,1001) " "
write(nout,1001) cpulbl,tcput,totlbl,tf(1)-ti(1),przlbl
write(nout,1001) elalbl,max(1.0d0,telat),totlbl,max(1.0d0,tf(2)-ti(2)),przlbl
!
tcpup = tf(1)
telap = tf(2)
!
return
 1001 format(1x,a,f20.4,1x,a,f20.4,1x,a)
end subroutine writime2
!---split



