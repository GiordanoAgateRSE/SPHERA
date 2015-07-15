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

