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

