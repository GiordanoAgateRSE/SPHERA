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

