!cfile BC_wall_elements.f90
!AA406 all the subroutines
!AA501 all the subroutine cycles are parallelized
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : BC_wall_elements
!
! Last updating : April 18, 2013
!
! Improvement traceback:
! 00 Amicarelli/Agate  30Nov11   First development
! 01 Amicarelli/Agate  18apr13   minor modifications
!
!************************************************************************************
! Module purpose : Estimating wall element density and pressure
!                  (Boundary SPH approach)
!
! Calling routine: Loop_Irre_2D,Loop_Irre_3D
!
!************************************************************************************
!
  subroutine BC_wall_elements
!
! Assigning modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
!
! Declarations
  implicit none
  character(len=lencard)  :: nomsub = "BC_wall_elements"
  integer(4)              :: contj,npartint,npi,npj
  double precision        :: rhoL,uL,uR,cL,rhostar,rhorif,c2,pstar,Ww_Shep
  double precision        :: uCartL(3)
  double precision, dimension(:), allocatable :: den
! Logic array to detect wall elements with fluid neighbours
  integer(4), dimension(:), allocatable :: neigh_w
!
! Statements
!
! Allocating and initializing the auxiliary variable
!AA601 sub start
  allocate(den(DBSPH%n_w++DBSPH%n_inlet++DBSPH%n_outlet))
  allocate(neigh_w(DBSPH%n_w++DBSPH%n_inlet++DBSPH%n_outlet))
!AA601 sub end  
  den = zero
  neigh_w = 0
!
!AA503 sub omp declarations
! free-surface condition (a sort of zeroing) for wall pressure and density
!$omp parallel do             &
!$omp default(none)           &
!$omp private(npi)            &
!$omp shared (DBSPH,pg_w,med)
  do npi=1,DBSPH%n_w
!
!AA406start
    pg_w(npi)%pres = zero
!AA503 sub 
    pg_w(npi)%dens = med(1)%den0
!    pg_w(npi)%wet = 0
!AA406end
!
  end do
!$omp end parallel do
!
!AA501test no parallelization of this cycle, because the critical sections are too slow 
!nomp parallel do   &
!nomp default(none) &
!nomp private(npi,contj,npj,rhoL,npartint,uCartL,uL,uR,cL,rhostar,rhorif,c2,pstar,Ww_Shep)   &
!nomp shared(nag,DBSPH,nPartIntorno_fw,NMAXPARTJ,neigh_w,pg,pg_w,rag_fw,Med,kernel_fw,den,PartIntorno_fw)
!
! Loop over computational fluid particles for density and pressure contributions to wall elements
  do npi=1,nag
!AA601
     if (pg(npi)%imed/=pg(1)%imed) cycle
! Loop over the wall neighbnours of the computational particle
     do contj=1,nPartIntorno_fw(npi)
        npartint = (npi-1)* NMAXPARTJ + contj
        npj = PartIntorno_fw(npartint)
! Zeroing wall density when a wall element has fluid neighbours (just the first time it is encountered)
        if (neigh_w(npj) == 0) then
!AA406test
           pg_w(npj)%dens = zero
!
           neigh_w(npj) = 1
        endif
! MUSCL reconstruction (first order): provides the initial left side conditions for the Riemann solver (subscripts L)
! The right state (R) is defined by the wall element properties
! ucartL represents the left state condition in Cartesian coordinates: uL is its projection over the wall normal
        rhoL=pg(npi)%dens-(pg(npi)%drho(1)*rag_fw(1,npartint)+pg(npi)%drho(2)*rag_fw(2,npartint)+ &
              pg(npi)%drho(3)*rag_fw(3,npartint))
        uCartL(1)=pg(npi)%var(1)-(pg(npi)%dvel(1,1)*rag_fw(1,npartint)+pg(npi)%dvel(1,2)*rag_fw(2,npartint)+ &
                   pg(npi)%dvel(1,3)*rag_fw(3,npartint))
        uCartL(2)=pg(npi)%var(2)-(pg(npi)%dvel(2,1)*rag_fw(1,npartint)+pg(npi)%dvel(2,2)*rag_fw(2,npartint)+ &
                   pg(npi)%dvel(2,3)*rag_fw(3,npartint))
        uCartL(3)=pg(npi)%var(3)-(pg(npi)%dvel(3,1)*rag_fw(1,npartint)+pg(npi)%dvel(3,2)*rag_fw(2,npartint)+ &
                   pg(npi)%dvel(3,3)*rag_fw(3,npartint))
        uL=uCartL(1)*pg_w(npj)%normal(1)+uCartL(2)*pg_w(npj)%normal(2)+uCartL(3)*pg_w(npj)%normal(3)
! Partial Riemann solver start
! The partial Riemann problem is here 1D, oriented along the wall normal 
! Solution (density) in the central zone (here the wall element): subscript star
        uR = pg_w(npj)%vel(1)*pg_w(npj)%normal(1)+pg_w(npj)%vel(2)*pg_w(npj)%normal(2)+pg_w(npj)%vel(3)*pg_w(npj)%normal(3)
        cL = Dsqrt(Med(pg(npi)%imed)%eps/rhoL)
!AA406 test start
!        rhoL = pg(npi)%dens
!        cL = Dsqrt(Med(pg(npi)%imed)%eps/rhoL) 
!        uL= pg(npi)%var(1)*pg_w(npj)%normal(1)+pg(npi)%var(2)*pg_w(npj)%normal(2)+pg(npi)%var(3)*pg_w(npj)%normal(3)
!        uR = zero
!AA406 test end
        rhostar = rhoL + (uL - uR) * rhoL / cL
!AA406 test        
!        if (rhoL == 1000.) rhostar = rhoL
!
        if (rhostar < 10.) then
           rhostar = 10.
        endif
! state equation for pressure
        rhorif = Med(pg(npi)%imed)%den0
        c2 = Med(pg(npi)%imed)%eps/rhorif
        pstar = c2 * (rhostar-rhorif)
        if (pstar < -99000.) then
           pstar = -99000.
        endif 
!
!AA406test 
!        kernel_fw(2,npartint) = rhostar
!        kernel_fw(3,npartint) = pstar
!
! partial Riemann solver end
!
! den(pnj) is an auxiliary vector to add contributions to the density and pressure denominators
! Ww_Shep: wall value weight 
!AA406test
        Ww_Shep = pg(npi)%mass / pg(npi)%dens * kernel_fw(1,npartint) 
!        Ww_Shep = pg(npi)%mass / pg(npi)%dens * kernel_fw(3,npartint) 
!         Ww_Shep = kernel_fw(1,npartint) / pg(npi)%uni 
!Ww_Shep = kernel_fw(3,npartint) / pg(npi)%uni 
!  Ww_Shep = kernel_fw(1,npartint) 
!
!AA501test
!nomp critical
!
!AA406test
        pg_w(npj)%dens = pg_w(npj)%dens + rhostar * Ww_Shep
        den(npj) = den(npj) + Ww_Shep
        pg_w(npj)%pres = pg_w(npj)%pres + pstar * Ww_Shep
!
!AA501test
!nomp end critical
!
     end do
  end do
!AA501test
!nomp end parallel do
!
! Updating wall element pressure and density
!$omp parallel do   &
!$omp default(none) &
!$omp private(npi)  &
!$omp shared(neigh_w,pg_w,den,DBSPH)
  do npi=1,DBSPH%n_w
     if (neigh_w(npi) == 1) then
!AA406test          
        pg_w(npi)%dens = pg_w(npi)%dens / den(npi)
        pg_w(npi)%pres = pg_w(npi)%pres / den(npi)
!
!AA501 !non negative pressure values
        if (pg_w(npi)%pres < 0.) pg_w(npi)%pres = zero 
!
     endif
  end do
!$omp end parallel do
!
! Deallocating the auxiliary variable
  deallocate(den)
  deallocate(neigh_w) 
  return
  end subroutine BC_wall_elements
!---split

!cfile viscomorris_wall_elements.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : viscomorris_wall_elements
!
! Last updating : Amicarelli/Agate  30 Nov 2011
!
!************************************************************************************
! Module purpose : Wall element contributions to Morris' viscosity term 
! 
! Calling routine: inter_EqMoto
!
! Called routines: /
!
!************************************************************************************
!
subroutine viscomorris_wall_elements(npi,npj,npartint,dervel,rvw)
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
! Declarations
! Implicit
implicit none
! Formal Arguments
integer(4),      intent(IN)  :: npi,npj,npartint
double precision,intent(IN)  :: dervel(3)
double precision,intent(OUT) :: rvw(3)
! Local Scalars
double precision :: rhoWw,rhotilde,anuitilde,factivis,dis2
!
! Statements
 rhoWw = pg_w(npj)%dens * pg_w(npj)%weight * kernel_fw(1,npartint) 
 rhotilde  = (pg(npi)%visc * pg(npi)%dens + pg(npi)%visc * pg_w(npj)%dens + 0.001d0)
 anuitilde = 4.0d0 * (pg(npi)%visc * pg(npi)%visc)      ! ATTENZIONE!! viscosita' cinematica
 factivis = rhoWw * anuitilde / rhotilde
 dis2 = (rag_fw(1,npartint)*rag_fw(1,npartint)+rag_fw(2,npartint)*rag_fw(2,npartint)+rag_fw(3,npartint)*rag_fw(3,npartint))
 rvw(:) = factivis * dervel(:) * (rag_fw(:,npartint)*pg_w(npj)%normal(:)) / dis2
!
return
end subroutine viscomorris_wall_elements
!---split

!cfile viscomon_wall_elements.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : viscomon_wall_elements
!
! Last updating : Amicarelli/Agate  30 Nov 2011
!
!************************************************************************************
! Module purpose : Wall element contributions of Monaghan's viscosity term 
! 
! Calling routine: inter_EqMoto
!
! Called routines: /
!
!************************************************************************************
!
subroutine viscomon_wall_elements(npi,npj,npartint,dervel,rvwalfa,rvwbeta)
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
integer(4),      intent(IN)  :: npi,npj,npartint
double precision,intent(IN)  :: dervel(3)
double precision,intent(OUT) :: rvwalfa(3), rvwbeta(3)
!
!.. Local Scalars ..
double precision :: rhoWw,rhotilde,celtilde,vrij,TermMon,dis2
!
!.. Executable Statements ..
!
 vrij = -dervel(1) * rag_fw(1,npartint) - dervel(2) * rag_fw(2,npartint) - dervel(3) * rag_fw(3,npartint)
 dis2 = (rag_fw(1,npartint)*rag_fw(1,npartint)+rag_fw(2,npartint)*rag_fw(2,npartint)+rag_fw(3,npartint)*rag_fw(3,npartint))
 if (vrij > zero) then
   rvwalfa = zero
   rvwbeta = zero
   else
      rhotilde = pg(npi)%dens + pg_w(npj)%dens
      celtilde = 2. * Med(pg(npi)%imed)%celerita 
      rhoWw = pg_w(npj)%dens * pg_w(npj)%weight * kernel_fw(1,npartint) 
      TermMon = Med(pg(npi)%imed)%alfaMon * celtilde * Domain%h / rhotilde
      rvwalfa(1:3) = rhoWw * TermMon * vrij * pg_w(npj)%normal(:) / dis2
      rvwbeta(1:3) = zero
 end if
!
return
end subroutine viscomon_wall_elements
!---split

!cfile Gradients_to_MUSCL.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : Gradients_to_MUSCL
!
! Last updating : November 30, 2011
!
! Creation : Amicarelli/Agate  30 Nov 2011
!
!************************************************************************************
! Module purpose : Pseudo consistent estimation, with boundary contributions and Shepard kernel,
!                  of the velocity and density gradients for the MUSCL 
!                  reconstruction (to feed the partial Riemann solver). 
!
!AA501
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!AA601 sub
! Called subroutines: DBSPH_inlet_outlet
!
!************************************************************************************
!
subroutine Gradients_to_MUSCL
!
!.. assign modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE
!
! Implicit Declarations
implicit none
!
! Local Parameters
integer(4)       :: npi,npj,npartint,contj
double precision :: vol_Shep,Ww_Shep
!
! Executable Statements
!AA601 sub
!$omp parallel do  &
!$omp default(none) &
!$omp private(npi,contj,npj,npartint,vol_Shep) &
!$omp shared (nag,pg,NMAXPARTJ,rag,nPartIntorno,Partintorno,PartKernel,nPartIntorno_fw,DBSPH)
 do npi = 1, nag ! loop over the fluid computational particles: start
!
!AA501
    if (nPartIntorno_fw(npi) > 0) then   
!    
! Initializing the gradients 
    pg(npi)%drho = zero
    pg(npi)%dvel = zero
!
    do contj = 1, nPartIntorno(npi)   !loop over the fluid neighbouring particles: start
       npartint = (npi-1)* NMAXPARTJ + contj
       npj = PartIntorno(npartint)
!AA601
       if (pg(npi)%imed==pg(npj)%imed) then
!AA406!!!test
          vol_Shep = pg(npj)%mass / pg(npj)%dens 
!       vol_Shep = pg(npj)%mass / pg(npj)%dens / pg(npi)%sigma
! Computation of the density gradient: fluid particle contributions
          pg(npi)%drho(:) = pg(npi)%drho(:) + rag(:,npartint) * (pg(npj)%dens-pg(npi)%dens) &
                            * Partkernel(1,npartint)  * vol_Shep
! Computation of the velocity gradient: fluid particle contributions
          pg(npi)%dvel(1,:) = pg(npi)%dvel(1,:) + rag(:,npartint) * (pg(npj)%var(1) - pg(npi)%var(1)) * &
                              Partkernel(1,npartint) * vol_Shep
          pg(npi)%dvel(2,:) = pg(npi)%dvel(2,:) + rag(:,npartint) * (pg(npj)%var(2) - pg(npi)%var(2)) * &
                              Partkernel(1,npartint) * vol_Shep
          pg(npi)%dvel(3,:) = pg(npi)%dvel(3,:) + rag(:,npartint) * (pg(npj)%var(3) - pg(npi)%var(3)) * &
                              Partkernel(1,npartint) * vol_Shep
!AA601 
       endif
    end do !loop over the fluid neighbouring particles: end
!AA601 sub
    if (DBSPH%MUSCL_boundary_flag==1) call Gradients_to_MUSCL_boundary(npi)
!AA501
    endif
! 
 end do ! loop over the fluid computational particles: end
!$omp end parallel do
!
return
!
return
end subroutine Gradients_to_MUSCL
!---split



!AA601 sub the whole subroutine
!cfile BC_wall_elements.f90
!************************************************************************************
!                          S P H E R A 5.0.3 
!
!                Smoothed Particle Hydrodinamics Code
!
!************************************************************************************
!
! File name     : Gradients_to_MUSCL_boundary
!
! Creation      : Amicarelli A., 30Nov11 (hard-coding)
! Updates       : Amicarelli A., 26Jan15 (integration with no hard-coding)
!
!************************************************************************************
! Module purpose : Boundary terms for the MUSCL reconstruction (DBSPH), if required in input
!
! Calling subroutines: Gradients_to_MUSCL
! Called subroutines:  /  
!
!************************************************************************************

subroutine Gradients_to_MUSCL_boundary(npi)

! Modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE

! Declarations ..
implicit none
integer(4), intent(in) :: npi
integer(4) :: contj,npartint,npj
double precision :: vol_Shep,Ww_Shep

! Statements
! Loop over the wall neighbouring elements: start (parallelization has no sense for this loop)
do contj=1,nPartIntorno_fw(npi)   
   npartint = (npi-1)* NMAXPARTJ + contj
   npj = PartIntorno_fw(npartint)
! Boundary term contributions for the density and velocity gradients
   Ww_Shep = pg_w(npj)%weight * kernel_fw(1,npartint) / pg(npi)%sigma
! Computation of the density gradient: wall contributions
   pg(npi)%drho(:) = pg(npi)%drho(:) +  (pg_w(npj)%dens-pg(npi)%dens) * pg_w(npj)%normal(:) * Ww_Shep
! Computation of the velocity gradient: wall contributions
   pg(npi)%dvel(1,:) = pg(npi)%dvel(1,:) + (pg_w(npj)%vel(1) - pg(npi)%var(1)) * pg_w(npj)%normal(:) * Ww_Shep
   pg(npi)%dvel(2,:) = pg(npi)%dvel(2,:) + (pg_w(npj)%vel(2) - pg(npi)%var(2)) * pg_w(npj)%normal(:) * Ww_Shep
   pg(npi)%dvel(3,:) = pg(npi)%dvel(3,:) + (pg_w(npj)%vel(3) - pg(npi)%var(3)) * pg_w(npj)%normal(:) * Ww_Shep
! Semi-particles
   vol_Shep = pg_w(npj)%mass / pg_w(npj)%dens / pg(npi)%sigma
! Computation of the density gradient: wall contributions
   pg(npi)%drho(:) = pg(npi)%drho(:) + rag_fw(:,npartint) * (pg_w(npj)%dens-pg(npi)%dens) &
                     * kernel_fw(2,npartint)  * vol_Shep
! Computation of the velocity gradient: wall contributions
   pg(npi)%dvel(1,:) = pg(npi)%dvel(1,:) + rag_fw(:,npartint) * (pg_w(npj)%vel(1) - pg(npi)%var(1)) * &
                       kernel_fw(2,npartint) * vol_Shep
   pg(npi)%dvel(2,:) = pg(npi)%dvel(2,:) + rag_fw(:,npartint) * (pg_w(npj)%vel(2) - pg(npi)%var(2)) * &
                       kernel_fw(2,npartint) * vol_Shep
   pg(npi)%dvel(3,:) = pg(npi)%dvel(3,:) + rag_fw(:,npartint) * (pg_w(npj)%vel(3) - pg(npi)%var(3)) * &
                       kernel_fw(2,npartint) * vol_Shep
end do   
! Loop over the wall neighbouring elements: end

return
end subroutine Gradients_to_MUSCL_boundary
!---split



!AA601 the whole subroutine
!cfile BC_wall_elements.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : DBSPH_inlet_outlet
!
! Creation      : Amicarelli, 26Jan15
!
!************************************************************************************
! Module purpose : Impose boundary conditions at the inlet and outlet sections (DBSPH) 
!
! Calling routines: inter_SmoothVelo_2D,inter_SmoothVelo_3D,Euler
!
! Called routines: /
!
!************************************************************************************
subroutine DBSPH_inlet_outlet(npi)

! Modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE

! Declarations 
implicit none
integer(4),intent(in)           :: npi

!Initializations

!Statements
! By hypothesis a fluid particle cannot be close to two or more inlet/outlet sections (otherwise it randomly takes an inlet velocity value)
! Inlet boundary conditions
if (pg(npi)%DBSPH_inlet_ID>0) pg(npi)%vel(:) = pg_w(pg(npi)%DBSPH_inlet_ID)%vel(:)
! Outlet boundary conditions
if (pg(npi)%DBSPH_outlet_ID>0) then
   pg(npi)%dens = Med(1)%den0
   pg(npi)%pres = pg_w(pg(npi)%DBSPH_outlet_ID)%pres
endif

!Deallocations

return
end subroutine DBSPH_inlet_outlet
!---split



!AA601 sub the whole subroutine
!cfile BC_wall_elements.f90
!************************************************************************************
!                          S P H E R A 5.0.3 
!
!                Smoothed Particle Hydrodinamics Code
!
!************************************************************************************
!
! File name     : wall_elements_pp
!
! Creation      : Amicarelli A.; 30Nov11
! Update        : Amicarelli A.; 26Jan15 
!                 Integration with no hard-coding; parallelization; 
!                 post-processing the wall surface element values (provided a selected region)
!                 post-processing the hydrodynamic normal force on DBSPH surface elements (provided a selected region)
!                 post-processing the wall surface element values (provided selected element IDs)
!
!************************************************************************************
! Module purpose : Smoothing wall element values for post-processing; initializing
!
! Calling routine: Memo_Ctl
!
!************************************************************************************

subroutine wall_elements_pp

! Modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE

! Declarations
implicit none
integer(4)       :: igridi,kgridi,jgridi,irang,krang,jrang,ncelj,jgrid1,jgrid2,i
integer(4)       :: npi,npj,npartint,index_rij_su_h,irestocell,ww
double precision :: rij_su_h,ke_coef,rij_su_h_quad,rijtemp,rijtemp2,wu,denom,gradmod
double precision :: ragtemp(3)
double precision, dimension(:), allocatable :: pres_wpp,den_wpp
character(255)   :: nomefilectl_wall
integer(4),external :: CellIndices, CellNumber
double precision :: Fx

! Initializations
Fx = 0.
ke_coef = Domain%coefke / Domain%h
allocate (pres_wpp(DBSPH%n_w))
allocate (den_wpp(DBSPH%n_w))
pres_wpp = zero
den_wpp = zero

! Statements
! Loop over the wall elements
!$omp parallel do default(none) &
!$omp shared(DBSPH,pg_w,ncord,Icont_w,NPartOrd_w,pres_wpp,Domain,den_wpp,doublesquareh,squareh) &
!$omp private(npi,irestocell,igridi,jgridi,jgrid1,jgrid2,kgridi,jrang,irang,krang,ncelj,ww,npj,rijtemp,ragtemp,rijtemp2,rij_su_h,rij_su_h_quad,index_rij_su_h,wu)
do npi = 1,DBSPH%n_w
   if (pg_w(npi)%cella == 0) cycle
   irestocell = CellIndices (pg_w(npi)%cella,igridi,jgridi,kgridi)
! Indici inizio-fine scansione Y per casi 2d o 3d
   jgrid1 = jgridi - (ncord-2)
   jgrid2 = jgridi + (ncord-2)
   loop_jrang: do jrang = jgrid1,jgrid2    ! ---- a  loop sulle 9 celle 
      loop_irang: do irang = igridi-1,igridi+1    ! ---- b  loop sulle 9 celle  
        loop_krang: do krang = kgridi-1,kgridi+1    ! ---- c  loop sulle 9 celle  
           ncelj = CellNumber (irang,jrang,krang)
           if (ncelj == 0) cycle    ! cella fuori campo
           if (Icont_w(ncelj+1) <= Icont_w(ncelj)) cycle
! Loop over the neighbouring wall particles in the cell
           do ww = Icont_w(ncelj),Icont_w(ncelj+1)-1
              npj = NPartOrd_w(ww)
! Relative positions and distances
              ragtemp(1:3) = pg_w(npi)%coord(1:3) - pg_w(npj)%coord(1:3)
              rijtemp = ragtemp(1)*ragtemp(1) + ragtemp(2)*ragtemp(2) + ragtemp(3)*ragtemp(3)
! Distance check
              if (rijtemp > doublesquareh) cycle
              rijtemp2 = rijtemp
              rij_su_h = Dsqrt(rijtemp) / Domain%h
              rij_su_h_quad = rijtemp2 / squareh
              index_rij_su_h = int(rij_su_h)
! Kernel computation
              if (index_rij_su_h >= 2) cycle
              wu = 0.666666667d0 + rij_su_h_quad * (rij_su_h * half - one)
              if (index_rij_su_h > 0) wu = (two - rij_su_h) * (two - rij_su_h) * (two - rij_su_h) * 0.166666666667d0
              if (pg_w(npj)%wet==1) then
                 pres_wpp(npi) = pres_wpp(npi) + wu * Domain%coefke * pg_w(npj)%pres * pg_w(npj)%weight
                 den_wpp(npi)  = den_wpp(npi) +  wu * Domain%coefke * pg_w(npj)%weight
              endif
           end do
        end do loop_krang   ! ---- c  loop sulle 9 celle    
      end do loop_irang   ! ---- b  loop sulle 9 celle    
   end do loop_jrang   ! ---- a  loop sulle 9 celle
end do 
!$omp end parallel do
! Loop over the wall particles to update the pressure values for post-processing
!$omp parallel do default(none) shared(DBSPH,den_wpp,pres_wpp) private(npi)
do npi=1,DBSPH%n_w
   if (den_wpp(npi) > 0.) pres_wpp(npi) = pres_wpp(npi) / den_wpp(npi) 
end do
!$omp end parallel do
! Writing the pressure force (post-processing) (provided a region)
if (DBSPH%n_monitor_regions==1) then
   write(nomefilectl_wall,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_wall_Fx_',it_corrente,".txt"
   open (unit_dbsph_Fx, file=nomefilectl_wall, status="unknown", form="formatted" )
   write (unit_dbsph_Fx,*) "Force "
   write (unit_dbsph_Fx,'(1x,2(a,10x))') " Time(s)","Fx(kgm/s^2)"
   call flush(unit_dbsph_Fx)
   write (unit_dbsph_Fx,*) "Force at boundaries (DBSPH)"
   do npi=1,DBSPH%n_w
      if ( (pg_w(npi)%coord(1)>=DBSPH%monitor_region(1)) .and. (pg_w(npi)%coord(1)<=DBSPH%monitor_region(2)) .and. &
           (pg_w(npi)%coord(2)>=DBSPH%monitor_region(3)) .and. (pg_w(npi)%coord(2)<=DBSPH%monitor_region(4)) .and. &
           (pg_w(npi)%coord(3)>=DBSPH%monitor_region(5)) .and. (pg_w(npi)%coord(3)<=DBSPH%monitor_region(6))        ) then
         if (pg_w(npi)%normal(1)/=0.) Fx=Fx+pres_wpp(npi)*pg_w(npi)%normal(1)*pg_w(npi)%weight
      endif
   end do
   write (unit_dbsph_Fx,'(2(1x,g14.7))') tempo,Fx
   close (unit_dbsph_Fx)
endif
! Writing the wall element pressure values derived from post-processing (provided a region)
if (DBSPH%n_monitor_regions==1) then
   write(nomefilectl_wall,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_wall_region_',it_corrente,".txt"
   open (unit_dbsph_se_reg, file=nomefilectl_wall, status="unknown", form="formatted" )
   write (unit_dbsph_se_reg,*) "Wall element values "
   write (unit_dbsph_se_reg,'(1x,2(a,10x),3(a,8x),3(a,5x),a,7x,a)') &
         " Time","Iter","ID","X Coord","Y Coord","Z Coord","X Velocity","Y Velocity","Z Velocity"," Pressure"," Pressure_pp","wet"
   call flush(unit_dbsph_se_reg)
   write (unit_dbsph_se_reg,*) "wall_elements"
   do npi=1,DBSPH%n_w
      if ( (pg_w(npi)%coord(1)>=DBSPH%monitor_region(1)) .and. (pg_w(npi)%coord(1)<=DBSPH%monitor_region(2)) .and. &
           (pg_w(npi)%coord(2)>=DBSPH%monitor_region(3)) .and. (pg_w(npi)%coord(2)<=DBSPH%monitor_region(4)) .and. &
           (pg_w(npi)%coord(3)>=DBSPH%monitor_region(5)) .and. (pg_w(npi)%coord(3)<=DBSPH%monitor_region(6))        ) then
         write (unit_dbsph_se_reg,'(g14.7,2(i14),8(1x,g14.7),i3)') tempo,it_corrente,i,pg_w(npi)%coord(:),pg_w(npi)%vel(:), &
                                                                   pg_w(npi)%pres,pres_wpp(npi),pg_w(npi)%wet 
      endif
   end do
   close (unit_dbsph_se_reg)
endif
! Writing the wall element pressure values derived from post-processing (provided the element IDs)
if (DBSPH%n_monitor_points>0) then
   write(nomefilectl_wall,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_wall_IDs_',it_corrente,".txt"
   open (unit_dbsph_se_ID, file=nomefilectl_wall, status="unknown", form="formatted" )
   write (unit_dbsph_se_ID,*) "Wall element values "
   write (unit_dbsph_se_ID,'(1x,2(a,10x),3(a,8x),3(a,5x),a,7x,a)') &
      " Time","Iter","ID","X Coord","Y Coord","Z Coord","X Velocity","Y Velocity","Z Velocity"," Pressure"," Pressure_pp","wet"
   call flush(unit_dbsph_se_ID)
   write (unit_dbsph_se_ID,*) "wall_elements"
   do i=1,DBSPH%n_monitor_points
      write (unit_dbsph_se_ID,'(g14.7,2(i14),8(1x,g14.7),i3)') tempo,it_corrente,i,pg_w(DBSPH%monitor_IDs(i))%coord(:),pg_w(DBSPH%monitor_IDs(i))%vel(:), &
                                                            pg_w(DBSPH%monitor_IDs(i))%pres,pres_wpp(DBSPH%monitor_IDs(i)),pg_w(DBSPH%monitor_IDs(i))%wet           
   end do   
   close (unit_dbsph_se_ID)
endif

! Deallocations
deallocate (pres_wpp)
deallocate (den_wpp)

return
end subroutine wall_elements_pp
!---split



!AA601 the whole subroutine
!cfile BC_wall_elements.f90
!************************************************************************************
!                          S P H E R A 5.0.3 
!
!                Smoothed Particle Hydrodinamics Code
!
!************************************************************************************
!
! Subroutine : wavy_inlet
!
! Creation   : Amicarelli A., 26Jan15; Main features (see Purpose)
!
!************************************************************************************
! Purpose : It provides a wavy flow at the inlet section: each layer is staggered 
!           by 0.5dx with respect to the previous and the following ones, which are instead aligned.
!           This is a numerical feature to reduce the SPH truncation error (necessary for DBSPH jets).
!
! Calling routine: GenerateSourceParticles_2D,GenerateSourceParticles_3D
! Called subroutines : / 
!
!************************************************************************************
subroutine wavy_inlet(i_inlet)

! Modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE

! Declarations
implicit none
!AA601!!! sub
integer(4),intent(in) :: i_inlet
integer(4) :: npart1,npart2,i
!AA601!!! sub
double precision :: Length1,Length2,rnd
double precision,dimension(3) :: cos_dir_1

!Initializations
Length1 = 0.d0
Length2 = 0.d0
npart1 = 0
npart2 = 0

! Statements
if (ncord==2) then
   npart1 = NumPartperLine(i_inlet)
   cos_dir_1(:) = BoundarySide(i_inlet)%T(:,1)
   else
! Length1 and Length 2 are the length scales of the inlet section: they are computed as the distance between the first and the last inlet vertices 
! and the third and the last inlet vertices, respectively.
! Particles are aligned with Plast-P1 and Plast-P3, where P1 is the first boundary vertex, ... and Plast the last boundary vertex. 
! In case of a triangular inlet, we have particles aligned with one direction: P3-P1.
! In case of a quadrilateral inlet, we have particles distributed along two directions: P4-P1 and P4-P3.
      do i=1,3
         Length1 = Length1 + (BoundaryFace(i_inlet)%Node(1)%GX(i) - BoundaryFace(i_inlet)%Node(BoundaryFace(i_inlet)%nodes)%GX(i))**2
         Length2 = Length2 + (BoundaryFace(i_inlet)%Node(3)%GX(i) - BoundaryFace(i_inlet)%Node(BoundaryFace(i_inlet)%nodes)%GX(i))**2
      end do
      Length1 = Dsqrt(Length1)
      Length2 = Dsqrt(Length2)
      npart1 = Int(Length1/Domain%dd+0.01d0)
      npart2 = Int(Length2/Domain%dd+0.01d0)
      npart1 = npart1*npart2 
      cos_dir_1(:) = BoundaryFace(i_inlet)%T(:,1)
   endif
! (ID particle = nag) indicates the last generated particle (of the on-going inlet section)
   select case (mod(itime_jet,4))
      case (1,3)  
         pg(nag)%coord(:) = pg(nag)%coord(:) - 0.25d0 * Domain%dd * cos_dir_1(:)
      case (2,0) 
         pg(nag)%coord(:) = pg(nag)%coord(:) + 0.25d0 * Domain%dd * cos_dir_1(:)
   end select
call random_number(rnd)
pg(nag)%coord(:) = pg(nag)%coord(:) + (two * rnd - one) * 0.1d0 * Domain%dd * cos_dir_1(:)

return
end subroutine wavy_inlet
!---split



!AA601 the whole subroutine
!cfile BC_wall_elements.f90
!************************************************************************************
!                          S P H E R A 5.0.3 
!
!                Smoothed Particle Hydrodinamics Code
!
!************************************************************************************
!
! Subroutine : Import_ply_surface_meshes
!
! Creation   : Amicarelli A., 26Jan15; Main features (see Purpose)
!
!************************************************************************************
! Purpose : It imports the surface meshes (generated by SnappyHexMesh -OpenFoam-) and 
!           then converted by Paraview into .ply. This subroutine is mandatory and 
!           activated only for DB-SPH.
! Calling procedures: gest_input
! Called procedures : area_triangle
!************************************************************************************

subroutine Import_ply_surface_meshes

! Modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE

! Declarations
implicit none
integer(4) :: file_stat,n_vertices,old_size_vert,old_size_face,new_size_vert,new_size_face,dealloc_stat,alloc_stat,n_faces,face_vert_num,i,j,k
integer(4) :: aux_face_vert(4)
character(80) :: file_name,aux_char_1,aux_char_2
type(vertex_der_type),dimension(:),allocatable :: aux_der_type_vert
type(face_der_type),dimension(:),allocatable :: aux_der_type_faces

! Interface blocks
interface
   subroutine area_triangle(P1,P2,P3,area,normal)
      implicit none
      double precision,intent(IN)    :: P1(3),P2(3),P3(3)
      double precision,intent(OUT)   :: area
      double precision,intent(OUT)   :: normal(3)
   end subroutine area_triangle
end interface
! Allocations
! Initializations
new_size_face = 0
! Statements
!Open the file name list (surface_mesh_list.inp)
open(unit_file_list,file="surface_mesh_list.inp",IOSTAT=file_stat)
if (file_stat/=0) then
   write(*,*) 'Error in opening surface_mesh_list.inp in Import_pl_surface_meshes; the program terminates here'
   stop
endif
read(unit_file_list,*,IOSTAT=file_stat) 
if (file_stat/=0) then
   write(*,*) 'Error in reading surface_mesh_list.inp in Import_pl_surface_meshes; the program terminates here'
   stop
endif
do
! Read the file name    
   read (unit_file_list,'(a)',IOSTAT=file_stat) file_name
   if (file_stat/=0) exit ! Exit the cicle at the end of file
   file_name = trim(file_name)
! Open the on-going mesh file
   open(unit_DBSPH_mesh,file=file_name,IOSTAT=file_stat)
   if (file_stat/=0) then
      write(*,*) 'Error in opening a surface mesh file in Import_pl_surface_meshes; the program terminates here'
      stop
   endif   
   read(unit_DBSPH_mesh,"(3/)",IOSTAT=file_stat)
   if (file_stat/=0) then
      write(*,*) 'Error in reading a surface mesh file in Import_pl_surface_meshes; the program terminates here'
      stop
   endif   
! Read the number of vertices in the file    
   read(unit_DBSPH_mesh,*,IOSTAT=file_stat) aux_char_1,aux_char_2,n_vertices
   if (file_stat/=0) then
      write(*,*) 'Error in reading a surface mesh file in Import_pl_surface_meshes; the program terminates here'
      stop
   endif
   read(unit_DBSPH_mesh,"(2/)",IOSTAT=file_stat)
   if (file_stat/=0) then
      write(*,*) 'Error in reading a surface mesh file in Import_pl_surface_meshes; the program terminates here'
      stop
   endif
! Read the number of faces in the file
   read(unit_DBSPH_mesh,*,IOSTAT=file_stat) aux_char_1,aux_char_2,n_faces
   if (file_stat/=0) then
      write(*,*) 'Error in reading a surface mesh file in Import_pl_surface_meshes; the program terminates here'
      stop
   endif   
   read(unit_DBSPH_mesh,"(1/)",IOSTAT=file_stat)
   if (file_stat/=0) then
      write(*,*) 'Error in reading a surface mesh file in Import_pl_surface_meshes; the program terminates here.'
      stop
   endif   
   if (.not.allocated(DBSPH%surf_mesh%vertices)) then
      allocate(DBSPH%surf_mesh%vertices(n_vertices),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(nout,*) 'Allocation of DBSPH%surf_mesh%vertices in Import_ply_surface_mesh failed; the program terminates here.'
         stop ! Stop the main program
      endif
      if (.not.allocated(aux_der_type_vert)) then
          allocate(aux_der_type_vert(n_vertices),STAT=alloc_stat)
          if (alloc_stat/=0) then
             write(nout,*) 'Allocation of aux_der_type_vert in Import_ply_surface_mesh failed; the program terminates here.'
             stop ! Stop the main program
          endif
      endif    
      old_size_vert = 0
      else
         old_size_vert = size(DBSPH%surf_mesh%vertices)
         new_size_vert = old_size_vert + n_vertices
         aux_der_type_vert(:) = DBSPH%surf_mesh%vertices(:)
         deallocate(DBSPH%surf_mesh%vertices,STAT=dealloc_stat)
         if (dealloc_stat/=0) then
            write(*,*) 'Deallocation of DBSPH%surf_mesh%vertices in Import_ply_surface_mesh failed; the program terminates here.'
            stop ! Stop the main program
         endif          
         allocate(DBSPH%surf_mesh%vertices(new_size_vert),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(*,*) 'Allocation of DBSPH%surf_mesh%vertices in Import_ply_surface_mesh failed; the program terminates here'
            stop ! Stop the main program
         endif         
         DBSPH%surf_mesh%vertices(1:old_size_vert) = aux_der_type_vert(:)
         if (allocated(aux_der_type_vert)) then
            deallocate(aux_der_type_vert,STAT=dealloc_stat)
            if (dealloc_stat/=0) then
               write(nout,*) 'Deallocation of aux_der_type_vert in Import_ply_surface_mesh failed; the program terminates here.'
               stop ! Stop the main program
            endif   
         endif         
         allocate(aux_der_type_vert(new_size_vert),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(nout,*) 'Allocation of aux_der_type_vert in Import_ply_surface_mesh failed; the program terminates here.'
            stop ! Stop the main program
         endif
   endif
! Read the vertex coordinates: start      
   do j=(old_size_vert+1),(old_size_vert+n_vertices)
      read (unit_DBSPH_mesh,*) DBSPH%surf_mesh%vertices(j)%pos(:)
   enddo
! Allocate or resize DBSPH%surf_mesh%faces on the maximum number of n_faces*2 (worst case with all quadrilateral faces)
   if (.not.allocated(DBSPH%surf_mesh%faces)) then
      if (ncord==3) then
         allocate(DBSPH%surf_mesh%faces(2*n_faces),STAT=alloc_stat)
         else
            allocate(DBSPH%surf_mesh%faces(n_faces),STAT=alloc_stat)
      endif
      if (alloc_stat/=0) then
         write(nout,*) 'Allocation of DBSPH%surf_mesh%faces in Import_ply_surface_mesh failed; the program terminates here.'
         stop ! Stop the main program
      endif
      if (.not.allocated(aux_der_type_faces)) then
         if (ncord==3) then
             allocate(aux_der_type_faces(2*n_faces),STAT=alloc_stat)
            else
             allocate(aux_der_type_faces(n_faces),STAT=alloc_stat)
         endif
         if (alloc_stat/=0) then
            write(nout,*) 'Allocation of aux_der_type_faces in Import_ply_surface_mesh failed; the program terminates here.'
            stop ! Stop the main program
         endif
      endif
      old_size_face = 0
      else
         old_size_face = size(DBSPH%surf_mesh%faces)
         if (ncord==3) then
            new_size_face = old_size_face + 2*n_faces
            else
               new_size_face = old_size_face + n_faces
         endif
         aux_der_type_faces(:) = DBSPH%surf_mesh%faces(:)
         deallocate(DBSPH%surf_mesh%faces,STAT=dealloc_stat)
         if (dealloc_stat/=0) then
            write(*,*) 'Deallocation of DBSPH%surf_mesh%faces in Import_ply_surface_mesh failed; the program terminates here'
            stop ! Stop the main program
         endif          
         allocate(DBSPH%surf_mesh%faces(new_size_face),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(*,*) 'Allocation of DBSPH%surf_mesh%faces in Import_ply_surface_mesh failed; the program terminates here'
            stop ! Stop the main program
         endif         
         DBSPH%surf_mesh%faces(1:old_size_face) = aux_der_type_faces(:)
         if (allocated(aux_der_type_faces)) then
            deallocate(aux_der_type_faces,STAT=dealloc_stat)
            if (dealloc_stat/=0) then
               write(nout,*) 'Deallocation of aux_der_type_faces in Import_ply_surface_mesh failed; the program terminates here.'
               stop ! Stop the main program
            endif   
         endif         
         allocate(aux_der_type_faces(new_size_face),STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(nout,*) 'Allocation of aux_der_type_faces in Import_ply_surface_mesh failed; the program terminates here.'
            stop ! Stop the main program
         endif
   endif 
! Read the face vertices: start
  k=old_size_face+1
   do j=1,n_faces
      read(unit_DBSPH_mesh,*) face_vert_num,aux_face_vert(:)
! Assignation of vertices with eventual conversion of any quadrilateral into 2 triangles; computation of area and normal 
      if (ncord==3) then
         DBSPH%surf_mesh%faces(k)%vert_list(1:3) = old_size_vert + aux_face_vert(1:3) + 1
         DBSPH%surf_mesh%faces(k+1)%vert_list(1) = old_size_vert + aux_face_vert(1) + 1
         DBSPH%surf_mesh%faces(k+1)%vert_list(2:3)= old_size_vert + aux_face_vert(3:4) + 1
         DBSPH%surf_mesh%faces(k)%vert_list(4) = 0
         k = k+2
! Computation of area and normal 
         call area_triangle(DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-2)%vert_list(1))%pos,DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-2)%vert_list(2))%pos, &
                            DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-2)%vert_list(3))%pos,DBSPH%surf_mesh%faces(k-2)%area,DBSPH%surf_mesh%faces(k-2)%normal)
         call area_triangle(DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(1))%pos,DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(2))%pos, &
                            DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(3))%pos,DBSPH%surf_mesh%faces(k-1)%area,DBSPH%surf_mesh%faces(k-1)%normal)
         else
            DBSPH%surf_mesh%faces(k)%vert_list(1:4) = old_size_vert + aux_face_vert(1:4) + 1
            k = k+1
! Computation of normal (area will be re-written)             
            call area_triangle(DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(1))%pos,DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(2))%pos, &
                               DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(3))%pos,DBSPH%surf_mesh%faces(k-1)%area,DBSPH%surf_mesh%faces(k-1)%normal)
! Computation of area in 2D (segment length)
            DBSPH%surf_mesh%faces(k-1)%area = dsqrt(dot_product(DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(1))%pos &
                                                               -DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(3))%pos, &
                                                                DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(1))%pos &
                                                               -DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(k-1)%vert_list(3))%pos))/dsqrt(2.d0)
      endif
   enddo  
   close(unit_DBSPH_mesh,IOSTAT=file_stat)
   if (file_stat/=0) then
      write(*,*) 'Error in closing a surface mesh file in Import_pl_surface_meshes; the program terminates here'
      stop
   endif
! Read the face vertices: end   
! Resize DBSPH%surf_mesh%faces on the actual number of faces 
   if (ncord==3) then
      new_size_face = size(DBSPH%surf_mesh%faces) - (k-old_size_face-1-2*n_faces)   
      else
         new_size_face = size(DBSPH%surf_mesh%faces) - (k-old_size_face-1-n_faces)   
   endif
   old_size_face = size(DBSPH%surf_mesh%faces)
   if (new_size_face>old_size_face) then
      aux_der_type_faces(:) = DBSPH%surf_mesh%faces(:)
      deallocate(DBSPH%surf_mesh%faces,STAT=dealloc_stat)
      if (dealloc_stat/=0) then
         write(*,*) 'Deallocation of DBSPH%surf_mesh%faces in Import_ply_surface_mesh failed; the program terminates here.'
         stop ! Stop the main program
      endif          
      allocate(DBSPH%surf_mesh%faces(new_size_face),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(*,*) 'Allocation of DBSPH%surf_mesh%faces in Import_ply_surface_mesh failed; the program terminates here.'
         stop ! Stop the main program
      endif         
      DBSPH%surf_mesh%faces(:) = aux_der_type_faces(1:old_size_face)
      if (allocated(aux_der_type_faces)) then
         deallocate(aux_der_type_faces,STAT=dealloc_stat)
         if (dealloc_stat/=0) then
            write(nout,*) 'Deallocation of aux_der_type_faces in Import_ply_surface_mesh failed; the program terminates here.'
            stop ! Stop the main program
         endif   
      endif         
      allocate(aux_der_type_faces(new_size_face),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(nout,*) 'Allocation of aux_der_type_faces in Import_ply_surface_mesh failed; the program terminates here.'
         stop ! Stop the main program
      endif
   endif
! Read the face vertices: end
enddo   
close(unit_file_list,IOSTAT=file_stat)
if (file_stat/=0) then
   write(*,*) 'Error in closing surface_mesh_list.inp in Import_pl_surface_meshes; the program terminates here.'
   stop
endif
! Initializing the number of surface elements
DBSPH%n_w = new_size_face 
! Deallocations
if (allocated(aux_der_type_vert)) then
   deallocate(aux_der_type_vert,STAT=dealloc_stat)
   if (dealloc_stat/=0) then
      write(nout,*) 'Deallocation of aux_der_type_vert in Import_ply_surface_mesh failed; the program terminates here.'
      stop ! Stop the main program
   endif   
endif 
if (allocated(aux_der_type_faces)) then
   deallocate(aux_der_type_faces,STAT=dealloc_stat)
   if (dealloc_stat/=0) then
      write(nout,*) 'Deallocation of aux_der_type_faces in Import_ply_surface_mesh failed; the program terminates here.'
      stop ! Stop the main program
   endif   
endif

return
end subroutine Import_ply_surface_meshes
! end of subroutine



!AA601 the whole subroutine
!cfile BC_wall_elements.f90
!************************************************************************************
!                          S P H E R A 5.0.3 
!
!                Smoothed Particle Hydrodinamics Code
!
!************************************************************************************
!
! Subroutine : DBSPH_IC_surface_elements
!
! Creation   : Amicarelli A., 16Jan15; Main features (see Purpose)
!
!************************************************************************************
! Purpose            : Initialization of wall surface elements for DBSPH
! Calling procedures : gest_input
! Called procedures  : /
!************************************************************************************

subroutine DBSPH_IC_surface_elements

! Modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE

! Declarations
implicit none
integer(4) :: alloc_stat,i,j
integer(4),external :: ParticleCellNumber

! Interface blocks

! Allocations

! Initializations

! Statements
if (.not.allocated(pg_w)) then
   allocate(pg_w(DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet),STAT=alloc_stat) 
   if (alloc_stat/=0) then
      write(nout,*) 'Allocation of pg_w in DBSPH_IC_surface_elements failed; the program terminates here.'
      call diagnostic (arg1=5,arg2=340)
      stop ! Stop the main program
      else
         pg_w(:)%cella = 0
         pg_w(:)%adjacent_faces(1) = 0
         pg_w(:)%adjacent_faces(2) = 0
         pg_w(:)%adjacent_faces(3) = 0         
         pg_w(:)%coord(1) = 0.d0
         pg_w(:)%coord(2) = 0.d0
         pg_w(:)%coord(3) = 0.d0
         pg_w(:)%vel(1) = 0.d0
         pg_w(:)%vel(2) = 0.d0         
         pg_w(:)%vel(3) = 0.d0
         pg_w(:)%dens = 0.d0
         pg_w(:)%pres = 0.d0
         pg_w(:)%normal(1) = 0.d0
         pg_w(:)%normal(2) = 0.d0
         pg_w(:)%normal(3) = 0.d0
         pg_w(:)%weight = 0.d0 
         pg_w(:)%wet = 0 
         pg_w(:)%mass = 0.d0 
         pg_w(:)%k_d = 0.d0
         pg_w(:)%volume = 0.d0
         write (nout,*) "Allocation of pg_w in DBSPH_IC_surface_elements successfully completed."
   endif   
endif 
!$omp parallel do default(none) shared(DBSPH,pg_w,Med,pg,ncord) private(i)
do i=1,DBSPH%n_w 
!AA601!!! sub start
   if (ncord==3) then
      pg_w(i)%coord(:) = ( DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(1))%pos(:) + &
                           DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(2))%pos(:) + &
                           DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(3))%pos(:) ) / 3.d0
      else
         pg_w(i)%coord(:) = ( DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(1))%pos(:) + &
                              DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(2))%pos(:) + &
                              DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(3))%pos(:) + &
                              DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(4))%pos(:) ) / 4.d0    
   endif
   pg_w(i)%cella = ParticleCellNumber(pg_w(i)%coord)
   pg_w(i)%vel(:) = 0.d0
   pg_w(i)%dens = Med(1)%den0                         
   pg_w(i)%pres = 0.d0                          
   pg_w(i)%normal(:) = DBSPH%surf_mesh%faces(i)%normal(:)                      
   pg_w(i)%weight = DBSPH%surf_mesh%faces(i)%area
   pg_w(i)%wet = 0       
   pg_w(i)%adjacent_faces(:) = 0
end do
!$omp end parallel do
!Initializing fictitious surface elements representing DBSPH inlet sections
!$omp parallel do default(none) shared(DBSPH,pg_w,Med,pg) private(i,j)
do i=(DBSPH%n_w+1),(DBSPH%n_w+DBSPH%n_inlet)
   j= i-DBSPH%n_w  
   pg_w(i)%coord(:) = DBSPH%inlet_sections(j,1:3)
   pg_w(i)%cella = ParticleCellNumber(pg_w(i)%coord)
   pg_w(i)%vel(:) = DBSPH%inlet_sections(j,7:9)
   pg_w(i)%dens = Med(1)%den0                         
   pg_w(i)%pres = 0.d0                          
   pg_w(i)%normal(:) = DBSPH%inlet_sections(j,4:6)                      
   pg_w(i)%weight = 0.d0
   pg_w(i)%wet = 0       
   pg_w(i)%adjacent_faces(:) = 0
end do
!$omp end parallel do
!Initializing fictitious surface elements representing DBSPH outlet sections
!$omp parallel do default(none) shared(DBSPH,pg_w,Med,pg) private(i,j)
do i=(DBSPH%n_w+DBSPH%n_inlet+1),(DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet)
   j= i-(DBSPH%n_w+DBSPH%n_inlet)  
   pg_w(i)%coord(:) = DBSPH%outlet_sections(j,1:3)
   pg_w(i)%cella = ParticleCellNumber(pg_w(i)%coord)
   pg_w(i)%vel(:) = 0.d0
   pg_w(i)%dens = Med(1)%den0                         
   pg_w(i)%pres = 0.d0                          
   pg_w(i)%normal(:) = DBSPH%outlet_sections(j,4:6)                      
   pg_w(i)%weight = 0.d0
   pg_w(i)%wet = 0       
   pg_w(i)%adjacent_faces(:) = 0
end do
!$omp end parallel do
! Deallocations

return
end subroutine DBSPH_IC_surface_elements
!end of the subroutine



!AA601 the whole subroutine
!cfile BC_wall_elements.f90
!************************************************************************************
!                          S P H E R A 5.0.3 
!
!                Smoothed Particle Hydrodinamics Code
!
!************************************************************************************
!
! Subroutine : DBSPH_find_close_faces
!
! Creation   : Amicarelli A., 26Jan15; Main features (see Purpose)
!
!************************************************************************************
! Purpose            : Find the adjacent surface elements of a generic 
!                      computational surface element, both in 3D -triangles- and 2D 
!                      -quadrilaterals- (DBSPH)
! Calling procedures : gest_input
! Called procedures  : /
!************************************************************************************

subroutine DBSPH_find_close_faces

! Modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE

! Declarations
implicit none
integer(4) :: npi,irestocell,i_cell_comp,j_cell_comp,k_cell_comp,j_cell_min,j_cell_max,ID_cell,i_cell,j_cell,k_cell,npj,i_vert,j_vert,aux_adjacent_faces,ww,n_vert
integer(4) :: aux_common_vertices
integer(4),external :: CellIndices,CellNumber

! Interface blocks

! Allocations

! Initializations

! Statements
!Loop over the DBSPH surface elements
!$omp parallel do default(none) &
!$omp shared(DBSPH,pg_w,Icont_w,ncord,NPartOrd_w) &
!$omp private(npi,irestocell,i_cell_comp,j_cell_comp,k_cell_comp,j_cell_min,j_cell_max,ID_cell,i_cell,j_cell,k_cell,npj,i_vert,j_vert,aux_adjacent_faces,ww,aux_common_vertices,n_vert)
loop_surface_elements: do npi = 1,DBSPH%n_w
   aux_adjacent_faces = 0
!Loop over the adjacent cells (background grid) 
   if (pg_w(npi)%cella == 0) cycle
   irestocell = CellIndices (pg_w(npi)%cella,i_cell_comp,j_cell_comp,k_cell_comp)
! Adjacent cell indices
   j_cell_min = j_cell_comp - (ncord-2)
   j_cell_max = j_cell_comp + (ncord-2)
   do j_cell = j_cell_min,j_cell_max    
      do i_cell = i_cell_comp-1,i_cell_comp+1    
         do k_cell = k_cell_comp-1,k_cell_comp+1  
            ID_cell = CellNumber(i_cell,j_cell,k_cell)
            if (ID_cell==0) cycle    
            if (Icont_w(ID_cell+1) <= Icont_w(ID_cell)) cycle
! Loop over the neighbouring wall particles in the cell
            loop_neighbour_surface_elements: do ww = Icont_w(ID_cell),Icont_w(ID_cell+1)-1
               npj = NPartOrd_w(ww)
               if (npi==npj) cycle
! Avoid considering inlet and outlet DBSPH elements    
               if (npj>DBSPH%n_w) cycle
               aux_common_vertices = 0
               if (ncord==3) then
                  n_vert = 3
                  else
                     n_vert = 4
               endif      
! Loop over the vertices of the computational element  
               do i_vert=1,n_vert
! Loop over the vertices of the neighbouring element
                  do j_vert=1,n_vert
                     if ( (DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(npi)%vert_list(i_vert))%pos(1) ==      &
                           DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(npj)%vert_list(j_vert))%pos(1) ) .and. & 
                          (DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(npi)%vert_list(i_vert))%pos(2) == &
                           DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(npj)%vert_list(j_vert))%pos(2) ) .and. & 
                          (DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(npi)%vert_list(i_vert))%pos(3) == &
                           DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(npj)%vert_list(j_vert))%pos(3) )  ) then
                        aux_common_vertices = aux_common_vertices + 1
                        if (aux_common_vertices==2) then
                           pg_w(npi)%adjacent_faces(aux_adjacent_faces+1) = npj 
                           aux_adjacent_faces = aux_adjacent_faces + 1
                           if (aux_adjacent_faces==ncord) cycle loop_surface_elements
                           cycle loop_neighbour_surface_elements
                        endif
                     endif
                  enddo
               enddo
            end do loop_neighbour_surface_elements
         end do    
      end do 
   end do 
end do loop_surface_elements
!$omp end parallel do

! Deallocations

return
end subroutine DBSPH_find_close_faces
!end of the subroutine



!AA601 the whole subroutine
!cfile BC_wall_elements.f90
!************************************************************************************
!                          S P H E R A 5.0.3 
!
!                Smoothed Particle Hydrodinamics Code
!
!************************************************************************************
!
! Subroutine : semi_particle_volumes
!
! Creation   : Amicarelli A., 26Jan15; Main features (see Purpose)
!
!************************************************************************************
! Purpose            : Compute the semi-particle shape coefficients and volumes 
! Calling procedures : gest_input
! Called procedures  : adjacent_faces_isolated_points
!************************************************************************************

subroutine semi_particle_volumes

! Modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE

! Declarations
implicit none
logical :: aux_false_hyp 
integer(4) :: j,i,ID_P_0_iso,ID_P_b_iso,dealloc_stat,aux_adjacent_faces,n_adj_faces_max
double precision :: aux_scalar,aux_scalar_2,k_parameter,alfa,alfa_sum
double precision,dimension(3) :: aux_vec
double precision,dimension(4,3) :: aux_face1,aux_face2

! Interface blocks
interface
   subroutine adjacent_faces_isolated_points(face1,face2,ID_face1_iso,ID_face2_iso,false_hyp)
      implicit none
      double precision, dimension(4,3), intent(in) :: face1,face2
      integer(4),intent(out) :: ID_face1_iso,ID_face2_iso
      logical,intent(out) :: false_hyp 
   end subroutine adjacent_faces_isolated_points
end interface

! Allocations

! Initializations
if (ncord==3) then
   n_adj_faces_max=3
   else
      n_adj_faces_max=2
endif

! Statements
! Loop over the DBSPH surface elements
!$omp parallel do default(none) shared(DBSPH,pg_w,ncord,Domain,Med,pg,nout,n_adj_faces_max) &
!$omp private(aux_scalar,i,j,k_parameter,alfa,alfa_sum,aux_adjacent_faces,aux_face1,aux_face2,aux_false_hyp,aux_scalar_2,aux_vec,ID_P_b_iso,ID_P_0_iso)
do i=1,DBSPH%n_w
   alfa_sum = 0.d0
   aux_adjacent_faces = 0
! Loop over the adjacent faces 
   do j=1,n_adj_faces_max
!   do j=1,3
      if (pg_w(i)%adjacent_faces(j)==0) cycle
! Provided 2 adjacent triangular faces, find the 2 vertices not in common, one per face
      aux_face1(1,:) = DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(1))%pos(:)
      aux_face1(2,:) = DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(2))%pos(:)
      aux_face1(3,:) = DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(3))%pos(:)
      aux_face2(1,:) = DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(pg_w(i)%adjacent_faces(j))%vert_list(1))%pos(:)
      aux_face2(2,:) = DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(pg_w(i)%adjacent_faces(j))%vert_list(2))%pos(:)
      aux_face2(3,:) = DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(pg_w(i)%adjacent_faces(j))%vert_list(3))%pos(:)
      if (ncord==2) then
         aux_face1(4,:) = DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(4))%pos(:)
         aux_face2(4,:) = DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(pg_w(i)%adjacent_faces(j))%vert_list(4))%pos(:)
      endif
      aux_adjacent_faces = aux_adjacent_faces + 1          
      call adjacent_faces_isolated_points(aux_face1,aux_face2,ID_P_0_iso,ID_P_b_iso,aux_false_hyp)
      if (aux_false_hyp==.true.) then
         write(nout,*) "Error! Two faces are not adjacent, but they should be (subroutine semi_particle_volume), Face1: ",i," Face2: ",pg_w(i)%adjacent_faces(j), &
                       "; the run terminates here."
         stop
      endif 
! Compute the spread angle between the normal vectors
      aux_scalar = dot_product(pg_w(i)%normal,pg_w(pg_w(i)%adjacent_faces(j))%normal)
      aux_vec(:) = aux_face2(ID_P_b_iso,:) - aux_face1(ID_P_0_iso,:) 
      aux_scalar_2 = dot_product(pg_w(i)%normal,aux_vec)
      if (aux_scalar_2>=0.000001d0) then
          alfa = PIGRECO + dacos(aux_scalar) 
          else if (aux_scalar_2<=-0.000001d0) then 
             alfa = PIGRECO - dacos(aux_scalar)
             else
                alfa = PIGRECO 
      endif       
! Update alfa_summation (algebric sum)
      alfa_sum = alfa_sum + alfa
   end do
! Compute k_d (shape coefficient) and semi-particle volume (area in 2D) and mass
   if (alfa_sum>(aux_adjacent_faces*PIGRECO*0.5d0)) then
      pg_w(i)%k_d = alfa_sum/(PIGRECO*0.5d0*aux_adjacent_faces)-1.d0
      else
         pg_w(i)%k_d = 0.d0
   endif
   pg_w(i)%volume = pg_w(i)%k_d * DBSPH%k_w * pg_w(i)%weight * Domain%dd / DBSPH%dx_dxw   
   pg_w(i)%mass = pg_w(i)%volume * Med(1)%den0  
end do
!$omp end parallel do

! Deallocations
if (allocated(DBSPH%surf_mesh%vertices)) then
   deallocate(DBSPH%surf_mesh%vertices,STAT=dealloc_stat) 
   if (dealloc_stat/=0) then
      write(nout,*) 'Deallocation of DBSPH%surf_mesh%vertices in DBSPH_find_close_faces failed; the program terminates here'
      stop ! Stop the main program
   endif   
endif 
if (allocated(DBSPH%surf_mesh%faces)) then
   deallocate(DBSPH%surf_mesh%faces,STAT=dealloc_stat) 
   if (dealloc_stat/=0) then
      write(nout,*) 'Deallocation of DBSPH%surf_mesh%faces in DBSPH_find_close_faces failed; the program terminates here'
      stop ! Stop the main program
   endif   
endif

return
end subroutine semi_particle_volumes
!end of the subroutine



!AA601 the whole subroutine
!cfile BC_wall_elements.f90
!************************************************************************************
!                          S P H E R A 5.0.3 
!
!                Smoothed Particle Hydrodinamics Code
!
!************************************************************************************
!
! File name     : DBSPH_kinematics
!
! Creation      : Amicarelli A., 26Jan16
!
!************************************************************************************
! Module purpose : Imposed kinematics for the DBSPH frontier system (linear interpolation of input data)
!                  (method derived from time_integration_body_dynamics)
!
! Calling routine: Euler
!
! Called routines: /
!
!************************************************************************************

subroutine DBSPH_kinematics

! Modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE

! Declarations
implicit none
integer(4) :: j

! Interface blocks

! Allocations

! Initializations

! Statements
! Loop over kinematics records to provide linear interpolation of imposed velocity to the first surface element
if (DBSPH%n_w>0) then
do j=1,DBSPH%n_kinematics_records
   if (DBSPH%kinematics(j,1)>=tempo) then
      if (DBSPH%kinematics(j,1)==tempo) then
         pg_w(1)%vel(:) = DBSPH%kinematics(j,2:4)
         else
            pg_w(1)%vel(:) = DBSPH%kinematics(j-1,2:4) + (DBSPH%kinematics(j,2:4)-DBSPH%kinematics(j-1,2:4))/ &
                             (DBSPH%kinematics(j,1)-DBSPH%kinematics(j-1,1)) * (tempo-DBSPH%kinematics(j-1,1))
      endif
      exit
   endif
enddo
endif

! Deallocations

return
end subroutine DBSPH_kinematics
!---split



!AA601 the whole subroutine
!cfile BC_wall_elements.f90
!************************************************************************************
!                          S P H E R A 5.0.3 
!
!                Smoothed Particle Hydrodinamics Code
!
!************************************************************************************
!
! File name     : adjacent_faces_isolated_points
! Creation      : Amicarelli A., 26Jan15
!
!************************************************************************************
! Module purpose : Provided 2 adjacent triangular/quadrilateral faces, find at least 2 vertices not in common, at least one per face.
!                  They are ID_face1_iso and ID_face2_iso. In case the faces are not adjacent, then false_hyp=.true.
!
! Calling routine: semi_particle_volumes
!
! Called routines: /
!
!************************************************************************************

subroutine adjacent_faces_isolated_points(face1,face2,ID_face1_iso,ID_face2_iso,false_hyp)

! Modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE

! Declarations
implicit none
double precision, dimension(4,3), intent(in) :: face1,face2
integer(4),intent(out) :: ID_face1_iso,ID_face2_iso
logical,intent(out) :: false_hyp 
integer(4) :: test_face1,test_face2,i,j,n_vert

! Interface blocks

! Allocations

! Initializations
test_face1 = 0
test_face2 = 0
false_hyp = .false.

! Statements
if (ncord==3) then
   n_vert = 3
   else
      n_vert = 4
endif
! do over the 3/4 vertices of the first face
do_vertices_face1: do i=1,n_vert
! do over the 3/4 vertices of the second face   
    do j=1,n_vert
       if ( (face1(i,1)==face2(j,1)) .and. (face1(i,2)==face2(j,2)) .and. (face1(i,3)==face2(j,3)) ) then
! In case the vertex is in common, the vertex of the first face is not anymore isolated and update test_face1/2 (the sum 
! of the 2 non-isolated vertex IDs of face 1/2; in 2D the ID squares are considered to avoid ambiguities in the following computations)
       if (ncord==3) then
          test_face1 = test_face1 + i 
          test_face2 = test_face2 + j 
          else
             test_face1 = test_face1 + i**2 
             test_face2 = test_face2 + j**2    
       endif
          cycle do_vertices_face1
       endif
   enddo
! end do over the vertices of the second face      
end do do_vertices_face1
! end do over the vertices of the first face
! The only/first vertex (in 3D/2D) of the 1st face, not contributing to test_face1, is finally found
if (ncord==3) then
   select case(test_face1)
      case(5)
      ID_face1_iso = 1
      case(4)
      ID_face1_iso = 2
      case(3)
      ID_face1_iso = 3
      case default
      ID_face1_iso = 0
   end select
   else
      select case(test_face1)
      case(13,20,25)
      ID_face1_iso = 1
      case(10,17)
      ID_face1_iso = 2
      case(5)
      ID_face1_iso = 3
      case default
      ID_face1_iso = 0
      end select   
endif
! The only/first vertex (in 3D/2D) of the 2nd face, not contributing to test_face2, is finally found
if (ncord==3) then
   select case(test_face2)
      case(5)
      ID_face2_iso = 1
      case(4)
      ID_face2_iso = 2
      case(3)
      ID_face2_iso = 3
      case default
      ID_face2_iso = 0
   end select  
   else
      select case(test_face2)
         case(13,20,25)
         ID_face2_iso = 1
         case(10,17)
         ID_face2_iso = 2
         case(5)
         ID_face2_iso = 3
         case default
         ID_face2_iso = 0
      end select 
endif
if ((ID_face1_iso==0).or.(ID_face2_iso==0)) false_hyp = .true.

! Deallocations

return
end subroutine adjacent_faces_isolated_points
!---split