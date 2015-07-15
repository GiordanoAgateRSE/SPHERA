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

