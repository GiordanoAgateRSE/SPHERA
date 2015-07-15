!cfile inter_EqCont_2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : inter_EqCont_2D
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Amicarelli/Agate  30Nov11        BSPH: wall element parameters and subroutines
!AA501b
! 04  Amicarelli-Agate  13nov12        Body dynamics
!AA504
! 06Amicarelli          08Apr14        (v5.04) Modifications for granular flows; second invariant of the strain rate tensor is computed only when necessary.
!
!************************************************************************************
! Module purpose : Module to accumulate the contributions of density variation in the
!                  continuity equation of the particles that are in the sphere of
!                  influence of the particle considered for the 2D case. The velocity 
!                  spatial derivatives and the second invariant of the deformation   
!                  velocity tensor are calculated.
! 
! Calling routine: Loop_Irre_2D
!
! Called routines: diffumorris
!
!************************************************************************************
!
subroutine inter_EqCont_2D 
! ex inter32d
!* implementa il meccanismo di ricerca delle particelle che agiscono
!* su quella i-esima.
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
!integer(4), parameter :: local_d = 500  ! num max part entro 2h
!
!.. Formal Arguments ..
!
!.. Local Scalars ..
integer(4) :: npi,npj,contj,npartint    !!!,npar2h
double precision :: rhoi,rhoj,amassj,moddervel,det,appo       !!! unity,pesoj,  
!modello bifluido
double precision :: factdiff, rvw, dervol   !, tdiff, cuei
!modello bifluido
!!!logical :: fix
!
!.. Local Arrays ..
double precision,dimension(3) :: pesogradj
! quantita' step eq continuita'
double precision,dimension(3) :: dvar      ! dervel, dvel per eq standard, dvar usa vel smoot di monaghan
double precision,dimension(4) :: derspa, aij
!
!.. Executable Statements ..
!
!!! !$omp private(i,unity,npar2h,derspa,dvar,aij,contj,npartint,j,fix,rhoi,rhoj,amassj,moddervel,pesoj,pesogradj) &
!
!AA504 sub omp directives
!AA406 sub
!$omp parallel do default(none) &
!$omp private(npi,derspa,dvar,aij,contj,npartint,npj,rhoi,rhoj,amassj,moddervel,pesogradj) &
!$omp private(dervol,factdiff,rvw,det,appo) &
!$omp shared(nag,pg,nPartIntorno,NMAXPARTJ,PartIntorno,PartKernel,rag,diffusione) &
!$omp shared(pg_w,nPartIntorno_fw,PartIntorno_fw,DBSPH,kernel_fw,rag_fw,Domain,Granular_flows_options)
!
!.. loops on all the active particles
!
 do npi = 1,nag
!
   if ( pg(npi)%vel_type /= "std" .or. pg(npi)%cella == 0) cycle
!
!* azzeramento quantita generali
!!!   npar2h    = 0
!!!   unity     = zero
   pg(npi)%dden  = zero
   pg(npi)%diffu = zero
   derspa(:) = zero
   dvar(:)   = zero
   aij(:)    = zero
!*_______________________________________________________________
!*prima passata per trovare celle interagenti e memorizzazione
!
   do contj = 1, nPartIntorno(npi)
!
     npartint = (npi-1)* NMAXPARTJ + contj
     npj = PartIntorno(npartint)
!
!!!     fix = .FALSE.
!
     rhoi    = pg(npi)%dens
     rhoj    = pg(npj)%dens
     amassj  = pg(npj)%mass
     dvar(:) = pg(npj)%var(:) - pg(npi)%var(:)
!
     if ( pg(npj)%vel_type /= "std" ) then          !non part fix o altro
       rhoj    = rhoi
       amassj  = pg(npi)%mass
       dvar(:) = pg(npj)%vel(:) - pg(npi)%var(:)        !!! resta definiti cosi' per particelle con moto fix ma non ferme
       if ( pg(npj)%vel(1) == zero .AND. pg(npj)%vel(2) == zero .AND. pg(npj)%vel(3) == zero ) then   !!! fix = .TRUE.
!!!       if ( pg(npj)%slip == "f" .AND. fix ) then 
         if ( pg(npj)%slip == "f" ) then 
           moddervel = -two * (pg(npi)%var(1) * pg(npj)%zer(1) + pg(npi)%var(3) * pg(npj)%zer(3))
           dvar(:) = moddervel * pg(npj)%zer(:)
         end if
       end if
     end if

!=============  EQUAZIONE DI CONTINUITA' ================
!
!* calcolo unita'
!!!     pesoj = amassj * PartKernel(4,npartint) / rhoj
     pesogradj(1:3) = amassj * rag(1:3,npartint) * PartKernel(1,npartint) / rhoj
!!!     unity = unity + pesoj
!!!     if ( pg(npj)%vel_type /= "fix" ) npar2h = npar2h + 1
!
!*eq. di continuita' calcolo drho ++++++++++++++++++++++++++++++++++++++++++++++++
!     pg(npi)%dden = pg(npi)%dden - amassj * PartKernel(1,npartint) * &
!                  (dervel(1)*rag(1,npartint) + dervel(2)*rag(2,npartint) + dervel(3)*rag(3,npartint)) 
!AA504 sub start
     if (Granular_flows_options%ID_erosion_criterion.ne.1) then
        appo = amassj * PartKernel(1,npartint) * &
               (dvar(1)*rag(1,npartint) + dvar(2)*rag(2,npartint) + dvar(3)*rag(3,npartint))
        else
             appo = rhoi * PartKernel(1,npartint) * (amassj/rhoj) * &
                    (dvar(1)*rag(1,npartint) + dvar(2)*rag(2,npartint) + dvar(3)*rag(3,npartint))
     endif
!AA504 sub end

!AA501btest
     pg(npi)%dden = pg(npi)%dden - appo

!*calcolo termine diffusivo
! modello bifluido
     if (diffusione) then
       dervol = pg(npj)%VolFra - pg(npi)%VolFra
       call diffumorris (npi,npj,npartint,dervol,factdiff,rvw)
       pg(npi)%diffu = pg(npi)%diffu + factdiff * rvw
     end if
! modello bifluido
!
!*calcolo derivate spaziali 
     if ( pg(npj)%vel_type /= "std" ) cycle
     
!AA504
     if ((Granular_flows_options%ID_erosion_criterion==1).or.(Granular_flows_options%ID_erosion_criterion==3)) then     
        derspa(1) = derspa(1) + pesogradj(1) * dvar(1)
        derspa(2) = derspa(2) + pesogradj(3) * dvar(1)
        derspa(3) = derspa(3) + pesogradj(1) * dvar(3)
        derspa(4) = derspa(4) + pesogradj(3) * dvar(3)
        aij(1) = aij(1) + (pg(npj)%coord(1) - pg(npi)%coord(1)) * pesogradj(1)
        aij(2) = aij(2) + (pg(npj)%coord(3) - pg(npi)%coord(3)) * pesogradj(1)
        aij(3) = aij(3) + (pg(npj)%coord(1) - pg(npi)%coord(1)) * pesogradj(3)
        aij(4) = aij(4) + (pg(npj)%coord(3) - pg(npi)%coord(3)) * pesogradj(3)
!AA504
     endif      
   end do
!
!AA406 start
! boundary contributions (BSPH)
    if (Domain%tipo == "bsph") then
       do contj = 1,nPartIntorno_fw(npi)
          npartint = (npi-1)* NMAXPARTJ + contj
          npj = PartIntorno_fw(npartint)
          dvar(:) = pg_w(npj)%vel(:) - pg(npi)%var(:)
!AA406test occhio nche al subtest Gamma
          appo = -pg_w(npj)%dens * pg_w(npj)%weight * kernel_fw(1,npartint) * pg(npi)%Gamma * &
                (dvar(1) * pg_w(npj)%normal(1) + dvar(2) * pg_w(npj)%normal(2) + dvar(3) * pg_w(npj)%normal(3))
          appo = appo - pg_w(npj)%mass * kernel_fw(2,npartint) * &
                (dvar(1) * rag_fw(1,npartint) + dvar(2) * rag_fw(2,npartint) + dvar(3) * rag_fw(3,npartint))
          pg(npi)%dden = pg(npi)%dden + appo
       end do
    endif
!AA406 end
! 
!prova!
!   pg(npi)%dden = floor(pg(npi)%dden * azzeramento) / azzeramento
!   if (dabs(pg(npi)%dden) < arrotondamento) pg(npi)%dden = zero 
!prova!
!
!AA504
   if ((Granular_flows_options%ID_erosion_criterion==1).or.(Granular_flows_options%ID_erosion_criterion==3)) then
      det = aij(1) * aij(4) - aij(2) * aij(3)
!
      if (det < 0.001d0) det = 0.001d0
!
      pg(npi)%dudx = ( derspa(1)*aij(4) - derspa(2)*aij(2)) / det       ! derivate norm accoppiate
      pg(npi)%dudy = (-derspa(1)*aij(3) + derspa(2)*aij(1)) / det       ! derivate norm accoppiate
      pg(npi)%dvdx = ( derspa(3)*aij(4) - derspa(4)*aij(2)) / det       ! derivate norm accoppiate
      pg(npi)%dvdy = (-derspa(3)*aij(3) + derspa(4)*aij(1)) / det       ! derivate norm accoppiate
!
!ç   pg(npi)%dudy = half * (pg(npi)%dudy + pg(npi)%dvdx)
!ç   pg(npi)%dvdx = pg(npi)%dudy
!
      pg(npi)%secinv = Dsqrt(half * pg(npi)%dudx*pg(npi)%dudx + half * pg(npi)%dvdy*pg(npi)%dvdy + quarter * &
                       (pg(npi)%dudy + pg(npi)%dvdx) * (pg(npi)%dudy + pg(npi)%dvdx)) !SaMa: forse perche' dudy=dvdx
!ç   pg(npi)%secinv = Dsqrt(half * pg(npi)%dudx*pg(npi)%dudx + half * pg(npi)%dvdy*pg(npi)%dvdy + pg(npi)%dudy*pg(npi)%dudy)
!AA504
   endif
!
!!!   pg(npi)%npar2h = npar2h
!
!!!   pg(npi)%uni = unity
!
 end do  ! nag
!
!$omp end parallel do

!AA501b
 call start_and_stop(3,12)
 call start_and_stop(2,19)
 if (n_bodies > 0) call body_particles_to_continuity
 call start_and_stop(3,19)
 call start_and_stop(2,12)

!
return
end subroutine inter_EqCont_2D
!---split

