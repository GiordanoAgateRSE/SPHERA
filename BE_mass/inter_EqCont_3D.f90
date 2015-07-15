!cfile inter_EqCont_3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : inter_EqCont_3D
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!AA501b
! 03  Amicarelli-Agate  13nov12        Body dynamics
!AA504
! 06Amicarelli          08Apr14        (v5.04) Modifications for granular flows; second invariant of the strain rate tensor is computed only when necessary.
!                                      Renormalization deactivated (because of errors). Correction of the formula for second invariant of the strain rate tensor. 
!
!************************************************************************************
!AA406 sub
! Module purpose : Module to accumulate the contributions of density variation in the
!
!                  continuity equation of the particles that are in the sphere of
!                  influence of the particle considered for the 3D case. The velocity 
!                  spatial derivatives and the second invariant of the deformation   
!                  velocity tensor are calculated.
!
! Calling routine: Loop_Irre_3D
!
! Called routines: diffumorris
!
!************************************************************************************
!
subroutine inter_EqCont_3D
! ex inter33d
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
!.. Local Scalars ..
integer(4) :: npi,npj,contj,npartint    !!!,npar2h
double precision :: rhoi,rhoj,amassj,moddervel,det,moddia,modout,appo    !!! unity,pesoj,   
! quantita' step2 eq del moto
!double precision pi,pj
! modello bifluido
double precision :: factdiff, rvw, dervol   !, tdiff, cuei
! modello bifluido
!!!logical :: fix
!
!.. Local Arrays ..
double precision,dimension (3) :: pesogradj
! quantita' step eq continuita'
double precision,dimension (3) :: dvar     ! dervel, dvel per eq standard, dvar usa vel smoot di monaghan
double precision,dimension (9) :: derspa,aij,bij,dvdi
!
!.. Executable Statements ..
!
!!! !$omp private(i,unity,npar2h,derspa,dvar,aij,bij,dvdi,contj,npartint,j,fix,rhoi,rhoj,amassj,pi,pj,moddervel) &
!!! !$omp private(pesoj,pesogradj,dervol,factdiff,rvw,det,moddia,modout) &.
!pi,pj,
!AA504 sub omp directives
!$omp parallel do default(none) &
!$omp private(npi,derspa,dvar,aij,bij,dvdi,contj,npartint,npj,rhoi,rhoj,amassj,moddervel) &
!$omp private(pesogradj,dervol,factdiff,rvw,det,moddia,modout,appo) &
!$omp shared(nag,pg,nPartIntorno,NMAXPARTJ,PartIntorno,PartKernel,rag,diffusione,Granular_flows_options)
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
   bij(:)    = zero
   dvdi(:)   = zero
!
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
     rhoi   = pg(npi)%dens
     rhoj   = pg(npj)%dens
     amassj = pg(npj)%mass
!     pj     = pg(npj)%pres
!     pi     = pg(npi)%pres
!
     dvar(:) = pg(npj)%var(:) - pg(npi)%var(:)
!
     if ( pg(npj)%vel_type /= "std" ) then          !non part fix o altro
!       pj    = pi
       rhoj    = rhoi
       amassj  = pg(npi)%mass
       dvar(:) = pg(npj)%vel(:) - pg(npi)%var(:)        !!! resta definiti così per particelle con moto fix ma non ferme
       if ( pg(npj)%vel(1) == zero .AND. pg(npj)%vel(2) == zero .AND. pg(npj)%vel(3) == zero ) then  !!! fix = .TRUE.
!!!       if ( pg(npj)%slip == "f" .AND. fix ) then 
         if ( pg(npj)%slip == "f") then 
           moddervel = -two * (pg(npi)%var(1) * pg(npj)%zer(1) + pg(npi)%var(3) * pg(npj)%zer(3))
           dvar(:) = moddervel * pg(npj)%zer(:)
         end if
       end if
     end if
!
!=============  EQUAZIONE DI CONTINUITA' ================
!      
!* calcolo unita'
!!!     pesoj = amassj * PartKernel(4,npartint) / rhoj
     pesogradj(1:3) = amassj * rag(1:3,npartint) * PartKernel(1,npartint) / rhoj
!!!     unity = unity + pesoj  
!!!     if ( pg(npj)%vel_type /= "fix" ) npar2h = npar2h + 1

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
     pg(npi)%dden = pg(npi)%dden - appo
!
!*calcolo termine diffusivo
! modello bifluido
     if (diffusione) then
       dervol = pg(npj)%VolFra - pg(npi)%VolFra
       call diffumorris (npi,npj,npartint,dervol,factdiff,rvw)
       pg(npi)%diffu = pg(npi)%diffu + factdiff * rvw
     end if
! modello bifluido

!*calcolo derivate spaziali 
     if ( pg(npj)%vel_type /= "std" ) cycle 
!AA504
     if ((Granular_flows_options%ID_erosion_criterion==1).or.(Granular_flows_options%ID_erosion_criterion==3)) then
!AA504 rm start     
!     derspa(1) = derspa(1) + pesogradj(1) * dvar(1)      !du/dx  !! LA MATRICE DERSPA
!     derspa(2) = derspa(2) + pesogradj(2) * dvar(1)      !du/dy  !! è ORDINATA 
!     derspa(3) = derspa(3) + pesogradj(3) * dvar(1)      !du/dz  !! PER COLONNE
!     derspa(4) = derspa(4) + pesogradj(1) * dvar(2)      !dv/dx
!     derspa(5) = derspa(5) + pesogradj(2) * dvar(2)      !dv/dy
!     derspa(6) = derspa(6) + pesogradj(3) * dvar(2)      !dv/dz
!     derspa(7) = derspa(7) + pesogradj(1) * dvar(3)      !dw/dx
!     derspa(8) = derspa(8) + pesogradj(2) * dvar(3)      !dw/dy
!     derspa(9) = derspa(9) + pesogradj(3) * dvar(3)      !dw/dz
!!.. normalizzazione del gradiente di velocita'
!     aij(1) = aij(1) + (pg(npj)%coord(1) - pg(npi)%coord(1)) * pesogradj(1)   !a11
!     aij(2) = aij(2) + (pg(npj)%coord(2) - pg(npi)%coord(2)) * pesogradj(1)   !a12
!     aij(3) = aij(3) + (pg(npj)%coord(3) - pg(npi)%coord(3)) * pesogradj(1)   !a13
!     aij(4) = aij(4) + (pg(npj)%coord(1) - pg(npi)%coord(1)) * pesogradj(2)   !a21
!     aij(5) = aij(5) + (pg(npj)%coord(2) - pg(npi)%coord(2)) * pesogradj(2)   !a22
!     aij(6) = aij(6) + (pg(npj)%coord(3) - pg(npi)%coord(3)) * pesogradj(2)   !a23
!     aij(7) = aij(7) + (pg(npj)%coord(1) - pg(npi)%coord(1)) * pesogradj(3)   !a31
!     aij(8) = aij(8) + (pg(npj)%coord(2) - pg(npi)%coord(2)) * pesogradj(3)   !a32
!     aij(9) = aij(9) + (pg(npj)%coord(3) - pg(npi)%coord(3)) * pesogradj(3)   !a33
!AA504 rm end
!AA504 start
     dvdi(1) = dvdi(1) + pesogradj(1) * dvar(1)      !du/dx  
     dvdi(2) = dvdi(2) + pesogradj(2) * dvar(1)      !du/dy  
     dvdi(3) = dvdi(3) + pesogradj(3) * dvar(1)      !du/dz  
     dvdi(4) = dvdi(4) + pesogradj(1) * dvar(2)      !dv/dx
     dvdi(5) = dvdi(5) + pesogradj(2) * dvar(2)      !dv/dy
     dvdi(6) = dvdi(6) + pesogradj(3) * dvar(2)      !dv/dz
     dvdi(7) = dvdi(7) + pesogradj(1) * dvar(3)      !dw/dx
     dvdi(8) = dvdi(8) + pesogradj(2) * dvar(3)      !dw/dy
     dvdi(9) = dvdi(9) + pesogradj(3) * dvar(3)      !dw/dz
!AA504 end
!AA504
     endif 
!============= FINE EQUAZIONE DI CONTINUITA' ================
!
   end do
!
!prova!
!   pg(npi)%dden = floor(pg(npi)%dden * azzeramento) / azzeramento
!   if (dabs(pg(npi)%dden) < arrotondamento) pg(npi)%dden = zero 
!prova!
!AA504
    if ((Granular_flows_options%ID_erosion_criterion==1).or.(Granular_flows_options%ID_erosion_criterion==3)) then
!AA504 rm start (this renormalization presents several errors; it should be re-written)
!   bij(1) = (aij(5) * aij(9) - aij(6) * aij(8))
!   bij(2) =-(aij(4) * aij(9) - aij(7) * aij(6))
!   bij(3) = (aij(4) * aij(8) - aij(5) * aij(7))
!   bij(4) =-(aij(2) * aij(9) - aij(3) * aij(8))
!   bij(5) = (aij(1) * aij(9) - aij(7) * aij(7))
!   bij(6) =-(aij(1) * aij(8) - aij(2) * aij(2))
!   bij(7) = (aij(2) * aij(6) - aij(3) * aij(5))
!   bij(8) =-(aij(1) * aij(6) - aij(3) * aij(4))
!   bij(9) = (aij(1) * aij(5) - aij(2) * aij(4))
!!
!   det = aij(1)*aij(5)*aij(9)+aij(2)*aij(6)*aij(7)+aij(3)*aij(4)*aij(8) - &
!         aij(1)*aij(6)*aij(8)-aij(4)*aij(2)*aij(9)-aij(7)*aij(5)*aij(3)
!!
!   if (det < 0.001d0) det = 0.001d0
!!
!   dvdi(1) = (bij(1)*derspa(1) - bij(2)*derspa(2) + bij(3)*derspa(3)) / det     !dudx
!   dvdi(2) = (bij(1)*derspa(4) - bij(2)*derspa(5) + bij(3)*derspa(6)) / det     !dudy
!   dvdi(3) = (bij(1)*derspa(7) - bij(2)*derspa(8) + bij(3)*derspa(9)) / det     !dudz
!   dvdi(4) = (bij(4)*derspa(1) - bij(5)*derspa(2) + bij(6)*derspa(3)) / det     !dvdx
!   dvdi(5) = (bij(4)*derspa(4) - bij(5)*derspa(5) + bij(6)*derspa(6)) / det     !dvdy
!   dvdi(6) = (bij(4)*derspa(7) - bij(5)*derspa(8) + bij(6)*derspa(9)) / det     !dvdz
!   dvdi(7) = (bij(7)*derspa(1) - bij(8)*derspa(2) + bij(9)*derspa(3)) / det     !dwdx
!   dvdi(8) = (bij(7)*derspa(4) - bij(8)*derspa(5) + bij(9)*derspa(6)) / det     !dwdy
!   dvdi(9) = (bij(7)*derspa(7) - bij(8)*derspa(8) + bij(9)*derspa(9)) / det     !dwdz
!AA504 rm end
!
! modifica solo dei termini utili per tensore vel def
!   dvdi(2) = half * (dvdi(2) + dvdi(4))
!   dvdi(3) = half * (dvdi(3) + dvdi(7))
!   dvdi(6) = half * (dvdi(6) + dvdi(8))
! 
!   moddia  = (dvdi(1)*dvdi(1) + dvdi(5)*dvdi(5) + dvdi(9)*dvdi(9))
!   modout  = (dvdi(2)*dvdi(2) + dvdi(3)*dvdi(3) + dvdi(6)*dvdi(6))
!   pg(npi)%secinv = Dsqrt( half*moddia + modout)
! 
   moddia  = (dvdi(1) * dvdi(1) + dvdi(5) * dvdi(5) + dvdi(9) * dvdi(9))
!AA504 sub 
   modout  = ( (dvdi(2)+dvdi(4))*(dvdi(2)+dvdi(4)) + (dvdi(3)+dvdi(7))*(dvdi(3)+dvdi(7)) + (dvdi(6)+dvdi(8))*(dvdi(6)+dvdi(8)) )
! old formula
! modout  = ( (dvdi(2)*dvdi(4))*(dvdi(2)*dvdi(4)) + (dvdi(3)*dvdi(7))*(dvdi(3)*dvdi(7)) + (dvdi(6)*dvdi(8))*(dvdi(6)*dvdi(8)) )
   
   pg(npi)%secinv = Dsqrt( half * moddia + quarter * modout)
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
!

!AA501b
 call start_and_stop(3,12)
 call start_and_stop(2,19)
 if (n_bodies > 0) call body_particles_to_continuity
 call start_and_stop(3,19)
 call start_and_stop(2,12)

return
end subroutine inter_EqCont_3D
!---split

