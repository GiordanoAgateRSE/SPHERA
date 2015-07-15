!cfile inter_SmoothVelo_2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : inter_SmoothVelo_2D
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
! 03  Amicarelli-Agate  13nov12        Body particle contributions 

!
!************************************************************************************
! Module purpose : Module to calculate the corrective term of the velocity
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D, Euler, Heun
!

!AA501b modified
!AA601 sub
! Called routines: body_to_smoothing_vel,DBSPH_inlet_outlet
 
!
!************************************************************************************
!
subroutine inter_SmoothVelo_2D 
! ex inter42d
!* implementa il meccanismo di ricerca delle particelle che agiscono su quella npi-esima.
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
  integer(4)       :: npi,i,ii,npj,contj,npartint
  integer(4)       :: icbs, Ncbs, IntNcbs, ibdt, ibdp, iside, sidestr   !Ncols, 
!AA406 sub
  double precision :: rhoi,rhoj,amassj,pesoj,  moddervel,unity
!
  double precision :: IntWdV  !xpmin,xpmax,interlen,xpi,ypi, deltai,sidelen, 
  character(4)     :: strtype
!
!.. Local Arrays ..
  integer(4),      dimension(1:PLANEDIM) :: acix
  double precision,dimension(1:PLANEDIM) :: sss, nnn, DVLoc, DVGlo, BCLoc, BCGlo
  double precision,dimension(1:PLANEDIM) :: IntLocXY
! quantita' step eq continuita'
  double precision,dimension(3) :: dervel     ! dvel per eq standard, dvar usa vel smoot di monaghan
  
!AA501b start
  double precision,dimension(:,:),allocatable :: dervel_mat
  double precision,dimension(:),allocatable :: unity_vec

  if (n_bodies > 0) then  
!Allocations
     allocate(dervel_mat(nag,3))
     allocate(unity_vec(nag))
!Initializations
     dervel_mat = 0.
     unity_vec = 0.
  endif
!AA501b end  
  
!.. Executable Statements ..
!
  acix(1) = 1
  acix(2) = 3
!

!AA406
  unity = zero
  
!AA501b
! Body particle contributions to pressure smoothing
  if (n_bodies > 0) then
     call start_and_stop(3,14)
     call start_and_stop(2,19)
     call body_to_smoothing_vel(dervel_mat,unity_vec)
     call start_and_stop(3,19)
     call start_and_stop(2,14)
  endif
!AA501b end   
  
!
!!! !$omp private(ii,npi,unity,contj,npartint,npj,rhoi,rhoj,amassj,dervel,moddervel,pesoj) &   !sidelen,
!AA501b modified
!$omp parallel do default(none) &
!$omp private(ii,npi,contj,npartint,npj,rhoi,rhoj,amassj,dervel,moddervel,pesoj) &
!$omp private(Ncbs,IntNcbs,ibdt,ibdp,icbs,IntLocXY,iside,sidestr,strtype,i) &
!$omp private(sss,nnn,DVGlo,DVLoc,IntWdV,BCLoc,BCGlo) &
!$omp shared(nag,Pg,Domain,Med,Tratto,acix,nPartIntorno,NMAXPARTJ,PartIntorno,PartKernel) &
!$omp shared(BoundaryDataPointer,BoundaryDataTab,BoundarySide,indarrayFlu,Array_Flu,esplosione,kernel_fw,unity,dervel_mat,unity_vec,n_bodies)
!
!.. loops on all the active particles in order to calculate the average velocity depending on 
!.. the smoothing coefficient of velocity ThetaV
!
!!!!!!  do npi = 1,nag
!!!!!!!
!!!!!!    if ( pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" .or. pg(npi)%state == "sol") cycle
!$$$$$
  do ii = 1,indarrayFlu
    npi = Array_Flu(ii)
!$$$$$
!
!* azzeramento quantita generali
    pg(npi)%var = zero
    pg(npi)%Envar = zero
!!!    unity   = zero
!
!.. contributo delle particelle interne
!
    do contj = 1, nPartIntorno(npi)
!
      npartint = (npi-1)* NMAXPARTJ + contj
      npj = PartIntorno(npartint)
!
      rhoi   = pg(npi)%dens
      rhoj   = pg(npj)%dens
      amassj = pg(npj)%mass
!
      dervel(:) = pg(npj)%vel(:) - pg(npi)%vel(:)
!
      if ( pg(npj)%vel_type /= "std" ) then          !non part fix o altro
        rhoj   = rhoi
        amassj = pg(npi)%mass
        moddervel = -two * (pg(npi)%vel(1)*pg(npj)%zer(1)+pg(npi)%vel(2)*pg(npj)%zer(2)+pg(npi)%vel(3)*pg(npj)%zer(3))
        dervel(:) = moddervel * pg(npj)%zer(:)    !+pg(npj)%vstart(:)
      end if
!!============= CORREZIONE VELOCITA' ===================
      if ( Med(pg(npj)%imed)%den0 /= Med(pg(npi)%imed)%den0 ) cycle
!
      pesoj = amassj * PartKernel(4,npartint) / rhoj
!
!AA406
      unity = unity + pesoj
!
!!!      unity = unity + pesoj  
      pg(npi)%var(:) = pg(npi)%var(:) + dervel(:) * pesoj   
! calcolo velocita' per modello diffusivo
!   pg(npj)%veldif(:) = pg(npj)%veldif(:) + pg(npj)%vel(:) * pesoj
!
!...................................... 2011 mar 08
!.. update of Specific Internal Energy
      if (esplosione) pg(npi)%Envar = pg(npi)%Envar + (pg(npj)%IntEn - pg(npi)%IntEn) * pesoj
!.. end update of Specific Internal Energy
!.........................................
!
    end do
!
!

!AA501b start
      if (n_bodies > 0) then
         pg(npi)%var(:) = pg(npi)%var(:) + dervel_mat(npi,:)
         unity = unity + unity_vec(npi)
      endif

!AA406 start
      if (Domain%tipo == "bsph") then
!AA406test
          pg(npi)%var(:) = pg(npi)%var(:)
!          pg(npi)%var(:) = pg(npi)%var(:) / unity
!AA601 sub start
! Impose boundary conditions at inlet and outlet sections (DB-SPH)
          call DBSPH_inlet_outlet(npi)
!AA601 sub end
     else
!AA406 end
!
!.. contributo dei contorni
!
    Ncbs = BoundaryDataPointer(1,npi)
    IntNcbs = BoundaryDataPointer(2,npi)
    ibdt = BoundaryDataPointer(3,npi)
!
    if (IntNcbs > 0) then
!
      do icbs = 1, IntNcbs
!
        ibdp = ibdt + icbs - 1
        IntLocXY(1:PLANEDIM) = BoundaryDataTab(ibdp)%LocXYZ(1:PLANEDIM)
        iside = BoundaryDataTab(ibdp)%CloBoNum
        sidestr = BoundarySide(iside)%stretch
        strtype = Tratto(sidestr)%tipo
!
        if (strtype == 'sour' .or. strtype == 'velo' .or. strtype == 'flow') then
          pg(npi)%var(:) = zero   
!!!          unity = one
          exit  
        end if
!
!        sidelen = BoundarySide(iside)%Length
        do i = 1, PLANEDIM
          sss(i) = BoundarySide(iside)%T(acix(i),1)
          nnn(i) = BoundarySide(iside)%T(acix(i),3)
          DVGlo(i) = two * (Tratto(sidestr)%velocity(acix(i)) - pg(npi)%vel(acix(i)))
        end do
        DVLoc(1) = sss(1) * DVGlo(1) + sss(2) * DVGlo(2)
        DVLoc(2) = nnn(1) * DVGlo(1) + nnn(2) * DVGlo(2)
!
!        IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(2)
        IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(3)
!
        if (strtype == 'fixe' .or. strtype == 'tapi') then
!
          BCLoc(1) = DVLoc(1) * IntWdV * Tratto(sidestr)%ShearCoeff
          BCLoc(2) = DVLoc(2) * IntWdV
          BCGlo(1) = sss(1) * BCLoc(1) + nnn(1) * BCLoc(2)
          BCGlo(2) = sss(2) * BCLoc(1) + nnn(2) * BCLoc(2)
!
          pg(npi)%var(1) = pg(npi)%var(1) + BCGlo(1)   
          pg(npi)%var(3) = pg(npi)%var(3) + BCGlo(2)   
!
!!!          unity = unity + IntWdV
!
        end if
      end do
    end if
!
!AA406
    endif
!
!
!!!    pg(npi)%uni = unity
!
!prova!
!    pg(npi)%var = floor(pg(npi)%var * azzeramento) / azzeramento
!    where (dabs(pg(npi)%var) < arrotondamento) pg(npi)%var = zero 
!prova!
!
  end do
!
!$omp end parallel do
!

!AA501b start
!Deallocations
  if (n_bodies > 0) then
     deallocate(dervel_mat)
     deallocate(unity_vec)
  endif
!AA501b end

return
end subroutine inter_SmoothVelo_2D
!---split

