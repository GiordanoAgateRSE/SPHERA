!cfile inter_SmoothVelo_3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : inter_SmoothVelo_3D
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
subroutine inter_SmoothVelo_3D
! ex inter43d
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
integer(4)       :: npi,i,j,npj,contj,npartint, ibdt, ibdp, ii
integer(4)       :: Ncbf, icbf, iface, facestr
double precision :: rhoi,rhoj,amassj,pesoj,  moddervel !unity
double precision :: IntWdV
character(4)     :: strtype
!
!.. Local Arrays ..
double precision,dimension(1:SPACEDIM) :: DVLoc, DVGlo, BCLoc, BCGlo
double precision,dimension(1:SPACEDIM) :: LocX
! quantita' step eq continuita'
double precision,dimension(3) :: dervel     ! dvel per eq standard, dvar usa vel smoot di monaghan

!AA501btest start
  double precision,dimension(:,:),allocatable :: dervel_mat
  double precision,dimension(:),allocatable :: unity_vec
  double precision :: unity

  if (n_bodies > 0) then  
!Allocations
     allocate(dervel_mat(nag,3))
     allocate(unity_vec(nag))
!Initializations
     dervel_mat = 0.
     unity_vec = 0.
     unity = 0.
  endif
!AA501btest end  

!
!.. Executable Statements ..
!!! !$omp private(ii,npi,unity,contj,npartint,npj,rhoi,rhoj,amassj,dervel,moddervel,pesoj) &
!

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

!AA501b modified
!AA501 sub
!$omp parallel do default(none) &
!$omp private(ii,npi,contj,npartint,npj,rhoi,rhoj,amassj,dervel,moddervel,pesoj) &
!$omp private(Ncbf,icbf,ibdt,ibdp,LocX,iface,facestr,strtype,DVGlo,DVLoc,i,j,IntWdV,BCLoc,BCGlo,unity) &
!$omp shared(nag,pg,Med,nPartIntorno,NMAXPARTJ,PartIntorno,PartKernel) &
!$omp shared(BoundaryDataPointer,BoundaryDataTab,BoundaryFace,Tratto,indarrayFlu,Array_Flu,esplosione,Domain,n_bodies,unity_vec,dervel_mat)
!
!.. loops on all the active particles in order to calculate the average velocity depending on 
!.. the smoothing coefficient of velocity ThetaV
!
!!!!!!  do npi = 1,nag
!!!!!!!
!!!!!!    if ( pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" .or. pg(npi)%state == "sol") cycle
!
!$$$$$
  do ii = 1,indarrayFlu
    npi = Array_Flu(ii)
!$$$$$
!
!* azzeramento quantita generali
    pg(npi)%var = zero
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
!!!      unity = unity + pesoj  
      pg(npi)%var(:) = pg(npi)%var(:) + dervel(:) * pesoj   
!
! calcolo velocita' per modello diffusivo
!      pg(npj)%veldif(:) = pg(npj)%veldif(:) + pg(npj)%vel(:) * pesoj
!
    end do
!
!...................................... 2011 mar 08
!.. update of Specific Internal Energy
      if (esplosione) pg(npi)%Envar = pg(npi)%Envar + (pg(npj)%IntEn - pg(npi)%IntEn) * pesoj
!.. end update of Specific Internal Energy
!.........................................
!
!AA406test
!

!AA501b start
      if (n_bodies > 0) then
         pg(npi)%var(:) = pg(npi)%var(:) + dervel_mat(npi,:)
         unity = unity + unity_vec(npi)
      endif
 
!AA601 sub start
! Impose boundary conditions at inlet and outlet sections (DB-SPH)
      if (Domain%tipo == "bsph") then
         call DBSPH_inlet_outlet(npi)
         else
!AA601 sub end
      
!.. contributo dei contorni
!
    ncbf = BoundaryDataPointer(1,npi)
    ibdt = BoundaryDataPointer(3,npi)
!
    if (Ncbf > 0) then        !Ci sono facce di contorno vicine
!
      do icbf = 1, Ncbf
!
        ibdp = ibdt + icbf - 1
        LocX(1:SPACEDIM) = BoundaryDataTab(ibdp)%LocXYZ(1:SPACEDIM)
        iface = BoundaryDataTab(ibdp)%CloBoNum
!
        facestr = BoundaryFace(iface)%stretch
        strtype = Tratto(facestr)%tipo
!
        if (strtype == 'sour' .or. strtype == 'velo' .or. strtype == 'flow') then
          pg(npi)%var(:) = zero   
!!!          unity = one
          exit  
        end if
!       
!       do i = 1, SPACEDIM
!         rrr(i) = BoundaryFace(iface)%T(i,1)
!         sss(i) = BoundaryFace(iface)%T(i,2)
!         nnn(i) = BoundaryFace(iface)%T(i,3)
!         DVGlo(i) = two * (Tratto(facestr)%velocity(i) - pg(npi)%vel(i))
!       end do
!       DVLoc(1) = rrr(1) * DVGlo(1) + rrr(2) * DVGlo(2) + rrr(3) * DVGlo(3)
!       DVLoc(2) = sss(1) * DVGlo(1) + sss(2) * DVGlo(2) + sss(3) * DVGlo(3)
!       DVLoc(3) = nnn(1) * DVGlo(1) + nnn(2) * DVGlo(2) + nnn(3) * DVGlo(3)
!       
        DVGlo(:) = two * (Tratto(facestr)%velocity(:) - pg(npi)%vel(:))
        do i = 1,SPACEDIM
          DVLoc(i) = zero
          do j = 1,SPACEDIM
            DVLoc(i) = DVLoc(i) + BoundaryFace(iface)%T(j,i) * DVGlo(j)
          end do
        end do
!
        IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(2)
!
        if (strtype == 'fixe' .or. strtype == 'tapi') then
!
          BCLoc(1) = DVLoc(1) * IntWdV * Tratto(facestr)%ShearCoeff
          BCLoc(2) = DVLoc(2) * IntWdV * Tratto(facestr)%ShearCoeff
          BCLoc(3) = DVLoc(3) * IntWdV
!          BCGlo(:) = rrr(:) * BCLoc(1) + sss(:) * BCLoc(2) + nnn(:) * BCLoc(3)
          do i = 1, SPACEDIM
            BCGlo(i) = zero
            do j = 1,SPACEDIM
              BCGlo(i) = BCGlo(i) + BoundaryFace(iface)%T(i,j) * BCLoc(j)
            end do
          end do
!
          pg(npi)%var(:) = pg(npi)%var(:) + BCGlo(:)   
!
        end if
      end do
    end if

!AA601
     endif   
    
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
end subroutine inter_SmoothVelo_3D
!---split

